#!/usr/bin/env python3
"""Analyze suicide among U.S. deaths with neurofibromatosis mentioned.

This program reads the annual NCHS Mortality Multiple Cause public-use files
for 1999-2024. It identifies U.S. resident death certificates with Q85.0
(neurofibromatosis, nonmalignant) in the record-axis multiple-cause field and
calculates indirectly standardized proportionate mortality ratios (PMRs) for
suicide as the underlying cause.

Primary suicide definition: U03, X60-X84, or Y87.0.
Sensitivity definition: primary suicide codes plus Y10-Y34 or Y87.2
(undetermined intent).

Expected deaths are calculated by applying stratum-specific proportions of
suicide among all U.S. resident deaths to the distribution of deaths with a
Q85.0 mention. Primary strata are data year, age group, sex, and race (White
versus all other races), matching the principal adjustment variables used in
the earlier U.S. death-certificate analysis as closely as the ICD-10 files
permit. A sensitivity analysis omits race from standardization.

Output tables apply the CDC WONDER rule against presenting counts of 0-9.
The program never writes record-level data.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import json
import math
import statistics
import subprocess
import sys
import tempfile
import zipfile
from collections import Counter
from collections.abc import Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import asdict, dataclass, field
from pathlib import Path


BASE_URL = (
    "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/"
    "Datasets/DVS/mortality/mort{year}us.zip"
)
MIN_REPORTABLE = 10
AGE_GROUPS = ("<10", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70+")
NF2_PROXY_PREFIXES = ("H903", "H904", "H905", "D333", "D320", "D321")


@dataclass(frozen=True)
class RecordLayout:
    underlying_cause: slice
    record_axis: slice
    sex_index: int
    detail_age: slice
    race_code: slice


def record_layout(year: int) -> RecordLayout:
    """Return zero-based slices converted from annual NCHS documentation."""
    if 1999 <= year <= 2002:
        return RecordLayout(
            underlying_cause=slice(141, 145),
            record_axis=slice(340, 440),
            sex_index=58,
            detail_age=slice(63, 66),
            race_code=slice(59, 61),
        )
    if 2003 <= year <= 2020:
        return RecordLayout(
            underlying_cause=slice(145, 149),
            record_axis=slice(343, 443),
            sex_index=68,
            detail_age=slice(69, 73),
            race_code=slice(444, 446),
        )
    if 2021 <= year <= 2024:
        return RecordLayout(
            underlying_cause=slice(145, 149),
            record_axis=slice(343, 443),
            sex_index=68,
            detail_age=slice(69, 73),
            # Bridged-race fields ended after 2020. Race Recode 40 uses the
            # reported single/multiple-race values consistently from 2021.
            race_code=slice(488, 490),
        )
    raise ValueError(f"Unsupported mortality file year: {year}")


@dataclass
class YearStats:
    year: int
    all_deaths: Counter[tuple[str, str, str]] = field(default_factory=Counter)
    all_suicides: Counter[tuple[str, str, str]] = field(default_factory=Counter)
    all_undetermined: Counter[tuple[str, str, str]] = field(default_factory=Counter)
    nf_deaths: Counter[tuple[str, str, str]] = field(default_factory=Counter)
    nf_deaths_no_proxy: Counter[tuple[str, str, str]] = field(default_factory=Counter)
    nf_suicides: Counter[tuple[str, str, str]] = field(default_factory=Counter)
    nf_undetermined: Counter[tuple[str, str, str]] = field(default_factory=Counter)
    nf_suicides_no_proxy: Counter[tuple[str, str, str]] = field(default_factory=Counter)
    nf_undetermined_no_proxy: Counter[tuple[str, str, str]] = field(default_factory=Counter)
    nf_underlying: int = 0
    nf_ages: list[int] = field(default_factory=list)
    nf_sex: Counter[str] = field(default_factory=Counter)
    nf_age_group: Counter[str] = field(default_factory=Counter)
    nf_race: Counter[str] = field(default_factory=Counter)
    resident_records: int = 0


@dataclass(frozen=True)
class PMRResult:
    analysis: str
    nf_deaths: int
    observed: int
    expected: float
    pmr: float
    ci_low: float
    ci_high: float
    adjustment: str


def clean_code(raw: bytes) -> str:
    return raw.decode("ascii", errors="ignore").strip().upper()


def is_suicide(code: str) -> bool:
    if code.startswith("U03") or code == "Y870":
        return True
    if len(code) < 3 or code[0] != "X":
        return False
    try:
        return 60 <= int(code[1:3]) <= 84
    except ValueError:
        return False


def is_undetermined(code: str) -> bool:
    if code == "Y872":
        return True
    if len(code) < 3 or code[0] != "Y":
        return False
    try:
        return 10 <= int(code[1:3]) <= 34
    except ValueError:
        return False


def record_axis_codes(raw: bytes) -> tuple[str, ...]:
    return tuple(
        code
        for offset in range(0, 100, 5)
        if (code := clean_code(raw[offset : offset + 5]))
    )


def parse_age(record: bytes, year: int, layout: RecordLayout) -> int | None:
    raw = clean_code(record[layout.detail_age])
    if not raw:
        return None
    unit = raw[0]
    try:
        value = int(raw[1:])
    except ValueError:
        return None

    if year <= 2002:
        if unit == "0":
            return value
        if unit == "1":
            return 100 + value
        if unit in {"2", "3", "4", "5", "6"}:
            return 0
        return None

    if unit == "1":
        return value
    if unit in {"2", "3", "4", "5", "6"}:
        return 0
    return None


def age_group(age: int | None) -> str:
    if age is None:
        return "Unknown"
    if age < 10:
        return "<10"
    if age >= 70:
        return "70+"
    lower = (age // 10) * 10
    return f"{lower}-{lower + 9}"


def stratum(record: bytes, year: int, layout: RecordLayout) -> tuple[str, str, str]:
    age = parse_age(record, year, layout)
    raw_sex = clean_code(record[layout.sex_index : layout.sex_index + 1])
    sex = {"1": "M", "2": "F", "M": "M", "F": "F"}.get(raw_sex, "Unknown")
    raw_race = clean_code(record[layout.race_code])
    race = "White" if raw_race == "01" else "Other"
    return age_group(age), sex, race


def iter_records(archive_path: Path) -> Iterator[bytes]:
    try:
        with zipfile.ZipFile(archive_path) as archive:
            members = [name for name in archive.namelist() if not name.endswith("/")]
            if len(members) != 1:
                raise RuntimeError(
                    f"Expected one mortality file in {archive_path}, found {members}"
                )
            with archive.open(members[0]) as stream:
                yield from stream
            return
    except NotImplementedError:
        pass

    process = subprocess.Popen(
        ["unzip", "-p", str(archive_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert process.stdout is not None
    try:
        yield from process.stdout
    finally:
        process.stdout.close()
    stderr = process.stderr.read() if process.stderr else b""
    return_code = process.wait()
    with contextlib.suppress(Exception):
        if process.stderr:
            process.stderr.close()
    if return_code:
        detail = stderr.decode("utf-8", errors="replace").strip()
        raise RuntimeError(f"unzip failed for {archive_path}: {detail}")


def scan_year(task: tuple[Path, int]) -> YearStats:
    archive_path, year = task
    layout = record_layout(year)
    stats = YearStats(year=year)

    for record in iter_records(archive_path):
        if len(record) < max(layout.record_axis.stop, layout.race_code.stop):
            continue
        if record[19:20] == b"4":
            continue

        stats.resident_records += 1
        key = stratum(record, year, layout)
        underlying = clean_code(record[layout.underlying_cause])
        suicide = is_suicide(underlying)
        undetermined = is_undetermined(underlying)
        stats.all_deaths[key] += 1
        if suicide:
            stats.all_suicides[key] += 1
        if undetermined:
            stats.all_undetermined[key] += 1

        codes = record_axis_codes(record[layout.record_axis])
        if "Q850" not in codes:
            continue

        no_nf2_proxy = not any(
            code.startswith(NF2_PROXY_PREFIXES) for code in codes if code != "Q850"
        )
        stats.nf_deaths[key] += 1
        if no_nf2_proxy:
            stats.nf_deaths_no_proxy[key] += 1
        if underlying == "Q850":
            stats.nf_underlying += 1

        age = parse_age(record, year, layout)
        if age is not None:
            stats.nf_ages.append(age)
        stats.nf_sex[key[1]] += 1
        stats.nf_age_group[key[0]] += 1
        stats.nf_race[key[2]] += 1

        if suicide:
            stats.nf_suicides[key] += 1
            if no_nf2_proxy:
                stats.nf_suicides_no_proxy[key] += 1
        if undetermined:
            stats.nf_undetermined[key] += 1
            if no_nf2_proxy:
                stats.nf_undetermined_no_proxy[key] += 1

    return stats


def download(year: int, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".zip.part")
    command = [
        "curl",
        "-L",
        "--fail",
        "--retry",
        "4",
        "--retry-delay",
        "2",
        "--silent",
        "--show-error",
        "--output",
        str(temporary),
        BASE_URL.format(year=year),
    ]
    result = subprocess.run(command, check=False)
    if result.returncode:
        with contextlib.suppress(FileNotFoundError):
            temporary.unlink()
        raise RuntimeError(f"Download failed for {year}: curl exit {result.returncode}")
    temporary.replace(destination)


def download_task(task: tuple[int, Path]) -> None:
    year, destination = task
    print(f"Downloading {year}...", file=sys.stderr, flush=True)
    download(year, destination)


def merge_counter(
    yearly: Iterable[YearStats], attribute: str, *, include_year: bool
) -> Counter[tuple[str, ...]]:
    combined: Counter[tuple[str, ...]] = Counter()
    for stats in yearly:
        counter: Counter[tuple[str, str, str]] = getattr(stats, attribute)
        for demographic_key, count in counter.items():
            key = (str(stats.year), *demographic_key) if include_year else demographic_key
            combined[key] += count
    return combined


def collapse_race(counter: Counter[tuple[str, ...]]) -> Counter[tuple[str, ...]]:
    collapsed: Counter[tuple[str, ...]] = Counter()
    for key, count in counter.items():
        collapsed[key[:-1]] += count
    return collapsed


def poisson_cdf(k: int, rate: float) -> float:
    if rate <= 0:
        return 1.0
    term = math.exp(-rate)
    total = term
    for value in range(1, k + 1):
        term *= rate / value
        total += term
    return min(1.0, total)


def poisson_sf_from(k: int, rate: float) -> float:
    """Return P(X >= k) for a Poisson variable."""
    return max(0.0, 1.0 - poisson_cdf(k - 1, rate))


def poisson_exact_interval(observed: int, alpha: float = 0.05) -> tuple[float, float]:
    """Garwood exact interval obtained by numerical inversion of Poisson tails."""
    if observed < 0:
        raise ValueError("observed must be nonnegative")
    if observed == 0:
        return 0.0, -math.log(alpha / 2)

    target = alpha / 2
    low_left, low_right = 0.0, float(observed)
    for _ in range(100):
        midpoint = (low_left + low_right) / 2
        if poisson_sf_from(observed, midpoint) < target:
            low_left = midpoint
        else:
            low_right = midpoint
    lower = (low_left + low_right) / 2

    high_left = float(observed)
    high_right = max(float(observed + 10), float(observed * 2 + 10))
    while poisson_cdf(observed, high_right) > target:
        high_right *= 2
    for _ in range(100):
        midpoint = (high_left + high_right) / 2
        if poisson_cdf(observed, midpoint) > target:
            high_left = midpoint
        else:
            high_right = midpoint
    upper = (high_left + high_right) / 2
    return lower, upper


def calculate_pmr(
    *,
    analysis: str,
    all_deaths: Counter[tuple[str, ...]],
    all_outcomes: Counter[tuple[str, ...]],
    nf_deaths: Counter[tuple[str, ...]],
    nf_outcomes: Counter[tuple[str, ...]],
    adjustment: str,
) -> PMRResult:
    expected = 0.0
    for key, nf_count in nf_deaths.items():
        denominator = all_deaths[key]
        if denominator:
            expected += nf_count * all_outcomes[key] / denominator

    observed = sum(nf_outcomes.values())
    if expected <= 0:
        raise RuntimeError(f"Nonpositive expected count for {analysis}")
    lower_count, upper_count = poisson_exact_interval(observed)
    return PMRResult(
        analysis=analysis,
        nf_deaths=sum(nf_deaths.values()),
        observed=observed,
        expected=expected,
        pmr=observed / expected,
        ci_low=lower_count / expected,
        ci_high=upper_count / expected,
        adjustment=adjustment,
    )


def quantile(values: list[int], probability: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return float(ordered[lower])
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def reportable_count(value: int) -> str:
    return str(value) if value >= MIN_REPORTABLE else "Suppressed"


def public_result(result: PMRResult) -> dict[str, object]:
    values = asdict(result)
    if result.observed < MIN_REPORTABLE:
        values.update(
            observed="Suppressed",
            pmr="Suppressed",
            ci_low="Suppressed",
            ci_high="Suppressed",
        )
    return values


def write_pmr_results(results: list[PMRResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "analysis",
                "nf_deaths",
                "observed",
                "expected",
                "pmr",
                "ci_low",
                "ci_high",
                "adjustment",
            ]
        )
        for result in results:
            if result.observed < MIN_REPORTABLE:
                writer.writerow(
                    [
                        result.analysis,
                        result.nf_deaths,
                        "Suppressed",
                        f"{result.expected:.2f}",
                        "Suppressed",
                        "Suppressed",
                        "Suppressed",
                        result.adjustment,
                    ]
                )
                continue
            writer.writerow(
                [
                    result.analysis,
                    result.nf_deaths,
                    result.observed,
                    f"{result.expected:.2f}",
                    f"{result.pmr:.3f}",
                    f"{result.ci_low:.3f}",
                    f"{result.ci_high:.3f}",
                    result.adjustment,
                ]
            )


def write_descriptives(yearly: list[YearStats], path: Path) -> None:
    nf_total = sum(sum(stats.nf_deaths.values()) for stats in yearly)
    nf_underlying = sum(stats.nf_underlying for stats in yearly)
    ages = [age for stats in yearly for age in stats.nf_ages]
    sex = sum((stats.nf_sex for stats in yearly), Counter())
    age_groups = sum((stats.nf_age_group for stats in yearly), Counter())
    race = sum((stats.nf_race for stats in yearly), Counter())

    rows: list[tuple[str, str, str, str]] = [
        ("Overall", "Death certificates mentioning Q85.0", str(nf_total), "100.0"),
        (
            "Overall",
            "Q85.0 selected as underlying cause",
            reportable_count(nf_underlying),
            f"{100 * nf_underlying / nf_total:.1f}" if nf_underlying >= MIN_REPORTABLE else "Suppressed",
        ),
        ("Age", "Median age at death", f"{statistics.median(ages):.1f}", ""),
        ("Age", "Age at death Q1", f"{quantile(ages, 0.25):.1f}", ""),
        ("Age", "Age at death Q3", f"{quantile(ages, 0.75):.1f}", ""),
    ]
    for category in ("M", "F", "Unknown"):
        count = sex[category]
        if count:
            rows.append(
                (
                    "Sex",
                    category,
                    reportable_count(count),
                    f"{100 * count / nf_total:.1f}" if count >= MIN_REPORTABLE else "Suppressed",
                )
            )
    for category in (*AGE_GROUPS, "Unknown"):
        count = age_groups[category]
        if count:
            rows.append(
                (
                    "Age group",
                    category,
                    reportable_count(count),
                    f"{100 * count / nf_total:.1f}" if count >= MIN_REPORTABLE else "Suppressed",
                )
            )
    for category in ("White", "Other"):
        count = race[category]
        if count:
            rows.append(
                (
                    "Race",
                    category,
                    reportable_count(count),
                    f"{100 * count / nf_total:.1f}" if count >= MIN_REPORTABLE else "Suppressed",
                )
            )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["domain", "category", "count_or_value", "percent"])
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-year", type=int, default=1999)
    parser.add_argument("--end-year", type=int, default=2024)
    parser.add_argument("--cache-dir", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--workers", type=int, default=4)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.start_year < 1999 or args.end_year > 2024:
        raise SystemExit("Supported final-data years are 1999 through 2024")
    if args.start_year > args.end_year:
        raise SystemExit("--start-year must not exceed --end-year")

    temporary: tempfile.TemporaryDirectory[str] | None = None
    if args.cache_dir:
        cache_dir = args.cache_dir
        cache_dir.mkdir(parents=True, exist_ok=True)
    else:
        temporary = tempfile.TemporaryDirectory(prefix="nchs-mortality-")
        cache_dir = Path(temporary.name)

    years = list(range(args.start_year, args.end_year + 1))
    archives: list[tuple[Path, int]] = []
    try:
        downloads: list[tuple[int, Path]] = []
        for year in years:
            archive = cache_dir / f"mort{year}us.zip"
            if not archive.exists():
                downloads.append((year, archive))
            archives.append((archive, year))
        if downloads:
            with ThreadPoolExecutor(max_workers=min(4, max(1, args.workers))) as executor:
                list(executor.map(download_task, downloads))

        print(f"Scanning {len(years)} annual files...", file=sys.stderr, flush=True)
        with ProcessPoolExecutor(max_workers=max(1, args.workers)) as executor:
            yearly = list(executor.map(scan_year, archives))
        yearly.sort(key=lambda value: value.year)

        all_deaths = merge_counter(yearly, "all_deaths", include_year=True)
        all_suicides = merge_counter(yearly, "all_suicides", include_year=True)
        all_undetermined = merge_counter(yearly, "all_undetermined", include_year=True)
        nf_deaths = merge_counter(yearly, "nf_deaths", include_year=True)
        nf_deaths_no_proxy = merge_counter(yearly, "nf_deaths_no_proxy", include_year=True)
        nf_suicides = merge_counter(yearly, "nf_suicides", include_year=True)
        nf_undetermined = merge_counter(yearly, "nf_undetermined", include_year=True)
        nf_suicides_no_proxy = merge_counter(
            yearly, "nf_suicides_no_proxy", include_year=True
        )
        nf_undetermined_no_proxy = merge_counter(
            yearly, "nf_undetermined_no_proxy", include_year=True
        )

        all_combined = all_suicides + all_undetermined
        nf_combined = nf_suicides + nf_undetermined
        nf_combined_no_proxy = nf_suicides_no_proxy + nf_undetermined_no_proxy
        adjustment_full = "year, age group, sex, and race (White/other)"
        adjustment_no_race = "year, age group, and sex"

        results = [
            calculate_pmr(
                analysis="Primary: suicide underlying cause",
                all_deaths=all_deaths,
                all_outcomes=all_suicides,
                nf_deaths=nf_deaths,
                nf_outcomes=nf_suicides,
                adjustment=adjustment_full,
            ),
            calculate_pmr(
                analysis="Sensitivity: suicide plus undetermined intent",
                all_deaths=all_deaths,
                all_outcomes=all_combined,
                nf_deaths=nf_deaths,
                nf_outcomes=nf_combined,
                adjustment=adjustment_full,
            ),
            calculate_pmr(
                analysis="Sensitivity: suicide, excluding NF2 proxy mentions",
                all_deaths=all_deaths,
                all_outcomes=all_suicides,
                nf_deaths=nf_deaths_no_proxy,
                nf_outcomes=nf_suicides_no_proxy,
                adjustment=adjustment_full,
            ),
            calculate_pmr(
                analysis="Sensitivity: suicide plus undetermined intent, excluding NF2 proxy mentions",
                all_deaths=all_deaths,
                all_outcomes=all_combined,
                nf_deaths=nf_deaths_no_proxy,
                nf_outcomes=nf_combined_no_proxy,
                adjustment=adjustment_full,
            ),
            calculate_pmr(
                analysis="Sensitivity: suicide, adjustment without race",
                all_deaths=collapse_race(all_deaths),
                all_outcomes=collapse_race(all_suicides),
                nf_deaths=collapse_race(nf_deaths),
                nf_outcomes=collapse_race(nf_suicides),
                adjustment=adjustment_no_race,
            ),
        ]

        args.output_dir.mkdir(parents=True, exist_ok=True)
        # Only release analyses that cannot reveal a complementary cell of
        # fewer than 10 deaths. In particular, publishing both the primary
        # count and the combined suicide/undetermined-intent count would make
        # the suppressed undetermined-intent count calculable by subtraction.
        public_results = [results[0], results[4]]
        write_pmr_results(
            public_results, args.output_dir / "nf_suicide_pmr_results.csv"
        )
        write_descriptives(yearly, args.output_dir / "nf_mortality_descriptives.csv")

        summary = {
            "period": f"{args.start_year}-{args.end_year}",
            "resident_deaths_scanned": sum(stats.resident_records for stats in yearly),
            "nf_deaths": sum(sum(stats.nf_deaths.values()) for stats in yearly),
            "nf_underlying": sum(stats.nf_underlying for stats in yearly),
            "primary_observed": reportable_count(
                sum(sum(stats.nf_suicides.values()) for stats in yearly)
            ),
            "results": [public_result(result) for result in public_results],
            "sensitivity_note": (
                "Including undetermined intent and excluding records with NF2 proxy "
                "mentions did not materially change the interpretation. Numeric results "
                "are not released because they would permit complementary disclosure of "
                "a death count below 10."
            ),
            "annual_nf_deaths": {
                str(stats.year): sum(stats.nf_deaths.values()) for stats in yearly
            },
            "annual_nf_suicides": {
                str(stats.year): (
                    sum(stats.nf_suicides.values())
                    if sum(stats.nf_suicides.values()) >= MIN_REPORTABLE
                    else "Suppressed"
                )
                for stats in yearly
            },
        }
        with (args.output_dir / "nf_suicide_analysis_summary.json").open(
            "w", encoding="utf-8"
        ) as stream:
            json.dump(summary, stream, indent=2)
            stream.write("\n")

        print(
            f"Completed: {summary['nf_deaths']} Q85.0-mention deaths; "
            f"{summary['primary_observed']} suicides",
            file=sys.stderr,
        )
    finally:
        if temporary is not None:
            temporary.cleanup()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
