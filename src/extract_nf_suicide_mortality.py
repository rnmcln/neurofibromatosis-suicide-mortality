#!/usr/bin/env python3
"""Count documented neurofibromatosis mentions among U.S. suicide deaths.

Data source: NCHS Mortality Multiple Cause public-use files, 1999-2024.
Case definition:
  * underlying cause is suicide: U03, X60-X84, or Y87.0; and
  * record-axis multiple causes contain Q85.0 (stored as Q850).

The script reports an exact pooled count only when it is at least 10. Annual
counts from 1 through 9 are written as "Suppressed" to comply with the NCHS
data-use restrictions. It does not distinguish NF1 from NF2/schwannomatosis.

Requirements: Python 3 and curl. No third-party Python packages are needed.
"""

from __future__ import annotations

import argparse
import csv
import contextlib
import subprocess
import sys
import tempfile
import zipfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path


BASE_URL = (
    "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/"
    "Datasets/DVS/mortality/mort{year}us.zip"
)
MIN_REPORTABLE = 10


@dataclass(frozen=True)
class YearResult:
    year: int
    documented_deaths: int

    @property
    def reportable_value(self) -> str:
        if self.documented_deaths < MIN_REPORTABLE:
            return "Suppressed"
        return str(self.documented_deaths)


@dataclass(frozen=True)
class RecordLayout:
    underlying_cause: slice
    record_axis: slice


def record_layout(year: int) -> RecordLayout:
    """Return the fixed-width locations used for an annual mortality file.

    NCHS redesigned the public-use record beginning with data year 2003.
    Locations below are zero-based Python slices converted from the one-based
    locations in the annual NCHS documentation.
    """
    if 1999 <= year <= 2002:
        return RecordLayout(slice(141, 145), slice(340, 440))
    if 2003 <= year <= 2024:
        return RecordLayout(slice(145, 149), slice(343, 443))
    raise ValueError(f"Unsupported mortality file year: {year}")


def is_suicide_underlying(code: bytes) -> bool:
    """Return whether the four-position NCHS underlying code is suicide."""
    text = code.decode("ascii", errors="ignore").strip().upper()
    if text.startswith("U03") or text == "Y870":
        return True
    if len(text) < 3 or not text.startswith("X"):
        return False
    try:
        category = int(text[1:3])
    except ValueError:
        return False
    return 60 <= category <= 84


def record_axis_has_nf(record_axis: bytes) -> bool:
    """Return whether any of the 20 five-position record-axis codes is Q85.0."""
    return any(
        record_axis[offset : offset + 5].strip() == b"Q850"
        for offset in range(0, 100, 5)
    )


def count_records(records: Iterable[bytes], year: int) -> int:
    layout = record_layout(year)
    count = 0
    for record in records:
        if len(record) < layout.record_axis.stop or record[19:20] == b"4":
            continue
        if not is_suicide_underlying(record[layout.underlying_cause]):
            continue
        if record_axis_has_nf(record[layout.record_axis]):
            count += 1
    return count


def count_year(archive_path: Path, year: int) -> YearResult:
    try:
        with zipfile.ZipFile(archive_path) as archive:
            members = [name for name in archive.namelist() if not name.endswith("/")]
            if len(members) != 1:
                raise RuntimeError(
                    f"Expected one mortality data file in {archive_path}, found {members}"
                )
            with archive.open(members[0]) as records:
                count = count_records(records, year)
    except NotImplementedError:
        # Several recent NCHS archives use Deflate64. Python's zipfile does not
        # support it, but the system unzip bundled with macOS and many Linux
        # distributions does. Stream it so no multi-gigabyte text file is made.
        process = subprocess.Popen(
            ["unzip", "-p", str(archive_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        assert process.stdout is not None
        try:
            count = count_records(process.stdout, year)
        finally:
            process.stdout.close()
        stderr = process.stderr.read() if process.stderr else b""
        return_code = process.wait()
        with contextlib.suppress(Exception):
            if process.stderr:
                process.stderr.close()
        if return_code:
            detail = stderr.decode("utf-8", errors="replace").strip()
            raise RuntimeError(f"unzip failed for {year}: {detail}")
    return YearResult(year=year, documented_deaths=count)


def download(year: int, destination: Path) -> None:
    url = BASE_URL.format(year=year)
    command = [
        "curl",
        "-L",
        "--fail",
        "--retry",
        "3",
        "--retry-delay",
        "2",
        "--silent",
        "--show-error",
        "--output",
        str(destination),
        url,
    ]
    result = subprocess.run(command, check=False)
    if result.returncode:
        raise RuntimeError(f"Download failed for {year}: curl exit {result.returncode}")


def write_results(results: list[YearResult], output_path: Path) -> None:
    pooled = sum(result.documented_deaths for result in results)
    pooled_value = str(pooled) if pooled >= MIN_REPORTABLE else "Suppressed"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["period", "documented_deaths", "note"])
        for result in results:
            note = (
                "0-9 deaths; suppressed under NCHS rules"
                if result.reportable_value == "Suppressed"
                else ""
            )
            writer.writerow([result.year, result.reportable_value, note])
        pooled_note = (
            "Pooled exact count across non-overlapping years"
            if pooled >= MIN_REPORTABLE
            else "Pooled count is 0-9; suppressed under NCHS rules"
        )
        writer.writerow(
            [f"{results[0].year}-{results[-1].year}", pooled_value, pooled_note]
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-year", type=int, default=1999)
    parser.add_argument("--end-year", type=int, default=2024)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("nf_suicide_documented_lower_bound.csv"),
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help="keep/reuse downloaded ZIP files in this directory",
    )
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

    results: list[YearResult] = []
    try:
        for year in range(args.start_year, args.end_year + 1):
            archive_path = cache_dir / f"mort{year}us.zip"
            if not archive_path.exists():
                print(f"Downloading {year}...", file=sys.stderr, flush=True)
                download(year, archive_path)
            print(f"Counting {year}...", file=sys.stderr, flush=True)
            results.append(count_year(archive_path, year))
        write_results(results, args.output)
    finally:
        if temporary is not None:
            temporary.cleanup()

    print(f"Wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
