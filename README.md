# Neurofibromatosis Mentions and Suicide Mortality in US Death Records

This repository contains reproducible analysis code and aggregate outputs for
an examination of suicide as the underlying cause of death among United States
death certificates that mention ICD-10 code Q85.0.

## Summary

The National Center for Health Statistics Mortality Multiple Cause public-use
files were examined for 1999-2024. Among 69,714,059 resident death records,
8,160 mentioned Q85.0. Suicide was the underlying cause in 14 records, compared
with 375.58 expected after indirect standardization by calendar year, age
group, sex, and race. The proportionate mortality ratio was 0.037 (95% exact
Poisson confidence interval, 0.020-0.063).

This result should not be interpreted as evidence that suicide risk is low
among people with neurofibromatosis. Inclusion depended on neurofibromatosis
being recorded on the death certificate. The marked deficit therefore
indicates that analyses beginning with a Q85.0 mention can differentially miss
suicide deaths. A clinically validated neurofibromatosis cohort linked to
mortality records is required to estimate suicide incidence or relative risk.

## Data source

National Center for Health Statistics Mortality Multiple Cause public-use
files, 1999-2024:

https://www.cdc.gov/nchs/data_access/vitalstatsonline.htm

The programs download the annual public-use archives directly from the
National Center for Health Statistics. No record-level data are stored in this
repository.

## Analysis

The primary outcome is suicide as the underlying cause of death, defined by
ICD-10 codes U03, X60-X84, or Y87.0. Neurofibromatosis-mention deaths are
identified by Q85.0 in the record-axis multiple-cause fields. Expected deaths
are calculated by applying stratum-specific suicide proportions among all
resident deaths to the distribution of Q85.0-mention deaths. Primary strata
are calendar year, age group, sex, and race (White or other).

The code also performs prespecified sensitivity analyses, including
standardization without race, inclusion of undetermined intent, and exclusion
of records containing diagnostic proxies that could indicate NF2-related
schwannomatosis. Small counts are suppressed, and complementary disclosure is
prevented in the public outputs.

## Repository contents

- `src/analyze_nf_suicide_pmr.py`: primary analysis and sensitivity analyses.
- `src/extract_nf_suicide_mortality.py`: documented lower-bound extraction.
- `results/nf_suicide_pmr_results.csv`: public PMR estimates.
- `results/nf_mortality_descriptives.csv`: aggregate descriptive statistics.
- `results/nf_suicide_analysis_summary.json`: machine-readable summary.
- `results/nf_suicide_documented_lower_bound.csv`: suppression-compliant
  annual output and pooled documented count.

## Reproduction

Python 3, `curl`, and `unzip` are required. No third-party Python packages are
needed.

```bash
python3 src/analyze_nf_suicide_pmr.py \
  --cache-dir mortality_cache \
  --output-dir results

python3 src/extract_nf_suicide_mortality.py \
  --cache-dir mortality_cache \
  --output results/nf_suicide_documented_lower_bound.csv
```

The complete analysis downloads and scans 26 large annual mortality files and
may require substantial time and storage.

## Author

Aaron Lawson McLean  
ORCID: https://orcid.org/0000-0001-5528-6905

## Citation

McLean, A. L. (2026). *Neurofibromatosis mentions and suicide mortality in US
death records* (Version 1.0.0) [Computer software]. Zenodo.
https://doi.org/10.5281/zenodo.24444191

## License

The software is released under the MIT License. The included aggregate outputs
may be reused with attribution.
