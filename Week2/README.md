# RBIF-111 Week 2 Homework

**Author:** Coddie Curry  
**Repository:** [github.com/bio-light-dev/RBIF111-Bioinformatics/tree/main/Week2](https://github.com/bio-light-dev/RBIF111-Bioinformatics/tree/main/Week2)  
**Files:** `HW1-Week-2.Rmd` (source), `HW1-Week-2.html` (knitted report)

Week 2 covers simulation-based statistics in R: how estimators behave across sample sizes, t-tests compared with brute-force resampling, the binomial distribution, normality testing, choosing a hypothesis test, and statistical power. Questions 5 and 6 apply these methods to a clinical table queried with SQL from the DuckDB/Parquet repository built in Week 1.

## What was completed

| Question | Task | Key result |
| --- | --- | --- |
| 1 | Simulated 1,000 samples from N(0, 1) for each n = 2–20; compared sd estimators with n − 1, n and n − 2 denominators; plotted their means and spreads | All three approach the true sd of 1 as n grows. The n − 1 estimator still slightly underestimates (≈ 0.8 at n = 2, ≈ 0.99 at n = 20), n underestimates most, n − 2 overestimates. Spread falls with n for all three |
| 2 | 10,000 pairs of samples (n = 5 and 7) from N(2, 3); difference in means, two-sided t-test p-values and brute-force rank-based p-values; plots a–f | Both samples come from the same distribution, so the null is true in every pair; t-test and brute-force p-values compared side by side |
| 3 | Two processes with the same mean but very different shapes: N(0, 1) vs an equal mix of N(−10, 1) and N(10, 1); 10,000 repeated Welch t-tests | p-values roughly uniform, about 5% below 0.05. The t-test only compares means, so it can't tell these processes apart |
| 4a | Binomial probability by hand: diagnostic test with sensitivity 0.8, 15 true-positive patients, P(exactly 10 positive) | C(15, 10) = 3003; 3003 × 0.8¹⁰ × 0.2⁵ ≈ 0.1032 (handwritten work included) |
| 4b | Confirmed with `dbinom(10, 15, 0.8)` | 0.1031823, matches the hand calculation |
| 4c | P(at least 10 positive) with `1 - pbinom(9, 15, 0.8)` and with `sum(dbinom(10:15, 15, 0.8))` | Both ≈ 0.9389; `all.equal()` returns `TRUE` |
| 4d | 10,000 simulated samples with `rbinom()`; proportion with ≥ 10 positives (and exactly 10) vs the analytic answers | Simulated proportions close to 0.9389 and 0.1032 |
| 4e | Smallest k with P(X ≥ k) ≤ 0.05, using `qbinom(0.95, 15, 0.8) + 1` | k = 15: all 15 testing positive happens only ≈ 3.5% of the time if sensitivity is 0.8. Because X is discrete, the actual tail is 0.035, not 0.05 |
| 5 setup | Reconnected to the Week 1 DuckDB repository and queried `UCI_HEALTH_UCI519_Heart_Failure_Clinical_Records` with SQL | 299 heart failure patients; outcome `DEATH_EVENT` (≈ 203 survived, ≈ 96 died) |
| 5a | Shapiro-Wilk test and QQ plot for six numeric features (age, CPK, ejection fraction, platelets, serum creatinine, serum sodium) within each outcome group | Most features are not normal; CPK and serum creatinine are strongly right-skewed |
| 5b | Boxplot of every feature, coloured by outcome (ggplot2, log10 y-axis because features differ by orders of magnitude) | Visible separation for serum creatinine and ejection fraction; heavy overlap for CPK and platelets |
| 5c | Compared each feature between groups, choosing Welch's t-test (both groups normal) or Wilcoxon rank-sum (otherwise) from the 5a results; Benjamini-Hochberg adjustment across the six tests | Results table with test used, group medians, raw and adjusted p-values |
| 6a | Highest-variance column among features that overlap between groups (not significant in 5c); 1,000 resamplings at n = 10, 15 and 20 per group, test chosen as in 5c, p-values saved in `Pvals` | Platelets: p-value histograms roughly flat at every sample size, with about 5% below 0.05 |
| 6b | Same process for the lowest-variance column among features that separate the groups (significant in 5c) | Serum creatinine: p-values pile up near 0, and the fraction below 0.05 (power) rises with n |
| 6c | Analytic power with `power.t.test()` using the observed difference in means and pooled sd, compared with simulated power; n needed for 80% power on the high-overlap column | Analytic and simulated power follow the same trend. Platelets would need far more patients per group than the dataset has |

## Data and setup notes

- **Data source:** the UCI Heart Failure Clinical Records table from the course repoquet repository (`profile = "quick"`, `~/repoquet_quick_example`). It was chosen over `NHANES_DEMO`, which has only a few truly numeric features (most of its numeric columns are category codes or survey weights).
- **Reconnecting instead of rebuilding:** the full build chunk from Week 1 is kept with `eval=FALSE`. The report reconnects by loading the existing `repository_config.R` with `load_repository_config()` and opening the database with `open_duckdb()`. `generate_real_world_repository(overwrite = FALSE)` can't be used for this because it stops when the repository already exists.
- **Outcome variable:** six 0/1 columns could serve as the outcome, and all have complete records (the same count). `DEATH_EVENT` was used because it is the clinical outcome the dataset was collected to study.
- **Excluded column:** `time` (follow-up period) was left out of the features because it depends on the outcome: patients who died were followed for less time.
- **Column types:** some columns came back from DuckDB as text or `integer64`, which `shapiro.test()` rejects. All features are converted with `as.numeric(as.character(x))` before testing.
- **Test choice:** "paired columns" was read as the same column compared across the two groups. The groups are different patients, so independent-sample tests were used, not paired tests.
- **Not in this repo:** the DuckDB database, Parquet files and downloaded data. The code regenerates them, and the source datasets have their own licence terms.

## AI reflection

I used Claude as a debugging partner and tutor throughout this assignment. I ran every piece of code myself and checked it against the actual output. Several suggestions were wrong, didn't fit my repository, or needed simplifying, and I had to verify or correct those.

**What I verified or corrected**

- **My own hand calculation.** In 4a I first substituted the number of successes (10) for p, which gave 1 − p = −9. That was the clue that p must be a probability between 0 and 1. I redid the calculation with p = 0.8 and confirmed it with `dbinom()`.
- **Guessed table names.** The AI assumed a breast cancer table called `UCI_BREAST_CANCER`. DuckDB said it didn't exist and suggested the heart failure table instead. I now check `DBI::dbListTables()` before writing any query.
- **A connection problem that took several attempts.** View registration failed for all 10 tables with a read-only error. The AI's first reconnect chunk opened the database read-only, and its next fixes (disconnecting, restarting R) didn't work. The real cause was that RStudio restored an old read-only `con` from `.RData` after every restart; `dbDisconnect()` failing with "Loading required namespace: duckdb" gave that away. I turned off saving and restoring `.RData`, and found the views already existed, so re-registering wasn't needed at all.
- **Data types from the database.** `shapiro.test()` failed with `is.numeric(x) is not TRUE` because some columns came back as text or `integer64`. I added a type check and conversion step.
- **Choosing a dataset.** I tried NHANES, a sickle cell GEO dataset (GSE35007) and breast cancer before settling on heart failure. The sickle cell data would have meant adding a new source and rebuilding the repository, so I stayed with data already in it. I also asked for shorter, simpler code when the first versions were longer than needed.

**What I understand differently now**

- "At least 10" means `1 - pbinom(9, ...)`, not `pbinom(10, ...)`, because the binomial is discrete. For the same reason `qbinom()` needs + 1 to find an upper-tail cut-off, and the real tail probability (0.035) is not exactly 0.05.
- With hundreds or thousands of observations, Shapiro-Wilk rejects normality for even tiny departures, so the QQ plot is needed to judge whether a departure matters.
- A test should be chosen from the data (Welch's t-test or Wilcoxon), not used by default. Wilcoxon compares distributions, roughly medians, which is why the results table reports medians.
- When there is no real difference, p-values are uniformly distributed and 5% fall below 0.05 by chance, as in Questions 2, 3 and 6a. When there is a real difference, they pile up near 0, and that fraction is the test's power.
- Power depends on effect size relative to noise. A small effect in noisy data (platelets) needs a sample far larger than the dataset has, which explains why the p-value histograms barely changed between n = 10 and n = 20.
