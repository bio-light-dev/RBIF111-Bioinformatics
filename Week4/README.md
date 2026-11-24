# RBIF-111 Week 4 Homework

**Author:** Coddie Curry  
**Repository:** [github.com/bio-light-dev/RBIF111-Bioinformatics/tree/main/Week4](https://github.com/bio-light-dev/RBIF111-Bioinformatics/tree/main/Week4)  
**Files:** `HW4-Week-4.Rmd` (source), `HW4-Week-4.html` (knitted report), `data/heart_failure_clinical_records.csv` (UCI dataset used in Question 4)

Week 4 covers the standard error of the mean, the Central Limit Theorem, confidence intervals and time series analysis. Part A (Questions 1–3) uses simulated data from distributions of known shape, so the behaviour of the SEM, CI coverage and CLT convergence can be checked against the truth. Part B (Questions 4–6) applies the same ideas to real data: the UCI Heart Failure Clinical Records dataset, a GEO microarray dataset with age as a continuous variable, and a repeated heart-rate series from the MIMIC-III demo.

## What was completed

| Question | Task | Key result |
| --- | --- | --- |
| 1a | SEM plot for the exponential distribution (`rexp`, rate λ = 2): 1,000 samples at n = 1–128, simulated vs theoretical SEM (1/λ)/√n | Mean of sample means stays at ≈ 0.50 = 1/λ (unbiased). Simulated SEM matches theory at every n (ratio 0.94–1.04), falling from 0.50 at n = 1 to 0.044 at n = 128, about √128 ≈ 11-fold |
| 1b | 95% CI of the mean for one exponential sample (n = 20) using `confint(lm(x ~ 1))` and manually with `mean`, `sd`, `qt` | Both methods give [0.172, 0.639] (x̄ = 0.406, SE = 0.111, t₀.₉₇₅,₁₉ = 2.093); `all.equal()` = TRUE; the interval contains the true mean 0.5 |
| 1c | CLT convergence for sums of n = 1–256 exponential variables, sampled explicitly and summed with `rowSums()`; histograms with normal overlay, Q-Q plots, skewness table | Strong right skew at n = 1 (skewness 1.94); the sum approaches normal as skewness falls as 2/√n (0.128 at n = 256). Simulated mean and variance match n/λ and n/λ² throughout. Looks close to normal by n ≈ 16–32 |
| 2a | SEM plot for the uniform distribution U(0, 6) (`runif`), theoretical SEM √((b−a)²/12)/√n | Mean of sample means ≈ 3 at every n; simulated SEM tracks √3/√n (ratio 0.98–1.04), from 1.72 at n = 1 to 0.151 at n = 128 |
| 2b | Same CI calculation for one U(0, 6) sample (n = 20) | Both methods give [2.425, 4.185] (x̄ = 3.305, SE = 0.421); identical; contains the true mean 3 |
| 2c | CLT convergence for sums of uniform variables; excess kurtosis tracked instead of skewness (the uniform is symmetric) | Flat at n = 1 (excess kurtosis −1.20), triangular at n = 2, visually normal by n ≈ 4–8. Converges much faster than the exponential because there is no skew to average out |
| 3 | CI coverage simulation: 2,000 samples at n = 2–1024 from N(0,1), Exp(1) and U(0,1); failure rate F = % of `confint()` intervals missing the true mean; noise band 4.04–5.96% | **Normal:** F within the noise band at every n (5.05% at n = 2). **Uniform:** elevated only at n = 2–4 (7.65%, 8.10%), in the band from n = 8. **Exponential:** 8.85–10.75% for n ≤ 8 and ≈ 8.3–8.5% at n = 16–32 (real confidence ≈ 90%), reaching the band at n = 64 (5.90%). The 95% holds exactly only when the t-interval's assumptions hold |
| 4a | UCI Heart Failure Clinical Records (299 patients) loaded into DuckDB; split by `death_event` (96 died, 203 survived); 1,000 replications of n = 30 per group, `t.test()` CI for the difference in mean ejection fraction, stored as a 2 × 1000 matrix via `replicate()` | Average CI length **12.09 percentage points** (range 8.62–14.99) |
| 4b | Population-level difference from the full dataset | Died 33.47% vs survived 40.27%: **−6.80 points**. Mean of the 1,000 sample differences −6.61 (SD 2.60) |
| 4c | Coverage of the 1,000 CIs | **975 of 1,000 (97.5%)** contained −6.80 (13 missed below, 12 above), above the 93.6–96.4% chance range: sampling 30 of only 96 patients without replacement makes the t-intervals conservative (finite population correction) |
| 4d | Repeated (a)–(c) for all seven continuous columns, plus the CI from the entire dataset | Clear differences for time (n = 30 detects it 99.9% of the time), serum creatinine (69%), ejection fraction (57%) and age (54%); a small one for serum sodium (32%); none for platelets and creatinine phosphokinase (p = 0.40, 0.37). Coverage 97.4–98.2% for symmetric columns but 94.0–94.5% for the skewed lab values. Full-dataset CIs are 1.86–2.14 times narrower than the average n = 30 interval, close to the √(1/30+1/30)/√(1/96+1/203) ≈ 2.1 predicted |
| 5 | GSE11882 (GPL570, 173 brain samples, quantile-normalized as in Week 3), `lm(age ~ expression)` for all 44,134 probes, `confint()` on the expression coefficient, genes ordered by p-value and all CIs plotted | 17,604 of 44,134 CIs exclude zero (≈ 2,200 expected by chance). Strongest: **CX3CL1** (823_at), slope −25.8 years per unit, 95% CI [−30.6, −21.0], entirely below zero. Weakest: **SRSF3** (232392_at), slope 0.0, CI [−7.4, 7.4], centred on zero and wider despite similar expression spread |
| 6a | MIMIC-III demo `CHARTEVENTS` (758,355 rows, 99 patients, 133 ICU stays) in DuckDB; hourly heart rate for ICU stay 249805 (1,115 readings, 35.2 days, 2.8% of hours interpolated); STL decomposition | A 24-hour seasonal term was not supported (seasonal strength 0.012, amplitude 2.4 bpm vs remainder SD 5.8), while the trend was strong (0.795). Decomposed as **trend + remainder**: multi-day swings between about 67 and 123 bpm, remainder SD 4.9 bpm with a few isolated spikes. Critical illness and ICU care blunt the normal circadian rhythm |
| 6b | Augmented Dickey-Fuller test (`adf.test()`) | ADF = −3.49, p = 0.043: stationary (d = 0), but only just; the fitted AR coefficients sum to ≈ 0.98, so deviations from the mean fade slowly |
| 6c | `auto.arima()` and a 24-hour forecast with 80%/95% prediction intervals | **ARIMA(4,0,1)**, mean 97.4 bpm, residual SD 4.8 bpm; Ljung-Box p = 0.70 (no autocorrelation left). Forecast 95.5 → 96.0 bpm, 95% interval widening from 86–105 bpm (1 h) to 74–118 bpm (24 h). Clinically plausible as a baseline, but it cannot anticipate abrupt clinical events |

## Data and setup notes

- **Course repository rebuild failed:** the repoquet repository from Weeks 1–2 (`~/repoquet_quick_example`) had been deleted. Rebuilding it failed because one source in the quick profile, UCI #891 (CDC Diabetes Health Indicators), repeatedly broke off mid-download ("connection reset by peer"). `FinalizeSchemaRegistry()` refuses to finish if any source has an error, so no tables were registered for any dataset. Filtering that row out of the setup table did not work either.
- **Question 4 workaround:** the Heart Failure file (UCI #519) had downloaded successfully into the repoquet cache, so it is loaded into an in-memory DuckDB database with `read_csv_auto()` and checked with `dbListTables()`/`DESCRIBE`. For reproducibility the code keeps a copy in `data/`; if it is missing, it takes the file from the repoquet cache or downloads it from UCI.
- **Question 4 sampling:** the 1,000 patient subsets are drawn once with `set.seed(123)` and reused for every column in 4d, in the same order as 4a, so ejection fraction reproduces parts (a)–(c) exactly. Binary 0/1 columns are excluded from 4d (their mean is a proportion). `time` is included but is partly defined by the outcome.
- **Question 5 data:** same dataset as Week 3 (GSE11882, age as the continuous variable), quantile-normalized because Week 3 showed an age-related intensity trend in the downloaded values. The series matrix and GPL570 annotation are downloaded into `GEOcache/` with a retry-and-resume helper (the NCBI download dropped at 9.9 of 39.4 MB). The prepared objects are saved to `GSE11882_Q5_data.rds`, so later knits reload them in seconds.
- **Question 6 data:** the quick repository's MIMIC-III tables (admissions, diagnoses, ICU stays, patients, procedures) contain no repeated measurements. `CHARTEVENTS.csv` (≈ 74 MB) was downloaded from the same open-access MIMIC-III demo on PhysioNet into `data/mimic_demo/`, with a `.complete` marker so it is only downloaded once. Heart rate uses item IDs 211 and 220045, limited to 20–250 bpm.
- **Plotting in the console:** running plot code line by line caused "plot.new has not been called yet" errors. Each figure is wrapped in `{ }` (and later in functions) so its parts always run together.
- **Not in this repo:** `GEOcache/`, `GSE11882_Q5_data.rds` and `data/mimic_demo/CHARTEVENTS.csv` (too large for GitHub). The code downloads or rebuilds them on the first run.

## AI reflection

I used Claude as a debugging partner and tutor for this assignment. I ran every piece of code myself and checked it against the actual output. Many problems came from my environment and from running code line by line, and several suggestions had to be corrected.

**What I verified or corrected**

- **Plots that wouldn't draw.** Running plot code line by line in the console gave "plot.new has not been called yet". Wrapping each figure in braces or a function fixed it. A later attempt to save a PNG came out empty because I pasted the placeholder comment instead of the plotting code.
- **Repeated or duplicated code.** In Question 3 the simulation loop ran twice without resetting the storage vectors, so the results doubled and the plot failed. A `stopifnot()` length check then caught the same problem when my chunk contained the loop twice.
- **The DuckDB repository.** The repository the course code expected no longer existed, and rebuilding it failed on an unrelated dataset download. Several suggested fixes (searching for the config file, searching for registered views or Parquet files, filtering the failing row) didn't solve it. Reading the cached CSV into DuckDB directly did.
- **Reproducibility.** The first working Question 4 version depended on a file that existed only on my machine. I asked what would happen if someone else ran it, which led to keeping a copy of the data in the project folder.
- **Slow and failed downloads.** Reading the Week 3 GEO files was slow because they were stored only in iCloud, and the code parsed the 41 MB series matrix twice. A direct NCBI download then dropped partway, so a retry-and-resume helper and an `.rds` cache were added.
- **Predictions that didn't match my results.** The expected write-ups said the uniform intervals would be close to 95% at every sample size and that the exponential would need "hundreds" of observations. My results showed elevated failure for the uniform at n = 2–4 and the exponential reaching the noise band at n = 64, so the interpretations were rewritten from my actual numbers.

**What I understand differently now**

- The SEM formula σ/√n needs only independent draws with finite variance, not normality. Ten times the precision costs one hundred times the data.
- The CLT is about shape. How fast a sum looks normal depends on the starting distribution's asymmetry: a symmetric uniform converges within a few terms, a skewed exponential needs dozens.
- "95% confidence" describes the procedure. A single interval either contains the true value or it doesn't. Coverage is exactly 95% only when the assumptions hold: lower for skewed data at small n, and higher when sampling without replacement from a small finite population.
- A confidence interval tells me more than a p-value: its position shows the direction and whether zero is plausible, and its width shows how precisely the effect is known.
- A narrower interval comes from more data, not a better method. The full Heart Failure dataset gave intervals about twice as narrow as samples of 30, as the 1/√n rule predicts.
- In ICU data an expected seasonal pattern may not be there. Checking seasonal strength and the ACF before keeping a seasonal term mattered, and a borderline ADF result means the series is stationary but very persistent.
- An ARIMA forecast extends past patterns. It gives a useful baseline and honest uncertainty, but it cannot predict a clinical event.
