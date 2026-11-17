# RBIF-111 Week 3 Homework

**Author:** Coddie Curry  
**Repository:** [github.com/bio-light-dev/RBIF111-Bioinformatics/tree/main/Week3](https://github.com/bio-light-dev/RBIF111-Bioinformatics/tree/main/Week3)  
**Files:** `HW3-Week-3.Rmd` (source), `HW3-Week-3.html` (knitted report)

Week 3 covers linear models and their extensions: fitting polynomial terms in a linear model, gene-by-gene regression on real expression data, non-linear smoothers (LOESS and GAM), array normalization, and checking model assumptions with residuals. Part A uses simulated data; Part B applies the methods to a GEO microarray dataset with age as the continuous variable.

## What was completed

| Question | Task | Key result |
| --- | --- | --- |
| 1 | Simulated `y = 20 + 10x + 5x² + N(0, 200)` for x = 1–20; fitted `y ~ 1`, `y ~ x` and `y ~ x + I(x^2)`; compared RSS, R², AIC and nested F-tests; plotted fitted values and residuals; repeated the fit 1,000 times at noise sd = 10, 50, 200, 500 and 1000 | The quadratic model fits best (R² 0.94, AIC 271.7 vs 280.2 linear and 324.5 intercept-only; adding x² is significant, p = 0.003). Linear-model residuals are U-shaped; quadratic residuals show no pattern. Only the x² coefficient is recovered well (4.95, CI 1.9–8.0, truth 5); intercept and slope CIs include 0. As noise rises, estimates stay unbiased but their spread grows, and power to detect x² falls from 100% (sd 10–50) to 87% (200), 25% (500) and 10% (1000) |
| 2 setup | Downloaded GSE11882 (Affymetrix GPL570, post-mortem brain tissue from cognitively normal adults) with `GEODataDownload()`; continuous variable = age | 44,134 probes × 173 samples. Ages 20–99 years (median 69, mean 60) from the `age (yrs):ch1` column, with no donors between about 52 and 64, so the samples form a younger and an older group |
| 2 | Looped `anova(lm(age ~ expression))` over every probe (21 s); plotted the p-value distribution; summaries, ANOVA tables and diagnostic plots for the most and least significant probes | Huge spike of small p-values: 25,696 probes (58%) with p < 0.05 and 23,243 with FDR < 0.05. Most significant: **CX3CL1** (823_at), slope −24.1, R² 0.36, F = 94.5, p < 2.2e-16 (higher expression in younger donors). Least significant: LOC101928111 (237764_at), p = 1, R² ≈ 0, a flat fit equivalent to the intercept-only model in Question 1 |
| 3 | Fitted LOESS (span 0.75) and a GAM (`mgcv`, REML) to CX3CL1; plotted all three fits on one scatterplot; compared AIC | LOESS AIC 1528.5 and GAM 1528.7 vs linear 1536.5 (ΔAIC ≈ 8). GAM edf = 3.97; F-test linear vs GAM p = 0.002; adjusted R² 0.40 vs 0.35. The curves show a threshold shape: age stays near 72 at low expression, falls steeply between expression ≈ 1 and 2, then levels off near 20. The non-linear fits capture real structure the line missed, and avoid the line's impossible predictions below age 20 |
| 4 | Checked whether the GEO data were already normalized; quantile-normalized with `limma::normalizeBetweenArrays()`; boxplots before and after; reran the Question 2 loop and compared | Submitters had already applied GC-RMA and GeneSpring normalization (per chip to the 50th percentile, per gene to the median), but sample medians still ranged 0.79–1.01 and **declined with age** (slope p = 2.9e-7). After quantile normalization: p < 0.05 fell from 25,696 to 17,604 and FDR < 0.05 from 23,243 to 13,962; Spearman correlation of p-values 0.62; 58 of the top 100 probes shared. CX3CL1 remained the top gene and got stronger (R² 0.36 → 0.40, F 94.5 → 112.9). Extra plots: per-gene −log10(p) before vs after, and sample median vs age before and after |
| 5 | For CX3CL1 after normalization: observed vs fitted plot, residuals vs fitted, histogram and Q-Q plot of residuals, Shapiro-Wilk test | R² = 0.40, model p = 1.4e-20. Shapiro-Wilk W = 0.966, p = 0.0003, so normality is rejected: the Q-Q plot shows a long negative tail (young donors predicted as old), and residual spread is much wider at fitted ages above 55. Age itself is non-normal (W = 0.92, p = 4.4e-8). The association is far too strong to disappear, but the exact p-value shouldn't be taken at face value |

## Data and setup notes

- **Dataset choice:** the Week 1 dataset (GSE266550) only had a treatment/control grouping, so it couldn't be used for a regression on a continuous variable. GSE11882 was chosen because it has age for every sample, and as an array dataset the metadata comes with the download.
- **gsm string:** `GEODataDownload()` runs a two-group limma contrast and needs a group string. Samples were coded 1 = at or above the median age and 0 = below, only to satisfy the function. All regressions use age as a continuous variable.
- **Expression scale:** the GEO values are GeneSpring-normalized ratios (each gene scaled to its median, range about 0–2.6), not log2 values. The plot axes labelled "log2 expression" should read "normalized expression".
- **Local cache:** the 41 MB series matrix timed out when `getGEO()` downloaded it (120-second default). It is now downloaded once with `options(timeout = 1200)` into a `GEOcache` folder next to the Rmd and read locally with `getGEO(filename = ..., getGPL = FALSE)`. The GPL570 annotation file is cached the same way.
- **Missing ID column:** reading the file with `getGPL = FALSE` leaves the feature data empty, so `topTable()` returned no `ID` column and `GEODataDownload()` failed with "undefined columns selected". A temporary `getGEO()` wrapper in Q2b adds the probe IDs back, and it is removed after the download.
- **AIC for LOESS:** `loess` objects have no `AIC()` method. AIC was computed from the Gaussian log-likelihood using the effective number of parameters (`trace.hat`), and the same formula was checked against `AIC()` for the linear model.
- **Non-independent samples:** GSE11882 includes several brain regions from the same donors, so samples are not independent and brain region is not in the model. This, together with the age-related intensity trend found in Question 4, explains why more than half of all probes look significant. All p-values in Part B are optimistic.
- **Not in this repo:** the `GEOcache` folder and downloaded GEO files. The code downloads them on the first run.

## AI reflection

I used Claude as a debugging partner and tutor for this assignment. I ran every piece of code myself and checked it against the actual output. Several suggestions were wrong, overcomplicated or didn't match my setup, and I had to correct or simplify them.

**What I verified or corrected**

- **Too many versions of the same script.** The AI gave me the Question 2 code as an Rmd, then as plain R, then as separate setup and analysis blocks, and later suggested a second approach reusing my Week 1 data. I ended up with pieces that didn't fit together. I asked for one clean Q2a (setup) and Q2b (analysis), which fixed the confusion.
- **A missing function.** The first script assumed `GEODataDownload()` was already defined and stopped at a `stopifnot()` check. I had to ask for the function from Week 1 Part 3 to be included.
- **Package and version warnings.** GEOquery was missing after an R update and had to be reinstalled with `BiocManager::install()`. The "limma was built under R version 4.6.1" message turned out to be a harmless warning, not the cause of the error.
- **Download failures.** `getGEO()` timed out at 90% of the file, then failed again with "cannot open the connection" because it contacts NCBI even when a cached copy exists. Caching the file and reading it locally fixed both.
- **A fix that caused a new error.** Reading the local file without the platform annotation removed the `ID` column the function needed. The error message ("undefined columns selected") pointed to the `subset()` call inside `GEODataDownload()`.
- **Wrong axis labels.** The code labelled expression as log2, but the boxplots showed values from 0 to about 2.6 centred near 0.9, which matches the submitters' per-gene median normalization, not a log scale.

**What I understand differently now**

- A model is "linear" when it is linear in its coefficients, not in x. That is why `lm()` can fit a curve by treating x² as another predictor.
- With noisy data, a term that is really there (x²) may not be detected in a single sample. Power depends on the effect size relative to the noise.
- Significance should be judged by the p-value, not the size of the slope, because slopes depend on each gene's expression scale.
- When more than half of all genes are "significant", that is a warning sign, not a discovery. Here it came from a technical intensity trend that tracked age and from repeated samples per donor.
- Data labelled as normalized on GEO still need checking. The sample medians in GSE11882 drifted downward with age, and quantile normalization removed about a third of the significant probes while making the strongest gene stronger.
- Extra flexibility (LOESS, GAM) always reduces the residual error. AIC, the GAM's effective degrees of freedom and an F-test show whether the improvement is real. For CX3CL1 it was: the relationship is a threshold, not a straight line.
- With 173 samples, Shapiro-Wilk rejects normality for modest departures, so the Q-Q plot and residual plot are needed to see what kind of problem it is. Here the residuals reflect the two age groups in the study design. Unequal variance and non-independent samples matter more for the p-value than non-normality.
