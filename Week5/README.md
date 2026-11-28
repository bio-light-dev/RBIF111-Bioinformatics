# RBIF-111 Week 5 Homework

**Author:** Coddie Curry  
**Repository:** [github.com/bio-light-dev/RBIF111-Bioinformatics/tree/main/Week5](https://github.com/bio-light-dev/RBIF111-Bioinformatics/tree/main/Week5)  
**Files:** `HW5-Week-5.Rmd` (source), `HW5-Week-5.html` (knitted report), `TCGA_BRCA_clinical_data.csv`, `TCGA_BRCA_mRNA_expression.csv`, `TCGA_BRCA_RPPA_protein_expression.csv` (TCGA-BRCA subset provided with the assignment)

Week 5 covers analysis of variance, contingency tables, generalized linear models, survival analysis, equivalence testing and clinical-trial effect measures. All questions use a 103-patient subset of the TCGA Breast Invasive Carcinoma (TCGA-BRCA) data, with matched tumor and normal mRNA expression, and focus on ERBB2 (HER2), the target of antibody-drug conjugates such as trastuzumab deruxtecan (T-DXd).

## What was completed

| Question | Task | Key result |
| --- | --- | --- |
| 1a | ERBB2 mRNA in matched tumor vs normal tissue (103 pairs): unpaired (Welch and Student) and paired t-tests on log2(x + 1) expression | Paired: mean difference **1.06 log2 units (≈ 2.1-fold higher in tumors)**, t = 6.19, df = 102, p = 1.2 × 10⁻⁸. Welch: t = 5.54, df = 196.7, p = 9.6 × 10⁻⁸. Both reject H₀. The paired test is the appropriate one (same patients, r = 0.20 between tumor and normal). Wilcoxon signed-rank (p = 7.6 × 10⁻⁸) and raw-scale tests agree |
| 1b | One-way ANOVA of tumor ERBB2 across Stage I–IV (n = 102; one "N" record excluded), equal-variance checks, post-hoc | **F(3, 98) = 0.487, p = 0.692**: fail to reject H₀, so no post-hoc needed (Tukey all p ≥ 0.70). Group sizes very unbalanced (14, 61, 25, 2). Variances equal on the log2 scale (Levene p = 0.587, Bartlett p = 0.183) but not on the raw scale (Bartlett p = 1.1 × 10⁻¹⁰). Welch ANOVA (p = 0.493) and Kruskal–Wallis (p = 0.728) agree |
| 2 | Logistic regression of ER-IHC status (71 Positive, 19 Negative) on log2 ERBB2; logistic vs linear predictions | **OR = 1.21 per doubling of ERBB2, 95% CI 0.85–1.83, p = 0.333**: not significant. The linear probability model predicts P > 1 beyond the observed range (1.08 at log2 ERBB2 = 25, 1.22 at 30); the logistic curve stays below 1 |
| 3 | Association between HER2-IHC and ER-IHC status (70 patients, 2 × 2 table) | ER-positive in 78% of HER2-negative and 75% of HER2-positive patients. One expected count < 5 (4.57), so Fisher's exact test was used: **p = 0.763**, OR = 0.85 (0.22–3.65). No association |
| 4a | Cox proportional hazards model of overall survival on HER2-IHC status (71 patients, 28 deaths), Kaplan–Meier curves | **HR = 1.83, 95% CI 0.80–4.20, p = 0.154**: HER2-positive patients had higher estimated risk of death, but not significantly. Log-rank p = 0.148 |
| 4b | Weibull accelerated failure time (AFT) model on the same outcome, compared with Cox | **Time ratio = 0.77, 95% CI 0.52–1.13, p = 0.184**; predicted median survival 82.4 vs 63.3 months. Weibull shape 1.99 (hazard rising over time), implied HR 1.69. Agrees with Cox on direction and non-significance; fitted curves track the Kaplan–Meier curves |
| 4c | Conceptual: Fine-Gray competing-risks model vs treating unrelated deaths as censoring | Censoring patients who died of other causes assumes they could still die of breast cancer, which overestimates cancer mortality. Fine-Gray gives a subdistribution HR describing the cumulative incidence (real-world probability) of breast cancer death |
| 5 | TOST equivalence test of ERBB2 in HER2-IHC Equivocal (n = 8) vs Negative (n = 50); margin = 0.5 × pooled SD = **±0.456 log2 units (≈ 1.37-fold)** | **TOST p = 0.708**, 90% CI [0.24, 0.88] extends past the upper margin: equivalence not shown. The ordinary Welch t-test was significant (p = 0.006): Equivocal tumors have **higher** ERBB2. With a 2-fold margin (±1 log2), equivalence would be supported (p = 0.014) |
| 6 | Absolute risk of death by HER2-IHC status, ARR/RRR/NNH, and effect modification by ER status (interaction logistic regression + forest plot) | Risk 36.0% (HER2−) vs 47.6% (HER2+): **absolute risk increase 11.6 points, 95% CI −13.5 to +36.8**, relative risk increase 32%, NNH ≈ 9 (Fisher p = 0.43). Stratified ORs: ER-negative 0.80 (0.09–6.85, n = 16), ER-positive 1.75 (0.52–5.89, n = 54). **Interaction p = 0.534** (LRT p = 0.532): no genuine subgroup effect |

## Data and setup notes

- **Data:** the three CSVs supplied with the assignment are in this folder. The mRNA file has genes as rows and samples as columns (`_01` = tumor, `_11` = normal). Clinical rows are matched to expression columns through `Sample_ID_Tumor` and `Sample_ID_Normal`. The RPPA protein file is included but not used in Questions 1–6.
- **log2 transform:** raw ERBB2 values are strongly right-skewed (694 to 587,741, with a few HER2-amplified tumors). All tests use log2(x + 1). Questions 1 and 2 also report raw-scale results, and conclusions are the same either way.
- **Exclusions:** each question keeps only the categories it asks for: Stage I–IV in 1b; ER Positive/Negative in Questions 2 and 6; HER2 Positive/Negative in Questions 3, 4 and 6; HER2 Equivocal/Negative in Question 5. NA, Indeterminate and Equivocal records are removed as stated in each question.
- **Survival coding:** `Overall_Survival_Status` "1:DECEASED" is the event and "0:LIVING" is censored. HER2-negative is the reference group.
- **Chunk dependencies:** the Q1a chunk creates `clin` and `erbb2`, which every later chunk uses, and Q4b uses the Cox and Kaplan–Meier objects from Q4a. Chunks must run in order.
- **Packages:** `tidyverse`, `car` (Levene's test) and `survival`.
- **Plotting in the notebook:** the base-R Kaplan–Meier plot (`plot()` + `legend()` + `mtext()`) failed with "plot.new has not been called yet" when run line by line, so it was rewritten in ggplot as a single statement.

## AI reflection

I used Claude to help write and debug the R code and to draft the interpretations. I ran every chunk myself, pasted the output back, and checked the write-ups against my actual results.

**What I verified or corrected**

- **Plot errors.** The Question 1b boxplot failed because the `geom_text()` call passed `aes()` twice, so ggplot read one as the `stat` argument. The Question 4a Kaplan–Meier plot failed in the notebook because base-R plotting split across lines. Both were rewritten.
- **Knitting.** The report wouldn't knit because I had pasted the Q1b chunk twice, giving a duplicate chunk label.
- **Missing output.** Some tables (the stage summary, the decision table, the Question 6 forest-plot data) didn't print when I ran code line by line, or I didn't copy all of the output. I reran those parts and confirmed the numbers.
- **A question's premise that didn't hold.** Question 5 asks why a non-significant t-test wouldn't prove equivalence, but here the t-test was significant (Equivocal > Negative). The write-up answers that part as a general principle and reports what my data actually showed.
- **Analysis choices.** I checked the choices that weren't spelled out in the questions: the log2 transform, using Fisher's test because of a small expected count, choosing the equivalence margin, and using Wald CIs in the forest plot so the plot matches the table.

**What I understand differently now**

- A paired design removes patient-to-patient variation. The same tumor/normal data gave a larger t and narrower CI when each patient was compared with themselves.
- ANOVA's equal-variance assumption matters most when group sizes are unbalanced, and it should be tested rather than assumed. A transformation can change whether it holds.
- Logistic regression is needed for a yes/no outcome because a straight line can predict probabilities outside 0–1.
- Cox and parametric survival models answer the same question in different units (hazard ratio vs time ratio). The Cox model makes no assumption about the hazard's shape; the Weibull model assumes one, but gives predicted survival times.
- Censoring assumes a patient could still have the event later. That is false for competing deaths, which is why Fine-Gray exists.
- "No significant difference" is not evidence of equivalence. TOST makes equivalence the thing that has to be proven, and the answer depends on the margin chosen.
- Subgroup estimates that look different are not a subgroup effect unless a formal interaction test supports it, especially with only 16 patients in a stratum.
