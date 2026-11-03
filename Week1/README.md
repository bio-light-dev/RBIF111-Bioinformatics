# RBIF-111 Week 1 Homework

**Author:** Coddie Curry  
**Repository:** [github.com/bio-light-dev/RBIF111-Bioinformatics/tree/main/Week1](https://github.com/bio-light-dev/RBIF111-Bioinformatics/tree/main/Week1)  
**Files:** `HW1-Week-1.Rmd` (source), `HW1-Week-1.html` (knitted report)

Week 1 covers programming fundamentals in R, getting open-source biomedical data (GEO, MIMIC-III, NHANES), and storing it in data containers (DuckDB, SummarizedExperiment).

## What was completed

| Question | Task | Key result |
| --- | --- | --- |
| 1a–1b | Loop-based mean function, with `NA` handling | Matches `mean()` / `mean(..., na.rm = TRUE)` |
| 1c | Built the course repoquet DuckDB repository, confirmed table and column with `DBI::dbListTables()` and `DESCRIBE`, pulled a numeric column with SQL | NHANES age (`RIDAGEYR`, `NHANES_DEMO`); own mean function agrees with `mean(na.rm = TRUE)` |
| 1d | Correlation test, age vs family income-to-poverty ratio (NHANES 2021–2023) | r = 0.157, t = 15.84, df = 9,890, p < 2.2e-16: significant but weak (r² ≈ 2.5%), not causal |
| 2a | Downloaded GSE266550 (human RNA-seq) with `GEODataDownload()`, saved as tab-delimited `.txt`, reloaded | 17,867 genes × 17 columns before and after reload |
| 2b | Added a `Subset_Comparison = "Treatment-Control"` column | Explained why labelling matters when combining DE results |
| 2c–2e | Per-sample and per-gene averages, with plots; top and bottom 5 genes | Per-gene values plotted on log2(x + 1) scale |
| 2f–2g | Volcano and MA plots (ggplot2) | 41 genes with adjusted p < 0.05 |
| 3a–3b | Built a SummarizedExperiment (raw counts + log2 CPM assays, DE results in rowData, group labels in colData) and checked all names line up | 17,867 × 10; all `identical()` checks TRUE |
| 3c–3e | Subset by treatment, by significance, and both; top 5 genes table; why use SummarizedExperiment | 41 genes × 6 treatment samples; top gene BICDL1 (padj ≈ 3e-7) |
| 4a | Why Bioconductor uses genes × samples but modelling tools expect samples × features | Wrong orientation runs silently in `lm()` and gives meaningless output |
| 4b | `gene_meta_cor()`: correlates one gene with a metadata variable, checking inputs first | BICDL1 vs treatment group: r = 0.934, p = 7.9e-05. A rowData variable is rejected with a length-mismatch error |

## Data and setup notes

- **GEO annotation fix:** NCBI returned a web page instead of the gene annotation file, so `GEODataDownload()` was changed to annotate Entrez IDs with Bioconductor's `org.Hs.eg.db`.
- **Sample count:** the counts file has 11 samples but the group string `"1011101100"` has 10, so GSM8250231 is left out of the comparison. The analysis uses the 10 samples in the comparison (6 treatment, 4 control).
- **DuckDB repository:** built with the instructor's [repoquet](https://github.com/brendangongol/repoquet) package (`profile = "quick"`) in `~/repoquet_quick_example`. The ClinVar source was dropped because its server blocked the download. The build chunk is set to `eval=FALSE` in the Rmd; the report connects to the existing database.
- **Not in this repo:** downloaded data, the DuckDB database and Parquet files. The code regenerates them, and MIMIC-III and NHANES have their own licence terms.

## AI reflection

I used Claude as a debugging partner and tutor throughout this assignment. I ran every piece of code myself and checked it against the actual output. Some of the AI's suggestions were wrong or didn't match the instructions, and I had to verify or correct those.

**What I verified or corrected**

- **Instructions over convenience.** For 1c the AI first used ICU length of stay from MIMIC-III. That is numeric, but it isn't the kind of variable the question gives as examples (a vital sign, a lab value, or an NHANES demographic or exam measure). I pushed back and switched to NHANES age. I also made sure the code used `DBI::dbListTables()` exactly as the question specified.
- **A misleading `TRUE`.** My first 1c comparison printed `TRUE`, but both values were missing. SQL accepted `SELECT los`, while DuckDB returned the column as `LOS`, so `$los` was `NULL`. `all.equal(NaN, NA)` returned `TRUE` even though nothing had been calculated. I now check `length()` of what I pull before trusting a comparison.
- **The first diagnosis was wrong.** The AI first blamed the failed GEO annotation download on a duplicated URL parameter. Fixing that didn't help. The 21 KB download (instead of about 4 MB) showed NCBI was returning a web page, which led to the `org.Hs.eg.db` workaround.
- **Silent data problems.** Building the SummarizedExperiment failed with "different row counts implied by arguments." Tracing it showed the data has 11 samples but the group string describes 10, so one sample had been silently dropped from the differential expression analysis. I wouldn't have noticed without that error.
- **Smaller corrections.** I removed a redundant `datatable()` call the AI added. I pointed out that the instruction in the SchemaReview workbook (D2) says to select a decision but the value actually goes in a different column (A2). I rewrote my 2b explanation after realising my first answer addressed why to edit data in R, not why to label combined results.

**What I understand differently now**

- Code that runs without errors can still be wrong. Mismatched row names in `merge()`, a transposed matrix in `lm()`, or lengths that happen to match can all produce confident but meaningless results. Checks with `dim()`, `length()`, `identical()` and `stopifnot()` are part of the analysis, not extras.
- How log2 CPM is calculated step by step (library size, counts per million, +1, log2) and why normalisation removes sequencing-depth effects. That is why BICDL1 correlated strongly with treatment (r = 0.93) but not with library size (r = −0.20).
- Why a gene's expression (one value per sample) can't be correlated with rowData (one value per gene), and why that check belongs inside the function.
- Knitting runs in a fresh session with the Rmd's folder as the working directory, which explains why installs, paths and database connections behaved differently there than in the console.
