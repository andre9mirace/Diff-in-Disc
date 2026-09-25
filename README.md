# Diff-in-Disc

**Difference-in-Discontinuities: a Method for Territorial Policy Evaluation — A Case Study on Rent Control in Paris**

André Miranda Acevedo — Master Statistiques et Économétrie, Université de Strasbourg (2023-24).
Internship at aivancity School for Technology, Business & Society.

The thesis formalizes the Difference-in-Discontinuities (Diff-in-Disc) estimator, a geographic regression
discontinuity design with a time dimension, and compares bandwidth selectors for it: the Imbens–Kalyanaraman
plug-in, the robust bias-corrected CCT bandwidth (`rdrobust`), and an adaptive stochastic gradient descent
(SGD) selector. It evaluates them on simulated data and then estimates the effect of the July 2019
*encadrement des loyers* (rent control) in Paris on property prices, using transactions near the Paris
city border (DVF 2019 to May 2021).

The full report is [`MIRANDA_Andre_M2_SE_2023-24_MEMOIRE.pdf`](MIRANDA_Andre_M2_SE_2023-24_MEMOIRE.pdf).

## Repository structure

```
Simulation/
  data_generation.Rmd      Simulated Paris / Issy-les-Moulineaux panel -> sim_data.csv
  Eval_methods_Sim.Rmd     DiD, RD and Diff-in-Disc estimates on the simulated data
  sim_data.csv
Bandwidth/
  Imbens/Imbens_BW.Rmd     Imbens–Kalyanaraman plug-in bandwidth
  SGD/BW_Functions.Rmd     Adaptive SGD bandwidth selector (final version)
  SGD/Functions/           Earlier drafts (Outdated/) and unfinished variants (Unfinished/)
Application/
  DVF_database.Rmd         Builds grouped_DVF.csv from raw DVF + commune shapes
  DVF_application.Rmd      Descriptive analysis, clustering and Diff-in-Disc estimation
  DVF/grouped_DVF.csv      Cleaned transaction data used by the analysis
  communes-dile-de-france-au-01-janvier/   Commune boundaries (shapefile)
  Graphs/                  Figures used in the report
Presentation/              Defense slides and images
Manuscript/                LaTeX source of the report body (V1.tex; pdfLaTeX + Biber, e.g. on Overleaf)
```

## Running the code

1. Open `Diff-in-Disc.Rproj` in RStudio. All paths are relative to the project root via `here::here()`,
   so the folder can live anywhere.
2. Restore the package versions with `renv::restore()`.
3. Knit the notebooks in this order:
   1. `Simulation/data_generation.Rmd` (seeded; regenerates `sim_data.csv`)
   2. `Simulation/Eval_methods_Sim.Rmd`
   3. `Bandwidth/Imbens/Imbens_BW.Rmd` and `Bandwidth/SGD/BW_Functions.Rmd`
   4. `Application/DVF_database.Rmd` (needs the raw DVF files, see below)
   5. `Application/DVF_application.Rmd` (runs directly from the versioned `grouped_DVF.csv`)

### Raw DVF data

`DVF_database.Rmd` reads the geolocated DVF extracts for 2019, 2020 and 2021. They are about 1 GB, so they
are not versioned. Download `full.csv.gz` for each year from
[files.data.gouv.fr/geo-dvf/latest/csv/](https://files.data.gouv.fr/geo-dvf/latest/csv/), unzip each one,
and save it as `Application/DVF/2019_full.csv`, `2020_full.csv` and `2021_full.csv`.
