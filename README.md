# Neuro Oligo Platform (NOP)

A SimBiology quantitative systems pharmacology (QSP) model of oligonucleotide PK/PD across
plasma, brain (superficial/deep), CSF, spinal cord and liver compartments. This repository
contains the model, data, and driver scripts used to **reproduce the manuscript figures**.

## Requirements

- **MATLAB** (R2024b was used)
- **SimBiology** — model construction and simulation
- **Statistics and Machine Learning Toolbox** — summary statistics used by the plotting helpers
  (and `partialcorr`/`makedist` for the sensitivity analysis)
- **Parallel Computing Toolbox** — *optional*. The scripts set `args.UseParallel = 1`; if you do
  not have this toolbox, set `args.UseParallel = 0` before running (or edit it in the script).

All required input files (`args_NOP.xlsx`, `parameters/`, and the curated datasets in `data/`) are
included in the repository.

## Setup

1. Clone/download the repository and open MATLAB.
2. In MATLAB, set the **Current Folder** to this project folder (`IT_NABM`).
3. When you press **Run** on a script for the first time, if MATLAB prompts, choose
   **"Change Folder"** so the current folder is this project folder. The scripts add the repository
   to the MATLAB path automatically (`addpath(genpath(pwd))`).

## Reproducing the manuscript figures

Open each script below in the MATLAB Editor and press **Run** (▶) or **F5**. Each script builds the
model, loads the data, and generates its figures. Figures are displayed on screen and exported to
the `figures/` folder (configurable in `args_NOP.xlsx`).

### 1. PK/PD figures — `simulate_NOP.m`
Builds the model and runs the BIIB080, APP siRNA (Mivelsiran), Tominersen, and Jacifusen scenarios
(two translation variants each), producing the model-vs-data PK/PD overlay figures.

### 2. Sensitivity analysis (one-at-a-time) — `SA_main.m`
Runs one-at-a-time parameter sweeps over the brain oligo half-life, brain protein half-life, IC50,
cellular uptake half-life, and mRNA half-life, producing the local sensitivity figures.

### 3. Allometric scaling — `scaling.m`
Standalone script (base MATLAB only). Produces the cross-species allometric scaling figures and
saves them as `IC50_final.png`, `oligo_hl_final.png`, and `prot_hl.png` in the project folder.

### 4. Global sensitivity analysis (PRCC) — `calibration/script_NOP_calibration.m`
Runs the global sensitivity analysis and produces the **PRCC** global sensitivity figure (and Sobol
figures). By default the parameter-estimation step is off (`args.switches.calibration = 0`) and
PRCC/Sobol are on, so simply pressing Run generates the sensitivity figures. This script locates the
project root automatically and can be run regardless of the current folder.

## Notes

- All exported figures go to the `figures/` folder by default; the location and file format are set
  in `args_NOP.xlsx`.
- The first run is slower because the model is compiled and (if enabled) a parallel pool is started.
- `calibration/script_NOP_calibration.m` (re)writes a large results file,
  `calibration/PSA_results_NOP.mat`, which is not tracked in version control.


## License
Please see Genentech Non-Commercial Software License.docx for license
