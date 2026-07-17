# **Quantification of Hepatic Fuel Usage with <sup>13</sup>C-Metabolic Flux Analysis**

## 🧭 Overview

This repository contains the computational framework used to quantify **hepatic energy metabolism and fuel utilization** through **in vivo <sup>13</sup>C metabolic flux analysis (MFA)**.

The framework integrates multiple isotope-tracing experiments with an Elementary Metabolite Unit (EMU)-based model of the liver, lumped peripheral tissues, and blood circulation. The models estimate hepatic metabolic fluxes under fasting and refeeding conditions and progressively incorporate additional tracer information, including hepatic triacylglycerol kinetics and dietary amino-acid tracing.

---

## 🗃️ Data Organization

The 🗂️ **`data`** folder serves as the central data hub directly accessed by the MFA models. It contains the circulating-nutrient tracing data reanalyzed from our previous work (Hui et al., *Cell Metabolism*, 2020). It also contains the compiled and finalized cleaned datasets used as direct model inputs.

The 🗂️ **`liver TAG tracing kinetics`** folder contains the linoleate pulse–chase data, and the 🗂️ **`Amino acid tracing`** folder contains the <sup>13</sup>C-protein tracing data. Cleaned datasets generated from both folders are exported to the central 🗂️ **`data`** folder before being read by the corresponding MFA models.

---

## ⚛️ Finalized MFA Models

The five finalized model folders are organized by nutritional state and increasing model complexity.

### Fasted State

🗂️ **`LvM_Fasted`**  Base model integrating circulating nutrient tracers to quantify hepatic metabolic fluxes and fuel utilization.

🗂️ **`LvM_Fasted_TAGkinetics`** Expanded model incorporating linoleate pulse–chase kinetics to deconvolute oxidation of circulating fatty acids and hepatic triacylglycerol.

### Refed State

🗂️ **`LvM_Refed`** Base model integrating circulating nutrient tracers to quantify hepatic metabolism during refeeding.

🗂️ **`LvM_Refed_TAGkinetics`** Expanded model incorporating linoleate pulse–chase kinetics to deconvolute oxidation of circulating fatty acids and hepatic triacylglycerol.

🗂️ **`LvM_Refed_TAGkinetics_portalAAs`** The full model further incorporating dietary protein tracing to quantify the contribution of dietary amino acids to hepatic energy metabolism.

---

## ⚙️ Computational Workflow

Each finalized model folder follows a numbered computational workflow. Files ending in `.R` are executable scripts. The script-generated objects are saved as `.RData` files with corresponding numerical prefixes, and these files are read as input in following scripts. The scripts should be run sequentially from Step 1 through Step 6.



Using the 🗂️ **`LvM_Fasted`** model as a representative example:



#### Step 1 — Supplementary Functions

📝 **`1_supplement_function.R`** defines utility functions required throughout the  MFA workflow, including numerical matrix operations, data-structure conversion, and replacement of flux notation with numerical flux values.

#### Step 2 — Stoichiometric and EMU Calculator

📝 **`2_calculator function.R`** models the stoichiometric network and carbon transition, performs EMU decomposition, and generates the EMU balance matrices and their derivatives. 

The metabolic network is specified in the Reaction and Carbon Atom Transition Table in **`LvM_fasted.xlsx`**. Each model folder contains its own table.  

#### Step 3 — Parallel EMU Decomposition

📝 **`3_decompose_parallel.R`** applies the `func.stoicEMU()` function defined in 📝`2_calculator function.R` to performs EMU decomposition for all tracer experiments in parallel, and compiles the tracer-specific and shared model objects required for flux optimization.

#### Step 4 — MFA Convergence and Optimization

📝 **`4_MFA_convergence.R`** simulates isotopologue labeling across all tracers, constructs flux-sensitivity matrices, and iteratively optimizes metabolic fluxes using constrained quadratic programming until convergence.

During each MFA run, temporary tracer-specific data are written to and read from the 🗂️ **`tracer_k_parallel`** folder. These intermediate files, including simulated labeling and sensitivity information, are updated during each optimization iteration.

#### Step 5 — Repeated Optimization and Optimal-Solution Selection

📝 **`5_1_repeat_to_confirm_optimal.R`** repeats the complete MFA optimization from multiple random initial flux sets to reduce the risk of convergence to a local minimum.

The output from each repeated MFA run is saved in the 🗂️ **`repeat_results`** folder. A graphical summary of each run, including convergence and model-fit diagnostics, is saved in the 🗂️ **`plots/converging_history`** folder.

📝 **`5_2_repeat_to_confirm_optimal_analysis.R`** compares the repeated runs and selects the repeat with the lowest cost as the optimal solution.

#### Step 6 — Confidence-Interval Analysis

📝 **`6_1_confidence_singleFlux.R`** calculates confidence intervals for individual reaction fluxes.

📝 **`6_2_confidence_netFlux.R`** calculates confidence intervals for net fluxes defined as the difference between paired forward and reverse reactions.

Confidence-profile results generated by both scripts are saved in the 🗂️ **`CI`** folder.

📝 **`6_3_confidence_analysis.R`** compiles the confidence-interval results and creates the confidence profiles.