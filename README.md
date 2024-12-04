**Combining seasonal malaria chemoprevention with novel therapeutics for malaria prevention: a mathematical modelling study**

Lydia Braunack-Mayer1,2, Josephine Malinga3,4, Narimane Nekkab1,2, Sherrie L Kelly1, Jörg J Möhrle1,2,4, Melissa A Penny1,2,3,5,*

1 Swiss Tropical and Public Health Institute, Allschwil, Switzerland
2 University of Basel, Basel, Switzerland
3 The Kids Research Institute Australia, Nedlands, WA, Australia
4 Medicines for Malaria Venture, Geneva, Switzerland
5 Centre for Child Health Research, The University of Western Australia, Crawley, WA, Australia
* Correspondence to:
Prof Melissa A Penny, melissa.penny@uwa.edu.au

In this study, we combined an individual-based malaria transmission model (https://github.com/SwissTPH/openmalaria/wiki) with explicit models of intervention dynamics for a range of hypothetical malaria prevention therapeutics. We used these models to estimate the impact of combining new therapeutics with seasonal malaria chemoprevention in children under five years old.


# Folders / Workflow Steps

## analysisworkflow

This workflow builds on the workflow presented in Golumbeanu (2021), Burgert (2021) and Braunack-Mayer (2024) to specify Target Product Profiles for new interventions against malaria. First, a set of simulated scenarios is defined. These are characterized by the delivery modality, tool specifications, and settings in which a malaria intervention is analysed. Second, health outcomes for these scenarios are simulated randomly over a large intervention parameter space. The resulting database of simulations is used to train a Gaussian process emulator, that predicts the health outcome given a set of input parameters. Third, the emulator is employed to perform sensitivity analysis of tool properties with respect to health outcomes. This analysis supports the identification of ideal product characteristics for new interventions to maximise their chance of achieving a desired health goal.

**Contributors (in chronological order): Melissa Penny, Guojing Yang, Monica Golumbeanu, Lydia Burgert, Mirjam Laager, Narimane Nekkab, Josephine Malinga, Lydia Braunack-Mayer**

### 0_scenarios
Contains exemplar XML files for each of the intervention models used to simulate data with OpenMalaria (https://github.com/SwissTPH/openmalaria/wiki).

### 1_OM_basic_workflow
Generates paramater table and XML scenarios from base scaffold.xml. Launches OM simulations with an output.txt file containing a table with four columns and no headers for survey measures. There is one line for each (5-day) time step.

### 2_postprocessing
Performs generalised post-processing of OM simulations by the settings specified in previous sets. For each setting, a split file is generated in “postprocessing/split” that specifies the parameters for this setting and based on that, a seeds file (for every simulation) and an agg file (aggregated over all seeds for one parameter set) is generated.

### 3_GP_train
Trains GPs using for a specified outcome calculated  in 2_postprocessing for each of the seeds files. 
- Predictors: therapeutic efficacy, half-life, and decay shape
- Predicted: health outcome 

### 4_sensitivity_analysis
Performs sensitivity analysis of health outcome using analysis of variance (sobol) for GPs trained in step 3 within the parameter bounds used for simulation (default).
- Predictors: therapeutic efficacy, half-life, and decay shape
- Predicted: health outcome 
Settings for the sobol analysis are defined in the function calc_sobol_idx in analysisworkflow/3_GP_train/GPtoolbox.R.

## data_and_visualisation

This folder contains the data generated during this study, along with the R scripts used to visualise data. There is a folder for each figure in the manuscript and supplement, containing:
- The .rds data file(s) corresponding to the figure,
- The Rscript(s) used to generate the figure, and
- A jpg version of the figure.
To reproduce a given figure, download the corresponding folder and update the file paths referenced in the corresponding Rscript.
