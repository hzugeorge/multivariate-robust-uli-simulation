Monte Carlo Simulation Operational Guide: Six-Step Workflow
This study aims to use Monte Carlo simulation to calculate various statistics and reliability probability tables for residual score reliability (ρ(U)) and difference score reliability (ρ(D)).
Prerequisites
	Environment: Ensure you have Python and the necessary libraries installed: numpy, pandas, os, datetime, itertools.
	Files: Place all six code files in the same working directory.
________________________________________
Step 1: Execute Monte Carlo Simulation and Data Generation (File: Step1_Monte_Carlo_Simulation_20251126)
Purpose: To run the core Monte Carlo simulation and generate sample statistics (ρ(U) and ρ(D)) under all parameter combinations.
	Check Parameter Settings: Open the file Step1_Monte_Carlo_Simulation_20251126 and confirm that the grid search parameters (e.g., VAR_X1_VALUES, RHO_X1X2_VALUES, RELIABILITY_X1_VALUES) and simulation parameters (e.g., N_SAMPLE_SIZE, N_ITERATIONS) meet your research requirements.
	Execute Code: Run Step1_Monte_Carlo_Simulation_20251126.
	Output: After the simulation is complete, the system will generate multiple batch statistics CSV files in the execution directory, formatted as:
	multivariate_simulation_stats_YYYYMMDDHHMM_Part_XX.csv
Step 2: Merge Simulation Statistics Files (File: Step2_Monte_Carlo_Simulation_20251126)
Purpose: To merge the numerous batch statistics files generated in Step 1 into larger files, facilitating subsequent analysis.
	File Movement (Recommended): Move all multivariate_simulation_stats_*.csv files generated in Step 1 into the DATA_FOLDER (default: multivariate_simulation_stats) directory.
	Check Parameter Settings: Open the file Step2_Monte_Carlo_Simulation_20251126 and confirm that DATA_FOLDER and BATCH_SIZE are correctly set.
	Execute Code: Run Step2_Monte_Carlo_Simulation_20251126.
	Output: The merged batch files are outputted, formatted as:
	multivariate_simulation_merged_Batch_XX_YYYYMMDDHHMM.csv
	Preparation for Analysis: Move all merged files to a dedicated analysis directory (e.g., multivariate_simulation_merged) for use in Steps 3 through 6.
Step 3: Analyze P(ρ(U) > ρ(D)) (File: Step3_Table2_20251126)
Purpose: To calculate the probability that residual reliability (ρ(U)) is greater than difference reliability (ρ(D)) and generate a summary table.
	Execute Code: Run Step3_Table2_20251126.
	Processing: The code reads the merged files from Step 2 and calculates the proportion of times ρ(U) > ρ(D) for each parameter combination.
	Output: A 9x9 summary table is generated (Index:ρ(x_1 ), Columns: ρ(x_2 )), including mean calculations in the margins.
	pivot_P_rhoU_gt_rhoD_YYYYMMDDHHMM.csv
Step 4: Analyze Mean ( ρ(U) > ρ(D)) (File: Step4_Table3_20251126)
Purpose: To calculate the mean difference between ρ(U) and ρ(D) and generate a summary table.
	Execute Code: Run Step4_Table3_20251126.
	Processing: The code reads the merged files from Step 2 and calculates the mean value of ρ(U) - ρ(D).
	Output: A 9x9 summary table is generated.
	pivot_Mean_RhoU_Minus_RhoD_YYYYMMDDHHMM.csv
Step 5: Analyze ρ(U) Reliability Probability 
(File: Step5_Table4_rhoU_20251126)
Purpose: To calculate the probability that ρ(U) meets or exceeds a defined threshold (default 0.7) and output detailed and summary tables.
	Check Parameter Settings: Open the file Step5_Table4_rhoU_20251126 and adjust the THRESHOLD (default $0.7$) as needed.
	Execute Code: Run Step5_Table4_rhoU_20251126.
	Output: Two tables are outputted:
	Full Detail Table: 1_Full729_rhoU_GE_0.7_YYYYMMDD_HHMM.csv (Results across three dimensions: ρ(x_1 ), ρ(x_2 ), ρ_x1x2)
	Averaged Summary Table: 2_Avg9x9_rhoU_GE_0.7_YYYYMMDD_HHMM.csv (Averaging out the effect of ρ_x1x2)
Step 6: Analyze ρ(D) Reliability 
Probability (File: Step6_Table5_rhoD_20251126)
Purpose: To calculate the probability that ρ(D) meets or exceeds a defined threshold (default 0.7) and output detailed and summary tables.
	Check Parameter Settings: Open the file Step6_Table5_rhoD_20251126 and adjust the THRESHOLD as needed.
	Execute Code: Run Step6_Table5_rhoD_20251126.
	Output: Two tables are outputted:
	Full Detail Table: 1_Full729_rhoD_GE_0.7_YYYYMMDD_HHMM.csv (Results across three dimensions: ρ(x_1 ), ρ(x_2 ), ρ_x1x2)
	Averaged Summary Table: 2_Avg9x9_rhoD_GE_0.7_YYYYMMDD_HHMM.csv (Averaging out the effect of ρ_x1x2)
