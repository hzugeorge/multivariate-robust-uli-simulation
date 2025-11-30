# multivariate-robust-uli-simulation

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
 
Step1_Monte_Carlo_Simulation_20251126
import numpy as np
import pandas as pd
import os
from datetime import datetime
from itertools import product 

# ======================================================================
# 🚀 Parameter Configuration Section (Grid Search) 🚀
# ======================================================================
# --- Varying Parameter Lists ---
VAR_X1_VALUES = [0.1, 0.5, 1.0] # Population Variance X1
VAR_X2_VALUES = [0.1, 0.5, 1.0] # Population Variance X2
RHO_X1X2_VALUES = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] # Population Correlation rho(X1, X2)
RELIABILITY_X1_VALUES = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] # Reliability rho(x1)
RELIABILITY_X2_VALUES = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] # Reliability rho(x2)

# --- Fixed Parameters ---
MEAN_X1 = 1.2
MEAN_X2 = 1.0

# Simulation Parameters
N_SAMPLE_SIZE = 50   # Sample size per draw (N)
N_ITERATIONS = 1000 # Number of sampling repetitions/iterations
# ======================================================================

# --- File Segmentation Settings (Dynamically Calculated) ---
MAX_RECORDS_TARGET = 1_000_000 # Target maximum number of records per output file
RECORDS_PER_COMBO = N_ITERATIONS * N_SAMPLE_SIZE # Total records generated per combination
# Calculate the maximum number of combinations per file to stay under the target limit
MAX_COMBOS_PER_FILE = int(MAX_RECORDS_TARGET / RECORDS_PER_COMBO)

# Ensure MAX_COMBOS_PER_FILE is at least 1
if MAX_COMBOS_PER_FILE < 1:
    MAX_COMBOS_PER_FILE = 1

# --- 1. Parameter Initialization and Valid Combination Generation ---
mu = np.array([MEAN_X1, MEAN_X2])

# Generate all combinations of the five factors
all_factor_values = (VAR_X1_VALUES, VAR_X2_VALUES, RHO_X1X2_VALUES, 
                     RELIABILITY_X1_VALUES, RELIABILITY_X2_VALUES)
all_combinations_full = list(product(*all_factor_values))

# Pre-filter all combinations that meet the constraint condition:
# |RHO_X1X2| must be <= sqrt(RELIABILITY_X1 * RELIABILITY_X2)
valid_combinations_list = []
combination_id_counter = 0

print("--- Pre-checking Constraint Conditions and Filtering Valid Combinations ---")
for combo in all_combinations_full:
    VAR_X1, VAR_X2, RHO_X1X2, RELIABILITY_X1, RELIABILITY_X2 = combo
    constraint_limit = np.sqrt(RELIABILITY_X1 * RELIABILITY_X2)
    
    # Check the constraint: the absolute value of the population correlation must be 
    # less than or equal to the geometric mean of the reliabilities.
    if np.abs(RHO_X1X2) <= constraint_limit:
        combination_id_counter += 1
        # Store the combination along with its unique Combination_ID
        valid_combinations_list.append((combination_id_counter, combo))

num_valid_combinations = len(valid_combinations_list)
N_total_combinations = len(all_combinations_full)

print(f"Initial total number of combinations: {N_total_combinations}")
print(f"Total valid combinations after filtering: {num_valid_combinations} types")
print("--- Simulation Parameter Settings ---")
print(f"Number of iterations per combination: {N_ITERATIONS}, Sample Size: {N_SAMPLE_SIZE}")
print(f"Total records generated per combination: {RECORDS_PER_COMBO} records")
print(f"**Based on the target limit of {MAX_RECORDS_TARGET:,} records, each raw data file will contain data for up to {MAX_COMBOS_PER_FILE} combinations**\n")

# --- 2. File Splitting and Batch Simulation Execution ---
timestamp = datetime.now().strftime("%Y%m%d%H%M")
Total_simulations_run = 0
Total_data_points_generated = 0

# Calculate how many files need to be split
num_files = int(np.ceil(num_valid_combinations / MAX_COMBOS_PER_FILE))

print(f"\n==================================================")
print(f"A total of {num_valid_combinations} valid combinations will be split into {num_files} files.")
print("==================================================")

# Execute simulation in batches and write to files
for i in range(num_files):
    # Reset lists to ensure clean collection for each batch and free up memory from the previous batch
    batch_samples_list = []  
    batch_stats_list = []    
    
    start_index = i * MAX_COMBOS_PER_FILE
    end_index = min((i + 1) * MAX_COMBOS_PER_FILE, num_valid_combinations)
    
    # Get the combinations to be processed in the current batch
    chunk_combinations = valid_combinations_list[start_index:end_index]
    
    print(f"--- Processing Batch {i+1}/{num_files} (Combination IDs: {start_index+1} to {end_index}) ---")
    
    for combination_id, combo in chunk_combinations:
        VAR_X1, VAR_X2, RHO_X1X2, RELIABILITY_X1, RELIABILITY_X2 = combo
        
        # Calculate the covariance matrix
        cov_x1x2 = RHO_X1X2 * np.sqrt(VAR_X1 * VAR_X2)
        Sigma = np.array([[VAR_X1, cov_x1x2],
                          [cov_x1x2, VAR_X2]])
        
        # Create population parameter dictionary
        population_params = {
            'Combination_ID': combination_id,
            'Population Mean X1': MEAN_X1,
            'Population Mean X2': MEAN_X2,
            'Population Var X1': VAR_X1,
            'Population Var X2': VAR_X2,
            'Population Corr X1X2': RHO_X1X2,
            'Population Reliability X1': RELIABILITY_X1,
            'Population Reliability X2': RELIABILITY_X2,
            'Sample Size N': N_SAMPLE_SIZE
        }
        
        # Inner loop: Execute N_ITERATIONS (1000 times) sampling
        for iteration in range(1, N_ITERATIONS + 1):
            Total_simulations_run += 1
            
            # a. Sampling (from a bivariate normal distribution)
            sample_data = np.random.multivariate_normal(mean=mu, cov=Sigma, size=N_SAMPLE_SIZE)
            df_sample = pd.DataFrame(sample_data, columns=['X1', 'X2'])
            
            # Attach ID tags to raw data and store in the batch list
            df_sample['Combination_ID'] = combination_id
            df_sample['Iteration'] = iteration
            batch_samples_list.append(df_sample)
            
            # b. Calculate sample statistics
            sample_var_x1 = df_sample['X1'].var(ddof=1) # Use ddof=1 for sample variance
            sample_var_x2 = df_sample['X2'].var(ddof=1)
            sample_std_x1 = np.sqrt(sample_var_x1)
            sample_std_x2 = np.sqrt(sample_var_x2)
            # Handle the case where df_sample[['X1', 'X2']].corr() might return a single value or a DataFrame
            corr_matrix = df_sample[['X1', 'X2']].corr()
            sample_corr_x1x2 = corr_matrix.iloc[0, 1] if not corr_matrix.empty and len(corr_matrix) > 1 else np.nan

            # c. Calculate rhoU (Reliability of Residual Scores) and rhoD (Reliability of Difference Scores)
            
            # rhoU (Reliability of Residual Scores Formula)
            denominator_U = 1 - sample_corr_x1x2**2
            if np.isclose(denominator_U, 0) or np.isnan(sample_corr_x1x2):
                rhoU = np.nan
            else:
                # This is the original formula from the Chinese code, kept for consistency:
                numerator_U = RELIABILITY_X1 + RELIABILITY_X2 * sample_corr_x1x2**2 - 2 * sample_corr_x1x2**2 
                rhoU = numerator_U / denominator_U

            # rhoD (Reliability of Difference Scores Formula)
            denominator_D = sample_var_x1 + sample_var_x2 - 2 * sample_std_x1 * sample_std_x2 * sample_corr_x1x2
            if np.isclose(denominator_D, 0) or np.isnan(sample_corr_x1x2):
                rhoD = np.nan
            else:
                # This is the original formula (approximation) from the Chinese code, kept for consistency:
                numerator_D_original = sample_var_x1 * RELIABILITY_X1 + sample_var_x2 * RELIABILITY_X2 - 2 * sample_std_x1 * sample_std_x2 * sample_corr_x1x2
                rhoD = numerator_D_original / denominator_D

            # d. Record statistics
            stats_dict = population_params.copy() 
            stats_dict.update({
                'Iteration': iteration,
                'Sample Mean X1': df_sample['X1'].mean(),
                'Sample Mean X2': df_sample['X2'].mean(),
                'Sample Var X1 (Sample)': sample_var_x1,
                'Sample Var X2 (Sample)': sample_var_x2,
                'Sample Corr (X1, X2)': sample_corr_x1x2,
                'rhoU': rhoU,
                'rhoD': rhoD
            })
            batch_stats_list.append(stats_dict)
    
    # Combine results for the current batch
    df_samples_chunk = pd.concat(batch_samples_list, ignore_index=True)
    df_stats_chunk = pd.DataFrame(batch_stats_list)
    
    Total_data_points_generated += len(df_samples_chunk)
    
    # Format file name (Part_01, Part_02, ...)
    part_number = f"Part_{i+1:02d}"
    
    raw_data_file_name = f'MCSdata_{timestamp}_{part_number}.csv'
    stats_file_name = f'multivariate_simulation_stats_{timestamp}_{part_number}.csv'
    
    # Save files
    df_samples_chunk.to_csv(raw_data_file_name, index=False, encoding='utf-8')
    df_stats_chunk.to_csv(stats_file_name, index=False, encoding='utf-8')
    
    print(f"✅ File Batch {i+1:02d} (containing {len(chunk_combinations)} combinations) saved successfully:")
    print(f"   - Raw Data: **{raw_data_file_name}** ({len(df_samples_chunk):,} records)")
    print(f"   - Statistics Summary: **{stats_file_name}** ({len(df_stats_chunk):,} records)")

# Final Report
print("\n==================================================")
print("             ✅ Monte Carlo Simulation Execution Complete ✅             ")
print("==================================================")
print(f"Total initial combinations: {N_total_combinations}")
print(f"Total valid combinations after constraint filtering: {num_valid_combinations} types")
print(f"Total simulations run ({N_ITERATIONS} iterations each): {Total_simulations_run:,}")
print(f"Total raw data points generated: {Total_data_points_generated:,} records")
print("==================================================")

 
Step2_Monte_Carlo_Simulation_20251126
import pandas as pd
import os
import glob
from datetime import datetime

# ======================================================================
# 🚀 PARAMETER CONFIGURATION 🚀
# ======================================================================
DATA_FOLDER = "multivariate_simulation_stats"
FILE_PATTERN = os.path.join(DATA_FOLDER, "multivariate_simulation_stats_*.csv") 
BATCH_SIZE = 20 
COMBINATION_ID_COLUMN = 'Combination_ID' 
# ======================================================================

def merge_files_in_batches_with_global_id(file_pattern, data_folder, batch_size, id_column):
    """
    Searches for all matching CSV files in the specified folder, reads them in batches,
    assigns a globally continuous sequential ID to the specified ID column, and
    outputs the result into multiple merged files.
    """
    
    if not os.path.isdir(data_folder):
        print(f"Error: Could not find folder '{data_folder}'. Please ensure it is in the same directory as this script.")
        return

    all_file_paths = glob.glob(file_pattern)
    
    if not all_file_paths:
        print(f"Error: No matching files found in the folder '{data_folder}' for pattern {os.path.basename(file_pattern)}.")
        return

    all_file_paths.sort()
    
    total_files = len(all_file_paths)
    print(f"🎉 Successfully found {total_files} files. Merging will proceed in batches of up to {batch_size} files.")
    
    global_id_counter = 0
    processed_count = 0
    batch_num = 1
    
    # ------------------------------------------------------------------
    # 🚀 LOOP AND PROCESS ALL FILES IN BATCHES 🚀
    # ------------------------------------------------------------------
    while processed_count < total_files:
        start_index = processed_count
        end_index = min(processed_count + batch_size, total_files)
        
        current_batch_paths = all_file_paths[start_index:end_index]
        current_batch_count_planned = len(current_batch_paths)

        print(f"\n--- 📦 Starting Batch {batch_num} (Planned Files: {current_batch_count_planned}) ---")
        
        df_list = []
        file_id_in_batch = 0
        files_successfully_read = 0
        
        # 2. Read, label, and concatenate files in the current batch
        for file_path in current_batch_paths:
            file_id_in_batch += 1 
            
            try:
                df_single = pd.read_csv(file_path)
                
                # --- Core Step: Generate Global Sequential ID ---
                row_count = len(df_single)
                
                new_ids = range(global_id_counter + 1, global_id_counter + row_count + 1)
                
                df_single[id_column] = list(new_ids)
                
                global_id_counter += row_count 
                
                # --- Original Labeling ---
                df_single['Source_File_ID'] = file_id_in_batch 
                df_single['Source_File_Name'] = os.path.basename(file_path)
                
                df_list.append(df_single)
                files_successfully_read += 1
                
            except Exception as e:
                print(f"   - ❌ CRITICAL WARNING: File read failed! {os.path.basename(file_path)} skipped. Error: {e}")
        
        if not df_list:
            print(f"Warning: All files in Batch {batch_num} failed to load or were empty. Skipping merge for this batch.")
            processed_count = end_index 
            batch_num += 1
            continue
            
        # 3. Concatenate all DataFrames vertically
        try:
            df_merged = pd.concat(df_list, ignore_index=True)
        except Exception as e:
            print(f"Error: An error occurred while merging Batch {batch_num} files. Please check for consistent file structure. Error: {e}")
            processed_count = end_index 
            batch_num += 1
            continue

        print(f"Batch {batch_num} files merged, total records: {len(df_merged)}")
        print(f"   - 📊 **Files Successfully Merged**: {files_successfully_read} / {current_batch_count_planned}")
        print(f"   - 🔑 {id_column} Range: {df_merged[id_column].min()} to {df_merged[id_column].max()}")

        # 4. Save results
        batch_str = str(batch_num).zfill(2)
        output_filename = f'multivariate_simulation_merged_Batch_{batch_str}_{datetime.now().strftime("%Y%m%d%H%M")}.csv'
        
        df_merged.to_csv(output_filename, index=False)
        print(f"Complete merged result saved to file: **{os.path.abspath(output_filename)}**")
        
        # Prepare for the next loop iteration
        processed_count = end_index
        batch_num += 1

# Execute the batch merge
merge_files_in_batches_with_global_id(
    FILE_PATTERN, 
    DATA_FOLDER, 
    BATCH_SIZE, 
    COMBINATION_ID_COLUMN
)
 
Step3_Table2_20251126

import pandas as pd
import numpy as np
import os
import glob
from datetime import datetime

# ======================================================================
# 🚀 PARAMETER AND PATH CONFIGURATION 🚀
# ======================================================================
# Folder name where merged files are stored (assumed to be in the same directory as this script)
MERGED_FOLDER = "multivariate_simulation_merged" 

# File pattern: matches all CSV files in the folder starting with 'multivariate_simulation_merged_Batch_'
FILE_PATTERN = os.path.join(MERGED_FOLDER, "multivariate_simulation_merged_Batch_*.csv") 

# Define the threshold needed for analysis (NOTE: Not used in the current P(rhoU > rhoD) calculation)
THRESHOLD = 0.7 

# Define the Population parameter columns required for analysis (must exist in the raw data)
# These columns define the pivot table's Index and Columns
POPULATION_INDEX_COL = 'Population Reliability X1'
POPULATION_COLUMNS_COL = 'Population Reliability X2'
# ======================================================================

def analyze_and_collect_stats(file_path, threshold):
    """
    Loads a single Batch CSV file, performs data cleaning, calculates the probability 
    statistic for each Combination_ID (P(rhoU > rhoD)), and returns the result for global merging.
    """
    file_name = os.path.basename(file_path)
    
    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        print(f"❌ Failed to read file {file_name}, skipping. Error: {e}")
        return None

    df.columns = df.columns.str.strip()
    # List of required columns for the calculation and grouping
    required_cols = ['Combination_ID', 'rhoU', 'rhoD', POPULATION_INDEX_COL, POPULATION_COLUMNS_COL]
    
    # Check if key columns exist
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"   ❌ Missing required columns: {missing_cols}, skipping.")
        return None

    # Data cleaning and type conversion
    df['rhoU'] = pd.to_numeric(df['rhoU'], errors='coerce')
    df['rhoD'] = pd.to_numeric(df['rhoD'], errors='coerce')
    # Remove rows where rhoU or rhoD calculation resulted in NaN
    df_cleaned = df.dropna(subset=['rhoU', 'rhoD'])

    if df_cleaned.empty:
        print("   Warning: Cleaned dataset is empty, skipping calculation.")
        return None

    # Create boolean column: True if rhoU > rhoD, False otherwise
    df_cleaned['Prob_rhoU_gt_rhoD'] = df_cleaned['rhoU'] > df_cleaned['rhoD']
    
    # Define columns for grouping
    group_cols = ['Combination_ID', POPULATION_INDEX_COL, POPULATION_COLUMNS_COL]

    # Perform group aggregation: Calculate the mean of the boolean column, which gives P(rhoU > rhoD)
    prob_stats = df_cleaned.groupby(group_cols).agg(
        Prob_rhoU_gt_rhoD=('Prob_rhoU_gt_rhoD', 'mean'),
    ).reset_index()

    return prob_stats

def generate_final_pivot_analysis():
    """
    Main function: Finds all Batch files, analyzes them in batches, merges results in memory, 
    and generates the final pivot table for P(rhoU > rhoD).
    """
    
    # 1. Find all files
    all_file_paths = glob.glob(FILE_PATTERN)
    all_file_paths.sort() 

    if not all_file_paths:
        print(f"❌ Error: No matching files found in the '{MERGED_FOLDER}' folder for pattern {os.path.basename(FILE_PATTERN)}.")
        return

    print(f"--- Starting Step-by-Step Analysis ---")
    print(f"Successfully found {len(all_file_paths)} Batch files ready for analysis.")

    # 2. Loop through all files and collect results
    all_results = []
    
    for i, file_path in enumerate(all_file_paths, 1):
        print(f"   [{i}/{len(all_file_paths)}] Analyzing {os.path.basename(file_path)}...")
        # Note: THRESHOLD is passed but not used in the current P(rhoU > rhoD) logic
        result = analyze_and_collect_stats(file_path, THRESHOLD) 
        if result is not None:
            all_results.append(result)
            
    if not all_results:
        print("\n❌ Warning: All found files could not be analyzed. Process terminated.")
        return

    # 3. Merge analysis results from all batches in memory
    print("\n--- Merging analysis statistics from all batches... ---")
    df_final_stats = pd.concat(all_results, ignore_index=True)
    
    # 4. Create the core pivot table
    print("--- Generating P(rhoU > rhoD) Pivot Table... ---")
    
    # Set index (Reliability X1), columns (Reliability X2), and use the calculated probability mean as values
    pivot_table = pd.pivot_table(
        df_final_stats, 
        values='Prob_rhoU_gt_rhoD', 
        index=POPULATION_INDEX_COL, 
        columns=POPULATION_COLUMNS_COL, 
        aggfunc='mean'
    )
    
    # --- 5. Calculate Marginal Aggregates (Sum and Mean) ---
    
    # Rename index and columns for better display
    pivot_table.index.name = 'Rhox1\\rhox2' 
    pivot_table.columns.name = 'RhoX2' 
    
    data_cols = [col for col in pivot_table.columns.tolist() if col not in ['Sum', 'Mean']]
    N_data_points = len(data_cols) # Number of columns (RhoX2 levels)

    # Calculate Row Sum and Row Mean (mean() automatically ignores NaN)
    pivot_table['Sum'] = pivot_table.sum(axis=1)
    pivot_table['Mean'] = pivot_table[data_cols].mean(axis=1)

    # Calculate Column Sum and Column Mean
    col_sum = pivot_table[data_cols].sum(axis=0)
    col_mean = pivot_table[data_cols].mean(axis=0)
    
    # Calculate Grand Mean 
    # NOTE: This calculation uses the total number of theoretical cells (len(pivot_table) * N_data_points), 
    # not the count of actual non-NaN cells, consistent with the original logic.
    total_sum = pivot_table['Sum'].sum()
    N_total_cells = len(pivot_table) * N_data_points
    grand_mean = total_sum / N_total_cells 

    # Construct Aggregate Rows
    sum_row = col_sum.to_dict()
    sum_row['Sum'] = total_sum
    sum_row['Mean'] = pivot_table['Mean'].sum() # Sum of means 
    
    mean_row = col_mean.to_dict()
    mean_row['Sum'] = grand_mean
    mean_row['Mean'] = grand_mean
    
    # Add aggregate rows to the pivot table
    pivot_table.loc['Sum'] = pd.Series(sum_row)
    pivot_table.loc['Mean'] = pd.Series(mean_row)

    # --- 6. Format Output and Save ---

    print("\n====================== P(rhoU > rhoD) Summary Pivot Table (Average Probability) ======================")
    print("Note: The data in the table are average probabilities aggregated across the fixed Var[X1], Var[X2], and RhoX1X2 conditions.")
    
    # Display results, formatted to 4 decimal places
    print(pivot_table.to_string(float_format='%.4f'))

    # Save results to CSV
    output_filename = f'pivot_P_rhoU_gt_rhoD_{datetime.now().strftime("%Y%m%d%H%M")}.csv'
    pivot_table.to_csv(output_filename, float_format='%.4f')
    
    print(f"\nComplete pivot table saved to file: **{os.path.abspath(output_filename)}**")
    print("=====================================================================================")

# Execute the analysis and pivot table generation
generate_final_pivot_analysis()

 
Step4_Table3_20251126

import pandas as pd
import numpy as np
import os
import glob
from datetime import datetime

# ======================================================================
# 🚀 PARAMETER AND PATH CONFIGURATION 🚀
# ======================================================================
# Folder name where merged files are stored (assumed to be in the same directory as this script)
MERGED_FOLDER = "multivariate_simulation_merged" 

# File pattern: matches all CSV files in the folder starting with 'multivariate_simulation_merged_Batch_'
FILE_PATTERN = os.path.join(MERGED_FOLDER, "multivariate_simulation_merged_Batch_*.csv") 

# Define the threshold needed for analysis (NOTE: Not used in this specific calculation)
THRESHOLD = 0.7 

# Define the required Population parameter columns
POPULATION_INDEX_COL = 'Population Reliability X1'
POPULATION_COLUMNS_COL = 'Population Reliability X2'
# ======================================================================

def analyze_and_collect_stats(file_path, threshold):
    """
    Loads a single Batch CSV file, performs data cleaning, calculates the Mean(rhoU - rhoD) 
    statistic for each Combination_ID, and returns the result for global merging.
    """
    file_name = os.path.basename(file_path)
    
    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        print(f"❌ Failed to read file {file_name}, skipping. Error: {e}")
        return None

    df.columns = df.columns.str.strip()
    required_cols = ['Combination_ID', 'rhoU', 'rhoD', POPULATION_INDEX_COL, POPULATION_COLUMNS_COL]
    
    # Check if key columns exist
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"   ❌ Missing required columns: {missing_cols}, skipping.")
        return None

    # Data cleaning and type conversion
    df['rhoU'] = pd.to_numeric(df['rhoU'], errors='coerce')
    df['rhoD'] = pd.to_numeric(df['rhoD'], errors='coerce')
    # Drop rows where rhoU or rhoD calculation resulted in NaN
    df_cleaned = df.dropna(subset=['rhoU', 'rhoD'])

    if df_cleaned.empty:
        print("   Warning: Cleaned dataset is empty, skipping calculation.")
        return None

    # --- Correction 1: Calculate the difference score (rhoU minus rhoD) ---
    DIFFERENCE_COL = 'RhoU_Minus_RhoD'
    TARGET_AGG_COL = 'Mean_RhoU_Minus_RhoD' # New output column name

    df_cleaned[DIFFERENCE_COL] = df_cleaned['rhoU'] - df_cleaned['rhoD']
    
    # Perform group aggregation: Calculate the mean of the difference score
    group_cols = ['Combination_ID', POPULATION_INDEX_COL, POPULATION_COLUMNS_COL]

    mean_stats = df_cleaned.groupby(group_cols).agg(
        **{TARGET_AGG_COL: (DIFFERENCE_COL, 'mean')} # Calculate the mean of the difference
    ).reset_index()

    return mean_stats

def generate_final_pivot_analysis():
    """
    Main function: Finds all Batch files, analyzes them in batches, merges results in memory, 
    and generates the final pivot table.
    """
    
    # 1. Find all files
    all_file_paths = glob.glob(FILE_PATTERN)
    all_file_paths.sort() 

    if not all_file_paths:
        print(f"❌ Error: No matching files found in the '{MERGED_FOLDER}' folder for pattern {os.path.basename(FILE_PATTERN)}.")
        return

    print(f"--- Starting Step-by-Step Analysis ---")
    print(f"Successfully found {len(all_file_paths)} Batch files ready for analysis.")

    # 2. Loop through all files and collect results
    all_results = []
    
    for i, file_path in enumerate(all_file_paths, 1):
        print(f"   [{i}/{len(all_file_paths)}] Analyzing {os.path.basename(file_path)}...")
        result = analyze_and_collect_stats(file_path, THRESHOLD)
        if result is not None:
            all_results.append(result)
            
    if not all_results:
        print("\n❌ Warning: All found files could not be analyzed. Process terminated.")
        return

    # 3. Merge analysis results from all batches in memory
    print("\n--- Merging analysis statistics from all batches... ---")
    df_final_stats = pd.concat(all_results, ignore_index=True)
    
    # 4. Create the core pivot table
    TARGET_VALUE_COLUMN = 'Mean_RhoU_Minus_RhoD'
    print(f"--- Generating {TARGET_VALUE_COLUMN} Pivot Table... ---")
    
    # Set index (Reliability X1), columns (Reliability X2), and use 'mean' as the aggregation function
    pivot_table = pd.pivot_table(
        df_final_stats, 
        values=TARGET_VALUE_COLUMN, # Use the mean difference value
        index=POPULATION_INDEX_COL, 
        columns=POPULATION_COLUMNS_COL, 
        aggfunc='mean'
    )
    
    # --- 5. Calculate Marginal Aggregates (Sum and Mean) ---
    
    # Rename index and columns for better display
    pivot_table.index.name = 'Rhox1\\rhox2' 
    pivot_table.columns.name = 'RhoX2' 
    
    # Ensure only data columns are included in marginal calculations
    data_cols = [col for col in pivot_table.columns.tolist() if col not in ['Sum', 'Mean']]
    N_data_points = len(data_cols)

    # Calculate Row Sum and Row Mean (sum and mean over RhoX2 levels)
    pivot_table['Sum'] = pivot_table[data_cols].sum(axis=1)
    pivot_table['Mean'] = pivot_table[data_cols].mean(axis=1)

    # Calculate Column Sum and Column Mean
    col_sum = pivot_table[data_cols].sum(axis=0)
    col_mean = pivot_table[data_cols].mean(axis=0)
    
    # Calculate Global Aggregates
    total_sum = pivot_table['Sum'].sum() # Sum of all Row Sums
    N_total_cells = len(pivot_table) * N_data_points
    grand_mean = total_sum / N_total_cells 

    # Construct Aggregate Rows
    sum_row = col_sum.to_dict()
    sum_row['Sum'] = total_sum
    sum_row['Mean'] = pivot_table['Mean'].sum() 
    
    mean_row = col_mean.to_dict()
    mean_row['Sum'] = grand_mean
    mean_row['Mean'] = grand_mean
    
    # Add aggregate rows to the pivot table
    pivot_table.loc['Sum'] = pd.Series(sum_row)
    pivot_table.loc['Mean'] = pd.Series(mean_row)

    # --- 6. Format Output and Save ---

    # --- Correction 2: Update Output Title and Filename ---
    print("\n====================== Mean(rhoU - rhoD) Summary Pivot Table (Average Difference) ======================")
    print("Note: The data in the table represents the average difference (RhoU minus RhoD) aggregated across the fixed Var[X1], Var[X2], and RhoX1X2 conditions.")
    
    # Display results, formatted to 4 decimal places
    print(pivot_table.to_string(float_format='%.4f'))

    # Save results to CSV
    output_filename = f'pivot_Mean_RhoU_Minus_RhoD_{datetime.now().strftime("%Y%m%d%H%M")}.csv'
    pivot_table.to_csv(output_filename, float_format='%.4f')
    
    print(f"\nComplete pivot table saved to file: **{os.path.abspath(output_filename)}**")
    print("=====================================================================================")

# Execute the analysis and pivot table generation
generate_final_pivot_analysis()
 
Step5_Table4_rhoU_20251126
# -*- coding: utf-8 -*-
"""
Multivariate Simulation Analysis Script
Calculates the probability of a target reliability/correlation measure 
exceeding a specified threshold (e.g., P(rhoU >= 0.7)) across simulation batches.
"""

import pandas as pd
import numpy as np
import os
import glob
from datetime import datetime

# ======================================================================
# PARAMETERS (EDIT ONLY THESE LINES!)
# ======================================================================
MERGED_FOLDER = "multivariate_simulation_merged"
FILE_PATTERN = os.path.join(MERGED_FOLDER, "multivariate_simulation_merged_Batch_*.csv")

ANALYSIS_TARGET = "rhoU"         # Change to "rhoU", "rhoD", or "rhoL"
THRESHOLD = 0.7

COL_RHO_X1      = 'Population Reliability X1'
COL_RHO_X2      = 'Population Reliability X2'
COL_CORR_X1X2   = 'Population Corr X1X2'
# ======================================================================

def analyze_batch(filepath):
    """
    Reads a batch file, calculates the hit rate (1 if TARGET >= THRESHOLD, 0 otherwise)
    and averages the hit rate across all combinations.
    """
    try:
        df = pd.read_csv(filepath)
    except:
        return None
    # Strip whitespace from column names for consistency
    df.columns = df.columns.str.strip()

    req = ['Combination_ID', ANALYSIS_TARGET, COL_RHO_X1, COL_RHO_X2, COL_CORR_X1X2]
    # Check if all required columns are present
    if not all(c in df.columns for c in req):
        print(f"Skipping {os.path.basename(filepath)}: Missing required columns.")
        return None

    # Drop rows where the analysis target is NaN
    df = df.dropna(subset=[ANALYSIS_TARGET])
    # Create the binary hit variable
    df['hit'] = (df[ANALYSIS_TARGET] >= THRESHOLD).astype(float)

    # Group by combination parameters and calculate the probability (mean of 'hit')
    grouped = df.groupby(['Combination_ID', COL_RHO_X1, COL_RHO_X2, COL_CORR_X1X2])['hit'].mean()
    # Reset index and name the resulting probability column 'prob'
    return grouped.reset_index(name='prob')

# ====================== MAIN PROGRAM ======================
# Find all batch files matching the pattern
files = sorted(glob.glob(FILE_PATTERN))
if not files:
    raise FileNotFoundError(f"Files not found: {FILE_PATTERN}")

print(f"Found {len(files)} batch files. Starting analysis of P({ANALYSIS_TARGET} ≥ {THRESHOLD}) ...\n")
all_data = []

# Process each file
for i, f in enumerate(files, 1):
    print(f"   [{i:>3}/{len(files)}] {os.path.basename(f)}")
    tmp = analyze_batch(f)
    if tmp is not None:
        all_data.append(tmp)

if not all_data:
    raise ValueError("No valid data found from any batch file.")

# Concatenate results from all batches
data = pd.concat(all_data, ignore_index=True)

# ====================== 1. FULL 9×9×9 DETAIL TABLE ======================
# Create a pivot table for the full breakdown
pivot_full = pd.pivot_table(
    data,
    values='prob',
    index=COL_RHO_X1,
    columns=[COL_RHO_X2, COL_CORR_X1X2],
    aggfunc='mean'
)

# Force completion of 81 columns (9 Rho_X2 levels * 9 Rho_X1X2 levels)
levels = np.round(np.linspace(0.1, 0.9, 9), 1)
full_cols = pd.MultiIndex.from_product([levels, levels], names=['Rho_X2', 'Rho_X1X2'])
pivot_full = pivot_full.reindex(columns=full_cols)

# Calculate the accurate Grand Mean (using only the 81 data cells)
data_grid = pivot_full.copy()
total_sum_full = data_grid.sum().sum(skipna=True)
n_cells_full   = data_grid.count().sum()
grand_mean = total_sum_full / n_cells_full if n_cells_full > 0 else np.nan

# Calculate Row Sum / Mean
data_grid[('Sum',   '')] = data_grid.sum(axis=1, skipna=True)
data_grid[('Mean', '')] = data_grid.mean(axis=1, skipna=True)

# Calculate Column Sum / Mean
col_sum  = data_grid.iloc[:, :81].sum(axis=0, skipna=True)
col_mean = data_grid.iloc[:, :81].mean(axis=0, skipna=True)

# Prepare Sum and Mean rows for concatenation
sum_row  = col_sum.copy()
mean_row = col_mean.copy()
sum_row.name  = 'Sum'
mean_row.name = 'Mean'
sum_row[('Sum',   '')] = data_grid[('Sum', '')].sum()
sum_row[('Mean', '')] = data_grid[('Mean', '')].sum() # Note: This sum is just for visual completeness, not meaningful
mean_row[('Sum',   '')] = grand_mean
mean_row[('Mean', '')] = grand_mean

# Combine the data grid with the calculated Sum and Mean rows
final_full = pd.concat([data_grid, pd.DataFrame([sum_row, mean_row])])
final_full.columns.names = ['Rho_X2', 'Rho_X1X2']
final_full.index.name = 'Rho_X1'

# Final sorting (ensure index is numeric for sorting)
data_part = final_full.iloc[:-2].copy()
# Convert index to numeric for sorting purposes
data_part.index = pd.to_numeric(data_part.index) 
# Sort by Rho_X1 (index) and then by the first level of the columns (Rho_X2)
data_part = data_part.sort_index().sort_index(axis=1, level=0)
final_full = pd.concat([data_part, final_full.iloc[-2:]])

# ====================== 2. CORRECT 9×9 AVERAGE TABLE (Key Fix!) ======================
# Method: Unstack the full table and average over the Rho_X1X2 dimension
temp = data_grid.iloc[:, :81].copy()  # Only take the 81 data columns

# Convert MultiIndex columns to regular columns (un-pivoting)
# This is a safe way to correctly average out Rho_X1X2
temp = temp.stack(['Rho_X2', 'Rho_X1X2'])
temp = temp.reset_index()
temp.columns = ['Rho_X1', 'Rho_X2', 'Rho_X1X2', 'value']

# Correctly average out the Rho_X1X2 effect and re-pivot to 9x9
table_9x9 = temp.groupby(['Rho_X1', 'Rho_X2'])['value'].mean().unstack('Rho_X2')

# Reindex to ensure all 0.1 to 0.9 levels are present
table_9x9 = table_9x9.reindex(index=levels, columns=levels)

# Add Marginals (Row Mean and Column Mean)
table_9x9['Row_Mean'] = table_9x9.mean(axis=1, skipna=True)
table_9x9.loc['Col_Mean'] = table_9x9.mean(axis=0, skipna=True)
# Set the Grand Mean at the intersection
table_9x9.loc['Col_Mean', 'Row_Mean'] = grand_mean

table_9x9.columns.name = 'Rho_X2'
table_9x9.index.name   = 'Rho_X1'
table_9x9 = table_9x9.round(4)

# ====================== OUTPUT ======================
timestamp = datetime.now().strftime("%Y%m%d_%H%M")

print("\n" + "="*120)
print(f"Analysis Complete! P({ANALYSIS_TARGET} ≥ {THRESHOLD})")
print(f"Actual Number of Cells Included: {int(n_cells_full)} / 729")
print(f"Grand Mean = {grand_mean:.6f}   (Consistent across both tables!)")
print("="*120)

# Full Table Output
print("\n1. Full 9×9×9 Table")
print(final_full.to_string(float_format=lambda x: "" if pd.isna(x) else f"{x:.4f}"))
file1 = f"1_Full729_P_{ANALYSIS_TARGET}_GE_{THRESHOLD}_{timestamp}.csv"
final_full.to_csv(file1, float_format='%.6f')

# 9x9 Average Table Output
print("\n\n2. Concise 9×9 Average Table (For Main Text, e.g., Table 4)")
print(table_9x9.to_string(float_format=lambda x: "" if pd.isna(x) else f"{x:.4f}"))
file2 = f"2_Table9x9_P_{ANALYSIS_TARGET}_GE_{THRESHOLD}_{timestamp}.csv"
table_9x9.to_csv(file2, float_format='%.4f')

print(f"\nFull table saved to → {os.path.abspath(file1)}")
print(f"9x9 table saved to → {os.path.abspath(file2)}")
print("\nDual tables complete!")


 
Step6_Table5_rhoD_20251126

# -*- coding: utf-8 -*-
"""
Multivariate Simulation Analysis Script (Fixed for P(rhoD >= 0.7))

Calculates the probability of the target discrimination reliability (rhoD) 
exceeding a specified threshold (0.7) across various combinations of 
population reliabilities (Rho_X1, Rho_X2) and correlation (Rho_X1X2).
"""

import pandas as pd
import numpy as np
import os
import glob
from datetime import datetime

# ======================================================================
# CONFIGURATION (Set for: rhoD ≥ 0.7)
# ======================================================================
MERGED_FOLDER = "multivariate_simulation_merged"
FILE_PATTERN = os.path.join(MERGED_FOLDER, "multivariate_simulation_merged_Batch_*.csv")

ANALYSIS_TARGET = "rhoD"         # Fixed analysis target: rhoD
THRESHOLD = 0.7                  # Fixed threshold: >= 0.7

COL_RHO_X1      = 'Population Reliability X1'
COL_RHO_X2      = 'Population Reliability X2'
COL_CORR_X1X2   = 'Population Corr X1X2'
# ======================================================================

def analyze_batch(filepath):
    """
    Reads a batch file, calculates the 'hit' variable (1 if TARGET >= THRESHOLD),
    and averages the hit rate across all combination parameters.
    """
    try:
        df = pd.read_csv(filepath)
    except:
        # Returns None if file reading fails
        return None
        
    # Standardize column names
    df.columns = df.columns.str.strip()

    req = ['Combination_ID', ANALYSIS_TARGET, COL_RHO_X1, COL_RHO_X2, COL_CORR_X1X2]
    # Check if all required columns are present
    if not all(c in df.columns for c in req):
        return None

    # Drop rows where the analysis target (rhoD) is NaN
    df = df.dropna(subset=[ANALYSIS_TARGET])
    
    # Key change: Define 'hit' as true if rhoD is >= THRESHOLD
    df['hit'] = (df[ANALYSIS_TARGET] >= THRESHOLD).astype(float)

    # Group by combination parameters and calculate the probability (mean of 'hit')
    grouped = df.groupby(['Combination_ID', COL_RHO_X1, COL_RHO_X2, COL_CORR_X1X2])['hit'].mean()
    # Reset index and name the resulting probability column 'prob'
    return grouped.reset_index(name='prob')

# ====================== MAIN PROGRAM ======================
files = sorted(glob.glob(FILE_PATTERN))
if not files:
    raise FileNotFoundError(f"No files found matching the pattern: {FILE_PATTERN}")

print(f"Analyzing P({ANALYSIS_TARGET} ≥ {THRESHOLD}). Found {len(files)} batch files...\n")
all_data = []

# Process each file
for i, f in enumerate(files, 1):
    print(f"   [{i:>3}/{len(files)}] Processing → {os.path.basename(f)}")
    tmp = analyze_batch(f)
    if tmp is not None:
        all_data.append(tmp)

if not all_data:
    raise ValueError("No valid data found across all files.")

# Concatenate results from all batches
data = pd.concat(all_data, ignore_index=True)

# ====================== 1. FULL 9×9×9 DETAIL TABLE ======================
# Create a pivot table for the full breakdown
pivot_full = pd.pivot_table(
    data,
    values='prob',
    index=COL_RHO_X1,
    columns=[COL_RHO_X2, COL_CORR_X1X2],
    aggfunc='mean'
)

# Define the nine levels (0.1 to 0.9)
levels = np.round(np.linspace(0.1, 0.9, 9), 1)
# Create MultiIndex for 9 Rho_X2 levels * 9 Rho_X1X2 levels
full_cols = pd.MultiIndex.from_product([levels, levels], names=['Rho_X2', 'Rho_X1X2'])
pivot_full = pivot_full.reindex(columns=full_cols)

# Calculate Grand Mean
data_grid = pivot_full.copy()
total_sum_full = data_grid.sum().sum(skipna=True)
n_cells_full   = data_grid.count().sum()
grand_mean = total_sum_full / n_cells_full if n_cells_full > 0 else np.nan

# Calculate Row Sum / Mean
data_grid[('Sum',   '')] = data_grid.sum(axis=1, skipna=True)
data_grid[('Mean', '')] = data_grid.mean(axis=1, skipna=True)

# Calculate Column Sum / Mean (Bottom rows)
col_sum  = data_grid.iloc[:, :81].sum(axis=0, skipna=True)
col_mean = data_grid.iloc[:, :81].mean(axis=0, skipna=True)

# Prepare Sum and Mean rows for concatenation
sum_row  = col_sum.copy()
mean_row = col_mean.copy()
sum_row.name  = 'Sum'
mean_row.name = 'Mean'
sum_row[('Sum',   '')] = data_grid[('Sum', '')].sum()
sum_row[('Mean', '')] = data_grid[('Mean', '')].sum()
mean_row[('Sum',   '')] = grand_mean
mean_row[('Mean', '')] = grand_mean

# Combine the data grid with the calculated Sum and Mean rows
final_full = pd.concat([data_grid, pd.DataFrame([sum_row, mean_row])])
final_full.columns.names = ['Rho_X2', 'Rho_X1X2']
final_full.index.name = 'Rho_X1'

# Final sorting
data_part = final_full.iloc[:-2].copy()
data_part.index = pd.to_numeric(data_part.index)
data_part = data_part.sort_index().sort_index(axis=1, level=0)
final_full = pd.concat([data_part, final_full.iloc[-2:]])

# ====================== 2. CORRECT 9×9 AVERAGE TABLE ======================
# Method: Unstack the full table and average over the Rho_X1X2 dimension
temp = data_grid.iloc[:, :81].copy()
# Un-pivot the MultiIndex columns to separate rows (future_stack=True avoids warnings)
temp = temp.stack(['Rho_X2', 'Rho_X1X2'], future_stack=True).reset_index()
temp.columns = ['Rho_X1', 'Rho_X2', 'Rho_X1X2', 'value']

# Correctly average out the Rho_X1X2 effect and re-pivot to 9x9
table_9x9 = temp.groupby(['Rho_X1', 'Rho_X2'])['value'].mean().unstack('Rho_X2')
table_9x9 = table_9x9.reindex(index=levels, columns=levels)

# Add Marginals (Row Mean and Column Mean)
table_9x9['Row_Mean'] = table_9x9.mean(axis=1, skipna=True)
table_9x9.loc['Col_Mean'] = table_9x9.mean(axis=0, skipna=True)
# Set the Grand Mean at the intersection
table_9x9.loc['Col_Mean', 'Row_Mean'] = grand_mean

table_9x9.columns.name = 'Rho_X2'
table_9x9.index.name   = 'Rho_X1'
table_9x9 = table_9x9.round(4)

# ====================== OUTPUT ======================
timestamp = datetime.now().strftime("%Y%m%d_%H%M")

print("\n" + "="*120)
print(f"   Analysis Complete! P({ANALYSIS_TARGET} ≥ {THRESHOLD})")
print(f"   Actual Number of Cells Included: {int(n_cells_full)} / 729")
print(f"   Grand Mean = {grand_mean:.6f}   (Consistent across both tables!)")
print("="*120)

# 1. Full Table Output (Suggested for Appendix)
print("\n1. Full 9×9×9 Detail Table (Suggested for Appendix)")
print(final_full.to_string(float_format=lambda x: "" if pd.isna(x) else f"{x:.4f}"))
file1 = f"1_Full729_{ANALYSIS_TARGET}_GE_{THRESHOLD}_{timestamp}.csv"
final_full.to_csv(file1, float_format='%.6f')

# 2. 9x9 Average Table Output (Suggested for Main Text Table)
print("\n\n2. Concise 9×9 Average Table (Suggested for Main Text Table)")
print(table_9x9.to_string(float_format=lambda x: "" if pd.isna(x) else f"{x:.4f}"))
file2 = f"2_Table9x9_{ANALYSIS_TARGET}_GE_{THRESHOLD}_{timestamp}.csv"
table_9x9.to_csv(file2, float_format='%.4f')

print(f"\nFull table saved to → {os.path.abspath(file1)}")
print(f"9x9 table saved to → {os.path.abspath(file2)}")
print("\nAnalysis finished!")

