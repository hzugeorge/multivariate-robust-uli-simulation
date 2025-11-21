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

