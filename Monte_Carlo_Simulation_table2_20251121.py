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