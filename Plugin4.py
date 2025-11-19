import os
import glob
import numpy as np
import pandas as pd

def process_csv_file(file_path, output_path):
    # Read the CSV file (assumes the first row is a header)
    df = pd.read_csv(file_path, header=0)
    # Note: the CSV is assumed to have 54 data rows (indices 0 to 53)
    
    # ==== STEP 1: Reorder the Rows ====
    # Identify the two rows (original indices 50 and 51, i.e. Excel rows 52 and 53)
    summary_indices = [50, 51]
    summary_df = df.loc[summary_indices]
    
    # Remove these rows from the original DataFrame
    remaining_df = df.drop(summary_indices)
    
    # Prepend the summary rows so that they become the new first two rows.
    # The final DataFrame 'new_df' will still have 54 rows.
    new_df = pd.concat([summary_df, remaining_df], ignore_index=True)
    
    # ==== STEP 2: Process Each FOV Group Over All 54 Data Rows ====
    # It is assumed that after the first column (which is an index/increasing number),
    # each FOV group is comprised of 6 columns in the order:
    # [Area, Mean, StdDev, Major, Minor, Angle].
    total_cols = new_df.shape[1]
    # Here we assume the CSV already contains Column A as the increasing numbers.
    # If not, adjust the starting index accordingly.
    num_fovs = (total_cols - 1) // 6  # Subtract 1 for the initial index column
    
    # Dictionary to collect processed FOV data
    selected_fov_data = {}
    
    # Use all 54 rows of new_df (indices 0 to 53)
    n_rows = new_df.shape[0]  # Should be 54
    
    for i in range(num_fovs):
        base = 1 + i * 6  # Starting index for this FOV group (assuming first column is already there)
        col_area   = new_df.columns[base]       # Area
        col_mean   = new_df.columns[base + 1]     # Mean
        col_stdev  = new_df.columns[base + 2]     # StdDev
        col_major  = new_df.columns[base + 3]     # Major
        col_minor  = new_df.columns[base + 4]     # Minor
        # The column for Angle (base + 5) is not used
        
        # Compute the ratio (Mean / StdDev) across ALL 54 rows
        ratio = new_df[col_mean] / new_df[col_stdev].replace(0, np.nan)
        avg_ratio = ratio.mean(skipna=True)
        
        # Process only if the average ratio is at least 0.8
        if avg_ratio >= 0.8:
            # Use all 54 rows for this FOV's data:
            mean_data = new_df[col_mean].iloc[:n_rows].reset_index(drop=True)
            stdev_data = new_df[col_stdev].iloc[:n_rows].reset_index(drop=True)
            
            # ==== Append Summary Values to the Mean Column ====
            # Use the summary values (from the first row of the moved rows)
            summary_area  = new_df[col_area].iloc[0]
            summary_major = new_df[col_major].iloc[0]
            summary_minor = new_df[col_minor].iloc[0]
            summary_series = pd.Series([summary_area, summary_major, summary_minor])
            
            # Append these summary values to the Mean column (making its length 54 + 3 = 57)
            mean_data = pd.concat([mean_data, summary_series], ignore_index=True)
            
            # For the StdDev column, append three blank cells so its length matches 57
            stdev_data = pd.concat([stdev_data, pd.Series([""] * 3)], ignore_index=True)
            
            # Save the processed data for this FOV group
            selected_fov_data[f"FOV{i+1}"] = pd.DataFrame({
                "Mean": mean_data,
                "StdDev": stdev_data
            })
    
    # ==== STEP 3: Create the Final Output File ====
    # The final output will have 57 rows per FOV group (54 data + 3 summary)
    final_nrows = n_rows + 3  # 54 + 3 = 57
    
    final_df = pd.DataFrame()
    
    # Create the extra first column with the following values:
    # First 54 rows: numbers 1 through 54; last 3 rows: "Area", "Major", "Minor"
    extra_col = pd.Series(list(range(1, n_rows + 1)) + ["Area", "Major", "Minor"])
    final_df["Index"] = extra_col
    
    # Add each selected FOV group (Mean and StdDev columns) side by side
    for fov, fov_df in selected_fov_data.items():
        final_df[f"{fov}_Mean"] = fov_df["Mean"]
        final_df[f"{fov}_StdDev"] = fov_df["StdDev"]
    
    # Write the final DataFrame to a new CSV file.
    final_df.to_csv(output_path, index=False)
    print(f"Processed file saved to: {output_path}")

# ==== Define Source and Destination Folders ====
input_folder = r"C:\Users\WerlandL\FRET"
output_folder = r"C:\Users\WerlandL\FRET2"
os.makedirs(output_folder, exist_ok=True)

# Process every CSV file in the source folder.
for file in glob.glob(os.path.join(input_folder, "*.csv")):
    filename = os.path.basename(file)
    output_file = os.path.join(output_folder, filename)
    print(f"Processing: {filename}")
    process_csv_file(file, output_file)
