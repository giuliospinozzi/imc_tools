import os
import pandas as pd

def load_csv(file_path):
    """Load the CSV file into a DataFrame."""
    df = pd.read_csv(file_path)
    return df

def calculate_mean_std(df):
    """Calculate mean and standard deviation for each marker."""
    marker_means = df.mean()
    marker_stds = df.std()
    return marker_means, marker_stds

def filter_markers(df, markers):
    """Filter the DataFrame to only include the specified markers."""
    return df[markers]

def calculate_coexpression_percentage(df, markers, marker_means):
    """Calculate the percentage of cells co-expressing the specified markers."""
    coexpressing_cells = df[(df[markers] >= marker_means[markers]).all(axis=1)]
    coexpression_percentage = (len(coexpressing_cells) / len(df)) * 100
    return coexpression_percentage

def process_roi_files(folder_path, marker_combinations, output_file):
    overall_means = []
    overall_stds = []
    coexpression_results = []
    
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".csv"):
            file_path = os.path.join(folder_path, file_name)
            df = load_csv(file_path)
            df = df.set_index("Object")  # Use the Object column as the index
            marker_means, marker_stds = calculate_mean_std(df)
            
            # Store overall means
            mean_dict = {"ROI": file_name}
            mean_dict.update(marker_means.to_dict())
            overall_means.append(mean_dict)
            
            # Store overall standard deviations
            std_dict = {"ROI": file_name}
            std_dict.update(marker_stds.to_dict())
            overall_stds.append(std_dict)
            
            for markers in marker_combinations:
                filtered_df = filter_markers(df, markers)
                coexpression_percentage = calculate_coexpression_percentage(filtered_df, markers, marker_means)
                
                coexpression_results.append({
                    "ROI": file_name,
                    "Markers": ", ".join(markers),
                    "Co-expression Percentage": coexpression_percentage
                })
    
    overall_means_df = pd.DataFrame(overall_means)
    overall_stds_df = pd.DataFrame(overall_stds)
    coexpression_results_df = pd.DataFrame(coexpression_results)
    
    with pd.ExcelWriter(output_file) as writer:
        overall_means_df.to_excel(writer, sheet_name="Means", index=False)
        overall_stds_df.to_excel(writer, sheet_name="Standard Deviations", index=False)
        coexpression_results_df.to_excel(writer, sheet_name="Co-expression Results", index=False)

def main(folder_path, marker_combinations, output_file):
    process_roi_files(folder_path, marker_combinations, output_file)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    # Path to the folder containing the CSV files
    folder_path = "C:\\Users\\39351\\OneDrive\\Desktop\\intensities"
    
    # Combinations of markers of interest
    marker_combinations = [
        ["PDGFRa", "Olig2","H3K27M"],
        ["BIIItubulin", "Olig2","H3K27M"],
        ["GFAP", "Tenascin","H3K27M"],
        ["Vimentin", "CD44","H3K27M"],
        ["EGFR", "Ki67", "Caspase3","H3K27M"]
    ]
    
    # Output file path
    output_file = "C:\\Users\\39351\\OneDrive\\Desktop\\intensities\\results.xlsx"
    
    main(folder_path, marker_combinations, output_file)
