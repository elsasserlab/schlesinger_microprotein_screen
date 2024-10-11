# Load necessary libraries
library(readr)
library(writexl)

# Specify the folder containing the TSV files and the output Excel file name
folder_path <- "./"
output_excel_file <- "output.xlsx"

# List all TSV files in the specified folder
tsv_files <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)

# Initialize an empty list to store data frames
sheets <- list()

# Loop through each TSV file and read it into a data frame
for (file in tsv_files) {
  # Read the TSV file into a data frame
  df <- read_tsv(file)
  # Use the file name (without extension) as the sheet name
  sheet_name <- tools::file_path_sans_ext(basename(file))
  # Add the data frame to the list of sheets
  sheets[[sheet_name]] <- df
}

# Write the list of data frames to an Excel file with multiple sheets
write_xlsx(sheets, output_excel_file)

cat("Excel file '", output_excel_file, "' has been created with sheets from TSV files in the folder '", folder_path, "'.\n", sep = "")
