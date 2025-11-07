#!/usr/bin/env bash
set -Eeuo pipefail

# Usage:
#   ./process_gz_files.sh <input_directory> <output_file>

# Description:
#   This script processes all `.gz` files in the specified directory:
#   1) Validates the integrity of each `.gz` file using `gzip -t`.
#   2) Skips corrupted files.
#   3) Decompresses and concatenates valid files into the specified output file.

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input_directory> <output_file>"
  exit 1
fi

INPUT_DIR="$1"
OUTPUT_FILE="$2"

# Check if input directory exists
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: Input directory '$INPUT_DIR' does not exist."
  exit 1
fi

# Create or empty the output file
> "$OUTPUT_FILE"

echo "Processing files in: $INPUT_DIR"
echo "Output will be written to: $OUTPUT_FILE"

# Process each .gz file in the input directory
for file in "$INPUT_DIR"/*.gz; do
  if [[ ! -f "$file" ]]; then
    echo "No .gz files found in the directory."
    break
  fi

  # echo "Validating: $file"

  # Check file integrity
  if ! gzip -t "$file" 2>/dev/null; then
    echo "Corrupted file detected: $file"
    continue
  fi

  #echo "Decompressing: $file"
  zcat "$file" >> "$OUTPUT_FILE"
done

echo "Processing completed."
