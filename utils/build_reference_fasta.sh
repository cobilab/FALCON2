#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 SOURCE_DIR OUTPUT_FILE"
    echo "  SOURCE_DIR   : Directory containing .fna files"
    echo "  OUTPUT_FILE  : Output file to store reference sequences"
    exit 1
}

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    usage
fi

# Assign arguments to variables
SOURCE_DIR="$1"
OUTPUT_FILE="$2"

# Check if the source directory exists
if [ ! -d "$SOURCE_DIR" ]; then
    echo "Error: Source directory '$SOURCE_DIR' does not exist."
    exit 1
fi

# Check if the source directory contains .fna files
if [ -z "$(find "$SOURCE_DIR" -maxdepth 1 -name '*.fna' -print -quit)" ]; then
    echo "Error: No .fna files found in '$SOURCE_DIR'."
    exit 1
fi

# Create or clear the output file
> "$OUTPUT_FILE"

# Process each .fna file in the source directory
echo "Building reference FASTA: $OUTPUT_FILE..."
for FILE in "$SOURCE_DIR"/*.fna; do
    # Append the contents of the .fna file to the output file
    cat "$FILE" >> "$OUTPUT_FILE"
done

echo "Reference sequences created: $OUTPUT_FILE"