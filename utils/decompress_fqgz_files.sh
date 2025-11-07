#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <source_directory> <destination_directory>"
  exit 1
fi

# Input arguments
SOURCE_DIR="$1"
DEST_DIR="$2"

# Check if the source directory exists
if [ ! -d "$SOURCE_DIR" ]; then
  echo "Source directory $SOURCE_DIR does not exist!"
  exit 1
fi

# Create the destination directory if it doesn't exist
if [ ! -d "$DEST_DIR" ]; then
  echo "Destination directory $DEST_DIR does not exist. Creating it..."
  mkdir -p "$DEST_DIR"
fi

# Loop through all .fq.gz files in the source directory
for FILE in "$SOURCE_DIR"/*.fq.gz; do
  if [ -f "$FILE" ]; then
    # Extract the file name without the .gz extension
    BASENAME=$(basename "$FILE" .gz)
    DEST_FILE="$DEST_DIR/$BASENAME"

    # Decompress the file into the destination directory
    echo "Decompressing $FILE to $DEST_FILE..."
    gunzip -c "$FILE" > "$DEST_FILE"
    
    if [ $? -eq 0 ]; then
      echo "Successfully decompressed $FILE to $DEST_FILE."
    else
      echo "Failed to decompress $FILE."
    fi
  else
    echo "No .fq.gz files found in $SOURCE_DIR."
  fi
done

echo "All decompression tasks completed."
