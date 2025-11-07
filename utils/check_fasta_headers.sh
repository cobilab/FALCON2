#!/usr/bin/env bash
set -Eeuo pipefail

###############################################################################
# check_fasta_headers.sh
#
# Usage:
#   ./check_fasta_headers.sh my_sequences.fna.gz
#   ./check_fasta_headers.sh uncompressed.fna
#
# Description:
#   - Detects if the input file is gzipped (ends with .gz)
#   - If gzipped, runs "gunzip -c" to uncompress on-the-fly
#   - Pipes resulting FASTA data into check_fasta_headers.awk
###############################################################################

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <fna or fna.gz file>"
  exit 1
fi

input_file="$1"
script_dir="$(dirname "$0")"

if [[ "$input_file" == *.gz ]]; then
  echo "Processing gzipped FASTA: $input_file"
  gunzip -c "$input_file" | "$script_dir/check_fasta_headers.awk"
else
  echo "Processing uncompressed FASTA: $input_file"
  "$script_dir/check_fasta_headers.awk" < "$input_file"
fi
