#!/usr/bin/env bash
set -Eeuo pipefail

# Usage: ./download_references_ncbi.sh [organism]
# Example: ./download_references_ncbi.sh bacteria
# Defaults to 'viruses' if no argument given

organism="${1:-viruses}"
query="${organism}[Organism] AND 'complete genome'[filter] AND latest_refseq[filter] AND (latest[filter] NOT anomalous[filter])"

echo "Querying NCBI Assembly for: $query"

mkdir -p "$organism"
cd "$organism"

# Step 1: ESearch
esearch -db assembly -query "$query" > "${organism}_search.xml"

# Step 2: ESummary
cat "${organism}_search.xml" | esummary -db assembly > "${organism}_summary.xml"

# Step 3: Extract metadata
cat "${organism}_summary.xml" | xtract \
    -pattern DocumentSummary \
    -element AssemblyAccession Organism FtpPath_RefSeq Taxid \
    -block 'Meta/Stats/Stat' \
    -if '@category' -equals 'total_length' \
    -element '.' \
| sed -E 's/\tStat +"?([0-9]+)"?/\t\1/g' \
| awk -F'\t' 'NF==5' \
> "${organism}.tsv"

# Step 4: Deduplicate and sort by genome size
sort -t $'\t' -k5,5n --stable "${organism}.tsv" \
| sort -t $'\t' -k4,4 -u --stable \
| sort -t $'\t' -k5,5n --stable \
> "${organism}_sorted.tsv"

# Step 5: Download via HTTPS
output_dir="reference_fastas"
mkdir -p "$output_dir"

echo "Starting HTTPS downloads..."
while IFS=$'\t' read -r acc org ftp taxid size; do
  if [[ -n "$ftp" && "$ftp" == ftp://* ]]; then
    # Convert FTP to HTTPS
    https_url="${ftp/ftp:/https:}"
    file_name="$(basename "$ftp")_genomic.fna.gz"
    full_url="${https_url}/${file_name}"

    echo "Downloading $file_name from $full_url ..."
    wget -q -P "$output_dir" "$full_url" || echo "Failed to download: $file_name"
  fi
done < "${organism}_sorted.tsv"

echo "All downloads complete."
echo "FASTA files saved in: ${organism}/${output_dir}"
