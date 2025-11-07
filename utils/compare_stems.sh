#!/usr/bin/env bash
# Compare two directories recursively and print ONLY missing files,
# matching by stem where stem = relative path without ".fna" or ".fna.gz".
# Only .fna and .fna.gz are considered.

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 DIR1 DIR2" >&2
  exit 1
fi

dir1="$(cd "$1" && pwd)"
dir2="$(cd "$2" && pwd)"

declare -A stems1 stems2

# Build stem sets (only .fna / .fna.gz)
while IFS= read -r -d '' rel; do
  s="$rel"
  if [[ "$s" == *.fna.gz ]]; then s="${s%.fna.gz}"; else s="${s%.fna}"; fi
  stems1["$s"]=1
done < <(cd "$dir1" && find . -type f \( -name "*.fna" -o -name "*.fna.gz" \) -printf '%P\0')

while IFS= read -r -d '' rel; do
  s="$rel"
  if [[ "$s" == *.fna.gz ]]; then s="${s%.fna.gz}"; else s="${s%.fna}"; fi
  stems2["$s"]=1
done < <(cd "$dir2" && find . -type f \( -name "*.fna" -o -name "*.fna.gz" \) -printf '%P\0')

# Files in DIR1 missing in DIR2 (by stem)
missing_in_dir2=()
while IFS= read -r -d '' rel; do
  s="$rel"
  if [[ "$s" == *.fna.gz ]]; then s="${s%.fna.gz}"; else s="${s%.fna}"; fi
  [[ -z "${stems2[$s]:-}" ]] && missing_in_dir2+=("$rel")
done < <(cd "$dir1" && find . -type f \( -name "*.fna" -o -name "*.fna.gz" \) -printf '%P\0')

# Files in DIR2 missing in DIR1 (by stem)
missing_in_dir1=()
while IFS= read -r -d '' rel; do
  s="$rel"
  if [[ "$s" == *.fna.gz ]]; then s="${s%.fna.gz}"; else s="${s%.fna}"; fi
  [[ -z "${stems1[$s]:-}" ]] && missing_in_dir1+=("$rel")
done < <(cd "$dir2" && find . -type f \( -name "*.fna" -o -name "*.fna.gz" \) -printf '%P\0')

echo "Missing in $dir2 (present in $dir1):"
((${#missing_in_dir2[@]})) && printf '%s\n' "${missing_in_dir2[@]}" | sort || echo "(none)"

echo
echo "Missing in $dir1 (present in $dir2):"
((${#missing_in_dir1[@]})) && printf '%s\n' "${missing_in_dir1[@]}" | sort || echo "(none)"
