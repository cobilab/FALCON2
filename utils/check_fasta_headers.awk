#!/usr/bin/awk -f
#
# check_fasta_headers.awk
#
# Purpose: Check multi-FASTA headers according to several rules:
#   1) "kraken:taxid:<digits>"
#   2) Entire ID is just digits
#   3) "accession-like" pattern (previously AB_1234.1, etc.)
#   4) NCBI typical headers containing NC_ or WP_ or "lcl|..."
#      (e.g. "lcl|NC_002162.1_cds_WP_006688755.1_1 ...")
#
# Print "OK:" if a header meets any rule, otherwise "ERROR:".

BEGIN {
  # 1) "kraken:taxid:\d+"
  krakenPat    = "kraken:taxid:[0-9]+"

  # 2) Entire ID is numeric (like ">12345")
  numericPat   = "^[0-9]+$"

  # 3) Simple "accession-like" pattern (letters/digits, maybe underscore, dot + version)
  accPat       = "^[A-Za-z0-9_]+\\.[0-9]+$"

  # 4) Typical NCBI headers:
  #    - Something with "NC_####" or "WP_####"
  #    - Or "lcl|..." at the start
  #    We make a broad pattern that checks for these substrings anywhere in the line
  #    (since many lines have "lcl|NC_002162.1_cds_WP_006688755.1_1 [gene=...]")
  ncbiPat      = "(NC_[0-9]+\\.[0-9]+)|(WP_[0-9]+\\.[0-9]+)|(lcl\\|)"
}

/^>/ {
  # Remove leading ">"
  header = substr($0, 2)
  # Trim leading/trailing spaces
  gsub(/^[ \t]+|[ \t]+$/, "", header)

  if (header ~ krakenPat) {
    print "OK:", $0
  }
  else if (header ~ numericPat) {
    print "OK:", $0
  }
  else if (header ~ accPat) {
    print "OK:", $0
  }
  else if (header ~ ncbiPat) {
    # Matches "lcl|" or "NC_xxx" or "WP_xxx"
    print "OK:", $0
  }
  else {
    print "ERROR: Unrecognized header =>", $0
    # If you want to stop on first error, uncomment:
    # exit 1
  }
  next
}

# Non-header lines are ignored
