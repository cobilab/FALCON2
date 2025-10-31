#!/bin/bash
# FALCON2 Comprehensive Test Suite
# Tests both synthetic and real data files
# Usage: ./test_falcon2.sh
#        NO_VALGRIND=1 ./test_falcon2.sh   # skip valgrind for speed
#        FALCON=./FALCON2 ./test_falcon2.sh # custom binary

# Configuration
FALCON="${FALCON:-./FALCON2}"
BASE_PARAMS="meta -v -F -t 15 -l 47"
BASE_PARAMS_NOF="${BASE_PARAMS/-F/}"  # without -F for permission tests
VALGRIND="valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes"
NO_VALGRIND="${NO_VALGRIND:-0}"

# Real test files (from provided test directory)
ROOT="${ROOT:-.}"
REAL_READS="$ROOT/reads.fq"
REAL_READS_GZ="$ROOT/reads.fq.gz"
REAL_DB1="$ROOT/VDB.fa"
REAL_DB1_GZ="$ROOT/VDB.fa.gz"
REAL_DB2="$ROOT/VDB2.fa"
REAL_DB2_GZ="$ROOT/VDB2.fa.gz"
REAL_FILTER="$ROOT/input_fasta.fa"
REAL_FILTER_GZ="$ROOT/input_fasta.fa.gz"
REAL_MODEL="$ROOT/falcon_model.fcm"

# Setup logs (persisted locally)
RUN_ID="$(date +'%Y%m%d_%H%M%S')"
LOG_DIR="./falcon_logs/$RUN_ID"
mkdir -p "$LOG_DIR"

WORK="$(mktemp -d)"
trap "chmod -f 644 $WORK/* 2>/dev/null; rm -rf $WORK" EXIT

echo "========================================="
echo "  FALCON2 Comprehensive Test Suite"
echo "========================================="
echo ""

# Check dependencies
if [ ! -x "$FALCON" ]; then
    echo "Error: FALCON2 not executable at $FALCON"
    exit 1
fi

if [[ "$NO_VALGRIND" != "1" ]]; then
    if ! command -v valgrind &> /dev/null; then
        echo "Error: valgrind not found"
        exit 1
    fi
    VG="$VALGRIND"
    echo "Valgrind: ENABLED"
else
    VG=""
    echo "Valgrind: DISABLED"
fi

GZIP_OK=0
if command -v gzip >/dev/null 2>&1; then
  GZIP_OK=1
fi

echo "Binary:   $FALCON"
echo "Version:  $($FALCON -V 2>&1 | head -n9 || echo 'unknown')"
echo "Logs:     $LOG_DIR"
echo "Root:     $ROOT"
echo ""

# Test counters
PASS=0
FAIL=0
SKIP=0
declare -a FAILURES=()

# Test runners
run_with_params() {
    local param_str="$1"; shift
    local name="$1"; local expect="$2"; shift 2
    local log="$LOG_DIR/${name// /_}.log"
    echo -n "  $name ... "
    
    # shellcheck disable=SC2086
    $VG $FALCON $param_str "$@" &> "$log"
    local rc=$?
    
    if [[ $rc -eq $expect ]]; then
        echo "PASS"
        ((PASS++))
    else
        echo "FAIL (exit=$rc, expected=$expect)"
        ((FAIL++))
        FAILURES+=("$name")
    fi
}

run_test()      { run_with_params "$BASE_PARAMS"      "$@"; }
run_test_nof()  { run_with_params "$BASE_PARAMS_NOF"  "$@"; }

skip_test() {
    echo "  $1 ... SKIP"
    ((SKIP++))
}

# Test data creators
mk_fq() {
    local n="${2:-2}"
    : > "$1"
    for i in $(seq 1 "$n"); do
        cat >> "$1" <<EOF
@read_$i
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
    done
}

mk_fa() {
    local n="${2:-2}"
    : > "$1"
    for i in $(seq 1 "$n"); do
        cat >> "$1" <<EOF
>seq_$i
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOF
    done
}

# Generate synthetic test data
echo "Generating synthetic test data..."
SMALL_FQ="$WORK/small.fq"; mk_fq "$SMALL_FQ" 5
SMALL_FA="$WORK/small.fa"; mk_fa "$SMALL_FA" 3
TINY_FQ="$WORK/tiny.fq"; mk_fq "$TINY_FQ" 1
TINY_FA="$WORK/tiny.fa"; mk_fa "$TINY_FA" 1
MEDIUM_FQ="$WORK/medium.fq"; mk_fq "$MEDIUM_FQ" 50
MEDIUM_FA="$WORK/medium.fa"; mk_fa "$MEDIUM_FA" 20
EMPTY_FQ="$WORK/empty.fq"; : > "$EMPTY_FQ"
EMPTY_FA="$WORK/empty.fa"; : > "$EMPTY_FA"
BAD_FQ="$WORK/bad.fq"; echo "not a valid fastq file" > "$BAD_FQ"
BAD_FA="$WORK/bad.fa"; echo "not a valid fasta file" > "$BAD_FA"
UNREADABLE_FQ="$WORK/unreadable.fq"; mk_fq "$UNREADABLE_FQ"; chmod 000 "$UNREADABLE_FQ" || true
MISSING_FQ="$WORK/missing.fq"
MISSING_FA="$WORK/missing.fa"
MISSING_MODEL="$WORK/missing.fcm"
MODEL_OUT="$WORK/trained_model.fcm"
MODEL_OUT2="$WORK/trained_model2.fcm"
OUT1="$WORK/output1.csv"
OUT2="$WORK/output2.csv"
OUT3="$WORK/output3.csv"
READONLY_OUT="$WORK/readonly.csv"; : > "$READONLY_OUT"; chmod 444 "$READONLY_OUT" || true

if [[ $GZIP_OK -eq 1 ]]; then
  SMALL_FQ_GZ="$WORK/small.fq.gz"; gzip -c "$SMALL_FQ" > "$SMALL_FQ_GZ"
  SMALL_FA_GZ="$WORK/small.fa.gz"; gzip -c "$SMALL_FA" > "$SMALL_FA_GZ"
  MEDIUM_FQ_GZ="$WORK/medium.fq.gz"; gzip -c "$MEDIUM_FQ" > "$MEDIUM_FQ_GZ"
fi
echo ""

# =========================
# Section 0: Early Exits
# =========================
echo "[0] Early Exit Tests"
run_test "help (-h)"        0 -h
run_test "version (-V)"     0 -V
run_test "show levels (-s)" 0 -s
echo ""

# =========================
# Section 1: Training Mode
# =========================
echo "[1] Training Mode (-T) - Synthetic Data"
run_test "train tiny reads"        0 -T "$TINY_FQ"
run_test "train small reads"       0 -T "$SMALL_FQ"
run_test "train medium reads"      0 -T "$MEDIUM_FQ"
if [[ $GZIP_OK -eq 1 ]]; then
  run_test "train small.fq.gz"     0 -T "$SMALL_FQ_GZ"
  run_test "train medium.fq.gz"    0 -T "$MEDIUM_FQ_GZ"
fi
run_test "train empty file"        0 -T "$EMPTY_FQ"
run_test "train bad format"        0 -T "$BAD_FQ"
run_test "train missing file"      1 -T "$MISSING_FQ"
run_test "train unreadable file"   1 -T "$UNREADABLE_FQ"
run_test "train no arguments"      1 -T
echo ""

# With real files if available
if [[ -f "$REAL_READS" ]]; then
    echo "[1b] Training Mode - Real Data"
    run_test "train real reads.fq"     0 -T "$REAL_READS"
    if [[ -f "$REAL_READS_GZ" ]]; then
        run_test "train real reads.fq.gz" 0 -T "$REAL_READS_GZ"
    fi
    echo ""
fi

# =========================
# Section 2: Model Operations
# =========================
echo "[2] Model Save/Load/Info"
run_test "save model (-S -M)"          0 -T -S -M "$MODEL_OUT" "$SMALL_FQ"
run_test "save model different path"   0 -T -S -M "$MODEL_OUT2" "$MEDIUM_FQ"

# Verify model was created
if [[ -f "$MODEL_OUT" ]]; then
    echo "  Model created: $MODEL_OUT ($(stat -f%z "$MODEL_OUT" 2>/dev/null || stat -c%s "$MODEL_OUT" 2>/dev/null || echo '?') bytes)"
    run_test "load saved model"        0 -T -L -M "$MODEL_OUT" "$SMALL_FQ"
    run_test "model info (-I)"         0 -T -I -L -M "$MODEL_OUT" "$SMALL_FQ"
fi

# Load real model if available
if [[ -f "$REAL_MODEL" ]]; then
    run_test "load real model"         0 -T -L -M "$REAL_MODEL" "$SMALL_FQ"
    run_test "real model info"         0 -T -I -L -M "$REAL_MODEL" "$SMALL_FQ"
fi

run_test "load without -M flag"        0 -T -L "$SMALL_FQ"
run_test "load missing model"          1 -T -L -M "$MISSING_MODEL" "$SMALL_FQ"
run_test "save without -M path"        0 -T -S "$SMALL_FQ"
echo ""

# =========================
# Section 3: Inference Mode
# =========================
echo "[3] Inference Mode - Synthetic Data"
run_test "infer tiny files"            0 "$TINY_FQ" "$TINY_FA"
run_test "infer small files"           0 "$SMALL_FQ" "$SMALL_FA"
run_test "infer medium files"          0 "$MEDIUM_FQ" "$MEDIUM_FA"
run_test "infer small reads + medium db" 0 "$SMALL_FQ" "$MEDIUM_FA"
if [[ $GZIP_OK -eq 1 ]]; then
  run_test "infer small.fq.gz + small.fa.gz" 0 "$SMALL_FQ_GZ" "$SMALL_FA_GZ"
  run_test "infer mixed: fq + fa.gz"   0 "$SMALL_FQ" "$SMALL_FA_GZ"
  run_test "infer mixed: fq.gz + fa"   0 "$SMALL_FQ_GZ" "$SMALL_FA"
fi

# Error cases
run_test "infer only 1 arg (reads)"    1 "$SMALL_FQ"
run_test "infer no arguments"          1
run_test "infer missing reads"         1 "$MISSING_FQ" "$SMALL_FA"
run_test "infer missing db"            1 "$SMALL_FQ" "$MISSING_FA"
run_test "infer unreadable reads"      1 "$UNREADABLE_FQ" "$SMALL_FA"
run_test "infer invalid reads format"  0 "$BAD_FQ" "$SMALL_FA"
run_test "infer invalid db format"     0 "$SMALL_FQ" "$BAD_FA"
run_test "infer empty reads"           0 "$EMPTY_FQ" "$SMALL_FA"
run_test "infer empty db"              0 "$SMALL_FQ" "$EMPTY_FA"
run_test "infer both empty"            0 "$EMPTY_FQ" "$EMPTY_FA"
echo ""

# Real data inference
if [[ -f "$REAL_READS" && -f "$REAL_DB1" ]]; then
    echo "[3b] Inference Mode - Real Data"
    run_test "infer reads.fq + VDB.fa"     0 "$REAL_READS" "$REAL_DB1"
    if [[ -f "$REAL_DB2" ]]; then
        run_test "infer reads.fq + VDB2.fa" 0 "$REAL_READS" "$REAL_DB2"
    fi
    if [[ -f "$REAL_READS_GZ" && -f "$REAL_DB1_GZ" ]]; then
        run_test "infer reads.fq.gz + VDB.fa.gz" 0 "$REAL_READS_GZ" "$REAL_DB1_GZ"
    fi
    if [[ -f "$REAL_DB2_GZ" ]]; then
        run_test "infer reads.fq + VDB2.fa.gz" 0 "$REAL_READS" "$REAL_DB2_GZ"
    fi
    echo ""
fi

# =========================
# Section 4: Magnet Integration
# =========================
echo "[4] Magnet Filter Tests"
run_test "magnet basic filter"         0 -mg -mf "$REAL_FILTER" "$REAL_READS" "$REAL_DB1"
run_test "magnet verbose (-mv)"        0 -mg -mv -mf "$REAL_FILTER" "$REAL_READS" "$REAL_DB1"
run_test "magnet with threshold"       0 -mg -mf "$REAL_FILTER" -mt 0.85 "$REAL_READS" "$REAL_DB1"
run_test "magnet with min length"      0 -mg -mf "$REAL_FILTER" -ml 40 "$REAL_READS" "$REAL_DB1"
run_test "magnet invert (-mi)"         0 -mg -mf "$REAL_FILTER" -mi "$REAL_READS" "$REAL_DB1"
run_test "magnet max positions"        0 -mg -mf "$REAL_FILTER" -mp 2 "$REAL_READS" "$REAL_DB1"
run_test "magnet all params"           0 -mg -mv -mf "$REAL_FILTER" -mt 0.85 -ml 40 -mi -mp 2 "$REAL_READS" "$REAL_DB1"
run_test "magnet gz filter"          0 -mg -mf "$REAL_FILTER_GZ" "$REAL_READS" "$REAL_DB1"

# Error cases
run_test "magnet missing filter"       1 -mg -mf "$MISSING_FA" "$REAL_READS" "$REAL_DB1"
run_test "magnet no -mf flag"          1 -mg "$REAL_READS" "$REAL_DB1"
run_test "magnet empty filter"         0 -mg -mf "$EMPTY_FA" "$REAL_READS" "$REAL_DB1"
echo ""

# =========================
# Section 5: Output Options
# =========================
echo "[5] Output File Tests"
run_test "output custom path (-x)"     0 -x "$OUT1" "$SMALL_FQ" "$SMALL_FA"

# Verify output was created
if [[ -f "$OUT1" ]]; then
    lines=$(wc -l < "$OUT1")
    echo "  Output verified: $OUT1 ($lines lines)"
fi

run_test "output different path"       0 -x "$OUT2" "$MEDIUM_FQ" "$MEDIUM_FA"
run_test "output overwrite (with -F)"  0 -x "$OUT1" "$SMALL_FQ" "$SMALL_FA"

# Permission test without -F
run_test_nof "output readonly no -F"   1 -x "$READONLY_OUT" "$SMALL_FQ" "$SMALL_FA"
chmod 644 "$READONLY_OUT" || true

# Output with different operations
run_test "output with training"        0 -T -x "$OUT3" "$SMALL_FQ"
echo ""

# =========================
# Section 6: Parameter Validation
# =========================
echo "[6] Thread Count (-n)"
run_test "threads n=1"                 0 -n 1 "$SMALL_FQ" "$SMALL_FA"
run_test "threads n=2"                 0 -n 2 "$SMALL_FQ" "$SMALL_FA"
run_test "threads n=8"                 0 -n 8 "$SMALL_FQ" "$SMALL_FA"
run_test "threads n=16"                0 -n 16 "$SMALL_FQ" "$SMALL_FA"
run_test "threads n=32"                0 -n 32 "$SMALL_FQ" "$SMALL_FA"
run_test "threads n=0 invalid"         1 -n 0 "$SMALL_FQ" "$SMALL_FA"
run_test "threads n=-1 invalid"        1 -n -1 "$SMALL_FQ" "$SMALL_FA"
echo ""

echo "[7] Top Results (-t)"
run_test "top t=1"                     0 -t 1 "$SMALL_FQ" "$SMALL_FA"
run_test "top t=5"                     0 -t 5 "$SMALL_FQ" "$SMALL_FA"
run_test "top t=10"                    0 -t 10 "$SMALL_FQ" "$SMALL_FA"
run_test "top t=50"                    0 -t 50 "$SMALL_FQ" "$SMALL_FA"
run_test "top t=100"                   0 -t 100 "$SMALL_FQ" "$SMALL_FA"
echo ""

echo "[8] Compression Level (-l)"
run_test "level l=20"                  0 -l 20 "$SMALL_FQ" "$SMALL_FA"
run_test "level l=30"                  0 -l 30 "$SMALL_FQ" "$SMALL_FA"
run_test "level l=36"                  0 -l 36 "$SMALL_FQ" "$SMALL_FA"
run_test "level l=40"                  0 -l 40 "$SMALL_FQ" "$SMALL_FA"
run_test "level l=47"                  0 -l 47 "$SMALL_FQ" "$SMALL_FA"
echo ""

echo "[9] Gamma Parameter (-g)"
run_test "gamma g=0.000015"            0 -g 0.000015 "$SMALL_FQ" "$SMALL_FA"
run_test "gamma g=0.0001"              0 -g 0.0001 "$SMALL_FQ" "$SMALL_FA"
run_test "gamma g=0.001"               0 -g 0.001 "$SMALL_FQ" "$SMALL_FA"
run_test "gamma g=0.01"                0 -g 0.01 "$SMALL_FQ" "$SMALL_FA"
run_test "gamma g=0.1"                 0 -g 0.1 "$SMALL_FQ" "$SMALL_FA"
run_test "gamma g=0.5"                 0 -g 0.5 "$SMALL_FQ" "$SMALL_FA"
run_test "gamma g=0.9"                 0 -g 0.9 "$SMALL_FQ" "$SMALL_FA"
echo ""

echo "[10] Collision Parameter (-c)"
run_test "collision c=1"               0 -c 1 "$SMALL_FQ" "$SMALL_FA"
run_test "collision c=10"              0 -c 10 "$SMALL_FQ" "$SMALL_FA"
run_test "collision c=50"              0 -c 50 "$SMALL_FQ" "$SMALL_FA"
run_test "collision c=100"             0 -c 100 "$SMALL_FQ" "$SMALL_FA"
echo ""

# =========================
# Section 11: Combined Parameters
# =========================
echo "[15] Parameter Combinations"
run_test "combo: threads + top"        0 -n 8 -t 10 "$SMALL_FQ" "$SMALL_FA"
run_test "combo: level + gamma"        0 -l 40 -g 0.01 "$SMALL_FQ" "$SMALL_FA"
run_test "combo: all inference params" 0 -n 4 -t 15 -l 47 -g 0.000015 -c 100 "$SMALL_FQ" "$SMALL_FA"
run_test "combo: training + save"      0 -T -S -M "$WORK/combo.fcm" "$SMALL_FQ"
run_test "combo: magnet + custom out"  0 -mg -mf "$SMALL_FA" -x "$WORK/magnet_out.csv" "$SMALL_FQ" "$SMALL_FA"
run_test "combo: local + profile"      0 -Z -y "$WORK/profile.txt" "$SMALL_FQ" "$SMALL_FA"
run_test "combo: sample + threads"     0 -p 10 -n 4 "$SMALL_FQ" "$SMALL_FA"
run_test "combo: multi-file + magnet"  0 -mg -mf "$SMALL_FA" "$SMALL_FQ:$MEDIUM_FQ" "$SMALL_FA:$MEDIUM_FA"
echo ""

# =========================
# Section 16: Advanced Parameters
# =========================
echo "[16] Advanced Parameters"
# Local similarity (-Z) with profile (-y)
run_test "local similarity (-Z)"       0 -Z "$SMALL_FQ" "$SMALL_FA"
run_test "profile output (-Z -y)"      0 -Z -y "$WORK/profile1.txt" "$SMALL_FQ" "$SMALL_FA"
if [[ -f "$WORK/profile1.txt" ]]; then
    echo "  Profile created: $WORK/profile1.txt ($(wc -l < "$WORK/profile1.txt") lines)"
fi

# Subsampling (-p)
run_test "subsample p=1 (all)"         0 -p 1 "$SMALL_FQ" "$SMALL_FA"
run_test "subsample p=2"               0 -p 2 "$SMALL_FQ" "$SMALL_FA"
run_test "subsample p=5"               0 -p 5 "$SMALL_FQ" "$SMALL_FA"
run_test "subsample p=10"              0 -p 10 "$SMALL_FQ" "$SMALL_FA"
run_test "subsample p=100"             0 -p 100 "$SMALL_FQ" "$SMALL_FA"

# =========================
# Section 17: Error Detection
# =========================
echo "[17] Invalid Flags & Arguments"
run_test "unknown flag --invalid"      0 --this-flag-does-not-exist "$SMALL_FQ" "$SMALL_FA"
run_test "typo flag -xyz"              0 -xyz123 "$SMALL_FQ" "$SMALL_FA"
run_test "double dash alone --"        0 -- "$SMALL_FQ" "$SMALL_FA"
run_test "flag without value -n"       1 -n "$SMALL_FQ" "$SMALL_FA"
echo ""

# =========================
# Section 9: Multi-File Input (colon-separated)
# =========================
echo "[13] Multi-File Input Tests"
# Multiple reads files
run_test "multi reads: tiny:small"     0 "$TINY_FQ:$SMALL_FQ" "$SMALL_FA"
run_test "multi reads: small:medium"   0 "$SMALL_FQ:$MEDIUM_FQ" "$SMALL_FA"
run_test "multi reads: 3 files"        0 "$TINY_FQ:$SMALL_FQ:$MEDIUM_FQ" "$SMALL_FA"

# Multiple database files
run_test "multi db: tiny:small"        0 "$SMALL_FQ" "$TINY_FA:$SMALL_FA"
run_test "multi db: small:medium"      0 "$SMALL_FQ" "$SMALL_FA:$MEDIUM_FA"
run_test "multi db: 3 files"           0 "$SMALL_FQ" "$TINY_FA:$SMALL_FA:$MEDIUM_FA"

# Both multi-file
run_test "multi both: 2:2"             0 "$TINY_FQ:$SMALL_FQ" "$TINY_FA:$SMALL_FA"
run_test "multi both: 3:2"             0 "$TINY_FQ:$SMALL_FQ:$MEDIUM_FQ" "$SMALL_FA:$MEDIUM_FA"

if [[ $GZIP_OK -eq 1 ]]; then
  # Mixed plain and gzipped in multi-file
  run_test "multi mixed: fq:fq.gz"     0 "$SMALL_FQ:$SMALL_FQ_GZ" "$SMALL_FA"
  run_test "multi mixed db: fa:fa.gz"  0 "$SMALL_FQ" "$SMALL_FA:$SMALL_FA_GZ"
  run_test "multi all mixed"           0 "$SMALL_FQ:$SMALL_FQ_GZ" "$SMALL_FA:$SMALL_FA_GZ"
fi

# Real files multi-file combinations
if [[ -f "$REAL_DB1" && -f "$REAL_DB2" ]]; then
    run_test "multi real: VDB:VDB2"    0 "$SMALL_FQ" "$REAL_DB1:$REAL_DB2"
    if [[ -f "$REAL_DB1_GZ" ]]; then
        run_test "multi real: VDB:VDB.gz" 0 "$SMALL_FQ" "$REAL_DB1:$REAL_DB1_GZ"
    fi
    if [[ -f "$REAL_DB2_GZ" ]]; then
        run_test "multi real mixed: VDB.fa:VDB2.fa.gz" 0 "$SMALL_FQ" "$REAL_DB1:$REAL_DB2_GZ"
    fi
fi

# Error cases with multi-file
run_test "multi missing first"         1 "$MISSING_FQ:$SMALL_FQ" "$SMALL_FA"
run_test "multi missing second"        1 "$SMALL_FQ:$MISSING_FQ" "$SMALL_FA"
run_test "multi missing db first"      1 "$SMALL_FQ" "$MISSING_FA:$SMALL_FA"
run_test "multi empty in list"         1 "$SMALL_FQ:" "$SMALL_FA"
run_test "multi only colons"           1 ":::" "$SMALL_FA"
echo ""

# =========================
# Section 10: Edge Cases
# =========================
echo "[14] Edge Cases"
run_test "very small reads (1)"        0 "$TINY_FQ" "$SMALL_FA"
run_test "very small db (1)"           0 "$SMALL_FQ" "$TINY_FA"
run_test "both very small"             0 "$TINY_FQ" "$TINY_FA"
run_test "thread count > cores"        0 -n 128 "$SMALL_FQ" "$SMALL_FA"
run_test "very high top value"         0 -t 1000 "$SMALL_FQ" "$SMALL_FA"
echo ""

# =========================
# Summary
# =========================
echo "========================================="
echo "  Test Summary"
echo "========================================="
echo "Total:  $((PASS + FAIL + SKIP))"
echo "Pass:   $PASS"
echo "Fail:   $FAIL"
echo "Skip:   $SKIP"
echo "Logs:   $LOG_DIR"
echo ""
echo "Test Coverage:"
echo "  [0]  Early Exit Tests (help, version, show)"
echo "  [1]  Training Mode - Synthetic"
echo "  [1b] Training Mode - Real Data"
echo "  [2]  Model Operations (save, load, info)"
echo "  [3]  Inference Mode - Synthetic"
echo "  [3b] Inference Mode - Real Data"
echo "  [4]  Magnet Filter Integration"
echo "  [5]  Output File Handling"
echo "  [6]  Thread Count Parameter"
echo "  [7]  Top Results Parameter"
echo "  [8]  Compression Level Parameter"
echo "  [9]  Gamma Parameter"
echo "  [10] Collision Parameter"
echo "  [13] Multi-File Input (colon-separated)"
echo "  [14] Edge Cases"
echo "  [15] Parameter Combinations"
echo "  [16] Advanced Parameters (-Z, -y, -p)"
echo "  [17] Invalid Flags & Arguments"
echo ""

if (( FAIL > 0 )); then
    echo "Failed tests:"
    for test in "${FAILURES[@]}"; do
        echo "  - $test"
    done
    echo ""
    echo "Check logs in: $LOG_DIR"
    echo "========================================="
    exit 1
else
    echo "All tests passed!"
    echo "========================================="
    exit 0
fi
