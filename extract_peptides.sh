#!/bin/bash

# This script runs the full pipeline consisting of two steps:
# Step 1: DNA sequence extraction based on HUMAnN2 output (array job)
# Step 2: Translation and BLAST validation of extracted sequences (single job)

# Example long usage:~~~~
# DATA_DIRECTORY="/home/armetcal/projects/rrg-bfinlay/armetcal/cmmi/pilot_g4h/mtx_extracted"
# ./extract_peptides.sh \
#    --working-dir $(pwd) \
#     --humann-out $DATA_DIRECTORY/humann_out/ur50/mtx \
#     --targets mtx_uniprot50_prev50_abun1_bug3_IDs.txt \
#     --fastq-dir $DATA_DIRECTORY/hostile/concatenated \
#     --output-dir extract_peptides_outputs \
#     --database /home/armetcal/scratch/humann_databases/uniref50/uniref50_201901b_full.dmnd \
#     --max-targets 5 \
#     --evalue 1e-3 \
#     --pident 90 \
#     --min-length 7 \
#     --sensitive \
#     --rerun

# Example short usage:~~~~
# DATA_DIRECTORY="/home/armetcal/projects/rrg-bfinlay/armetcal/cmmi/pilot_g4h/mtx_extracted"
# ./extract_peptides.sh \
#     -H $DATA_DIRECTORY/humann_out/ur50/mtx \
#     -t mtx_uniprot50_prev50_abun1_bug3_IDs.txt \
#     -f $DATA_DIRECTORY/hostile/concatenated \
#     -d /home/armetcal/scratch/humann_databases/uniref50/uniref50_201901b_full.dmnd \
#     -s \
#     -r
#~~~~~~~~~~~~~~~~~~~

# Usage function
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Runs the full pipeline: Step 1 (array job) → Step 2 (single job)"
    echo "Number of samples for Step 1 is automatically determined from HUMAnN2 output directory"
    echo ""
    echo "GLOBAL & STEP 1 OPTIONS:"
    echo "  -w, --working-dir DIR      Working directory (default: current directory)"
    echo "  -H, --humann-out DIR      HUMAnN2 output directory, where each sample has a temp folder (required)"
    echo "  -t, --targets FILE         Line-separated list of UniRef50 ID targets (required)"
    echo "  -f, --fastq-dir DIR        Directory with cleaned, unannotated FASTQ files (required)"
    echo "  -O, --output-dir DIR       Output directory (required)"
    echo ""
    echo "STEP 2 OPTIONS:"
    echo "  -d, --database FILE        Diamond database file (required)"
    echo "  -m, --max-targets INT      Max targets per query (default: 5)"
    echo "  -e, --evalue FLOAT         E-value threshold (default: 1e-3)"
    echo "  -p, --pident INT           Percent identity cutoff (default: 90)"
    echo "  -l, --min-length INT       Minimum amino acid length (default: 7)"
    echo "  -s, --sensitive            Use sensitive mode (default: false)"
    echo "  -r, --rerun                Rerun filtering only using existing BLAST output"
    echo ""
    echo "GENERAL OPTIONS:"
    echo "  -h, --help                 Show this help message"
    echo ""
    echo "EXAMPLE:"
    echo "  $0 --working-dir /scratch/user/project --humann-out /path/to/humann_out --targets targets.txt --fastq-dir /path/to/fastq --output-dir /path/to/output --database /path/to/database.dmnd --pident 95 --min-length 10 --sensitive"
}

# Wrapper function to run the full pipeline
run_full_pipeline() {
    local WORKING_DIR="$1"
    local HUMANN_OUT="$2"
    local TARGETS="$3"
    local FASTQ_DIR="$4"
    local STEP1_OUTPUT="$5"
    local STEP2_OUTPUT="$6"
    local ARRAY_SIZE="$7"
    local DATABASE="$8"
    local MAX_TARGETS="$9"
    local EVALUE="${10}"
    local PIDENT_CUTOFF="${11}"
    local MIN_LENGTH="${12}"
    local SENSITIVE_MODE="${13}"
    local RERUN_MODE="${14}"
    
    local STEP1_SCRIPT="./step1_function.sh"
    local STEP2_SCRIPT="./step2_function.sh"
    
    echo "=== Starting Full Pipeline ==="
    date
    echo "Step 1 Parameters:"
    echo "  Working dir: $WORKING_DIR"
    echo "  HUMAnN output folder: $HUMANN_OUT"
    echo "  Targets: $TARGETS"
    echo "  FASTQ dir: $FASTQ_DIR"
    echo "  Output dir: $STEP1_OUTPUT"
    echo "  Array size: $ARRAY_SIZE"
    echo "Step 2 Parameters:"
    echo "  Working dir: $WORKING_DIR"
    echo "  Output dir: $STEP2_OUTPUT"
    echo "  Database: $DATABASE"
    echo "  Max targets: $MAX_TARGETS"
    echo "  E-value: $EVALUE"
    echo "  Pident: $PIDENT_CUTOFF"
    echo "  Min length: $MIN_LENGTH"
    echo "  Sensitive: $SENSITIVE_MODE"
    echo "  Rerun: $RERUN_MODE"
    
    # Validate scripts
    if [[ ! -x "$STEP1_SCRIPT" ]]; then
        chmod +x "$STEP1_SCRIPT"
    fi
    if [[ ! -x "$STEP2_SCRIPT" ]]; then
        chmod +x "$STEP2_SCRIPT"
    fi

    # Prep output folders
    mkdir -p "$STEP1_OUTPUT/logs"
    mkdir -p "$STEP2_OUTPUT/logs"

    # Step 1: Submit array job
    echo "=== Step 1: DNA Sequence Extraction ==="
    echo "Submitting Step 1 array job ($ARRAY_SIZE samples)..."
    
    local step1_cmd=("$STEP1_SCRIPT"
        --working-dir "$WORKING_DIR"
        --humann-out "$HUMANN_OUT"
        --targets "$TARGETS"
        --output-dir "$STEP1_OUTPUT"
        --fastq-dir "$FASTQ_DIR"
    )
    
    # Usually runs in ~3 mins, should be safe with 30 mins
    # Use a temporary script file - necessary to properly handle array job submission with multiple parameters
    local temp_script=$(mktemp)
    cat > "$temp_script" << EOF
#!/bin/bash
#SBATCH --job-name=step1_extract_relevant_reads
#SBATCH --array=1-$ARRAY_SIZE
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --output=$STEP1_OUTPUT/logs/step1_%A_%a.out

source /home/armetcal/miniconda3/etc/profile.d/conda.sh
conda activate humann
module load seqkit/2.5.1 StdEnv/2023

exec ${step1_cmd[@]}
EOF

    local step1_job_id=$(sbatch --parsable "$temp_script")
    rm "$temp_script"

    echo "Step 1 job ID: $step1_job_id"


    # Wait for Step 1 completion
    echo "Waiting for Step 1 to complete..."
    local wait_count=0
    while squeue -j "$step1_job_id" 2>/dev/null | grep -q "$step1_job_id"; do
        sleep 60
        ((wait_count++))
        if ((wait_count % 5 == 0)); then
            echo "Step 1 still running... ($wait_count minutes elapsed)"
        fi
    done
    echo "Step 1 completed successfully"
    
    # Step 2: Submit single job
    echo "=== Step 2: Translation and BLAST Validation ==="
    echo "Submitting Step 2..."
    
    local step2_cmd=("$STEP2_SCRIPT"
        --working-dir "$WORKING_DIR"
        --database "$DATABASE"
        --targets "$TARGETS"
        --input-dir "$STEP1_OUTPUT"
        --output-dir "$STEP2_OUTPUT"
        --max-targets "$MAX_TARGETS"
        --evalue "$EVALUE"
        --pident "$PIDENT_CUTOFF"
        --min-length "$MIN_LENGTH"
    )
    
    # Add boolean flags if enabled
    [[ "$SENSITIVE_MODE" == "true" ]] && step2_cmd+=(--sensitive)
    [[ "$RERUN_MODE" == "true" ]] && step2_cmd+=(--rerun)
    
    local step2_job_id=$(sbatch --parsable \
        --job-name=step2_translate_and_blast \
        --cpus-per-task=10 \
        --mem=64G \
        --time=1:00:00 \
        --output="$STEP2_OUTPUT/logs/step2.out" \
        --dependency=afterok:$step1_job_id \
        "${step2_cmd[@]}")
    
    echo "Step 2 job ID: $step2_job_id"
    echo "=== Step 2 Submitted to SLURM ==="
    echo "NOTE: pipeline is NOT complete until Step 2 finishes!!"
    echo "Step 2 will run automatically after Step 1 finishes"
    echo "Monitor: squeue -j $step2_job_id"
    echo "Step 1 results: $STEP1_OUTPUT"
    echo "Step 2 results: $STEP2_OUTPUT"
    
    return 0
}

# Parse command-line arguments with defaults
WORKING_DIR=$(pwd)
HUMANN_OUT=""
TARGETS=""
FASTQ_DIR=""
OUTPUT_DIR="$(pwd)/extract_peptides_outputs"
DATABASE=""
MAX_TARGETS=5
EVALUE="1e-3"
PIDENT_CUTOFF=90
MIN_LENGTH=7
SENSITIVE_MODE="false"
RERUN_MODE="false"

while [[ $# -gt 0 ]]; do
    case $1 in
        # Step 1 options
        -w|--working-dir)
            WORKING_DIR="$2"
            shift 2
            ;;
        -H|--humann-out)
            HUMANN_OUT="$2"
            shift 2
            ;;
        -t|--targets)
            TARGETS="$2"
            shift 2
            ;;
        -f|--fastq-dir)
            FASTQ_DIR="$2"
            shift 2
            ;;
        -O|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        # Step 2 options
        -d|--database)
            DATABASE="$2"
            shift 2
            ;;
        -m|--max-targets)
            MAX_TARGETS="$2"
            shift 2
            ;;
        -e|--evalue)
            EVALUE="$2"
            shift 2
            ;;
        -p|--pident)
            PIDENT_CUTOFF="$2"
            shift 2
            ;;
        -l|--min-length)
            MIN_LENGTH="$2"
            shift 2
            ;;
        -s|--sensitive)
            SENSITIVE_MODE="true"
            shift 1
            ;;
        -r|--rerun)
            RERUN_MODE="true"
            shift 1
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Set output directories. Currently pools both steps into same output dir, feel free to modify as needed.
STEP1_OUTPUT=$OUTPUT_DIR
STEP2_OUTPUT=$OUTPUT_DIR

# Calculate array_size (1 per sample)
samples=("$HUMANN_OUT"/*_humann_temp)
if [[ ${#samples[@]} -eq 0 ]]; then
    echo "ERROR: No humann_temp directories found in $HUMANN_OUT"
    return 1
fi
ARRAY_SIZE=${#samples[@]}
echo "Automatically detected $ARRAY_SIZE samples"

# Data size warnings
if [[ $ARRAY_SIZE -gt 50 ]]; then
    echo "WARNING: Large number of samples ($ARRAY_SIZE). Consider adjusting time limits, especially for Step 2."
fi
TARGET_NUM=$(wc -l < "$TARGETS")
if [[ $TARGET_NUM -gt 5000 ]]; then
    echo "WARNING: Large number of target IDs ($TARGET_NUM). Consider increasing memory/time limits."
fi

# Validate required arguments
if [[ -z "$WORKING_DIR" || -z "$HUMANN_OUT" || -z "$TARGETS" || -z "$FASTQ_DIR" || -z "$DATABASE" ]]; then
    echo "ERROR: Missing required arguments"
    usage
    exit 1
fi

# Validate numeric parameters
if ! [[ "$ARRAY_SIZE" =~ ^[0-9]+$ ]] || [[ "$ARRAY_SIZE" -lt 1 ]]; then
    echo "ERROR: Array size must be positive integer: $ARRAY_SIZE"
    exit 1
fi

# Run the pipeline
run_full_pipeline "$WORKING_DIR" "$HUMANN_OUT" "$TARGETS" "$FASTQ_DIR" "$STEP1_OUTPUT" "$STEP2_OUTPUT" \
    "$ARRAY_SIZE" "$DATABASE" "$MAX_TARGETS" "$EVALUE" "$PIDENT_CUTOFF" "$MIN_LENGTH" "$SENSITIVE_MODE" "$RERUN_MODE"

EXIT_CODE=$?
if [[ $EXIT_CODE -eq 0 ]]; then
    echo "Full pipeline submitted successfully"
else
    echo "Pipeline submission failed"
fi

date
exit $EXIT_CODE
