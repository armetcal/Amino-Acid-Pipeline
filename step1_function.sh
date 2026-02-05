#!/bin/bash

# This script extracts DNA sequences from HUMAnN2 alignment results by identifying reads 
# that were confidently assigned to target UniRef50 IDs, then retrieves the original DNA 
# sequences from the input FASTQ files for subsequent six-frame translation and validation.

# Example usage:~~~~
#./step1_function.sh \
#  --working-dir /scratch/armetcal/cmmi/pilot_g4h/extract_relevant_aa_seqs \
#  --humann-root /home/armetcal/projects/rrg-bfinlay/armetcal/cmmi/pilot_g4h/mtx_extracted/humann_out/ur50/mtx \
#  --targets mtx_uniprot50_prev50_abun1_bug3_IDs.txt \
#  --output-dir temp_pipeline/step1 \
#  --fastq-dir /home/armetcal/projects/rrg-bfinlay/armetcal/cmmi/pilot_g4h/mtx_extracted/hostile/concatenated \
#  --sample-index 1
#~~~~~~~~~~~~~~~~~~~

# Load environment
source /home/armetcal/miniconda3/etc/profile.d/conda.sh
conda activate humann
module load seqkit/2.5.1 StdEnv/2023

# Usage function
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Extracts DNA sequences from HUMAnN2 results for target UniRef50 IDs"
    echo ""
    echo "OPTIONS:"
    echo "  -w, --working-dir DIR      Working directory (required)"
    echo "  -r, --humann-root DIR      HUMAnN2 output root directory (required)"
    echo "  -t, --targets FILE         Targets file with UniRef50 IDs (required)"
    echo "  -o, --output-dir DIR       Output directory (required)"
    echo "  -f, --fastq-dir DIR        Directory with original FASTQ files (required)"
    echo "  -i, --sample-index NUM     Sample index (1-based, defaults to SLURM_ARRAY_TASK_ID)"
    echo "  -h, --help                 Show this help message"
    echo ""
    echo "EXAMPLE:"
    echo "  $0 --working-dir /scratch/user/project --humann-root /path/to/humann_out --targets targets.txt --output-dir /path/to/output --fastq-dir /path/to/fastq_files --sample-index 1"
}

# Function to extract DNA sequences from HUMAnN2 results
extract_target_dna_sequences() {
    local WORKING_DIRECTORY="$1"
    local ROOT="$2"
    local TARGETS="$3"
    local OUTROOT="$4"
    local FASTQ_DIR="$5"
    local SAMPLE_INDEX="$6"
    
    # Record start time
    local START_TIME=$(date)
    
    # Create output directory
    mkdir -p "$OUTROOT"
    
    # Get sample directory
    local samples=("$ROOT"/*_humann_temp)
    if [[ ${#samples[@]} -eq 0 ]]; then
        echo "ERROR: No humann_temp directories found in $ROOT"
        return 1
    fi
    
    if [[ $SAMPLE_INDEX -lt 1 || $SAMPLE_INDEX -gt ${#samples[@]} ]]; then
        echo "ERROR: SAMPLE_INDEX $SAMPLE_INDEX out of range (1-${#samples[@]})"
        return 1
    fi
    
    local tempdir="${samples[$SAMPLE_INDEX-1]}"
    local sample_prefix=$(basename "$tempdir" | sed 's/_humann_temp$//')
    
    # Sample-specific output directory
    local OUTDIR="$OUTROOT/${sample_prefix}_dna_seqs"
    mkdir -p "$OUTDIR"
    
    echo "=== Step 1: DNA Sequence Extraction ==="
    echo "Start time: $START_TIME"
    echo "Sample: $sample_prefix"
    echo "Working directory: $WORKING_DIRECTORY"
    echo "HUMAnN2 root: $ROOT"
    echo "Targets file: $TARGETS"
    echo "Output directory: $OUTROOT"
    echo "FASTQ directory: $FASTQ_DIR"
    echo "Sample index: $SAMPLE_INDEX"
    
    # Validate input files
    if [[ ! -f "$TARGETS" ]]; then
        echo "ERROR: Targets file not found: $TARGETS"
        return 1
    fi
    
    if [[ ! -d "$FASTQ_DIR" ]]; then
        echo "ERROR: FASTQ directory not found: $FASTQ_DIR"
        return 1
    fi
    
    # 1) Build target lookup
    declare -A TARGET_SET
    while IFS= read -r target; do
        target_id="${target%%|*}"
        TARGET_SET["$target_id"]=1
    done < "$TARGETS"
    
    echo "Processing $sample_prefix with ${#TARGET_SET[@]} target IDs"
    
    # 2) Extract reads assigned to targets from TSV
    local DIAMOND_TSV="$tempdir/${sample_prefix}_diamond_aligned.tsv"
    local read_count=0
    local extracted_count=0
    
    if [[ -s "$DIAMOND_TSV" ]]; then
        awk '
            BEGIN {
                while ((getline < "'"$TARGETS"'") > 0) {
                    target = $0
                    sub(/\|.*/, "", target)
                    target_set[target] = 1
                }
            }
            {
                uniref_id = $2
                sub(/\|.*/, "", uniref_id)
                if (uniref_id in target_set) {
                    print $1
                }
            }
        ' "$DIAMOND_TSV" > "$OUTDIR/target_read_ids.txt"
        
        read_count=$(wc -l < "$OUTDIR/target_read_ids.txt")
        echo "Found $read_count reads assigned to target IDs"
        
        if [[ $read_count -eq 0 ]]; then
            echo "No target-assigned reads for $sample_prefix"
            rm -rf "$OUTDIR"
            # Still create flag file for completeness
            local FLAG_FILE="$OUTROOT/${sample_prefix}_step1_completed.flag"
            cat > "$FLAG_FILE" << EOF
STEP1_COMPLETED: $(date)
SAMPLE: $sample_prefix
STATUS: NO_TARGET_READS
TARGET_IDS_PROCESSED: ${#TARGET_SET[@]}
READS_ASSIGNED: 0
SEQUENCES_EXTRACTED: 0
EOF
            return 0
        fi
    else
        echo "No DIAMOND TSV found for $sample_prefix"
        return 0
    fi
    
    # 3) Extract original reads from FASTQ (pattern matching)
    local FASTQ_FILE="$FASTQ_DIR/${sample_prefix}.fastq.gz"
    if [[ -s "$FASTQ_FILE" ]]; then
        # Remove HUMAnN2's |151 suffix to match original FASTQ headers
        sed 's/|.*$//' "$OUTDIR/target_read_ids.txt" > "$OUTDIR/target_read_ids_clean.txt"
        
        # Fast extraction using seqkit - output is fastq format
        seqkit grep -f "$OUTDIR/target_read_ids_clean.txt" "$FASTQ_FILE" -o "$OUTDIR/target_dna_sequences.fq"
        
        # Convert to FASTA
        seqkit fq2fa "$OUTDIR/target_dna_sequences.fq" -o "$OUTDIR/target_dna_sequences.fa"
        
        # Count and clean up
        extracted_count=$(grep -c '^>' "$OUTDIR/target_dna_sequences.fa" 2>/dev/null || echo 0)
        echo "Extracted $extracted_count DNA sequences"
        
        # Cleanup intermediate files
        rm -f "$OUTDIR/target_read_ids.txt" "$OUTDIR/target_read_ids_clean.txt" "$OUTDIR/target_dna_sequences.fq"
    else
        echo "ERROR: FASTQ file not found: $FASTQ_FILE"
        return 1
    fi
    
    # Record end time and calculate duration
    local END_TIME=$(date)
    local START_SEC=$(date -d "$START_TIME" +%s)
    local END_SEC=$(date -d "$END_TIME" +%s)
    local DURATION=$((END_SEC - START_SEC))
    
    # Create flag file with summary information
    local FLAG_FILE="$OUTROOT/${sample_prefix}_step1_completed.flag"
    cat > "$FLAG_FILE" << EOF
STEP1_COMPLETED: $(date)
SAMPLE: $sample_prefix
STATUS: SUCCESS
TARGET_IDS_PROCESSED: ${#TARGET_SET[@]}
READS_ASSIGNED: $read_count
SEQUENCES_EXTRACTED: $extracted_count
DURATION_SECONDS: $DURATION
DURATION_HUMAN: $(date -u -d @$DURATION +'%H:%M:%S')
OUTPUT_FILE: ${sample_prefix}_dna_seqs/target_dna_sequences.fa
EOF
    
    echo "Flag file created: $FLAG_FILE"
    echo "Step 1 complete for $sample_prefix"
    echo "Duration: $(date -u -d @$DURATION +'%H:%M:%S')"
    
    return 0
}

# Parse command-line arguments
WORKING_DIRECTORY=""
ROOT=""
TARGETS=""
OUTROOT=""
FASTQ_DIR=""
SAMPLE_INDEX="$SLURM_ARRAY_TASK_ID"

while [[ $# -gt 0 ]]; do
    case $1 in
        -w|--working-dir)
            WORKING_DIRECTORY="$2"
            shift 2
            ;;
        -r|--humann-root)
            ROOT="$2"
            shift 2
            ;;
        -t|--targets)
            TARGETS="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTROOT="$2"
            shift 2
            ;;
        -f|--fastq-dir)
            FASTQ_DIR="$2"
            shift 2
            ;;
        -i|--sample-index)
            SAMPLE_INDEX="$2"
            shift 2
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

# Validate required arguments
if [[ -z "$WORKING_DIRECTORY" || -z "$ROOT" || -z "$TARGETS" || -z "$OUTROOT" || -z "$FASTQ_DIR" ]]; then
    echo "ERROR: Missing required arguments"
    usage
    exit 1
fi

# Validate sample index is a number
if ! [[ "$SAMPLE_INDEX" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Sample index must be a positive integer: $SAMPLE_INDEX"
    exit 1
fi

# Validate directories and files exist
if [[ ! -d "$WORKING_DIRECTORY" ]]; then
    echo "ERROR: Working directory does not exist: $WORKING_DIRECTORY"
    exit 1
fi

if [[ ! -d "$ROOT" ]]; then
    echo "ERROR: HUMAnN2 root directory does not exist: $ROOT"
    exit 1
fi

if [[ ! -f "$TARGETS" ]]; then
    echo "ERROR: Targets file does not exist: $TARGETS"
    exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: FASTQ directory does not exist: $FASTQ_DIR"
    exit 1
fi

# Call the function
extract_target_dna_sequences "$WORKING_DIRECTORY" "$ROOT" "$TARGETS" "$OUTROOT" "$FASTQ_DIR" "$SAMPLE_INDEX"

# Function to create summary table from flag files
create_step1_summary_table() {
    local INPUT_DIR="$1"
    local OUTPUT_FILE="$2"
    
    echo "Creating Step 1 summary table..."
    
    # Get all flag files
    local flag_files=("$INPUT_DIR"/*_step1_completed.flag)
    if [[ ${#flag_files[@]} -eq 0 ]]; then
        echo "No flag files found for summary"
        return 0
    fi
    
    # Extract all unique field names from all flag files
    declare -A all_fields
    for flag_file in "${flag_files[@]}"; do
        while IFS=': ' read -r field value; do
            if [[ -n "$field" ]]; then
                all_fields["$field"]=1
            fi
        done < <(grep ':' "$flag_file")
    done
    
    # Create header row
    echo -n "SAMPLE" > "$OUTPUT_FILE"
    for field in "${!all_fields[@]}"; do
        echo -ne "\t$field" >> "$OUTPUT_FILE"
    done
    echo >> "$OUTPUT_FILE"  # Newline
    
    # Process each flag file
    for flag_file in "${flag_files[@]}"; do
        local sample_prefix=$(basename "$flag_file" | sed 's/_step1_completed.flag$//')
        echo -n "$sample_prefix" >> "$OUTPUT_FILE"
        
        # Create associative array for this sample's data
        declare -A sample_data
        while IFS=': ' read -r field value; do
            if [[ -n "$field" ]]; then
                # Remove leading/trailing whitespace from value
                value=$(echo "$value" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')
                sample_data["$field"]="$value"
            fi
        done < "$flag_file"
        
        # Output data in the same order as header
        for field in "${!all_fields[@]}"; do
            if [[ -n "${sample_data[$field]}" ]]; then
                echo -ne "\t${sample_data[$field]}" >> "$OUTPUT_FILE"
            else
                echo -ne "\tNA" >> "$OUTPUT_FILE"
            fi
        done
        echo >> "$OUTPUT_FILE"  # Newline
    done
    
    echo "Summary table created: $OUTPUT_FILE"
    echo "Rows: ${#flag_files[@]}, Columns: $((${#all_fields[@]} + 1))"
}

# Call the table function
SUMMARY_TABLE="$OUTPUT_DIR/step1_samples_summary.tsv"
create_step1_summary_table "$OUTPUT_DIR" "$SUMMARY_TABLE"

# Final status
EXIT_CODE=$?
if [[ $EXIT_CODE -eq 0 ]]; then
    echo "Step 1 completed successfully for sample index $SAMPLE_INDEX"
else
    echo "Step 1 failed for sample index $SAMPLE_INDEX with exit code $EXIT_CODE"
fi

date
exit $EXIT_CODE
