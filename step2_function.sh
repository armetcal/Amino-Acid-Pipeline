#!/bin/bash

# This script performs six-frame translation of DNA sequences using bacterial genetic code, then validates the translations 
# against the UniRef50 database via Diamond BLAST to identify high-quality protein sequences matching target IDs with 
# ≥90% identity and ≥7 amino acid length.
# The script combines sequences from all samples, deduplicates them, and processes them through a quality-controlled 
# validation pipeline to produce final amino acid sequences for downstream peptide matching.

# Example usage:~~~~
# ./step2_function.sh \
#  --working-dir /scratch/armetcal/cmmi/pilot_g4h/extract_relevant_aa_seqs \
#  --database /home/armetcal/scratch/humann_databases/uniref50/uniref50_201901b_full.dmnd \
#  --targets mtx_uniprot50_prev50_abun1_bug3_IDs.txt \
#  --input-dir temp_pipeline/step1 \
#  --output-dir temp_pipeline/step2 \
#  --max-targets 5 \
#  --evalue 1e-3 \
#  --pident 90 \
#  --min-length 7 \
#  --sensitive \ # Remove if you want fast mode
#  --rerun # Optional: Rerun filtering using existing BLAST output
#~~~~~~~~~~~~~~~~~~~

# Load environment
source /home/armetcal/miniconda3/etc/profile.d/conda.sh
conda activate humann
module load seqkit/2.5.1 StdEnv/2023

# Usage function
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Performs six-frame translation and Diamond BLAST validation of DNA sequences"
    echo ""
    echo "OPTIONS:"
    echo "  -w, --working-dir DIR      Working directory (required)"
    echo "  -d, --database FILE        Diamond database file (required)"
    echo "  -t, --targets FILE         Targets file with UniRef50 IDs (required)"
    echo "  -i, --input-dir DIR        Input directory from Step 1 (required)"
    echo "  -o, --output-dir DIR       Output directory (required)"
    echo "  -m, --max-targets INT      Max targets per query (default: 5)"
    echo "  -e, --evalue FLOAT         E-value threshold (default: 1e-3)"
    echo "  -p, --pident INT           Percent identity cutoff (default: 90)"
    echo "  -l, --min-length INT       Minimum amino acid length (default: 7)"
    echo "  -s, --sensitive            Use sensitive mode (default: true)"
    echo "  -r, --rerun                Rerun filtering only using existing BLAST output (default: false)"
    echo "  -h, --help                 Show this help message"
    echo ""
    echo "EXAMPLE:"
    echo "  $0 --working-dir /scratch/user/project --database /path/to/uniref50.dmnd --targets targets.txt --input-dir /path/to/step1 --output-dir /path/to/step2 --pident 95 --min-length 10 --sensitive"
}

# Function to perform translation and BLAST validation
translate_and_validate_sequences() {
    local WORKING_DIRECTORY="$1"
    local DATABASE="$2"
    local TARGETS="$3"
    local INPUT_DIR="$4"
    local OUTPUT_DIR="$5"
    local MAX_TARGETS="$6"
    local EVALUE="$7"
    local PIDENT_CUTOFF="$8"
    local MIN_LENGTH="$9"
    local SENSITIVE_MODE="${10}"
    local RERUN_MODE="${11}"
    
    # Record start time
    local START_TIME=$(date)
    
    # Create output directory
    mkdir -p "$OUTPUT_DIR"
    
    echo "=== Step 2: Translation and BLAST Validation ==="
    echo "Start time: $START_TIME"
    echo "Working directory: $WORKING_DIRECTORY"
    echo "Database: $DATABASE"
    echo "Targets file: $TARGETS"
    echo "Input directory: $INPUT_DIR"
    echo "Output directory: $OUTPUT_DIR"
    echo "Parameters: max-targets=$MAX_TARGETS, evalue=$EVALUE, pident=$PIDENT_CUTOFF, min-length=$MIN_LENGTH, sensitive=$SENSITIVE_MODE, rerun=$RERUN_MODE"
    
    # File paths
    local COMBINED_DNA="$OUTPUT_DIR/combined_dna_sequences.fa"
    local DEDUP_DNA="$OUTPUT_DIR/deduplicated_dna_sequences.fa"
    local TRANSLATED_AA="$OUTPUT_DIR/six_frame_translated.fa"
    local BLAST_OUTPUT="$OUTPUT_DIR/diamond_blast_results.tsv"
    local FILTERED_HITS="$OUTPUT_DIR/filtered_blast_hits.tsv"
    local HIGH_QUALITY_AA="$OUTPUT_DIR/high_quality_aa_sequences.fa"
    local REFORMATTED_AA="$OUTPUT_DIR/final_format_aa_sequences.faa"
    
    if [[ "$RERUN_MODE" == "true"  && -s "$BLAST_OUTPUT" ]]; then        
        local actual_lines=$(wc -l < "$BLAST_OUTPUT" 2>/dev/null || echo 0)
        echo "Rerun mode: using existing BLAST output for filtering"
    else
        echo "Rerun mode not enabled or BLAST output missing; proceeding with full pipeline execution"
        RERUN_MODE="false"
    fi
    
    if [[ "$RERUN_MODE" == "false" ]]; then
        # Full pipeline execution
        # Check if Step 1 completed successfully
        local step1_flags=("$INPUT_DIR"/*_step1_completed.flag)
        if [[ ${#step1_flags[@]} -eq 0 ]]; then
            echo "ERROR: No Step 1 completion flags found in $INPUT_DIR"
            return 1
        fi

        echo "Found ${#step1_flags[@]} Step 1 completion flags"

        # Verify each sample has DNA sequences available
        local samples_with_dna=0
        for flag_file in "${step1_flags[@]}"; do
            local sample_prefix=$(basename "$flag_file" | sed 's/_step1_completed.flag$//')
            local sample_dna_file="$INPUT_DIR/${sample_prefix}_dna_seqs/target_dna_sequences.fa"
            
            if [[ -s "$sample_dna_file" ]] && (grep -q "STATUS: SUCCESS" "$flag_file" || grep -q "SEQUENCES_EXTRACTED: [1-9]" "$flag_file"); then
                ((samples_with_dna++))
            else
                echo "WARNING: Step 1 incomplete or no DNA sequences for sample: $sample_prefix"
            fi
        done

        if [[ $samples_with_dna -eq 0 ]]; then
            echo "ERROR: No samples have DNA sequences available for processing"
            return 1
        fi

        echo "Verified $samples_with_dna samples have DNA sequences available"
        
        # 1) Combine all DNA sequences from Step 1
        echo "Combining DNA sequences from Step 1..."
        > "$COMBINED_DNA"  # Clear file
        for flag_file in "${step1_flags[@]}"; do
            local sample_prefix=$(basename "$flag_file" | sed 's/_step1_completed.flag$//')
            local sample_dna_file="$INPUT_DIR/${sample_prefix}_dna_seqs/target_dna_sequences.fa"
            
            if [[ -s "$sample_dna_file" ]]; then
                cat "$sample_dna_file" >> "$COMBINED_DNA"
            fi
        done

        local total_dna_seqs=$(grep -c '^>' "$COMBINED_DNA" 2>/dev/null || echo 0)
        echo "Combined $total_dna_seqs DNA sequences"

        if [[ $total_dna_seqs -eq 0 ]]; then
            echo "ERROR: No DNA sequences found to process"
            return 1
        fi

        # 2) Deduplicate DNA sequences (identical only)
        echo "Deduplicating DNA sequences..."
        seqkit rmdup -s "$COMBINED_DNA" -o "$DEDUP_DNA"

        local dedup_count=$(grep -c '^>' "$DEDUP_DNA" 2>/dev/null || echo 0)
        echo "After deduplication: $dedup_count unique DNA sequences"

        # 3) Six-frame translation
        echo "Performing six-frame translation..."
        seqkit translate -f 6 -F -T 11 "$DEDUP_DNA" -o "$TRANSLATED_AA"

        local translated_count=$(grep -c '^>' "$TRANSLATED_AA" 2>/dev/null || echo 0)
        echo "Translated to $translated_count amino acid sequences"

        # 4) Diamond BLASTp against UniRef50 database
        echo "Running Diamond BLASTp..."
        local diamond_cmd="diamond blastp \
            --db \"$DATABASE\" \
            --query \"$TRANSLATED_AA\" \
            --out \"$BLAST_OUTPUT\" \
            --outfmt 6 qseqid sseqid pident length evalue bitscore \
            --max-target-seqs \"$MAX_TARGETS\" \
            --evalue \"$EVALUE\" \
            --threads 8"

        if [[ "$SENSITIVE_MODE" == "true" ]]; then
            diamond_cmd="$diamond_cmd --sensitive"
        else
            diamond_cmd="$diamond_cmd --fast"
        fi

        eval "$diamond_cmd"

        local blast_hits=$(wc -l < "$BLAST_OUTPUT" 2>/dev/null || echo 0)
        echo "Diamond BLAST completed: $blast_hits hits"
    else
        # Rerun mode: load existing counts
        local total_dna_seqs=0
        local dedup_count=0
        local translated_count=$(grep -c '^>' "$TRANSLATED_AA" 2>/dev/null || echo 0)
        local blast_hits=$(wc -l < "$BLAST_OUTPUT" 2>/dev/null || echo 0)
        echo "Rerun mode: Using existing data (translated: $translated_count, BLAST hits: $blast_hits)"
    fi
    
    # 5) Filter BLAST results
    echo "Filtering BLAST hits (${PIDENT_CUTOFF}% identity, ${MIN_LENGTH}+ AA length, target ID match)..."
    
    # Build target set for filtering
    declare -A TARGET_SET
    while IFS= read -r target; do
        target_id="${target%%|*}"
        TARGET_SET["$target_id"]=1
    done < "$TARGETS"
    
    # Filter BLAST results
    awk -v pident_cutoff="$PIDENT_CUTOFF" -v min_length="$MIN_LENGTH" -v targets="$(printf "%s\n" "${!TARGET_SET[@]}")" '
        BEGIN {
            split(targets, arr, "\n")
            for (i in arr) target_set[arr[i]] = 1
        }
        {
            # Extract UniRef50 ID from subject (column 2)
            subject = $2
            sub(/\|.*/, "", subject)
            
            # Check filters: target ID match, identity cutoff, length cutoff
            if (subject in target_set && $3 >= pident_cutoff && $4 >= min_length) {
                print $0
            }
        }
    ' "$BLAST_OUTPUT" > "$FILTERED_HITS"
    
    local filtered_count=$(wc -l < "$FILTERED_HITS" 2>/dev/null || echo 0)
    echo "Filtered to $filtered_count high-quality hits"
    
    # 6) Statistics
    echo "Calculating statistics for flag file..."
    
    # Target statistics
    cut -f2 "$FILTERED_HITS" | sed 's/|.*//' | sort -u > "$OUTPUT_DIR/matched_targets.txt"
    local unique_targets_count=$(wc -l < "$OUTPUT_DIR/matched_targets.txt" 2>/dev/null || echo 0)
    local original_targets_matched=$(grep -Fxf "$TARGETS" "$OUTPUT_DIR/matched_targets.txt" | wc -l 2>/dev/null || echo 0)
    local new_targets_discovered=$((unique_targets_count - original_targets_matched))
    
    # Frame statistics
    local frames_with_hits=$(cut -f1 "$FILTERED_HITS" | sed 's/_frame=.*//' | sort -u | wc -l 2>/dev/null || echo 0)
    local total_frames_searched=$((dedup_count * 6))
    
    # Quality distribution
    local perfect_hits=$(awk -v min_len="$MIN_LENGTH" '$3 == 100 && $4 >= min_len' "$FILTERED_HITS" | wc -l 2>/dev/null || echo 0)
    local high_quality_hits=$(awk -v min_len="$MIN_LENGTH" '$3 >= 95 && $4 >= min_len' "$FILTERED_HITS" | wc -l 2>/dev/null || echo 0)
    
    # 7) Extract high-quality AA sequences
    local final_count=0
    if [[ $filtered_count -gt 0 ]]; then
        # Extract unique query sequence IDs from filtered hits
        cut -f1 "$FILTERED_HITS" | sort -u > "$OUTPUT_DIR/high_quality_query_ids.txt"
        
        # Extract the corresponding AA sequences
        seqkit grep -f "$OUTPUT_DIR/high_quality_query_ids.txt" "$TRANSLATED_AA" -o "$HIGH_QUALITY_AA"
        
        final_count=$(grep -c '^>' "$HIGH_QUALITY_AA" 2>/dev/null || echo 0)
        echo "Extracted $final_count high-quality amino acid sequences"
        
        # Cleanup
        rm -f "$OUTPUT_DIR/high_quality_query_ids.txt"
    else
        echo "No high-quality hits found"
        touch "$HIGH_QUALITY_AA"  # Create empty file
    fi
    
    echo "Reformatting headers to proteomics format..."
    if [[ $final_count -gt 0 ]]; then
        # First, map query IDs to UniRef50 IDs and create temporary FASTA with UniRef headers
        local TEMP_FASTA="$OUTPUT_DIR/temp_uniref_mapped.faa"

        echo "Mapping FASTA headers to UniRef50 IDs..."
        echo "Input files:"
        echo "  BLAST: $FILTERED_HITS"
        echo "  FASTA: $HIGH_QUALITY_AA" 
        echo "  Output: $TEMP_FASTA"

        awk -F'\t' '
            # First pass: Build mapping from BLAST results
            NR == FNR {
                query_id = $1
                gsub(/[[:space:]]+$/, "", query_id)
                uniref_full = $2
                sub(/\|.*/, "", uniref_full)
                best_hit[query_id] = uniref_full
                next
            }
            
            # Second pass: Process FASTA file
            /^>/ {
                fasta_header = substr($0, 2)
                gsub(/[[:space:]]+$/, "", fasta_header)
                getline
                sequence = $0
                
                if (fasta_header in best_hit) {
                    uniref_id = best_hit[fasta_header]
                    printf ">%s\n%s\n", uniref_id, sequence
                }
            }
        ' "$FILTERED_HITS" "$HIGH_QUALITY_AA" > "$TEMP_FASTA"

        local temp_count=$(grep -c '^>' "$TEMP_FASTA" 2>/dev/null || echo 0)
        echo "Temporary FASTA created with $temp_count sequences"

        if [[ $temp_count -eq 0 ]]; then
            echo "ERROR: Mapping failed - check file paths and permissions"
            return 1
        fi

        # Now add unique numbering for identical UniRef IDs with different sequences
        echo "Adding unique numbering for identical UniRef IDs..."
        awk '
            /^>/ {
                if (sequence) {
                    # Store the previous sequence
                    gsub(/^>/, "", header)
                    uniref_base = header
                    
                    if (!(uniref_base in count)) {
                        count[uniref_base] = 0
                    }
                    count[uniref_base]++
                    
                    # Output with unique numbering
                    print ">" uniref_base "_" count[uniref_base] " GN=" uniref_base "_" count[uniref_base]
                    print sequence
                }
                header = $0
                sequence = ""
                next
            }
            {
                sequence = sequence $0
            }
            END {
                if (sequence) {
                    gsub(/^>/, "", header)
                    uniref_base = header
                    
                    if (!(uniref_base in count)) {
                        count[uniref_base] = 0
                    }
                    count[uniref_base]++
                    
                    print ">" uniref_base "_" count[uniref_base] " GN=" uniref_base "_" count[uniref_base]
                    print sequence
                }
            }
        ' "$TEMP_FASTA" > "$REFORMATTED_AA"
        
        local reformatted_count=$(grep -c '^>' "$REFORMATTED_AA" 2>/dev/null || echo 0)
        echo "Reformatted $reformatted_count sequences"
        
        # Cleanup temporary file
        rm -f "$TEMP_FASTA"
    else
        touch "$REFORMATTED_AA"
        echo "No sequences to reformat"
    fi

    # Record end time and calculate duration
    local END_TIME=$(date)
    local START_SEC=$(date -d "$START_TIME" +%s)
    local END_SEC=$(date -d "$END_TIME" +%s)
    local DURATION=$((END_SEC - START_SEC))
    
    # Create flag file with comprehensive summary information
    local FLAG_FILE="$OUTPUT_DIR/step2_completed.flag"
    cat > "$FLAG_FILE" << EOF
STEP2_COMPLETED: $(date)
MODE: $([[ "$RERUN_MODE" == "true" ]] && echo "RERUN" || echo "FULL")
PARAMETERS: max-targets=$MAX_TARGETS, evalue=$EVALUE, pident=$PIDENT_CUTOFF, min-length=$MIN_LENGTH, sensitive=$SENSITIVE_MODE
STATUS: SUCCESS
TARGET_IDS_PROCESSED: ${#TARGET_SET[@]}
ORIGINAL_TARGETS_MATCHED: $original_targets_matched
NEW_TARGETS_DISCOVERED: $new_targets_discovered
TOTAL_TARGETS_WITH_HITS: $unique_targets_count
DNA_SEQUENCES_COMBINED: $total_dna_seqs
DNA_SEQUENCES_DEDUP: $dedup_count
FRAMES_SEARCHED: $total_frames_searched
FRAMES_WITH_HITS: $frames_with_hits
AA_SEQUENCES_TRANSLATED: $translated_count
BLAST_HITS: $blast_hits
HIGH_QUALITY_HITS: $filtered_count
PERFECT_HITS_100PCT: $perfect_hits
HIGH_QUALITY_HITS_95PCT: $high_quality_hits
FINAL_AA_SEQUENCES: $final_count
DURATION_SECONDS: $DURATION
DURATION_HUMAN: $(date -u -d @$DURATION +'%H:%M:%S')
UNFORMATTED_OUTPUT_FILE: high_quality_aa_sequences.fa
FORMATTED_OUTPUT_FILE: final_format_aa_sequences.faa
BLAST_OUTPUT: diamond_blast_results.tsv
EOF
    
    echo "Flag file created: $FLAG_FILE"
    echo "Step 2 complete"
    echo "Duration: $(date -u -d @$DURATION +'%H:%M:%S')"
    
    return 0
}

# Parse command-line arguments
WORKING_DIRECTORY=""
DATABASE=""
TARGETS=""
INPUT_DIR=""
OUTPUT_DIR=""
MAX_TARGETS=5
EVALUE="1e-3"
PIDENT_CUTOFF=90
MIN_LENGTH=7
SENSITIVE_MODE="true"
RERUN_MODE="false"

while [[ $# -gt 0 ]]; do
    case $1 in
        -w|--working-dir)
            WORKING_DIRECTORY="$2"
            shift 2
            ;;
        -d|--database)
            DATABASE="$2"
            shift 2
            ;;
        -t|--targets)
            TARGETS="$2"
            shift 2
            ;;
        -i|--input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
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
            return 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$WORKING_DIRECTORY" || -z "$DATABASE" || -z "$TARGETS" || -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: Missing required arguments"
    usage
    return 1
fi

# Validate numeric parameters
if ! [[ "$MAX_TARGETS" =~ ^[0-9]+$ ]] || [[ "$MAX_TARGETS" -lt 1 ]]; then
    echo "ERROR: Max targets must be a positive integer: $MAX_TARGETS"
    return 1
fi

if ! [[ "$EVALUE" =~ ^[0-9.eE-]+$ ]]; then
    echo "ERROR: E-value must be a numeric value: $EVALUE"
    return 1
fi

if ! [[ "$PIDENT_CUTOFF" =~ ^[0-9]+$ ]] || [[ "$PIDENT_CUTOFF" -lt 0 || "$PIDENT_CUTOFF" -gt 100 ]]; then
    echo "ERROR: Percent identity must be between 0 and 100: $PIDENT_CUTOFF"
    return 1
fi

if ! [[ "$MIN_LENGTH" =~ ^[0-9]+$ ]] || [[ "$MIN_LENGTH" -lt 1 ]]; then
    echo "ERROR: Minimum length must be a positive integer: $MIN_LENGTH"
    return 1
fi

# Validate files and directories
if [[ ! -d "$WORKING_DIRECTORY" ]]; then
    echo "ERROR: Working directory does not exist: $WORKING_DIRECTORY"
    return 1
fi

if [[ ! -f "$DATABASE" ]]; then
    echo "ERROR: Database file does not exist: $DATABASE"
    return 1
fi

if [[ ! -f "$TARGETS" ]]; then
    echo "ERROR: Targets file does not exist: $TARGETS"
    return 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory does not exist: $INPUT_DIR"
    return 1
fi

# Call the function
echo "Starting Step 2: Translation and BLAST Validation"
translate_and_validate_sequences "$WORKING_DIRECTORY" "$DATABASE" "$TARGETS" "$INPUT_DIR" "$OUTPUT_DIR" "$MAX_TARGETS" "$EVALUE" "$PIDENT_CUTOFF" "$MIN_LENGTH" "$SENSITIVE_MODE" "$RERUN_MODE"

# Final status
EXIT_CODE=$?
if [[ $EXIT_CODE -eq 0 ]]; then
    echo "Step 2 completed successfully"
else
    echo "Step 2 failed with exit code $EXIT_CODE"
fi

date
exit $EXIT_CODE
