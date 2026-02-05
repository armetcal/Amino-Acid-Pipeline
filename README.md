# Metagenomic Protein Sequence Extraction Pipeline

## Overview
This pipeline extracts protein sequences that correspond to UniRef IDs of interest from metagenomic data that has been processed using HUMAnN2. It identifies DNA sequences assigned to target UniRef50 IDs, performs six-frame translation, and validates the resulting protein sequences via Diamond BLAST against the UniRef50 database. The output protein sequences are combined into a single file, ready for comparison with proteomic data.

## Pipeline Architecture

### Two-Step Process:
1. **Step 1**: DNA sequence extraction from HUMAnN2 results (parallel array job)
2. **Step 2**: Six-frame translation + Diamond BLAST validation (single batch job)

## Quick Start

### Basic Usage:
```bash
# Run the full pipeline with automatic sample detection
./extract_peptides.sh \
    --working-dir $(pwd) \
    --humann-out path/to/humann/temp/folders \
    --targets list_of_UniRef_IDs.txt \
    --fastq-dir path/to/clean/unannotated/reads \
    --output-dir extract_peptides_outputs \
    --database path/to/diamond/database.dmnd \
    --max-targets 5 \
    --evalue 1e-3 \
    --pident 90 \
    --min-length 7 \
    --sensitive \
    --rerun
```

## File Structure

```
working_directory/
├── extract_peptides_outputs/
│   ├── sample1_dna_seqs/               # Step 1: Extracted DNA sequences
│   │   └── target_dna_sequences.fa
│   ├── step1_samples_summary.tsv       # Step 1: Sample statistics
│   ├── *_step1_completed.flag          # Step 1: Completion markers
│   ├── diamond_blast_results.tsv       # Step 2: BLAST results (reusable)
│   ├── filtered_blast_hits.tsv         # Step 2: Quality-filtered hits
│   ├── all_samples_deduplicated.faa    # Step 2: Final protein sequences
│   └── logs/                           # SLURM job logs                             
└── targets.txt  # Target UniRef IDs
```

## Input Requirements

### Required Files:
1. **HUMAnN2 Output**: Directory containing `*_humann_temp` folders with Diamond alignment files
2. **Target IDs**: Text file with UniRef50 accession numbers (one per line)
3. **Original FASTQ**: Quality-controlled reads used as HUMAnN2 input
4. **Diamond Database**: UniRef50 database in `.dmnd` format

### File Formats:
- **Target IDs**: `UniRef50_A0A010ZXN6` (one per line)
- **FASTQ**: Original sequencing reads (gzipped supported)
- **HUMAnN2 files**: Standard HUMAnN2 output structure (folder of temp folders)
- - Only the diamond_aligned temp files are currently supported - others can be deleted to reduce storage requirements

## Output Files

### Step 1 Outputs:
- `target_dna_sequences.fa`: DNA sequences assigned to target IDs
- `step1_samples_summary.tsv`: Tab-separated summary statistics
- Completion flag files for pipeline tracking

### Step 2 Outputs:
- `all_samples_deduplicated.faa`: Final protein sequences in proteomics format
- `diamond_blast_results.tsv`: Raw BLAST results (preserved for reruns)
- `filtered_blast_hits.tsv`: Quality-filtered BLAST hits

### Final Sequence Format:
```
>UniRef50_Q2NFV2_1 GN=UniRef50_Q2NFV2_1
MPVVYIMFPLVSFIVSIYFKYYSACVYWISIYVVYMCF*VDSISYLCIYC
>UniRef50_Q2NFV2_2 GN=UniRef50_Q2NFV2_2
AV*KACCLVIWILSYKATLPAPPH
```

## Quality Control Parameters

### Default Thresholds:
- **Identity**: ≥90% sequence identity to UniRef50
- **Length**: ≥7 amino acids
- **E-value**: 1e-3 (matches HUMAnN2 default)
- **Max targets**: 5 hits per query

### adjustable Parameters:
```bash
--pident 90      # Percent identity cutoff (0-100)
--min-length 7   # Minimum protein length
--evalue 1e-3    # BLAST e-value threshold
--max-targets 5  # Maximum hits per query
--sensitive      # Use sensitive Diamond mode
```

## Rerun Capabilities

The pipeline supports efficient reruns using existing BLAST results:

```bash
# Rerun with different thresholds (no BLAST recomputation)
sbatch run_full_pipeline.sh [...] --rerun
```

## SLURM Configuration

### Resource Requirements:
- **Step 1**: 4 CPUs, 8GB RAM per sample (array job)
- **Step 2**: 10 CPUs, 64GB RAM (single job)
- **Time**: 30 minutes (Step 1) + 1 hour (Step 2) (May need to adjust for very large jobs)

### Job Dependencies:
Step 2 automatically waits for Step 1 completion via SLURM job dependencies.

## Validation Features

### Quality Controls:
1. **Six-frame translation** with bacterial genetic code
2. **Diamond BLAST validation** against full UniRef50 database
3. **Stringent filtering** (identity + length thresholds)
4. **Duplicate removal** at DNA and protein levels

### Statistics Collected:
- Targets matched vs. new targets discovered
- Frame translation success rates
- Hit quality distributions (perfect hits, high-quality hits)
- Processing duration and sample counts

## Troubleshooting

### Common Issues:
1. **Missing HUMAnN2 files**: Ensure `*_humann_temp` directories exist
2. **Permission errors**: Check output directory write permissions
3. **Empty results**: Verify target IDs exist in UniRef50 database
4. **SLURM errors**: Check job resource requirements match cluster limits

### Debugging:
Check intermediate files, logs, and flag file contents.

## Dependencies

### Software:
- Diamond alignment tool
- seqkit sequence processing
- HUMAnN2 metagenomic pipeline
- SLURM workload manager

### Environments:
- Conda environment with `humann` package
- Standard Environment modules

## Citation

If you use this pipeline in your research, please cite the relevant tools:
- HUMAnN v3
- Diamond v2.1.11
- seqkit v2.5.1

## Support

For issues or questions, check the SLURM job logs in the `logs/` directory and ensure all input files meet the specified format requirements.
