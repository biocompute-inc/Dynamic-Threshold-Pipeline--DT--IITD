# Methylation-Based DNA Data Decoder

Complete automated pipeline for decoding binary messages encoded in DNA methylation patterns.

## Overview

This pipeline takes BAM files from nanopore sequencing with methylation calling and automatically:
1. Indexes BAM files (creates .bai files)
2. Generates methylation BED files using modkit
3. Detects optimal thresholds using conservative dynamic detection
4. Decodes binary sequences from methylation patterns
5. Converts binary to ASCII messages

## Features

- **Fully Automated**: Single command processes entire pipeline
- **Conservative Dynamic Thresholding**: Automatically optimizes per-sample thresholds (bias=0.8)
- **Statistical Analysis**: Provides fold-change, Cohen's d, signal quality metrics
- **Batch Processing**: Handle multiple BAM files simultaneously
- **Comprehensive Reports**: HTML and text summaries with heatmaps

## Requirements

### System Dependencies
- Python 3.8+
- samtools (for BAM indexing)
- modkit (for methylation calling)

### Python Dependencies
```bash
pip install -r requirements.txt
```

Required packages:
- numpy
- scipy
- python-docx
- matplotlib
- seaborn
- pandas

## Quick Start

### 1. Basic Usage

```bash
# Process all BAM files in a folder
python master_pipeline.py \
    --bam-dir /path/to/bam/files \
    --reference reference.fasta \
    --oligo-doc "methylated oligo sequences.docx" \
    --output-dir results/

# Process single BAM file
python master_pipeline.py \
    --bam sample.bam \
    --reference reference.fasta \
    --oligo-doc "methylated oligo sequences.docx" \
    --output-dir results/
```

### 2. With Custom Parameters

```bash
python master_pipeline.py \
    --bam-dir bam_files/ \
    --reference reference.fasta \
    --oligo-doc "methylated oligo sequences.docx" \
    --output-dir results/ \
    --conservative-bias 0.8 \
    --modkit-confidence 0.7 \
    --expected-message "EpiC" \
    --threads 4
```

## Pipeline Steps

### Step 1: BAM Indexing
Creates .bai index files for all BAM files using samtools.

### Step 2: BED Generation
Runs modkit pileup to extract methylation data:
- Confidence threshold: 0.7 (default)
- CpG-only filtering
- bedMethyl format output

### Step 3: Threshold Detection
Uses conservative dynamic threshold detector:
- Separates expected 0s and 1s
- Calculates fold-change and Cohen's d
- Applies conservative bias (0.8 default)
- Detects optimal threshold per sample

### Step 4: Binary Decoding
Extracts methylation at 36 CpG positions:
- Applies detected threshold (methylation% > threshold → bit 1)
- Generates 36-bit binary sequence

### Step 5: ASCII Conversion
Converts first 32 bits to ASCII:
- Groups 8 bits → 1 byte
- Converts to decimal → ASCII character
- Returns 4-character message

## Output Structure

```
results/
├── bam_indices/          # .bai index files
├── bed_files/            # bedMethyl files from modkit
├── thresholds/           # JSON threshold detection results
├── decoded/              # Final decoded sequences
│   ├── summary.txt
│   ├── detailed_report.html
│   └── sample_heatmaps/
└── pipeline.log          # Complete pipeline log
```

## Output Files

### summary.txt
Quick overview of all samples:
```
Sample          Threshold  Signal    Decoded  Match%  Status
sample1.bam     10.0%      WEAK      'EpiC'   100.0%  ✓ PERFECT
sample2.bam     9.5%       WEAK      'EpiC'   100.0%  ✓ PERFECT
```

### detailed_report.html
Interactive HTML report with:
- Methylation heatmaps
- Bit-by-bit comparison
- Statistical metrics
- Signal quality analysis

### threshold JSON
Per-sample threshold detection:
```json
{
  "sample1": {
    "threshold": 10.0,
    "confidence": "HIGH_CONSERVATIVE",
    "fold_change": 13.62,
    "cohens_d": 5.80
  }
}
```

## Advanced Usage

### Custom Conservative Bias

```bash
# More conservative (avoid false positives)
python master_pipeline.py --bam-dir bam/ --conservative-bias 1.0

# Balanced (original dynamic)
python master_pipeline.py --bam-dir bam/ --conservative-bias 0.0

# Recommended (best performance)
python master_pipeline.py --bam-dir bam/ --conservative-bias 0.8
```

### Decode Different Message

```bash
# Expecting "AHEAD" instead of "EpiC"
python master_pipeline.py \
    --bam-dir bam/ \
    --expected-message "AHEAD"
```

### Skip Steps (Resume Pipeline)

```bash
# If BAM files already indexed
python master_pipeline.py --bam-dir bam/ --skip-indexing

# If BED files already exist
python master_pipeline.py --bam-dir bam/ --skip-modkit --bed-dir existing_beds/
```

## Command-Line Options

```
Required:
  --bam-dir DIR              Directory containing BAM files
  --bam FILE                 Single BAM file to process
  --reference FILE           Reference FASTA file
  --oligo-doc FILE           Methylated oligo sequences (DOCX)
  --output-dir DIR           Output directory

Optional:
  --conservative-bias FLOAT  Conservative bias 0.0-1.0 (default: 0.8)
  --modkit-confidence FLOAT  Modkit filter threshold (default: 0.7)
  --expected-message TEXT    Expected encoded message (default: "EpiC")
  --threads INT              Number of parallel threads (default: 4)
  --skip-indexing            Skip BAM indexing step
  --skip-modkit              Skip modkit pileup step
  --bed-dir DIR              Use existing BED files
  --verbose                  Verbose output
  --quiet                    Minimal output
```

## Example Workflows

### 1. Process New Sequencing Run

```bash
# Complete pipeline from scratch
python master_pipeline.py \
    --bam-dir /sequencing_run/bam_files/ \
    --reference reference.fasta \
    --oligo-doc "methylated oligo sequences.docx" \
    --output-dir run_2026_05_17/ \
    --threads 8
```

### 2. Re-analyze with Different Threshold

```bash
# Use existing BED files, re-run threshold detection
python master_pipeline.py \
    --skip-modkit \
    --bed-dir previous_run/bed_files/ \
    --reference reference.fasta \
    --oligo-doc "methylated oligo sequences.docx" \
    --conservative-bias 0.9 \
    --output-dir reanalysis/
```

### 3. Quality Control Check

```bash
# Run with verbose output and maximum conservative bias
python master_pipeline.py \
    --bam sample.bam \
    --reference reference.fasta \
    --oligo-doc "methylated oligo sequences.docx" \
    --output-dir qc_check/ \
    --conservative-bias 1.0 \
    --verbose
```

## Interpreting Results

### Signal Quality

- **Fold-change > 10×**: Excellent separation between 0s and 1s
- **Cohen's d > 2.0**: Statistically significant separation
- **Practical viability**: Mean 1s > 10%, absolute difference > 5%

### Decoding Status

- **✓ PERFECT (100%)**: All 32 bits match expected pattern
- **★ EXCELLENT (≥95%)**: 1-2 bit errors
- **◐ GOOD (≥90%)**: 3-4 bit errors
- **○ PARTIAL (≥75%)**: 5-8 bit errors
- **✗ FAILED (<75%)**: Poor signal or wrong encoding

### Conservative Bias Impact

| Bias | Behavior | Use Case |
|------|----------|----------|
| 0.0  | Balanced accuracy | Clean data, trust signal |
| 0.5  | Moderate conservative | Some background noise |
| **0.8**  | **Recommended** | **Standard methylation data** |
| 1.0  | Maximum specificity | High noise, critical applications |

## Troubleshooting

### "samtools not found"
```bash
conda install -c bioconda samtools
```

### "modkit not found"
```bash
wget https://github.com/nanoporetech/modkit/releases/latest/download/modkit-x86_64-unknown-linux-gnu.tar.gz
tar -xzf modkit-*.tar.gz
export PATH=$PATH:$(pwd)
```

### "No methylation signal detected"
- Check BAM file has methylation tags (MM, ML)
- Verify correct reference FASTA
- Try lower conservative bias (0.5)
- Check modkit confidence threshold

### "All samples decode to wrong message"
- Verify expected message is correct
- Check oligo sequences DOCX file
- Ensure reference FASTA matches experimental design

## Citation

If you use this pipeline in your research, please cite:

```
Methylation-Based DNA Data Decoder
BioCompute × IIT-Delhi, 2026
https://github.com/your-repo/methylation-decoder
```

## License

MIT License - See LICENSE file

## Contact

For questions or issues:
- GitHub Issues: https://github.com/your-repo/methylation-decoder/issues
- Email: support@example.com

## Version History

- **v1.0.0** (2026-05-17): Initial release
  - Conservative dynamic thresholding
  - Batch processing support
  - Comprehensive reporting
