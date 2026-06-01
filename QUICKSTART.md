# Quick Start Guide

Get started with the Methylation Decoder in 5 minutes!

## 1. Install Dependencies

```bash
# Install Python packages
pip install -r requirements.txt

# Install samtools (choose one method)
conda install -c bioconda samtools
# OR
sudo apt-get install samtools

# Install modkit
wget https://github.com/nanoporetech/modkit/releases/latest/download/modkit-x86_64-unknown-linux-gnu.tar.gz
tar -xzf modkit-*.tar.gz
export PATH=$PATH:$(pwd)
```

## 2. Prepare Your Data

You need:
- ✓ BAM file(s) with methylation tags (from dorado basecaller)
- ✓ Reference FASTA file
- ✓ Methylated oligo sequences DOCX file

## 3. Run the Pipeline

### Option A: Use the Example Script (Easiest)

```bash
# Edit the paths in example_run.sh
nano example_run.sh

# Run it
./example_run.sh
```

### Option B: Run Directly

```bash
python3 master_pipeline.py \
    --bam-dir /path/to/bam/files \
    --reference reference.fasta \
    --oligo-doc "methylated oligo sequences.docx" \
    --output-dir results/
```

### Option C: Single BAM File

```bash
python3 master_pipeline.py \
    --bam sample.bam \
    --reference reference.fasta \
    --oligo-doc "methylated oligo sequences.docx" \
    --output-dir results/
```

## 4. View Results

```bash
# Quick summary
cat results/decoded/summary.txt

# Detailed JSON results
cat results/decoded/decoded_results.json

# View heatmaps
ls results/decoded/sample_heatmaps/*.png
```

## Example Output

```
======================================================================================================
  METHYLATION-BASED DNA DATA DECODER - PIPELINE RESULTS
======================================================================================================

Sample                         Threshold    Signal     Decoded      Match%     Status
------------------------------------------------------------------------------------------------------
0_0.7                                10.0% WEAK       'EpiC'        100.0% ✓
1_0.7                                 9.5% WEAK       'EpiC'        100.0% ✓
2_0.7                                 9.5% WEAK       'EpiC'        100.0% ✓
======================================================================================================

SUCCESS RATE: 3/3 (100.0%) perfect matches
```

## Common Issues

### "samtools not found"
```bash
conda install -c bioconda samtools
```

### "modkit not found"  
Make sure modkit is in your PATH:
```bash
export PATH=$PATH:/path/to/modkit
```

### "No methylation signal"
- Check BAM has MM/ML tags: `samtools view sample.bam | head | grep MM`
- Try lower conservative bias: `--conservative-bias 0.5`

## Advanced Usage

### Custom Conservative Bias

```bash
# More conservative (avoid false positives)
python3 master_pipeline.py ... --conservative-bias 1.0

# Balanced
python3 master_pipeline.py ... --conservative-bias 0.5

# Recommended (default)
python3 master_pipeline.py ... --conservative-bias 0.8
```

### Decode Different Message

```bash
python3 master_pipeline.py ... --expected-message "AHEAD"
```

### Use Existing BED Files

```bash
# Skip modkit step, use pre-generated BEDs
python3 master_pipeline.py \
    --skip-modkit \
    --bed-dir existing_beds/ \
    --reference reference.fasta \
    --oligo-doc "methylated oligo sequences.docx" \
    --output-dir reanalysis/
```

## What's Next?

- Read full [README.md](README.md) for detailed documentation
- Check the [docs/](docs/) folder for pipeline details
- Review example outputs in `test_data/`

## Need Help?

- Check pipeline.log for detailed error messages
- Review the troubleshooting section in README.md
- Open an issue on GitHub
