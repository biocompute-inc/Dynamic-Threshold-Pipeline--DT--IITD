#!/bin/bash
#
# Example: Run methylation decoder pipeline
#
# This script demonstrates how to use the pipeline with your data.
# Modify the paths below to match your file locations.
#

# ============================================================================
# CONFIGURATION - EDIT THESE PATHS
# ============================================================================

# Directory containing your BAM files
BAM_DIR="../epic_new_bam"

# Reference FASTA file
REFERENCE="../reference.fasta"

# Methylated oligo sequences document
OLIGO_DOC="../methylated oligo sequences.docx"

# Output directory for results
OUTPUT_DIR="./results_$(date +%Y%m%d_%H%M%S)"

# ============================================================================
# OPTIONAL PARAMETERS
# ============================================================================

# Conservative bias (0.0-1.0, recommended: 0.8)
CONSERVATIVE_BIAS=0.8

# Modkit confidence threshold (0.0-1.0, recommended: 0.7)
MODKIT_CONFIDENCE=0.7

# Expected encoded message
EXPECTED_MESSAGE="EpiC"

# Number of parallel threads
THREADS=4

# ============================================================================
# RUN PIPELINE
# ============================================================================

echo "========================================================================================================"
echo "  METHYLATION DECODER PIPELINE"
echo "========================================================================================================"
echo ""
echo "BAM directory:        $BAM_DIR"
echo "Reference:            $REFERENCE"
echo "Oligo document:       $OLIGO_DOC"
echo "Output directory:     $OUTPUT_DIR"
echo ""
echo "Conservative bias:    $CONSERVATIVE_BIAS"
echo "Modkit confidence:    $MODKIT_CONFIDENCE"
echo "Expected message:     '$EXPECTED_MESSAGE'"
echo "Threads:              $THREADS"
echo ""
echo "========================================================================================================"
echo ""

python3 master_pipeline.py \
    --bam-dir "$BAM_DIR" \
    --reference "$REFERENCE" \
    --oligo-doc "$OLIGO_DOC" \
    --output-dir "$OUTPUT_DIR" \
    --conservative-bias "$CONSERVATIVE_BIAS" \
    --modkit-confidence "$MODKIT_CONFIDENCE" \
    --expected-message "$EXPECTED_MESSAGE" \
    --threads "$THREADS"

# Check if pipeline succeeded
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================================================================================"
    echo "  PIPELINE COMPLETED SUCCESSFULLY!"
    echo "========================================================================================================"
    echo ""
    echo "Results saved to: $OUTPUT_DIR"
    echo ""
    echo "View summary:"
    echo "  cat $OUTPUT_DIR/decoded/summary.txt"
    echo ""
    echo "View detailed results:"
    echo "  ls -lh $OUTPUT_DIR/decoded/"
    echo ""
else
    echo ""
    echo "========================================================================================================"
    echo "  PIPELINE FAILED - CHECK LOGS"
    echo "========================================================================================================"
    echo ""
    echo "Check log file:"
    echo "  cat $OUTPUT_DIR/pipeline.log"
    echo ""
    exit 1
fi
