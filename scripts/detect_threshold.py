#!/usr/bin/env python3
"""
Step 3: Detect optimal thresholds using conservative dynamic detection

Analyzes signal separation and determines optimal threshold per sample.
"""

import json
import logging
from pathlib import Path
import sys

# Add lib directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))

from dynamic_threshold_detector import DynamicThresholdDetector


logger = logging.getLogger(__name__)


def message_to_binary(message):
    """Convert ASCII message to binary representation"""
    binary = ''
    for char in message:
        byte = format(ord(char), '08b')
        binary += byte
    return binary


def detect_thresholds(bed_files, reference, oligo_doc, expected_message='EpiC',
                     conservative_bias=0.8, output_dir=None, verbose=False):
    """
    Detect optimal thresholds for all BED files.

    Args:
        bed_files: List of BED file paths
        reference: Path to reference FASTA
        oligo_doc: Path to methylated oligo sequences DOCX
        expected_message: Expected encoded message
        conservative_bias: Conservative bias (0.0-1.0)
        output_dir: Output directory for JSON results
        verbose: Verbose output

    Returns:
        Dictionary of threshold results per sample
    """
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    # Convert message to binary
    expected_binary = message_to_binary(expected_message)

    logger.info(f"Expected message: '{expected_message}' = {expected_binary[:32]}")
    logger.info(f"Conservative bias: {conservative_bias}")

    # Initialize detector
    detector = DynamicThresholdDetector(
        verbose=verbose,
        conservative_bias=conservative_bias
    )

    results = {}

    logger.info(f"\nAnalyzing {len(bed_files)} sample(s)...")

    for bed_file in bed_files:
        bed_path = Path(bed_file)
        sample_name = bed_path.stem

        logger.info(f"\n{'='*80}")
        logger.info(f"  {sample_name}")
        logger.info(f"{'='*80}")

        try:
            result = detector.detect_threshold(str(bed_file), expected_binary)

            results[sample_name] = {
                'threshold': float(result.threshold),
                'confidence': result.confidence,
                'signal_type': result.signal_type,
                'recommendation': result.recommendation,
                'practical_warning': result.practical_warning,
                'statistics': {
                    'mean_zeros': float(result.statistics.mean_zeros) if result.statistics else 0.0,
                    'mean_ones': float(result.statistics.mean_ones) if result.statistics else 0.0,
                    'fold_change': float(result.statistics.fold_change_mean) if result.statistics else 0.0,
                    'cohens_d': float(result.statistics.separation_score) if result.statistics else 0.0,
                    'is_practically_viable': bool(result.statistics.is_practically_viable) if result.statistics else False
                }
            }

            logger.info(f"\n✓ Threshold: {result.threshold}% ({result.confidence})")

        except Exception as e:
            logger.error(f"✗ Failed to detect threshold for {sample_name}: {e}")
            results[sample_name] = {
                'threshold': 10.0,
                'confidence': 'ERROR',
                'signal_type': 'FAILED',
                'error': str(e)
            }

    # Save results
    if output_dir:
        output_file = output_dir / 'threshold_results.json'
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"\n✓ Threshold results saved to {output_file}")

    # Print summary
    logger.info(f"\n{'='*80}")
    logger.info("  THRESHOLD DETECTION SUMMARY")
    logger.info(f"{'='*80}")
    logger.info(f"{'Sample':<30} {'Threshold':<12} {'Confidence':<30} {'Fold-Change':<12}")
    logger.info("-"*80)

    for name, res in results.items():
        fold = res.get('statistics', {}).get('fold_change', 0)
        logger.info(f"{name:<30} {res['threshold']:>10.1f}% "
                   f"{res['confidence']:<30} {fold:>10.2f}x")

    logger.info("="*80)

    return results


if __name__ == '__main__':
    # Standalone usage
    import argparse

    parser = argparse.ArgumentParser(description='Detect optimal thresholds')
    parser.add_argument('bed_files', nargs='+', help='BED files to analyze')
    parser.add_argument('--reference', required=True, help='Reference FASTA')
    parser.add_argument('--oligo-doc', required=True, help='Oligo sequences DOCX')
    parser.add_argument('--expected-message', default='EpiC',
                       help='Expected message (default: EpiC)')
    parser.add_argument('--conservative-bias', type=float, default=0.8,
                       help='Conservative bias 0.0-1.0 (default: 0.8)')
    parser.add_argument('--output-dir', help='Output directory')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO if not args.verbose else logging.DEBUG,
                       format='%(asctime)s [%(levelname)s] %(message)s')

    detect_thresholds(
        args.bed_files,
        args.reference,
        args.oligo_doc,
        args.expected_message,
        args.conservative_bias,
        args.output_dir,
        args.verbose
    )
