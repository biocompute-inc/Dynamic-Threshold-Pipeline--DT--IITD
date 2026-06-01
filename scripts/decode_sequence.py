#!/usr/bin/env python3
"""
Step 4: Decode binary sequences from methylation patterns

Uses detected thresholds to convert methylation data to binary and ASCII.
"""

import json
import logging
from pathlib import Path
import sys

# Add lib directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))

from brick_position_extractor import extract_brick_positions, load_bed_methylation


logger = logging.getLogger(__name__)


def message_to_binary(message):
    """Convert ASCII message to binary"""
    return ''.join(format(ord(char), '08b') for char in message)


def decode_single_sample(bed_file, threshold, brick_data, expected_binary):
    """
    Decode a single sample.

    Args:
        bed_file: Path to BED file
        threshold: Methylation threshold (%)
        brick_data: List of CpG position data
        expected_binary: Expected binary pattern

    Returns:
        Dictionary with decoding results
    """
    # Load methylation data
    meth_data = load_bed_methylation(bed_file)

    # Extract methylation at 36 positions and convert to bits
    pattern_bits = []
    meth_values = []

    for brick in brick_data[:36]:
        pos = brick['cpg_pos']
        if pos in meth_data:
            pct = meth_data[pos]
            meth_values.append(pct)
            if pct > threshold:
                pattern_bits.append(1)
            else:
                pattern_bits.append(0)
        else:
            meth_values.append(0.0)
            pattern_bits.append(0)

    pattern_str = ''.join(map(str, pattern_bits))

    # Decode first 32 bits to ASCII
    decoded_chars = []
    for i in range(0, 32, 8):
        byte = pattern_str[i:i+8]
        decimal = int(byte, 2)
        char = chr(decimal) if 32 <= decimal <= 126 else f'[{decimal}]'
        decoded_chars.append(char)
    decoded = ''.join(decoded_chars)

    # Compare to ground truth
    matches = sum(1 for i in range(min(32, len(expected_binary)))
                  if i < len(pattern_str) and pattern_str[i] == expected_binary[i])
    match_pct = 100 * matches / min(32, len(expected_binary))

    return {
        'binary_sequence': pattern_str,
        'decoded_message': decoded,
        'methylation_values': meth_values,
        'match_percentage': match_pct,
        'matches': matches,
        'total_bits': min(32, len(expected_binary))
    }


def decode_sequences(bed_files, threshold_results, reference, oligo_doc,
                    expected_message='EpiC', output_dir=None, heatmap_dir=None):
    """
    Decode all samples using detected thresholds.

    Args:
        bed_files: List of BED file paths
        threshold_results: Dictionary of threshold results
        reference: Path to reference FASTA
        oligo_doc: Path to oligo sequences DOCX
        expected_message: Expected encoded message
        output_dir: Output directory for results
        heatmap_dir: Output directory for heatmaps

    Returns:
        Dictionary of decoding results per sample
    """
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    if heatmap_dir:
        heatmap_dir = Path(heatmap_dir)
        heatmap_dir.mkdir(parents=True, exist_ok=True)

    # Extract brick positions once
    logger.info("Extracting CpG positions from reference...")
    brick_data = extract_brick_positions(reference, oligo_doc)
    logger.info(f"✓ Extracted {len(brick_data)} CpG positions")

    # Convert expected message to binary
    expected_binary = message_to_binary(expected_message)

    results = {}

    logger.info(f"\nDecoding {len(bed_files)} sample(s)...")

    for bed_file in bed_files:
        bed_path = Path(bed_file)
        sample_name = bed_path.stem

        # Get threshold for this sample
        if sample_name not in threshold_results:
            logger.warning(f"No threshold found for {sample_name}, using default 10%")
            threshold = 10.0
            signal_type = 'UNKNOWN'
            confidence = 'DEFAULT'
        else:
            threshold = threshold_results[sample_name]['threshold']
            signal_type = threshold_results[sample_name]['signal_type']
            confidence = threshold_results[sample_name]['confidence']

        logger.info(f"\n{sample_name}:")
        logger.info(f"  Threshold: {threshold}% ({signal_type})")

        try:
            # Decode
            result = decode_single_sample(
                bed_file,
                threshold,
                brick_data,
                expected_binary
            )

            # Determine status
            match_pct = result['match_percentage']
            if match_pct == 100.0:
                status = "✓ PERFECT"
            elif match_pct >= 95.0:
                status = "★ EXCELLENT"
            elif match_pct >= 90.0:
                status = "◐ GOOD"
            elif match_pct >= 75.0:
                status = "○ PARTIAL"
            else:
                status = "✗ FAILED"

            logger.info(f"  Decoded: '{result['decoded_message']}'")
            logger.info(f"  Match: {match_pct:.1f}% {status}")

            # Store results
            results[sample_name] = {
                'threshold': threshold,
                'signal_type': signal_type,
                'confidence': confidence,
                'binary': result['binary_sequence'],
                'decoded': result['decoded_message'],
                'match_pct': match_pct,
                'status': status
            }

            # Generate heatmap if directory provided
            if heatmap_dir:
                try:
                    generate_heatmap(
                        sample_name,
                        result['methylation_values'],
                        result['binary_sequence'],
                        threshold,
                        expected_binary,
                        expected_message,
                        heatmap_dir
                    )
                except Exception as e:
                    logger.warning(f"Could not generate heatmap: {e}")

        except Exception as e:
            logger.error(f"✗ Failed to decode {sample_name}: {e}")
            results[sample_name] = {
                'threshold': threshold,
                'decoded': 'ERROR',
                'match_pct': 0.0,
                'error': str(e)
            }

    # Save results
    if output_dir:
        output_file = output_dir / 'decoded_results.json'
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"\n✓ Decoding results saved to {output_file}")

    return results


def generate_heatmap(sample_name, meth_values, binary_str, threshold,
                     expected_binary, expected_message, output_dir):
    """Generate heatmap visualization (simplified version)"""
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    from matplotlib.patches import Rectangle

    # Create figure
    fig, ax = plt.subplots(figsize=(20, 4))

    # Create single-row heatmap
    data_matrix = np.array([meth_values[:36]])

    # Create x-axis labels from expected binary (extend to 36 with ?)
    if len(expected_binary) < 36:
        x_labels = list(expected_binary) + ['?'] * (36 - len(expected_binary))
    else:
        x_labels = list(expected_binary[:36])

    # Create heatmap
    sns.heatmap(
        data_matrix,
        annot=True,
        fmt='.1f',
        cmap='YlOrRd',
        vmin=0,
        vmax=100,
        center=50,
        xticklabels=x_labels,
        yticklabels=[sample_name],
        cbar_kws={'label': 'Methylation %'},
        linewidths=0.5,
        linecolor='lightgray',
        ax=ax
    )

    # Highlight detected 1s
    for i in range(min(36, len(binary_str))):
        if binary_str[i] == '1':
            rect = Rectangle((i, 0), 1, 1, fill=False, edgecolor='lime', linewidth=3)
            ax.add_patch(rect)

    # Add character labels
    colors = ['#d62728', '#2ca02c', '#1f77b4', '#9467bd', '#7f7f7f']
    for i, char in enumerate(expected_message[:4]):
        x_pos = i * 8 + 4
        ax.text(x_pos, -0.15, char, fontsize=14, fontweight='bold',
                ha='center', va='bottom', color=colors[i % 5])

    plt.title(f"{sample_name} - Threshold: {threshold}%\n"
             f"Green boxes = detected 1s (methylation >{threshold}%)",
             fontsize=14, fontweight='bold', pad=20)

    ax.set_xlabel("Expected Bit", fontsize=10)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    # Save
    output_file = Path(output_dir) / f"{sample_name}_heatmap.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"  → Heatmap: {output_file.name}")


if __name__ == '__main__':
    # Standalone usage
    import argparse

    parser = argparse.ArgumentParser(description='Decode methylation sequences')
    parser.add_argument('bed_files', nargs='+', help='BED files to decode')
    parser.add_argument('--thresholds', required=True,
                       help='JSON file with threshold results')
    parser.add_argument('--reference', required=True, help='Reference FASTA')
    parser.add_argument('--oligo-doc', required=True, help='Oligo sequences DOCX')
    parser.add_argument('--expected-message', default='EpiC',
                       help='Expected message (default: EpiC)')
    parser.add_argument('--output-dir', help='Output directory')
    parser.add_argument('--heatmap-dir', help='Heatmap output directory')

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                       format='%(asctime)s [%(levelname)s] %(message)s')

    # Load thresholds
    with open(args.thresholds) as f:
        threshold_results = json.load(f)

    decode_sequences(
        args.bed_files,
        threshold_results,
        args.reference,
        args.oligo_doc,
        args.expected_message,
        args.output_dir,
        args.heatmap_dir
    )
