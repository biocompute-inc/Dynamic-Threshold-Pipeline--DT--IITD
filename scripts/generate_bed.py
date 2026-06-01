#!/usr/bin/env python3
"""
Step 2: Generate BED files using modkit pileup

Extracts methylation data from BAM files into bedMethyl format.
"""

import subprocess
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed


logger = logging.getLogger(__name__)


def generate_single_bed(bam_file, reference, output_dir, confidence_threshold=0.7):
    """
    Generate BED file for a single BAM file.

    Args:
        bam_file: Path to BAM file
        reference: Path to reference FASTA
        output_dir: Output directory for BED file
        confidence_threshold: Modkit filter-threshold (0.0-1.0)

    Returns:
        Path to generated BED file
    """
    bam_path = Path(bam_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Output BED file name
    bed_file = output_dir / f"{bam_path.stem}_{confidence_threshold}.bed"

    logger.info(f"Processing {bam_path.name} → {bed_file.name}")

    try:
        # Run modkit pileup
        cmd = [
            'modkit', 'pileup',
            str(bam_path),
            str(bed_file),
            '--ref', str(reference),
            '--cpg',
            '--filter-threshold', str(confidence_threshold)
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )

        if not bed_file.exists():
            raise FileNotFoundError(f"Expected BED file not created: {bed_file}")

        # Get basic stats
        line_count = sum(1 for _ in open(bed_file))
        logger.info(f"  ✓ {bed_file.name} ({line_count} CpG sites)")

        return bed_file

    except subprocess.CalledProcessError as e:
        logger.error(f"  ✗ modkit failed for {bam_path.name}")
        logger.error(f"    Error: {e.stderr}")
        raise
    except Exception as e:
        logger.error(f"  ✗ Error: {e}")
        raise


def generate_bed_files(bam_files, reference, output_dir, confidence_threshold=0.7, threads=4):
    """
    Generate BED files for multiple BAM files in parallel.

    Args:
        bam_files: List of BAM file paths
        reference: Path to reference FASTA
        output_dir: Output directory
        confidence_threshold: Modkit filter-threshold
        threads: Number of parallel threads

    Returns:
        List of generated BED file paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    bed_files = []

    logger.info(f"Generating BED files for {len(bam_files)} BAM file(s)...")
    logger.info(f"  Confidence threshold: {confidence_threshold}")
    logger.info(f"  Parallel threads: {threads}")

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(
                generate_single_bed,
                bam, reference, output_dir, confidence_threshold
            ): bam
            for bam in bam_files
        }

        for future in as_completed(futures):
            try:
                bed_file = future.result()
                bed_files.append(bed_file)
            except Exception as e:
                bam = futures[future]
                logger.error(f"Failed to generate BED for {bam}: {e}")

    logger.info(f"Successfully generated {len(bed_files)}/{len(bam_files)} BED files")

    return bed_files


if __name__ == '__main__':
    # Standalone usage
    import argparse

    parser = argparse.ArgumentParser(description='Generate BED files with modkit')
    parser.add_argument('bam_files', nargs='+', help='BAM files to process')
    parser.add_argument('--reference', required=True, help='Reference FASTA file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--confidence', type=float, default=0.7,
                       help='Modkit filter-threshold (default: 0.7)')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads')

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                       format='%(asctime)s [%(levelname)s] %(message)s')

    generate_bed_files(
        args.bam_files,
        args.reference,
        args.output_dir,
        args.confidence,
        args.threads
    )
