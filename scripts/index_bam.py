#!/usr/bin/env python3
"""
Step 1: Index BAM files using samtools

Creates .bai index files for all input BAM files.
"""

import subprocess
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed


logger = logging.getLogger(__name__)


def index_single_bam(bam_file, output_dir=None):
    """
    Index a single BAM file.

    Args:
        bam_file: Path to BAM file
        output_dir: Optional directory to copy .bai file to

    Returns:
        Path to .bai file
    """
    bam_path = Path(bam_file)
    bai_path = bam_path.parent / f"{bam_path.name}.bai"

    logger.info(f"Indexing {bam_path.name}...")

    try:
        # Run samtools index
        result = subprocess.run(
            ['samtools', 'index', str(bam_path)],
            capture_output=True,
            text=True,
            check=True
        )

        if not bai_path.exists():
            raise FileNotFoundError(f"Expected .bai file not created: {bai_path}")

        # Copy to output directory if specified
        if output_dir:
            output_dir = Path(output_dir)
            output_bai = output_dir / bai_path.name
            import shutil
            shutil.copy2(bai_path, output_bai)
            logger.info(f"  → Copied index to {output_bai}")
            return output_bai

        logger.info(f"  ✓ {bai_path.name}")
        return bai_path

    except subprocess.CalledProcessError as e:
        logger.error(f"  ✗ Failed to index {bam_path.name}")
        logger.error(f"    Error: {e.stderr}")
        raise
    except Exception as e:
        logger.error(f"  ✗ Error: {e}")
        raise


def index_bam_files(bam_files, output_dir=None, threads=4):
    """
    Index multiple BAM files in parallel.

    Args:
        bam_files: List of BAM file paths
        output_dir: Optional directory to save .bai files
        threads: Number of parallel threads

    Returns:
        List of .bai file paths
    """
    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    indexed_files = []

    logger.info(f"Indexing {len(bam_files)} BAM file(s) using {threads} threads...")

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(index_single_bam, bam, output_dir): bam
            for bam in bam_files
        }

        for future in as_completed(futures):
            try:
                bai_file = future.result()
                indexed_files.append(bai_file)
            except Exception as e:
                bam = futures[future]
                logger.error(f"Failed to index {bam}: {e}")

    logger.info(f"Successfully indexed {len(indexed_files)}/{len(bam_files)} files")

    return indexed_files


if __name__ == '__main__':
    # Standalone usage
    import argparse

    parser = argparse.ArgumentParser(description='Index BAM files with samtools')
    parser.add_argument('bam_files', nargs='+', help='BAM files to index')
    parser.add_argument('--output-dir', help='Directory to save .bai files')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads')

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                       format='%(asctime)s [%(levelname)s] %(message)s')

    index_bam_files(args.bam_files, args.output_dir, args.threads)
