#!/usr/bin/env python3
"""
Methylation-Based DNA Data Decoder - Master Pipeline

Complete automated pipeline from BAM files to decoded ASCII messages.

Usage:
    python master_pipeline.py --bam-dir bam_files/ --reference ref.fasta \
        --oligo-doc oligos.docx --output-dir results/
"""

import os
import sys
import argparse
import json
import logging
from pathlib import Path
from datetime import datetime
import subprocess

# Add lib directory to path
sys.path.insert(0, str(Path(__file__).parent / 'lib'))
sys.path.insert(0, str(Path(__file__).parent / 'scripts'))


class MethylationPipeline:
    """Master pipeline orchestrator"""

    def __init__(self, args):
        self.args = args
        self.setup_logging()
        self.setup_output_dirs()
        self.bam_files = self.find_bam_files()

    def setup_logging(self):
        """Setup logging configuration"""
        log_file = Path(self.args.output_dir) / 'pipeline.log'
        log_file.parent.mkdir(parents=True, exist_ok=True)

        level = logging.DEBUG if self.args.verbose else \
                logging.WARNING if self.args.quiet else logging.INFO

        logging.basicConfig(
            level=level,
            format='%(asctime)s [%(levelname)s] %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)

    def setup_output_dirs(self):
        """Create output directory structure"""
        base = Path(self.args.output_dir)
        self.dirs = {
            'base': base,
            'indices': base / 'bam_indices',
            'beds': base / 'bed_files',
            'thresholds': base / 'thresholds',
            'decoded': base / 'decoded',
            'heatmaps': base / 'decoded' / 'sample_heatmaps'
        }

        for dir_path in self.dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"Output directory: {base}")

    def find_bam_files(self):
        """Find BAM files to process"""
        if self.args.bam:
            bam_files = [Path(self.args.bam)]
        elif self.args.bam_dir:
            bam_dir = Path(self.args.bam_dir)
            bam_files = list(bam_dir.glob('*.bam'))
        else:
            self.logger.error("Must provide --bam or --bam-dir")
            sys.exit(1)

        if not bam_files:
            self.logger.error(f"No BAM files found in {self.args.bam_dir}")
            sys.exit(1)

        self.logger.info(f"Found {len(bam_files)} BAM file(s)")
        for bam in bam_files:
            self.logger.info(f"  - {bam.name}")

        return bam_files

    def check_dependencies(self):
        """Check required tools are available"""
        self.logger.info("Checking dependencies...")

        tools = {
            'samtools': 'conda install -c bioconda samtools',
            'modkit': 'https://github.com/nanoporetech/modkit/releases'
        }

        missing = []
        for tool, install_cmd in tools.items():
            if subprocess.run(['which', tool], capture_output=True).returncode != 0:
                missing.append((tool, install_cmd))

        if missing:
            self.logger.error("Missing dependencies:")
            for tool, cmd in missing:
                self.logger.error(f"  - {tool}: {cmd}")
            sys.exit(1)

        # Check Python packages
        try:
            import numpy, scipy, docx, matplotlib, seaborn, pandas
        except ImportError as e:
            self.logger.error(f"Missing Python package: {e}")
            self.logger.error("Run: pip install -r requirements.txt")
            sys.exit(1)

        self.logger.info("✓ All dependencies found")

    def step1_index_bams(self):
        """Step 1: Index BAM files"""
        if self.args.skip_indexing:
            self.logger.info("STEP 1: Skipping BAM indexing (--skip-indexing)")
            return

        self.logger.info("="*80)
        self.logger.info("STEP 1: INDEXING BAM FILES")
        self.logger.info("="*80)

        from index_bam import index_bam_files

        indexed = index_bam_files(
            self.bam_files,
            output_dir=self.dirs['indices'],
            threads=self.args.threads
        )

        self.logger.info(f"✓ Indexed {len(indexed)} BAM files")

    def step2_generate_beds(self):
        """Step 2: Generate BED files with modkit"""
        if self.args.skip_modkit:
            if self.args.bed_dir:
                self.logger.info("STEP 2: Using existing BED files from --bed-dir")
                self.bed_files = list(Path(self.args.bed_dir).glob('*.bed'))
            else:
                self.logger.error("--skip-modkit requires --bed-dir")
                sys.exit(1)
            return

        self.logger.info("="*80)
        self.logger.info("STEP 2: GENERATING BED FILES (modkit pileup)")
        self.logger.info("="*80)

        from generate_bed import generate_bed_files

        self.bed_files = generate_bed_files(
            bam_files=self.bam_files,
            reference=self.args.reference,
            output_dir=self.dirs['beds'],
            confidence_threshold=self.args.modkit_confidence,
            threads=self.args.threads
        )

        self.logger.info(f"✓ Generated {len(self.bed_files)} BED files")

    def step3_detect_thresholds(self):
        """Step 3: Detect optimal thresholds"""
        self.logger.info("="*80)
        self.logger.info("STEP 3: DETECTING OPTIMAL THRESHOLDS")
        self.logger.info("="*80)

        from detect_threshold import detect_thresholds

        self.threshold_results = detect_thresholds(
            bed_files=self.bed_files,
            reference=self.args.reference,
            oligo_doc=self.args.oligo_doc,
            expected_message=self.args.expected_message,
            conservative_bias=self.args.conservative_bias,
            output_dir=self.dirs['thresholds'],
            verbose=not self.args.quiet
        )

        self.logger.info(f"✓ Detected thresholds for {len(self.threshold_results)} samples")

    def step4_decode_sequences(self):
        """Step 4: Decode binary sequences"""
        self.logger.info("="*80)
        self.logger.info("STEP 4: DECODING BINARY SEQUENCES")
        self.logger.info("="*80)

        from decode_sequence import decode_sequences

        self.decode_results = decode_sequences(
            bed_files=self.bed_files,
            threshold_results=self.threshold_results,
            reference=self.args.reference,
            oligo_doc=self.args.oligo_doc,
            expected_message=self.args.expected_message,
            output_dir=self.dirs['decoded'],
            heatmap_dir=self.dirs['heatmaps']
        )

        self.logger.info(f"✓ Decoded {len(self.decode_results)} samples")

    def step5_generate_report(self):
        """Step 5: Generate final report"""
        self.logger.info("="*80)
        self.logger.info("STEP 5: GENERATING FINAL REPORT")
        self.logger.info("="*80)

        self.generate_summary_report()
        self.logger.info(f"✓ Report saved to {self.dirs['decoded']}")

    def generate_summary_report(self):
        """Generate summary report"""
        summary_file = self.dirs['decoded'] / 'summary.txt'

        with open(summary_file, 'w') as f:
            f.write("="*100 + "\n")
            f.write("  METHYLATION-BASED DNA DATA DECODER - PIPELINE RESULTS\n")
            f.write("="*100 + "\n\n")

            f.write(f"Run Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Samples Processed: {len(self.decode_results)}\n")
            f.write(f"Expected Message: '{self.args.expected_message}'\n")
            f.write(f"Conservative Bias: {self.args.conservative_bias}\n")
            f.write(f"Modkit Confidence: {self.args.modkit_confidence}\n\n")

            f.write("="*100 + "\n")
            f.write(f"{'Sample':<30} {'Threshold':<12} {'Signal':<10} {'Decoded':<12} {'Match%':<10} {'Status'}\n")
            f.write("-"*100 + "\n")

            perfect = 0
            for sample_name, result in self.decode_results.items():
                if result['match_pct'] == 100.0:
                    perfect += 1

                status_symbol = "✓" if result['match_pct'] == 100 else \
                               "★" if result['match_pct'] >= 95 else \
                               "◐" if result['match_pct'] >= 90 else \
                               "○" if result['match_pct'] >= 75 else "✗"

                f.write(f"{sample_name:<30} {result['threshold']:>10.1f}% "
                       f"{result['signal_type']:<10} '{result['decoded']}'      "
                       f"{result['match_pct']:>7.1f}% {status_symbol}\n")

            f.write("="*100 + "\n\n")

            f.write(f"SUCCESS RATE: {perfect}/{len(self.decode_results)} "
                   f"({100*perfect/len(self.decode_results):.1f}%) perfect matches\n\n")

            f.write("Output Files:\n")
            f.write(f"  - Detailed results: {self.dirs['decoded'] / 'decoded_results.json'}\n")
            f.write(f"  - Heatmaps: {self.dirs['heatmaps']}/\n")
            f.write(f"  - Threshold data: {self.dirs['thresholds'] / 'threshold_results.json'}\n")
            f.write(f"  - Pipeline log: {self.dirs['base'] / 'pipeline.log'}\n")

        self.logger.info(f"\n{open(summary_file).read()}")

    def run(self):
        """Run complete pipeline"""
        start_time = datetime.now()

        self.logger.info("="*80)
        self.logger.info("  METHYLATION-BASED DNA DATA DECODER")
        self.logger.info("="*80)
        self.logger.info(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        self.logger.info("")

        try:
            self.check_dependencies()
            self.step1_index_bams()
            self.step2_generate_beds()
            self.step3_detect_thresholds()
            self.step4_decode_sequences()
            self.step5_generate_report()

            elapsed = datetime.now() - start_time
            self.logger.info("")
            self.logger.info("="*80)
            self.logger.info(f"  PIPELINE COMPLETED SUCCESSFULLY")
            self.logger.info("="*80)
            self.logger.info(f"Total time: {elapsed}")
            self.logger.info(f"Results: {self.dirs['decoded']}")

        except Exception as e:
            self.logger.error(f"Pipeline failed: {e}", exc_info=True)
            sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Methylation-Based DNA Data Decoder Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Process all BAM files in a directory
    python master_pipeline.py --bam-dir bam_files/ --reference ref.fasta \\
        --oligo-doc oligos.docx --output-dir results/

    # Process single BAM file with custom bias
    python master_pipeline.py --bam sample.bam --reference ref.fasta \\
        --oligo-doc oligos.docx --output-dir results/ --conservative-bias 0.9

    # Resume from existing BED files
    python master_pipeline.py --skip-modkit --bed-dir beds/ \\
        --reference ref.fasta --oligo-doc oligos.docx --output-dir results/
        """
    )

    # Input files
    input_group = parser.add_argument_group('Input Files')
    input_group.add_argument('--bam', help='Single BAM file to process')
    input_group.add_argument('--bam-dir', help='Directory containing BAM files')
    input_group.add_argument('--reference', required=True, help='Reference FASTA file')
    input_group.add_argument('--oligo-doc', required=True,
                           help='Methylated oligo sequences DOCX file')

    # Output
    output_group = parser.add_argument_group('Output')
    output_group.add_argument('--output-dir', required=True,
                            help='Output directory for all results')

    # Parameters
    param_group = parser.add_argument_group('Pipeline Parameters')
    param_group.add_argument('--conservative-bias', type=float, default=0.8,
                           help='Conservative bias for threshold detection (0.0-1.0, default: 0.8)')
    param_group.add_argument('--modkit-confidence', type=float, default=0.7,
                           help='Modkit filter-threshold (default: 0.7)')
    param_group.add_argument('--expected-message', default='EpiC',
                           help='Expected encoded message (default: "EpiC")')
    param_group.add_argument('--threads', type=int, default=4,
                           help='Number of parallel threads (default: 4)')

    # Control flow
    control_group = parser.add_argument_group('Pipeline Control')
    control_group.add_argument('--skip-indexing', action='store_true',
                             help='Skip BAM indexing step')
    control_group.add_argument('--skip-modkit', action='store_true',
                             help='Skip modkit pileup step')
    control_group.add_argument('--bed-dir', help='Use existing BED files from directory')

    # Verbosity
    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('--verbose', action='store_true',
                               help='Verbose output')
    verbosity_group.add_argument('--quiet', action='store_true',
                               help='Minimal output')

    args = parser.parse_args()

    # Validation
    if not args.bam and not args.bam_dir:
        parser.error("Must provide either --bam or --bam-dir")

    if args.skip_modkit and not args.bed_dir:
        parser.error("--skip-modkit requires --bed-dir")

    if args.conservative_bias < 0 or args.conservative_bias > 1:
        parser.error("--conservative-bias must be between 0.0 and 1.0")

    # Run pipeline
    pipeline = MethylationPipeline(args)
    pipeline.run()


if __name__ == '__main__':
    main()
