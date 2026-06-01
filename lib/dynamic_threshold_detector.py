#!/usr/bin/env python3
"""
Dynamic Threshold Detector for Methylation Signal Analysis

This module automatically detects optimal thresholds for methylation-based binary decoding
by analyzing the signal separation between expected 0s and 1s.

Author: BioCompute x IIT-Delhi
Date: May 16, 2026
"""

import numpy as np
from scipy import stats
from docx import Document
import argparse
import json
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, asdict


@dataclass
class SignalStatistics:
    """Statistics for signal separation analysis"""
    mean_zeros: float
    mean_ones: float
    median_zeros: float
    median_ones: float
    std_zeros: float
    std_ones: float
    n_zeros: int
    n_ones: int
    fold_change_mean: float
    fold_change_median: float
    separation_score: float
    is_significant: bool
    is_practically_viable: bool
    viability_reason: str


@dataclass
class ThresholdResult:
    """Result of dynamic threshold detection"""
    threshold: float
    confidence: str  # 'HIGH', 'MEDIUM', 'LOW', 'STATISTICALLY_SIGNIFICANT_BUT_WEAK'
    signal_type: str  # 'STRONG', 'WEAK', 'FAILED', 'TOO_WEAK'
    statistics: SignalStatistics
    recommendation: str
    practical_warning: str


class DynamicThresholdDetector:
    """
    Automatically detect optimal methylation threshold based on signal separation.

    The detector analyzes the distribution of methylation values for expected 0s and 1s,
    calculates separation metrics, and determines an appropriate threshold.
    """

    # Threshold detection parameters
    MIN_FOLD_CHANGE = 3.0  # Minimum fold-change to consider signals separable
    HIGH_FOLD_CHANGE = 5.0  # Fold-change for high confidence
    MIN_SEPARATION_SCORE = 2.0  # Minimum Cohen's d for significance

    # Practical viability thresholds
    MIN_PRACTICAL_SIGNAL = 10.0  # Minimum mean methylation for "1s" to be practically viable (%)
    MIN_ABSOLUTE_DIFFERENCE = 5.0  # Minimum absolute difference between 1s and 0s (%)

    def __init__(self, verbose: bool = True, conservative_bias: float = 0.0):
        """
        Initialize the detector.

        Args:
            verbose: If True, print detailed analysis information
            conservative_bias: Bias threshold towards 1s to avoid false positives on 0s.
                             0.0 = midpoint (default)
                             0.5 = 75% towards 1s, 25% towards 0s
                             0.8 = 90% towards 1s, 10% towards 0s
                             1.0 = use lower bound of 1s distribution
        """
        self.verbose = verbose
        self.conservative_bias = conservative_bias
        self.brick_data = None

    def reverse_complement(self, seq: str) -> str:
        """Return reverse complement of DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement.get(base, base) for base in reversed(seq))

    def extract_brick_positions(self) -> List[Dict]:
        """Extract CpG positions from methylated oligo sequences"""
        with open('reference.fasta') as f:
            ref_seq = ''.join(line.strip() for line in f if not line.startswith('>'))

        doc = Document("methylated oligo sequences.docx")
        sequences_with_markers = []

        for para in doc.paragraphs:
            text = para.text.strip()
            if '/i5MedC/' in text:
                parts = text.split()
                for part in parts:
                    if '/i5MedC/' in part:
                        cleaned = part
                        marker_end = cleaned.find('/i5MedC/') + len('/i5MedC/')
                        after_marker = cleaned[marker_end:]
                        if '/' in after_marker:
                            extra_slash = after_marker.find('/')
                            cleaned = cleaned[:marker_end + extra_slash]
                        sequences_with_markers.append(cleaned)

        sequences_with_markers = list(dict.fromkeys(sequences_with_markers))

        brick_data = []
        for idx, seq_with_marker in enumerate(sequences_with_markers[:36]):
            parts = seq_with_marker.split('/i5MedC/')
            if len(parts) != 2:
                continue
            before_marker = parts[0]
            after_marker = parts[1]
            full_oligo = before_marker + 'C' + after_marker
            rc_oligo = self.reverse_complement(full_oligo)
            match_pos = ref_seq.find(rc_oligo)
            if match_pos == -1:
                continue
            meth_c_in_oligo = len(before_marker)
            meth_c_in_rc = len(full_oligo) - 1 - meth_c_in_oligo
            meth_g_in_ref = match_pos + meth_c_in_rc
            if meth_g_in_ref > 0 and ref_seq[meth_g_in_ref-1] == 'C' and ref_seq[meth_g_in_ref] == 'G':
                cpg_pos = meth_g_in_ref - 1
            else:
                cpg_pos = meth_g_in_ref
            brick_data.append({'brick': idx + 1, 'cpg_pos': cpg_pos})

        return brick_data

    def load_bed_file(self, bed_file: str) -> Dict[int, float]:
        """
        Load methylation data from BED file.

        Args:
            bed_file: Path to bedMethyl file

        Returns:
            Dictionary mapping position to methylation percentage
        """
        meth_data = {}
        with open(bed_file) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 10 and fields[3] == 'm':
                    pos = int(fields[1])
                    pct_meth = float(fields[9].split()[1])
                    meth_data[pos] = pct_meth
        return meth_data

    def calculate_signal_statistics(
        self,
        zeros_meth: List[float],
        ones_meth: List[float]
    ) -> SignalStatistics:
        """
        Calculate comprehensive statistics for signal separation.

        Args:
            zeros_meth: Methylation values for expected 0 positions
            ones_meth: Methylation values for expected 1 positions

        Returns:
            SignalStatistics object with all metrics
        """
        # Basic statistics
        mean_zeros = np.mean(zeros_meth)
        mean_ones = np.mean(ones_meth)
        median_zeros = np.median(zeros_meth)
        median_ones = np.median(ones_meth)
        std_zeros = np.std(zeros_meth)
        std_ones = np.std(ones_meth)

        # Fold changes
        fold_change_mean = mean_ones / mean_zeros if mean_zeros > 0 else np.inf
        fold_change_median = median_ones / median_zeros if median_zeros > 0 else np.inf

        # Cohen's d - effect size for separation
        pooled_std = np.sqrt((std_zeros**2 + std_ones**2) / 2)
        cohens_d = (mean_ones - mean_zeros) / pooled_std if pooled_std > 0 else 0

        # Determine if separation is significant
        is_significant = (
            fold_change_mean >= self.MIN_FOLD_CHANGE and
            cohens_d >= self.MIN_SEPARATION_SCORE
        )

        # Determine if signal is practically viable
        absolute_difference = mean_ones - mean_zeros
        is_practically_viable = (
            mean_ones >= self.MIN_PRACTICAL_SIGNAL and
            absolute_difference >= self.MIN_ABSOLUTE_DIFFERENCE
        )

        # Generate viability reason
        if is_practically_viable:
            viability_reason = "Signal strength adequate for practical use"
        else:
            reasons = []
            if mean_ones < self.MIN_PRACTICAL_SIGNAL:
                reasons.append(f"'1' signal too weak ({mean_ones:.1f}% < {self.MIN_PRACTICAL_SIGNAL}%)")
            if absolute_difference < self.MIN_ABSOLUTE_DIFFERENCE:
                reasons.append(f"Absolute separation too small ({absolute_difference:.1f}% < {self.MIN_ABSOLUTE_DIFFERENCE}%)")
            viability_reason = "; ".join(reasons)

        return SignalStatistics(
            mean_zeros=mean_zeros,
            mean_ones=mean_ones,
            median_zeros=median_zeros,
            median_ones=median_ones,
            std_zeros=std_zeros,
            std_ones=std_ones,
            n_zeros=len(zeros_meth),
            n_ones=len(ones_meth),
            fold_change_mean=fold_change_mean,
            fold_change_median=fold_change_median,
            separation_score=cohens_d,
            is_significant=is_significant,
            is_practically_viable=is_practically_viable,
            viability_reason=viability_reason
        )

    def determine_optimal_threshold(
        self,
        zeros_meth: List[float],
        ones_meth: List[float],
        statistics: SignalStatistics
    ) -> Tuple[float, str, str]:
        """
        Determine optimal threshold based on signal characteristics.

        Args:
            zeros_meth: Methylation values for expected 0s
            ones_meth: Methylation values for expected 1s
            statistics: Pre-calculated signal statistics

        Returns:
            Tuple of (threshold, confidence_level, signal_type)
        """
        # Check if signals are separable
        if not statistics.is_significant:
            return 10.0, 'LOW', 'FAILED'

        # Calculate base threshold using multiple methods

        # Method 1: Midpoint between means
        threshold_mean = (statistics.mean_zeros + statistics.mean_ones) / 2

        # Method 2: Midpoint between medians
        threshold_median = (statistics.median_zeros + statistics.median_ones) / 2

        # Method 3: Conservative approaches to avoid false positives
        # 3a. Upper bound of 0s distribution (mean + 2*std)
        threshold_zeros_upper = statistics.mean_zeros + 2 * statistics.std_zeros

        # 3b. Lower bound of 1s distribution (mean - 1*std)
        threshold_ones_lower = max(5.0, statistics.mean_ones - 1 * statistics.std_ones)

        # 3c. Maximum observed 0 + safety margin
        max_zero = max(zeros_meth)
        min_one = min(ones_meth)
        threshold_max_zero_margin = max_zero + (min_one - max_zero) * 0.3  # 30% above max zero

        # 3d. 95th percentile of 0s
        threshold_zeros_p95 = np.percentile(zeros_meth, 95)

        # Method 4: Maximum separation with specificity bias
        all_values = list(zeros_meth) + list(ones_meth)
        all_values_sorted = sorted(all_values)

        best_threshold = threshold_mean
        best_score = 0

        # Weight specificity higher to avoid false positives
        specificity_weight = 2.0 if self.conservative_bias > 0 else 1.0
        sensitivity_weight = 1.0

        # Test candidate thresholds
        for candidate in all_values_sorted:
            if statistics.mean_zeros < candidate < statistics.mean_ones:
                tp = sum(1 for v in ones_meth if v > candidate)  # True positives
                tn = sum(1 for v in zeros_meth if v <= candidate)  # True negatives
                fp = sum(1 for v in zeros_meth if v > candidate)  # False positives
                fn = sum(1 for v in ones_meth if v <= candidate)  # False negatives

                # Calculate specificity and sensitivity
                specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
                sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0

                # Weighted score (prioritize specificity if conservative)
                score = (specificity_weight * specificity + sensitivity_weight * sensitivity) / (specificity_weight + sensitivity_weight)

                if score > best_score:
                    best_score = score
                    best_threshold = candidate

        # Apply conservative bias to shift threshold towards 1s
        if self.conservative_bias > 0:
            # Blend between balanced threshold and conservative threshold
            conservative_threshold = max(
                threshold_zeros_upper,
                threshold_zeros_p95,
                threshold_max_zero_margin
            )

            # Weighted average based on bias
            # bias=0.0 → 100% best_threshold
            # bias=0.5 → 50% best_threshold, 50% conservative
            # bias=1.0 → 100% conservative (or ones_lower bound)
            if self.conservative_bias >= 1.0:
                final_threshold = min(conservative_threshold, threshold_ones_lower)
            else:
                final_threshold = (1 - self.conservative_bias) * best_threshold + \
                                self.conservative_bias * conservative_threshold
        else:
            final_threshold = best_threshold

        # Check practical viability first
        if not statistics.is_practically_viable:
            # Statistically significant but practically weak
            confidence = 'STATISTICALLY_SIGNIFICANT_BUT_WEAK'
            signal_type = 'TOO_WEAK'
            threshold = round(final_threshold * 2) / 2
            threshold = max(5.0, min(threshold, 95.0))
            return threshold, confidence, signal_type

        # Determine confidence and signal type for practically viable signals
        if statistics.fold_change_mean >= self.HIGH_FOLD_CHANGE:
            confidence = 'HIGH'
            signal_type = 'STRONG' if statistics.mean_ones > 50 else 'WEAK'
        elif statistics.fold_change_mean >= self.MIN_FOLD_CHANGE:
            confidence = 'MEDIUM'
            signal_type = 'WEAK'
        else:
            confidence = 'LOW'
            signal_type = 'WEAK'

        # Add bias indicator to confidence if conservative
        if self.conservative_bias > 0:
            confidence = f"{confidence}_CONSERVATIVE"

        # Round threshold to nearest 0.5%
        threshold = round(final_threshold * 2) / 2

        # Apply constraints
        threshold = max(5.0, min(threshold, 95.0))

        return threshold, confidence, signal_type

    def detect_threshold(
        self,
        bed_file: str,
        expected_binary: str = "01000101011100000110100101000011"
    ) -> ThresholdResult:
        """
        Main method: Detect optimal threshold for a BED file.

        Args:
            bed_file: Path to bedMethyl file
            expected_binary: Expected binary pattern (default: "EpiC")

        Returns:
            ThresholdResult object with threshold and statistics
        """
        if self.verbose:
            print(f"\n{'='*100}")
            print(f"  DYNAMIC THRESHOLD DETECTION: {bed_file}")
            print(f"{'='*100}\n")

        # Load data
        if self.brick_data is None:
            self.brick_data = self.extract_brick_positions()
            if self.verbose:
                print(f"✓ Extracted {len(self.brick_data)} CpG positions")

        meth_data = self.load_bed_file(bed_file)
        if self.verbose:
            print(f"✓ Loaded methylation data from BED file")

        # Separate 0s and 1s
        zeros_meth = []
        ones_meth = []

        n_positions = min(len(expected_binary), len(self.brick_data))

        for i in range(n_positions):
            brick = self.brick_data[i]
            pos = brick['cpg_pos']

            if pos in meth_data:
                pct = meth_data[pos]
                if expected_binary[i] == '0':
                    zeros_meth.append(pct)
                elif expected_binary[i] == '1':
                    ones_meth.append(pct)

        if len(zeros_meth) == 0 or len(ones_meth) == 0:
            if self.verbose:
                print("✗ ERROR: Insufficient data for threshold detection")
            return ThresholdResult(
                threshold=10.0,
                confidence='LOW',
                signal_type='FAILED',
                statistics=None,
                recommendation="Insufficient data - using default threshold of 10%",
                practical_warning=None
            )

        # Calculate statistics
        statistics = self.calculate_signal_statistics(zeros_meth, ones_meth)

        # Determine threshold
        threshold, confidence, signal_type = self.determine_optimal_threshold(
            zeros_meth, ones_meth, statistics
        )

        # Generate recommendation and practical warning
        practical_warning = ""

        if not statistics.is_practically_viable:
            recommendation = (
                f"⚠ CAUTION: While statistically significant ({statistics.fold_change_mean:.1f}x fold-change, "
                f"Cohen's d={statistics.separation_score:.2f}), the absolute signal levels are too low for "
                f"practical use. {statistics.viability_reason}. "
                f"This sample may represent background noise rather than real methylation encoding."
            )
            practical_warning = (
                f"⚠ PRACTICAL VIABILITY WARNING: Mean '1' signal = {statistics.mean_ones:.2f}%, "
                f"Mean '0' signal = {statistics.mean_zeros:.2f}%. "
                f"Absolute difference = {statistics.mean_ones - statistics.mean_zeros:.2f}%. "
                f"Signal may be indistinguishable from technical noise. "
                f"Recommend treating this sample as FAILED for practical applications."
            )
        elif statistics.fold_change_mean >= self.HIGH_FOLD_CHANGE:
            recommendation = (
                f"Excellent signal separation detected ({statistics.fold_change_mean:.1f}x fold-change). "
                f"Dynamic threshold of {threshold}% provides optimal discrimination. "
                f"Signal is practically viable (Mean 1s: {statistics.mean_ones:.1f}%)."
            )
        elif statistics.is_significant:
            recommendation = (
                f"Moderate signal separation detected ({statistics.fold_change_mean:.1f}x fold-change). "
                f"Dynamic threshold of {threshold}% recommended with {confidence} confidence. "
                f"Signal is practically viable (Mean 1s: {statistics.mean_ones:.1f}%)."
            )
        else:
            recommendation = (
                f"Poor signal separation ({statistics.fold_change_mean:.1f}x fold-change). "
                f"Decoding may be unreliable. Default threshold of {threshold}% applied."
            )

        result = ThresholdResult(
            threshold=threshold,
            confidence=confidence,
            signal_type=signal_type,
            statistics=statistics,
            recommendation=recommendation,
            practical_warning=practical_warning
        )

        if self.verbose:
            self._print_results(result)

        return result

    def _print_results(self, result: ThresholdResult):
        """Print detailed results"""
        stats = result.statistics

        print(f"{'='*100}")
        print(f"  SIGNAL ANALYSIS RESULTS")
        print(f"{'='*100}")
        print(f"Expected 0s:")
        print(f"  Mean:   {stats.mean_zeros:>8.2f}%  ±{stats.std_zeros:.2f}%")
        print(f"  Median: {stats.median_zeros:>8.2f}%")
        print(f"  Count:  {stats.n_zeros:>8}")
        print()
        print(f"Expected 1s:")
        print(f"  Mean:   {stats.mean_ones:>8.2f}%  ±{stats.std_ones:.2f}%")
        print(f"  Median: {stats.median_ones:>8.2f}%")
        print(f"  Count:  {stats.n_ones:>8}")
        print()
        print(f"Separation Metrics:")
        print(f"  Fold-change (mean):   {stats.fold_change_mean:>6.2f}x")
        print(f"  Fold-change (median): {stats.fold_change_median:>6.2f}x")
        print(f"  Effect size (Cohen's d): {stats.separation_score:>6.2f}")
        print(f"  Statistical significance: {'✓ YES' if stats.is_significant else '✗ NO'}")
        print(f"  Practical viability: {'✓ YES' if stats.is_practically_viable else '⚠ NO'}")
        if not stats.is_practically_viable:
            print(f"    Reason: {stats.viability_reason}")
        print()
        print(f"{'='*100}")
        print(f"  THRESHOLD DETERMINATION")
        print(f"{'='*100}")
        print(f"Optimal Threshold: {result.threshold}%")
        print(f"Confidence Level:  {result.confidence}")
        print(f"Signal Type:       {result.signal_type}")
        print()
        print(f"Recommendation:")
        print(f"  {result.recommendation}")
        if result.practical_warning:
            print()
            print(f"{'!'*100}")
            print(f"  PRACTICAL WARNING")
            print(f"{'!'*100}")
            print(f"  {result.practical_warning}")
            print(f"{'!'*100}")
        print(f"{'='*100}\n")

    def batch_process(
        self,
        bed_files: List[str],
        output_json: Optional[str] = None
    ) -> Dict[str, ThresholdResult]:
        """
        Process multiple BED files and detect thresholds for each.

        Args:
            bed_files: List of BED file paths
            output_json: Optional path to save results as JSON

        Returns:
            Dictionary mapping file names to ThresholdResult objects
        """
        results = {}

        for bed_file in bed_files:
            sample_name = bed_file.split('/')[-1].replace('.bed', '')
            result = self.detect_threshold(bed_file)
            results[sample_name] = result

        if output_json:
            # Convert to JSON-serializable format
            json_data = {}
            for name, r in results.items():
                stats_dict = None
                if r.statistics:
                    stats_dict = asdict(r.statistics)
                    # Convert bool to int for JSON serialization
                    stats_dict['is_significant'] = int(stats_dict['is_significant'])
                    stats_dict['is_practically_viable'] = int(stats_dict['is_practically_viable'])

                json_data[name] = {
                    'threshold': r.threshold,
                    'confidence': r.confidence,
                    'signal_type': r.signal_type,
                    'statistics': stats_dict,
                    'recommendation': r.recommendation,
                    'practical_warning': r.practical_warning
                }

            with open(output_json, 'w') as f:
                json.dump(json_data, f, indent=2)

            print(f"✓ Results saved to {output_json}")

        return results


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(
        description='Dynamic Threshold Detector for Methylation Signal Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Detect threshold for a single file
    python dynamic_threshold_detector.py --bed sample.bed

    # Detect threshold with custom expected pattern
    python dynamic_threshold_detector.py --bed sample.bed --pattern "01001000"

    # Batch process multiple files
    python dynamic_threshold_detector.py --batch epic_analysis/*.bed --output results.json

    # Quiet mode
    python dynamic_threshold_detector.py --bed sample.bed --quiet
        """
    )

    parser.add_argument('--bed', help='Path to single BED file')
    parser.add_argument('--batch', nargs='+', help='Process multiple BED files')
    parser.add_argument('--pattern', default="01000101011100000110100101000011",
                       help='Expected binary pattern (default: "EpiC")')
    parser.add_argument('--output', help='Output JSON file for batch results')
    parser.add_argument('--quiet', action='store_true', help='Suppress detailed output')
    parser.add_argument('--conservative', type=float, default=0.0, metavar='BIAS',
                       help='Conservative bias (0.0-1.0). 0.0=balanced (default), '
                            '0.5=moderate bias towards 1s, 1.0=maximum specificity. '
                            'Higher values avoid false positives on 0s.')

    args = parser.parse_args()

    detector = DynamicThresholdDetector(verbose=not args.quiet, conservative_bias=args.conservative)

    if args.bed:
        result = detector.detect_threshold(args.bed, args.pattern)
        print(f"\n✓ Detected threshold: {result.threshold}% ({result.confidence} confidence)")

    elif args.batch:
        results = detector.batch_process(args.batch, args.output)

        print(f"\n{'='*100}")
        print(f"  BATCH PROCESSING SUMMARY")
        print(f"{'='*100}")
        print(f"{'Sample':<30} {'Threshold':<12} {'Confidence':<38} {'Signal Type':<12} {'Fold-Change':<12} {'Viable':<8}")
        print(f"{'-'*130}")

        for name, result in results.items():
            fold = result.statistics.fold_change_mean if result.statistics else 0
            viable = '✓' if (result.statistics and result.statistics.is_practically_viable) else '⚠'
            conf_display = result.confidence if len(result.confidence) < 38 else result.confidence[:35] + '...'
            print(f"{name:<30} {result.threshold:>10.1f}% {conf_display:<38} {result.signal_type:<12} {fold:>10.2f}x {viable:>6}")

        print(f"{'='*100}\n")

    else:
        parser.print_help()


if __name__ == '__main__':
    main()
