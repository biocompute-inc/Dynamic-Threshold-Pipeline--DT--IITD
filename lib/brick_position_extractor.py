#!/usr/bin/env python3
"""
Utility functions for extracting CpG positions and loading methylation data.

These utilities are shared across the pipeline scripts.
"""

from pathlib import Path
from docx import Document


def reverse_complement(seq):
    """Return reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


def extract_brick_positions(reference_file, oligo_doc_file):
    """
    Extract CpG positions from reference and methylated oligo sequences.

    Args:
        reference_file: Path to reference FASTA file
        oligo_doc_file: Path to methylated oligo sequences DOCX file

    Returns:
        List of dictionaries with brick number and CpG position
    """
    # Read reference sequence
    with open(reference_file) as f:
        ref_seq = ''.join(line.strip() for line in f if not line.startswith('>'))

    # Parse DOCX to get sequences with methylation markers
    doc = Document(oligo_doc_file)
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

    # Remove duplicates while preserving order
    sequences_with_markers = list(dict.fromkeys(sequences_with_markers))

    # Extract CpG positions for first 36 sequences
    brick_data = []
    for idx, seq_with_marker in enumerate(sequences_with_markers[:36]):
        parts = seq_with_marker.split('/i5MedC/')
        if len(parts) != 2:
            continue

        before_marker = parts[0]
        after_marker = parts[1]
        full_oligo = before_marker + 'C' + after_marker

        # Find position in reference (using reverse complement)
        rc_oligo = reverse_complement(full_oligo)
        match_pos = ref_seq.find(rc_oligo)

        if match_pos == -1:
            continue

        # Calculate CpG position
        meth_c_in_oligo = len(before_marker)
        meth_c_in_rc = len(full_oligo) - 1 - meth_c_in_oligo
        meth_g_in_ref = match_pos + meth_c_in_rc

        # Check if it's a CpG dinucleotide
        if meth_g_in_ref > 0 and ref_seq[meth_g_in_ref-1] == 'C' and ref_seq[meth_g_in_ref] == 'G':
            cpg_pos = meth_g_in_ref - 1
        else:
            cpg_pos = meth_g_in_ref

        brick_data.append({
            'brick': idx + 1,
            'cpg_pos': cpg_pos
        })

    return brick_data


def load_bed_methylation(bed_file):
    """
    Load methylation data from bedMethyl file.

    Args:
        bed_file: Path to bedMethyl file

    Returns:
        Dictionary mapping position to methylation percentage
    """
    meth_data = {}

    with open(bed_file) as f:
        for line in f:
            fields = line.strip().split('\t')

            # Check if this is a methylation row (field 4 = 'm')
            if len(fields) >= 10 and fields[3] == 'm':
                pos = int(fields[1])  # Column 2: position (0-based)

                # Column 10 contains space-separated values
                # Format: "coverage methylation_pct ..."
                # We want the 2nd value (index 1)
                column10_values = fields[10].split()
                if len(column10_values) >= 2:
                    pct_meth = float(fields[10])
                    meth_data[pos] = pct_meth

    return meth_data


def get_methylation_at_positions(bed_file, positions):
    """
    Get methylation values at specific positions.

    Args:
        bed_file: Path to bedMethyl file
        positions: List of positions to extract

    Returns:
        Dictionary mapping position to methylation percentage
    """
    all_meth = load_bed_methylation(bed_file)

    result = {}
    for pos in positions:
        if pos in all_meth:
            result[pos] = all_meth[pos]
        else:
            result[pos] = 0.0

    return result
