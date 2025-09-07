"""Utility functions for codoff analysis."""

import gzip
import sys
import traceback
from typing import Dict, List, Set, Tuple
import pyrodigal
from Bio import SeqIO
from Bio.Seq import Seq


def print_progress_bar(
    iteration: int,
    total: int,
    prefix: str = '',
    suffix: str = '',
    decimals: int = 1,
    length: int = 100,
    fill: str = '█',
    print_end: str = "\r"
) -> None:
    """
    Display a progress bar in the terminal.

    Function from: https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters
    Call in a loop to create terminal progress bar.

    Args:
        iteration: Current iteration number
        total: Total number of iterations
        prefix: Prefix string to display before the bar
        suffix: Suffix string to display after the bar
        decimals: Number of decimal places in percent complete
        length: Character length of the progress bar
        fill: Character used to fill the progress bar
        print_end: End character (e.g. "\\r", "\\r\\n")
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    sys.stderr.write(f'\r{prefix} |{bar}| {percent}% {suffix}')
    # Print New Line on Complete
    if iteration == total:
        sys.stderr.write('\n')


def parse_cds_coord(str_gbk_loc: str) -> Tuple[List[List[int]], int, int, str, bool]:
    """
    Parse GenBank CDS coordinate string into structured format.

    Args:
        str_gbk_loc: GenBank location string (e.g., "[1:100](+)")

    Returns:
        Tuple containing:
            - all_coords: List of coordinate lists [start, end, direction]
            - start: Overall start coordinate
            - end: Overall end coordinate
            - direction: Strand direction ('+' or '-')
            - is_multi_part: Whether the feature spans multiple coordinates
    """
    try:
        start = None
        end = None
        direction = None
        all_coords = []
        is_multi_part = False
        if 'join' not in str(str_gbk_loc) and 'order' not in str(str_gbk_loc):
            start = min([int(x.strip('>').strip('<')) for x in
                        str(str_gbk_loc)[1:].split(']')[0].split(':')]) + 1
            end = max([int(x.strip('>').strip('<')) for x in
                      str(str_gbk_loc)[1:].split(']')[0].split(':')])
            direction = str(str_gbk_loc).split('(')[1].split(')')[0]
            all_coords.append([start, end, direction])
        elif 'order' in str(str_gbk_loc):
            is_multi_part = True
            all_starts = []
            all_ends = []
            all_directions = []
            for exon_coord in str(str_gbk_loc)[6:-1].split(', '):
                ec_start = min(
                    [int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
                ec_end = max(
                    [int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
                ec_direction = exon_coord.split('(')[1].split(')')[0]
                all_starts.append(ec_start)
                all_ends.append(ec_end)
                all_directions.append(ec_direction)
                all_coords.append([ec_start, ec_end, ec_direction])
            assert (len(set(all_directions)) == 1)
            start = min(all_starts)
            end = max(all_ends)
            direction = all_directions[0]
        else:
            is_multi_part = True
            all_starts = []
            all_ends = []
            all_directions = []
            for exon_coord in str(str_gbk_loc)[5:-1].split(', '):
                ec_start = min(
                    [int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
                ec_end = max(
                    [int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
                ec_direction = exon_coord.split('(')[1].split(')')[0]
                all_starts.append(ec_start)
                all_ends.append(ec_end)
                all_directions.append(ec_direction)
                all_coords.append([ec_start, ec_end, ec_direction])
            assert (len(set(all_directions)) == 1)
            start = min(all_starts)
            end = max(all_ends)
            direction = all_directions[0]
        return (all_coords, start, end, direction, is_multi_part)
    except Exception:
        sys.stderr.write("Issues with GenBank CDS processing.\n")
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)


def gene_call_using_pyrodigal(
    genome_fasta_file: str,
    focal_scaffold: str,
    focal_start_coord: int,
    focal_end_coord: int
) -> Tuple[Dict[str, str], Set[str]]:
    """
    Call genes using pyrodigal for FASTA files.

    Uses pyrodigal to predict open reading frames in FASTA sequences and identifies
    which genes fall within the specified focal region coordinates.

    Args:
        genome_fasta_file: Path to FASTA file containing genome sequences
        focal_scaffold: Scaffold identifier for the focal region
        focal_start_coord: Start coordinate of the focal region
        focal_end_coord: End coordinate of the focal region

    Returns:
        Tuple containing:
            - locus_tag_sequences: Dictionary mapping gene IDs to DNA sequences
            - focal_lts: Set of gene IDs that fall within the focal region
    """
    locus_tag_sequences = {}
    focal_lts = set([])
    try:
        full_genome_seqs = []
        with open(genome_fasta_file) as off:
            for rec in SeqIO.parse(off, 'fasta'):
                full_genome_seqs.append(str(rec.seq))
        spacer = 'N' * 10000
        full_genome = Seq(spacer.join(full_genome_seqs))
        orf_finder = pyrodigal.GeneFinder()
        orf_finder.train(bytes(full_genome))
        locus_tag_id = 1
        with open(genome_fasta_file) as off:
            for rec in SeqIO.parse(off, "fasta"):
                scaffold = rec.id
                scaff_seq = str(rec.seq)
                for pred in orf_finder.find_genes(bytes(rec.seq)):
                    start_coord = pred.begin
                    end_coord = pred.end
                    strand = pred.strand
                    gene_seq = Seq(scaff_seq[start_coord - 1:end_coord])
                    if strand == -1:
                        gene_seq = str(Seq(scaff_seq[start_coord - 1:end_coord]).reverse_complement())
                    locus_tag_name = 'CDS_' + str(locus_tag_id)
                    locus_tag_sequences[locus_tag_name] = gene_seq
                    if scaffold == focal_scaffold and start_coord >= focal_start_coord and end_coord <= focal_end_coord:
                        focal_lts.add(locus_tag_name)
                    locus_tag_id += 1
    except Exception:
        sys.stderr.write("Issues with pyrodigal gene calling.\n")
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)
    return (locus_tag_sequences, focal_lts)


def check_is_genbank_with_cds(gbk_file: str) -> bool:
    """
    Check if GenBank file contains CDS features.

    Validates that the specified GenBank file contains at least one CDS (coding sequence)
    feature, which is required for codon usage analysis.

    Args:
        gbk_file: Path to GenBank file to check

    Returns:
        True if file contains CDS features, False otherwise
    """
    try:
        gbk_with_cds = False

        ogf = None
        if gbk_file.endswith('.gz'):
            ogf = gzip.open(gbk_file, 'rt')
        else:
            ogf = open(gbk_file)
        for rec in SeqIO.parse(ogf, 'genbank'):
            for feat in rec.features:
                if feat.type == 'CDS':
                    gbk_with_cds = True
                    break
        ogf.close()

        return gbk_with_cds
    except Exception:
        sys.stderr.write("Issues with attemping to parse suspected GenBank file %s" % gbk_file)
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)


def confirm_fasta(fasta_file: str) -> bool:
    """
    Check if file is a valid FASTA file.

    Validates that the specified file can be parsed as FASTA format and contains
    at least one sequence record.

    Args:
        fasta_file: Path to FASTA file to validate

    Returns:
        True if file contains FASTA records, False otherwise
    """
    try:
        fasta_recs = 0
        with open(fasta_file) as off:
            for rec in SeqIO.parse(off, 'fasta'):
                fasta_recs += 1
        return fasta_recs > 0
    except Exception:
        sys.stderr.write("Issues with attemping to parse suspected FASTA file %s" % fasta_file)
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)
