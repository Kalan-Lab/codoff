"""
Utility functions for codoff.

This module provides utility functions for file handling, gene calling,
and coordinate parsing used by the codoff package.
"""

import sys
import traceback
import gzip
from typing import List, Tuple, Set, Dict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
import pyrodigal


def print_progress_bar(iteration: int, total: int, prefix: str = '', 
                      suffix: str = '', decimals: int = 1, length: int = 100, 
                      fill: str = '█', print_end: str = "\r") -> None:
    """
    Create a terminal progress bar.
    
    Function from: https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters
    
    Parameters
    ----------
    iteration : int
        Current iteration number
    total : int
        Total number of iterations
    prefix : str, optional
        Prefix string to display before the bar, by default ''
    suffix : str, optional
        Suffix string to display after the bar, by default ''
    decimals : int, optional
        Number of decimal places in percent complete, by default 1
    length : int, optional
        Character length of the progress bar, by default 100
    fill : str, optional
        Character used to fill the progress bar, by default '█'
    print_end : str, optional
        End character (e.g., "\\r", "\\r\\n"), by default "\\r"
    """
    percent = ("{0:." + str(decimals) + "f}").format(
        100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    sys.stderr.write(f'\r{prefix} |{bar}| {percent}% {suffix}')
    
    # Print new line on complete
    if iteration == total:
        sys.stderr.write('\n')


def parse_cds_coord(location: FeatureLocation) -> Tuple[List[List[int]], int, int, str, bool]:
    """
    Parse a Biopython FeatureLocation object to extract start, end, and strand.
    
    This function replaces manual string parsing with direct use of Biopython's
    object attributes for robustness, especially with multi-exon eukaryotic features.
    
    Parameters
    ----------
    location : FeatureLocation
        Biopython FeatureLocation object from a SeqFeature.
    
    Returns
    -------
    Tuple[List[List[int]], int, int, str, bool]
        Tuple containing:
        - all_coords: List of [start, end, direction] for each exon
        - start: Minimum start coordinate (1-based)
        - end: Maximum end coordinate
        - direction: Strand direction ('+' or '-')
        - is_multi_part: Whether the CDS spans multiple exons
    """
    try:
        all_coords = []
        is_multi_part = len(location.parts) > 1

        # Overall direction and bounds
        direction = '+' if location.strand >= 0 else '-'
        start = int(location.start) + 1  # Convert from 0-based to 1-based
        end = int(location.end)
        
        # Handle compound locations (e.g., joins for exons)
        if is_multi_part:
            for part in location.parts:
                part_direction = '+' if part.strand >= 0 else '-'
                all_coords.append([int(part.start) + 1, int(part.end), part_direction])
        else:
            all_coords.append([start, end, direction])
            
        return all_coords, start, end, direction, is_multi_part
        
    except Exception as e:
        sys.stderr.write("Issues with GenBank CDS processing.\n")
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)


def gene_call_using_pyrodigal(genome_fasta_file: str, focal_scaffold: str, 
                             focal_start_coord: int, focal_end_coord: int) -> Tuple[Dict[str, str], Set[str], Dict[str, Tuple[str, int, int]]]:
    """
    Perform gene calling using pyrodigal on a FASTA file.
    
    Parameters
    ----------
    genome_fasta_file : str
        Path to the input FASTA file
    focal_scaffold : str
        Scaffold identifier for the focal region
    focal_start_coord : int
        Start coordinate of the focal region
    focal_end_coord : int
        End coordinate of the focal region
    
    Returns
    -------
    Tuple[Dict[str, str], Set[str], Dict[str, Tuple[str, int, int]]]
        Tuple containing:
        - locus_tag_sequences: Dictionary mapping locus tags to sequences
        - focal_lts: Set of locus tags in the focal region
        - gene_coords: Dictionary mapping locus tags to (scaffold, start, end)
    """
    locus_tag_sequences = {}
    focal_lts = set()
    gene_coords = {}
    
    try:
        # Read all sequences and concatenate with spacers
        full_genome_seqs = []
        with open(genome_fasta_file) as off:
            for rec in SeqIO.parse(off, 'fasta'):
                full_genome_seqs.append(str(rec.seq))
        
        spacer = 'N' * 10000
        full_genome = Seq(spacer.join(full_genome_seqs))
        
        # Train pyrodigal on the full genome
        orf_finder = pyrodigal.GeneFinder()
        orf_finder.train(bytes(full_genome))
        
        # Find genes in each scaffold
        locus_tag_id = 1
        with open(genome_fasta_file) as off:
            for rec in SeqIO.parse(off, "fasta"):
                scaffold = rec.id
                scaff_seq = str(rec.seq)
                
                for pred in orf_finder.find_genes(bytes(rec.seq)):
                    start_coord = pred.begin
                    end_coord = pred.end
                    strand = pred.strand
                    
                    gene_seq = Seq(scaff_seq[start_coord-1:end_coord])
                    if strand == -1:
                        gene_seq = str(Seq(scaff_seq[start_coord-1:end_coord])
                                     .reverse_complement())
                    
                    locus_tag_name = 'CDS_' + str(locus_tag_id)
                    locus_tag_sequences[locus_tag_name] = gene_seq
                    gene_coords[locus_tag_name] = (scaffold, start_coord, end_coord)
                    
                    # Check if gene is in focal region
                    if (scaffold == focal_scaffold and 
                        start_coord >= focal_start_coord and 
                        end_coord <= focal_end_coord):
                        focal_lts.add(locus_tag_name)
                    
                    locus_tag_id += 1
                    
    except Exception as e:
        sys.stderr.write("Issues with pyrodigal gene calling.\n")
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)
    
    return locus_tag_sequences, focal_lts, gene_coords


def check_is_genbank_with_cds(gbk_file: str) -> bool:
    """
    Check if a file is a valid GenBank file containing CDS features.
    
    Parameters
    ----------
    gbk_file : str
        Path to the file to check
    
    Returns
    -------
    bool
        True if the file is a valid GenBank file with CDS features, False otherwise
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
        
    except Exception as e:
        sys.stderr.write(f"Issues with attempting to parse suspected GenBank file {gbk_file}")
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)


def confirm_fasta(fasta_file: str) -> bool:
    """
    Check if a file is a valid FASTA file.
    
    Parameters
    ----------
    fasta_file : str
        Path to the file to check
    
    Returns
    -------
    bool
        True if the file is a valid FASTA file, False otherwise
    """
    try:
        fasta_recs = 0
        with open(fasta_file) as off:
            for rec in SeqIO.parse(off, 'fasta'):
                fasta_recs += 1
        return fasta_recs > 0
        
    except Exception as e:
        sys.stderr.write(f"Issues with attempting to parse suspected FASTA file {fasta_file}")
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)


def calculate_n50(scaffold_lengths: List[int]) -> int:
    """
    Calculate the N50 statistic for an assembly.
    
    N50 is defined as the sequence length of the shortest contig at 50% of the 
    total genome length.
    
    Parameters
    ----------
    scaffold_lengths : List[int]
        List of scaffold/contig lengths
    
    Returns
    -------
    int
        The N50 value
    """
    if not scaffold_lengths:
        return 0
    
    sorted_lengths = sorted(scaffold_lengths, reverse=True)
    total_length = sum(sorted_lengths)
    target_length = total_length / 2
    
    cumulative_length = 0
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= target_length:
            return length
    
    return sorted_lengths[-1]


def check_assembly_quality(genome_file: str, focal_region_size: int, 
                           verbose: bool = True) -> None:
    """
    Check if the assembly N50 is larger than the focal region size.
    
    If the N50 is smaller than or equal to the focal region size, this indicates
    the assembly is too fragmented for reliable analysis and an error is raised.
    
    Parameters
    ----------
    genome_file : str
        Path to the genome file (GenBank or FASTA)
    focal_region_size : int
        Size of the focal region in base pairs
    verbose : bool, optional
        Whether to print progress messages, by default True
    
    Raises
    ------
    SystemExit
        If N50 is less than or equal to the focal region size
    """
    try:
        scaffold_lengths = []
        
        # Determine file type and parse
        is_genbank = check_is_genbank_with_cds(genome_file)
        
        if is_genbank:
            ogf = None
            if genome_file.endswith('.gz'):
                ogf = gzip.open(genome_file, 'rt')
            else:
                ogf = open(genome_file)
            
            for rec in SeqIO.parse(ogf, 'genbank'):
                scaffold_lengths.append(len(rec.seq))
            ogf.close()
        else:
            # Try as FASTA
            with open(genome_file) as off:
                for rec in SeqIO.parse(off, 'fasta'):
                    scaffold_lengths.append(len(rec.seq))
        
        if not scaffold_lengths:
            if verbose:
                sys.stderr.write('Error: No sequences found in genome file.\n')
            sys.exit(1)
        
        n50 = calculate_n50(scaffold_lengths)
        
        if n50 <= focal_region_size:
            if verbose:
                sys.stderr.write(f'Error: Assembly N50 ({n50} bp) is less than or equal to the focal region size ({focal_region_size} bp).\n')
                sys.stderr.write('The assembly is too fragmented for reliable analysis of this region size.\n')
            sys.exit(1)
        
        if verbose:
            sys.stderr.write(f'Assembly N50: {n50} bp (sufficient for focal region size: {focal_region_size} bp)\n')
    
    except SystemExit:
        raise
    except Exception as e:
        sys.stderr.write('Issues with checking assembly quality.\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)


# Backward compatibility aliases
printProgressBar = print_progress_bar
parseCDSCoord = parse_cds_coord
geneCallUsingPyrodigal = gene_call_using_pyrodigal
checkIsGenBankWithCDS = check_is_genbank_with_cds
confirmFasta = confirm_fasta
calculateN50 = calculate_n50
checkAssemblyQuality = check_assembly_quality