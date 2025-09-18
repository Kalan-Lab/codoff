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


def parse_cds_coord(str_gbk_loc: str) -> Tuple[List[List[int]], int, int, str, bool]:
    """
    Parse GenBank CDS coordinate string to extract start, end, direction.
    
    Parameters
    ----------
    str_gbk_loc : str
        GenBank location string (e.g., "complement(123..456)")
    
    Returns
    -------
    Tuple[List[List[int]], int, int, str, bool]
        Tuple containing:
        - all_coords: List of [start, end, direction] for each exon
        - start: Minimum start coordinate
        - end: Maximum end coordinate  
        - direction: Strand direction ('+' or '-')
        - is_multi_part: Whether the CDS spans multiple exons
    """
    try:
        start = None
        end = None
        direction = None
        all_coords = []
        is_multi_part = False
        
        if not ('join' in str(str_gbk_loc) or 'order' in str(str_gbk_loc)):
            # Single exon case
            coords_str = str(str_gbk_loc)[1:].split(']')[0]
            coords = [int(x.strip('>').strip('<')) for x in coords_str.split(':')]
            start = min(coords) + 1
            end = max(coords)
            direction = str(str_gbk_loc).split('(')[1].split(')')[0]
            all_coords.append([start, end, direction])
            
        elif 'order' in str(str_gbk_loc):
            # Multi-exon case with order
            is_multi_part = True
            all_starts = []
            all_ends = []
            all_directions = []
            
            coords_str = str(str_gbk_loc)[6:-1]  # Remove 'order(' and ')'
            for exon_coord in coords_str.split(', '):
                ec_coords_str = exon_coord[1:].split(']')[0]
                ec_coords = [int(x.strip('>').strip('<')) 
                            for x in ec_coords_str.split(':')]
                ec_start = min(ec_coords) + 1
                ec_end = max(ec_coords)
                ec_direction = exon_coord.split('(')[1].split(')')[0]
                
                all_starts.append(ec_start)
                all_ends.append(ec_end)
                all_directions.append(ec_direction)
                all_coords.append([ec_start, ec_end, ec_direction])
            
            assert len(set(all_directions)) == 1
            start = min(all_starts)
            end = max(all_ends)
            direction = all_directions[0]
            
        else:
            # Multi-exon case with join
            is_multi_part = True
            all_starts = []
            all_ends = []
            all_directions = []
            
            coords_str = str(str_gbk_loc)[5:-1]  # Remove 'join(' and ')'
            for exon_coord in coords_str.split(', '):
                ec_coords_str = exon_coord[1:].split(']')[0]
                ec_coords = [int(x.strip('>').strip('<')) 
                            for x in ec_coords_str.split(':')]
                ec_start = min(ec_coords) + 1
                ec_end = max(ec_coords)
                ec_direction = exon_coord.split('(')[1].split(')')[0]
                
                all_starts.append(ec_start)
                all_ends.append(ec_end)
                all_directions.append(ec_direction)
                all_coords.append([ec_start, ec_end, ec_direction])
            
            assert len(set(all_directions)) == 1
            start = min(all_starts)
            end = max(all_ends)
            direction = all_directions[0]
            
        return all_coords, start, end, direction, is_multi_part
        
    except Exception as e:
        sys.stderr.write("Issues with GenBank CDS processing.\n")
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)


def gene_call_using_pyrodigal(genome_fasta_file: str, focal_scaffold: str, 
                             focal_start_coord: int, focal_end_coord: int) -> Tuple[Dict[str, str], Set[str]]:
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
    Tuple[Dict[str, str], Set[str]]
        Tuple containing:
        - locus_tag_sequences: Dictionary mapping locus tags to sequences
        - focal_lts: Set of locus tags in the focal region
    """
    locus_tag_sequences = {}
    focal_lts = set()
    
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
    
    return locus_tag_sequences, focal_lts


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


# Backward compatibility aliases
printProgressBar = print_progress_bar
parseCDSCoord = parse_cds_coord
geneCallUsingPyrodigal = gene_call_using_pyrodigal
checkIsGenBankWithCDS = check_is_genbank_with_cds
confirmFasta = confirm_fasta