"""
Codon usage analysis and simulation functions for codoff.

This module provides the main functionality for analyzing codon usage
bias in genomic regions using Monte Carlo simulations.
"""

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import os
import random
import sys
import gzip
import math
import traceback
from collections import defaultdict
from operator import itemgetter
from typing import Dict, Any, Optional, List, Set

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from scipy import stats, spatial

from codoff import util

# Constants
VALID_BASES = {'A', 'C', 'G', 'T'}

def check_data_type(variable: Any, expected_type: type) -> bool:
    """
    Check if the variable is of the expected type.
    
    Parameters
    ----------
    variable : Any
        The variable to check
    expected_type : type
        The expected type of the variable
    
    Returns
    -------
    bool
        True if the variable is of the expected type, False otherwise
    """
    return isinstance(variable, expected_type)


def extract_genome_codon_data(full_genome_file: str, 
                             verbose: bool = True) -> Dict[str, Any]:
    """
    Extract codon usage data from the full genome file once to avoid redundant computation.
    
    Parameters
    ----------
    full_genome_file : str
        The path to the full genome file in GenBank format
    verbose : bool, optional
        Whether to print progress messages to stderr, by default True
    
    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - 'locus_tag_sequences': dict mapping locus_tag to nucleotide sequence
        - 'gene_codons': dict mapping locus_tag to codon counts
        - 'gene_list': list of all locus_tags
        - 'all_cods': set of all codons found
        - 'total_cds_length': total CDS length in the genome
    """
    try:
        assert check_data_type(full_genome_file, str)
        assert os.path.isfile(full_genome_file)
        assert util.checkIsGenBankWithCDS(full_genome_file)
    except Exception:
        sys.stderr.write('The full genome file must be a string to a file.\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    locus_tag_sequences = {}
    gene_codons = defaultdict(lambda: defaultdict(int))
    gene_list = []
    all_cods = set()
    all_codon_counts = defaultdict(int)
    total_cds_length = 0
    
    try:
        ofgbk = None
        if full_genome_file.endswith('.gz'):
            ofgbk = gzip.open(full_genome_file, 'rt')
        else:
            ofgbk = open(full_genome_file)
        
        for rec in SeqIO.parse(ofgbk, 'genbank'):
            full_sequence = str(rec.seq)
            for feature in rec.features:
                if not feature.type == 'CDS': 
                    continue
                locus_tag = None
                try:
                    locus_tag = feature.qualifiers.get('locus_tag')[0]
                except Exception:
                    try:
                        locus_tag = feature.qualifiers.get('gene')[0]
                    except Exception:
                        locus_tag = feature.qualifiers.get('protein_id')[0]
                assert locus_tag is not None
                
                all_coords, start, end, direction, is_multi_part = util.parseCDSCoord(str(feature.location))
                if end >= len(full_sequence): 
                    end = len(full_sequence)
                nucl_seq = ''
                for sc, ec, dc in sorted(all_coords, key=itemgetter(0), reverse=False):
                    if ec >= len(full_sequence):
                        nucl_seq += full_sequence[sc - 1:]
                    else:
                        nucl_seq += full_sequence[sc - 1:ec]

                if direction == '-':
                    nucl_seq = str(Seq(nucl_seq).reverse_complement())
                locus_tag_sequences[locus_tag] = nucl_seq

                if len(str(nucl_seq)) % 3 == 0:
                    total_cds_length += len(nucl_seq)
                    gene_list.append(locus_tag)
                    seq = locus_tag_sequences[locus_tag]
                    codon_seq = [str(seq)[i:i + 3] for i in range(0, len(str(seq)), 3)]
                    for cod in list(codon_seq):
                        if not(len(cod) == 3 and cod[0] in VALID_BASES and cod[1] in VALID_BASES and cod[2] in VALID_BASES): 
                            continue
                        gene_codons[locus_tag][cod] += 1
                        all_codon_counts[cod] += 1
                        all_cods.add(cod)
        ofgbk.close()
        
        if total_cds_length == 0:
            if verbose:
                sys.stderr.write('Error: The genome appears to have no CDS features. Please check the input.\n')
            sys.exit(1)
            
    except Exception:
        sys.stderr.write('Issues with extracting codon usage info from genome GenBank file.\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    # Convert defaultdict to regular dict for better performance and pickling
    gene_codons_dict = {gene: dict(codons) for gene, codons in gene_codons.items()}
    
    return {
        'locus_tag_sequences': locus_tag_sequences,
        'gene_codons': gene_codons_dict,
        'gene_list': gene_list,
        'all_cods': all_cods,
        'all_codon_counts': dict(all_codon_counts),
        'total_cds_length': total_cds_length
    }


def process_bgc_with_cached_data(genome_data: Dict[str, Any], focal_genbank_files: List[str], outfile: Optional[str] = None, plot_outfile: Optional[str] = None, verbose: bool = True, num_sims: int = 10000) -> Dict[str, Any]:
    """
    Process BGC regions using pre-computed genome data to avoid redundant computation.
    
    Parameters
    ----------
    genome_data : dict
        Pre-computed genome data from extract_genome_codon_data()
    focal_genbank_files : list
        The list of paths to the GenBank files for the focal region.
    outfile : str
        The path to the output file. If not provided, the output will be 
        printed to stdout.
    plot_outfile : str
        The path to the plot output file. If not provided, no plot will 
        be made.
    verbose : bool
        Whether to print progress messages to stderr. Default is True.
    
    Returns
    -------
    dict
        Result dictionary from _stat_calc_and_simulation()
    """
    try:
        assert check_data_type(focal_genbank_files, list)
        for reg_gbk in focal_genbank_files:
            assert check_data_type(reg_gbk, str)
            assert os.path.isfile(reg_gbk)
            assert util.checkIsGenBankWithCDS(reg_gbk)
    except Exception:
        sys.stderr.write('The focal region GenBank files must be a list of string paths to existing files.\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    if outfile != 'stdout' and outfile is not None and os.path.isfile(outfile):
        sys.stderr.write('The outfile must be a string path to a file which does not already exist.\n')
        sys.stderr.write('Outfile: %s\n' % str(outfile))
        sys.exit(1)

    if plot_outfile != 'stdout' and plot_outfile is not None and os.path.isfile(plot_outfile):
        sys.stderr.write('The plot outfile must be a string path to a file which does not already exist.\n')
        sys.stderr.write('Plot outfile: %s\n' % str(plot_outfile))
        sys.exit(1)    

    try:
        assert check_data_type(verbose, bool)
    except Exception:
        sys.stderr.write('The verbose flag must be a boolean.\n')
        sys.stderr.write('Verbose: %s\n' % str(verbose))
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    # Extract pre-computed data
    locus_tag_sequences = genome_data['locus_tag_sequences']
    gene_codons = genome_data['gene_codons']
    gene_list = genome_data['gene_list']
    all_cods = genome_data['all_cods']
    all_codon_counts = genome_data['all_codon_counts']
    total_cds_length = genome_data['total_cds_length']

    cod_freq_dict_focal = defaultdict(int)
    cod_freq_dict_background = defaultdict(int)
    foc_codon_count = 0
    focal_cds_length = 0

    try:
        # Parse GenBank files of focal regions for locus_tags
        focal_lts = set()
        for foc_gbk in focal_genbank_files:
            ogbk = None
            if foc_gbk.endswith('.gz'):
                ogbk = gzip.open(foc_gbk, 'rt')
            else:
                ogbk = open(foc_gbk)
            for rec in SeqIO.parse(ogbk, 'genbank'):
                for feature in rec.features:
                    if not feature.type == 'CDS': 
                        continue
                    locus_tag = None
                    try:
                        locus_tag = feature.qualifiers.get('locus_tag')[0]
                    except Exception:
                        try:
                            locus_tag = feature.qualifiers.get('gene')[0]
                        except Exception:
                            locus_tag = feature.qualifiers.get('protein_id')[0]
                    assert locus_tag is not None
                    focal_lts.add(locus_tag)
            ogbk.close()
        
        # Check for missing focal locus tags and issue warnings
        missing_focal_lts = focal_lts - set(locus_tag_sequences.keys())
        if missing_focal_lts and verbose:
            sys.stderr.write('Warning: The following focal region locus tags were not found in the full genome GenBank file:\n')
            for missing_lt in sorted(missing_focal_lts):
                sys.stderr.write('  %s\n' % missing_lt)
            sys.stderr.write('These locus tags will be ignored in the analysis.\n')
        
        # Calculate focal region statistics using pre-computed per-gene codon counts
        # First, compute focal locus total nucleotide length (for size checks)
        for locus_tag in focal_lts:
            if locus_tag in locus_tag_sequences:
                seq = locus_tag_sequences[locus_tag]
                if len(str(seq)) % 3 == 0:
                    focal_cds_length += len(seq)
        # Build focal codon counts by summing cached gene codon counts
        for locus_tag in focal_lts:
            if locus_tag in gene_codons:
                for cod, count in gene_codons[locus_tag].items():
                    cod_freq_dict_focal[cod] += count
                    foc_codon_count += count
        # Background is genome total minus focal (v1.2.1 behavior)
        for cod in all_cods:
            cod_freq_dict_background[cod] = all_codon_counts.get(cod, 0) - cod_freq_dict_focal[cod]

        if focal_cds_length == 0:
            if verbose:
                sys.stderr.write('Error: The focal region appears to have no CDS features. This might be because locus_tags in the focal region(s) GenBank file do not match locus_tags in the full genome GenBank file. Please check the input.\n')
            sys.exit(1)

        size_comparison = focal_cds_length/total_cds_length
        if size_comparison >= 0.05:
            sys.stderr.write('Error: The size of the focal region is >5%% of the full genome. This is not be appropriate for codoff analysis.\n')
            sys.stderr.write('Focal region(s): %s\n' % str(focal_genbank_files))
            sys.exit(1)

    except Exception:
        sys.stderr.write('Issues with determining codon usage info for inputs with focal region and genome provided as GenBank.\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    # Run rest of codoff using pre-computed data
    result = _stat_calc_and_simulation(all_cods, cod_freq_dict_focal, cod_freq_dict_background, gene_list, gene_codons, foc_codon_count, all_codon_counts, outfile=outfile, plot_outfile=plot_outfile, verbose=verbose,  num_sims=num_sims)
    return result


def codoff_main_gbk(full_genome_file: str, focal_genbank_files: List[str], outfile: Optional[str] = None, plot_outfile: Optional[str] = None, verbose: bool = True, genome_data: Optional[Dict[str, Any]] = None, num_sims: int = 10000) -> Dict[str, Any]:
    """
    A full genome and a specific region must each be provided in 
    GenBank format, with locus_tags overlapping. locus_tags in the
    focal region GenBank that are not in the full genome GenBank will 
    be ignored.

    It calls the private function _stat_calc_and_simulation() to 
    perform the main statistical calculations and simulations used 
    for inference of the empirical P-value.
    
    The function takes the following arguments:
    * full_genome_file: str
        The path to the full genome file in GenBank format.
    * focal_genbank_files: list
        The list of paths to the GenBank files for the focal region.
    * outfile: str
        The path to the output file. If not provided, the output will be 
        printed to stdout.
    * plot_outfile: str
        The path to the plot output file. If not provided, no plot will 
        be made.
    * verbose: bool
        Whether to print progress messages to stderr. Default is True.
    * genome_data: dict, optional
        Pre-computed genome data from extract_genome_codon_data() to avoid 
        redundant computation. If provided, full_genome_file will be ignored.
    """
    
    # If genome_data is provided, use it instead of processing the full_genome_file
    if genome_data is not None:
        return process_bgc_with_cached_data(genome_data, focal_genbank_files, outfile, plot_outfile, verbose, num_sims)
    
    try:
        assert check_data_type(full_genome_file, str)
        assert os.path.isfile(full_genome_file)
        assert util.checkIsGenBankWithCDS(full_genome_file)
    except Exception:
        sys.stderr.write('The full genome file must be a string to a file.\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    try:
        assert check_data_type(focal_genbank_files, list)
        for reg_gbk in focal_genbank_files:
            assert check_data_type(reg_gbk, str)
            assert os.path.isfile(reg_gbk)
            assert util.checkIsGenBankWithCDS(reg_gbk)
    except Exception:
        sys.stderr.write('The focal region GenBank files must be a list of string paths to existing files.\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)
 
    if outfile != 'stdout' and outfile is not None and os.path.isfile(outfile):
        sys.stderr.write('The outfile must be a string path to a file which does not already exist.\n')
        sys.stderr.write('Outfile: %s\n' % str(outfile))
        sys.exit(1)

    if plot_outfile != 'stdout' and plot_outfile is not None and os.path.isfile(plot_outfile):
        sys.stderr.write('The plot outfile must be a string path to a file which does not already exist.\n')
        sys.stderr.write('Plot outfile: %s\n' % str(plot_outfile))
        sys.exit(1)    

    try:
        assert check_data_type(verbose, bool)
    except Exception:
        sys.stderr.write('The verbose flag must be a boolean.\n')
        sys.stderr.write('Verbose: %s\n' % str(verbose))
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    cod_freq_dict_focal = defaultdict(int)
    cod_freq_dict_background = defaultdict(int)
    all_cods = set()
    all_codon_counts = defaultdict(int)
    gene_codons = defaultdict(lambda: defaultdict(int))
    gene_list = []
    foc_codon_count = 0
    try:
        # parse GenBank files of focal regions for locus_tags
        focal_lts = set()
        for foc_gbk in focal_genbank_files:
            ogbk = None
            if foc_gbk.endswith('.gz'):
                ogbk = gzip.open(foc_gbk, 'rt')
            else:
                ogbk = open(foc_gbk)
            for rec in SeqIO.parse(ogbk, 'genbank'):
                for feature in rec.features:
                    if not feature.type == 'CDS': continue
                    locus_tag = None
                    try:
                        locus_tag = feature.qualifiers.get('locus_tag')[0]
                    except Exception:
                        try:
                            locus_tag = feature.qualifiers.get('gene')[0]
                        except Exception:
                            locus_tag = feature.qualifiers.get('protein_id')[0]
                    assert locus_tag is not None
                    focal_lts.add(locus_tag)
            ogbk.close()
        
        
        # parse nucleotide coding sequences from full-genome Genbank
        locus_tag_sequences = {}
        focal_cds_length = 0
        total_cds_length = 0 
        if not os.path.isfile(full_genome_file):
            sys.stderr.write('Unable to find full-genome GenBank file %s.\n' % foc_gbk)
            sys.exit(1)

        ofgbk = None
        if full_genome_file.endswith('.gz'):
            ofgbk = gzip.open(full_genome_file, 'rt')
        else:
            ofgbk = open(full_genome_file)
        for rec in SeqIO.parse(ofgbk, 'genbank'):
            full_sequence = str(rec.seq)
            for feature in rec.features:
                if not feature.type == 'CDS': continue
                locus_tag = None
                try:
                    locus_tag = feature.qualifiers.get('locus_tag')[0]
                except Exception:
                    try:
                        locus_tag = feature.qualifiers.get('gene')[0]
                    except Exception:
                        locus_tag = feature.qualifiers.get('protein_id')[0]
                assert locus_tag is not None
                all_coords, start, end, direction, is_multi_part = util.parseCDSCoord(str(feature.location))
                if end >= len(full_sequence): end = len(full_sequence)
                nucl_seq = ''
                for sc, ec, dc in sorted(all_coords, key=itemgetter(0), reverse=False):
                    if ec >= len(full_sequence):
                        nucl_seq += full_sequence[sc - 1:]
                    else:
                        nucl_seq += full_sequence[sc - 1:ec]

                if direction == '-':
                    nucl_seq = str(Seq(nucl_seq).reverse_complement())
                locus_tag_sequences[locus_tag] = nucl_seq

                if len(str(nucl_seq)) % 3 == 0:
                    if locus_tag in focal_lts:
                        focal_cds_length += len(nucl_seq)
                    total_cds_length += len(nucl_seq)
        ofgbk.close()
        
        # Check for missing focal locus tags and issue warnings
        missing_focal_lts = focal_lts - set(locus_tag_sequences.keys())
        if missing_focal_lts and verbose:
            sys.stderr.write('Warning: The following focal region locus tags were not found in the full genome GenBank file:\n')
            for missing_lt in sorted(missing_focal_lts):
                sys.stderr.write('  %s\n' % missing_lt)
            sys.stderr.write('These locus tags will be ignored in the analysis.\n')
        
        if total_cds_length == 0:
            if verbose:
                sys.stderr.write('Error: The genome appears to have no CDS features. Please check the input.\n')
            sys.exit(1)
        
        if focal_cds_length == 0:
            if verbose:
                sys.stderr.write('Error: The focal region appears to have no CDS features. This might be because locus_tags in the focal region(s) GenBank file do not match locus_tags in the full genome GenBank file. Please check the input.\n')
            sys.exit(1)

        size_comparison = focal_cds_length/total_cds_length
        if size_comparison >= 0.05:
            sys.stderr.write('Error: The size of the focal region is >5%% of the full genome. This is not be appropriate for codoff analysis.\n')
            sys.stderr.write('Focal region(s): %s\n' % str(focal_genbank_files))
            sys.exit(1)

        # get codon frequencies for CDS in BGC and background genome
        for locus_tag in locus_tag_sequences:
            gene_list.append(locus_tag)
            seq = locus_tag_sequences[locus_tag]
            if not len(str(seq)) % 3 == 0:
                if verbose:
                    sys.stderr.write("The locus tag %s is ignored because it was not of length 3.\n" % locus_tag)
                continue
            codon_seq = [str(seq)[i:i + 3] for i in range(0, len(str(seq)), 3)]
            for cod in list(codon_seq):
                if not(len(cod) == 3 and cod[0] in VALID_BASES and cod[1] in VALID_BASES and cod[2] in VALID_BASES): continue
                gene_codons[locus_tag][cod] += 1
                all_codon_counts[cod] += 1
                if locus_tag in focal_lts:
                    foc_codon_count += 1
                    cod_freq_dict_focal[cod] += 1
                else:
                    cod_freq_dict_background[cod] += 1
                all_cods.add(cod)
    except Exception:
        sys.stderr.write('Issues with determinining codon usage info for inputs with focal region and genome provided as GenBank.\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    # run rest of codoff (separate function to avoid redundancy between codoff_main_gbk() and codoff_main_coords()
    # Use all genes for simulation pool to match v1.2.1 behavior
    result = _stat_calc_and_simulation(all_cods, cod_freq_dict_focal, cod_freq_dict_background, gene_list, gene_codons, foc_codon_count, all_codon_counts, outfile=outfile, plot_outfile=plot_outfile, verbose=verbose,  num_sims=num_sims)
    return result


def codoff_main_coords(full_genome_file: str, focal_scaffold: str, focal_start_coord: int, focal_end_coord: int, outfile: Optional[str] = None, plot_outfile: Optional[str] = None, verbose: bool = True, num_sims: int = 10000) -> Dict[str, Any]:
    """
    A full genome file can be provided in either GenBank or FASTA 
    format. If the latter, pyrodigal is used for gene calling, 
    so it only works for bacteria.
    
    Afterwards, coordinates provided by users for the focal 
    region of interest are used to partition which locus tags for
    CDS features belong to the focal region and which belong to the
    background genome.

    It calls the private function _stat_calc_and_simulation() to 
    perform the main statistical calculations and simulations used 
    for inference of the empirical P-value.
    
    The function takes the following arguments:
    * full_genome_file: str
        The path to the full genome file in GenBank format.
    * focal_scaffold: str
        The scaffold identifier for the focal region.
    * focal_start_coord: int
        The start coordinate for the focal region.  
    * focal_end_coord: int
        The end coordinate for the focal region.
    * outfile: str
        The path to the output file. If not provided, the output will be 
        printed to stdout.
    * plot_outfile: str
        The path to the plot output file. If not provided, no plot will
        be made.
    * verbose: bool
        Whether to print progress messages to stderr. Default is True.
    """
    
    try:
        assert check_data_type(full_genome_file, str)
        assert os.path.isfile(full_genome_file)
    except Exception:
        sys.stderr.write('The full genome file %s must be a string path to an existing file.\n' % full_genome_file)
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    try:
        assert check_data_type(focal_scaffold, str)
        assert check_data_type(focal_start_coord, int)
        assert check_data_type(focal_end_coord, int)
        assert focal_start_coord < focal_end_coord
    except Exception:
        sys.stderr.write('The focal region must be specified via a string for the scaffold and two integers for the start and end coordinates.\n')
        sys.stderr.write('Scaffold: %s\n' % focal_scaffold)
        sys.stderr.write('Start coord: %s\n' % str(focal_start_coord))
        sys.stderr.write('End coord: %s\n' % str(focal_end_coord))
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)
    
    if outfile != 'stdout' and outfile is not None and os.path.isfile(outfile):
        sys.stderr.write('The outfile must be a string path to a file which does not already exist.\n')
        sys.stderr.write('Outfile: %s\n' % str(outfile))
        sys.exit(1)

    if plot_outfile != 'stdout' and plot_outfile is not None and os.path.isfile(plot_outfile):
        sys.stderr.write('The plot outfile must be a string path to a file which does not already exist.\n')
        sys.stderr.write('Plot outfile: %s\n' % str(plot_outfile))
        sys.exit(1)    

    try:
        assert check_data_type(verbose, bool)
    except Exception:
        sys.stderr.write('The verbose flag must be a boolean.\n')
        sys.stderr.write('Verbose: %s\n' % str(verbose))
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    cod_freq_dict_focal = defaultdict(int)
    cod_freq_dict_background = defaultdict(int)
    all_cods = set()
    all_codon_counts = defaultdict(int)
    gene_codons = defaultdict(lambda: defaultdict(int))
    gene_list = []
    foc_codon_count = 0

    if util.checkIsGenBankWithCDS(full_genome_file):
        try:
            # parse nucleotide coding sequences from full-genome Genbank
            locus_tag_sequences = {}
            focal_lts = set()
            focal_cds_length = 0
            total_cds_length = 0
                
            ofgbk = None
            if full_genome_file.endswith('.gz'):
                ofgbk = gzip.open(full_genome_file, 'rt')
            else:
                ofgbk = open(full_genome_file)
            for rec in SeqIO.parse(ofgbk, 'genbank'):
                full_sequence = str(rec.seq)
                for feature in rec.features:
                    if not feature.type == 'CDS': continue
                    locus_tag = None
                    try:
                        locus_tag = feature.qualifiers.get('locus_tag')[0]
                    except Exception:
                        try:
                            locus_tag = feature.qualifiers.get('gene')[0]
                        except Exception:
                            locus_tag = feature.qualifiers.get('protein_id')[0]
                    assert locus_tag is not None
                    all_coords, start, end, direction, is_multi_part = util.parseCDSCoord(str(feature.location))
                    if end >= len(full_sequence): end = len(full_sequence)
                    nucl_seq = ''
                    for sc, ec, dc in sorted(all_coords, key=itemgetter(0), reverse=False):
                        if ec >= len(full_sequence):
                            nucl_seq += full_sequence[sc - 1:]
                        else:
                            nucl_seq += full_sequence[sc - 1:ec]

                    if direction == '-':
                        nucl_seq = str(Seq(nucl_seq).reverse_complement())
                    locus_tag_sequences[locus_tag] = nucl_seq
                    if rec.id == focal_scaffold and start >= focal_start_coord and end <= focal_end_coord:
                        focal_lts.add(locus_tag)

                    if len(str(nucl_seq)) % 3 == 0:
                        if locus_tag in focal_lts:
                            focal_cds_length += len(nucl_seq)
                        total_cds_length += len(nucl_seq)
            ofgbk.close()
            
            # Check if any focal locus tags were found and issue warning if none
            if not focal_lts and verbose:
                sys.stderr.write('Warning: No CDS features were found in the specified focal region coordinates.\n')
                sys.stderr.write('Focal scaffold: %s, Start: %d, End: %d\n' % (focal_scaffold, focal_start_coord, focal_end_coord))
                sys.stderr.write('Please check that the coordinates are correct and contain CDS features.\n')
            
            if total_cds_length == 0:
                sys.stderr.write('Error: The genome appears to have no CDS features. Please check the input.\n')
                sys.exit(1)

            if focal_cds_length == 0:
                sys.stderr.write('Error: The focal region appears to have no CDS features. Please check the input.\n')
                sys.exit(1)
            
            size_comparison = focal_cds_length/total_cds_length
            if size_comparison >= 0.05:
                sys.stderr.write('Error: The size of the focal region is >5%% of the full genome. This is not be appropriate for codoff analysis\n')
                sys.stderr.write('Focal region size: %d\n' % focal_cds_length)
                sys.stderr.write('Genome size: %d\n' % total_cds_length)
                sys.exit(1)
            
            # get codon frequencies for CDS in focal and background genome
            for locus_tag in locus_tag_sequences:
                gene_list.append(locus_tag)
                seq = locus_tag_sequences[locus_tag]
                if not len(str(seq)) % 3 == 0:
                    sys.stderr.write("The locus tag %s is ignored because it was not of length 3.\n" % locus_tag)
                    continue
                codon_seq = [str(seq)[i:i + 3] for i in range(0, len(str(seq)), 3)]
                for cod in list(codon_seq):
                    if not(len(cod) == 3 and cod[0] in VALID_BASES and cod[1] in VALID_BASES and cod[2] in VALID_BASES): continue
                    gene_codons[locus_tag][cod] += 1
                    all_codon_counts[cod] += 1
                    if locus_tag in focal_lts:
                        foc_codon_count += 1
                        cod_freq_dict_focal[cod] += 1
                    else:
                        cod_freq_dict_background[cod] += 1
                    all_cods.add(cod)
        except Exception:
            sys.stderr.write('Issues with determinining codon usage info for inputs with full genome provided as GenBank.\n')
            sys.stderr.write(traceback.format_exc() + '\n')
            sys.exit(1)

    else:
        try:
            try:
                assert util.confirmFasta(full_genome_file)
            except Exception:
                sys.stderr.write('Input genome does not appear to be either GenBank with CDS features nor a FASTA file.\n')
                sys.stderr.write(traceback.format_exc() + '\n')
                sys.exit(1)            
        
            # parse nucleotide coding sequences from full-genome Genbank
            locus_tag_sequences, focal_lts = util.geneCallUsingPyrodigal(full_genome_file, focal_scaffold, focal_start_coord, focal_end_coord)
            
            total_cds_length = 0
            focal_cds_length = 0
            for locus_tag in locus_tag_sequences:
                seq = locus_tag_sequences[locus_tag]
                if len(str(seq)) % 3 == 0:
                    total_cds_length += len(str(seq))
                    if locus_tag in focal_lts:
                        focal_cds_length += len(str(seq))

            if total_cds_length == 0:
                if verbose:
                    sys.stderr.write('Error: The genome appears to have no CDS features. Please check the input.\n')
                sys.exit(1)
        
            if focal_cds_length == 0:
                if verbose:
                    sys.stderr.write('Error: The focal region appears to have no CDS features. Perhaps the region coordinates are too small?\n')
                sys.exit(1)
            
            size_comparison = focal_cds_length/total_cds_length
            if size_comparison >= 0.05:
                sys.stderr.write('Error: The size of the focal region is >5%% of the full genome. This is not be appropriate for codoff analysis.\n')
                sys.stderr.write('Focal region size: %s\n' % str(focal_cds_length))
                sys.stderr.write('Genome size: %s\n' % str(total_cds_length))
                sys.exit(1)    

            # get codon frequencies for CDS in focal and background genome
            for locus_tag in locus_tag_sequences:
                gene_list.append(locus_tag)
                seq = locus_tag_sequences[locus_tag]
                if not len(str(seq)) % 3 == 0:
                    if verbose:
                        sys.stderr.write("The locus tag %s is ignored because it was not of length 3.\n" % locus_tag)
                    continue
                codon_seq = [str(seq)[i:i + 3] for i in range(0, len(str(seq)), 3)]
                for cod in list(codon_seq):
                    if not(len(cod) == 3 and cod[0] in VALID_BASES and cod[1] in VALID_BASES and cod[2] in VALID_BASES): continue
                    gene_codons[locus_tag][cod] += 1
                    all_codon_counts[cod] += 1
                    if locus_tag in focal_lts:
                        foc_codon_count += 1
                        cod_freq_dict_focal[cod] += 1
                    else:
                        cod_freq_dict_background[cod] += 1
                    all_cods.add(cod)
        except Exception:
            sys.stderr.write('Issues with determinining codon usage info for inputs with full genome provided as FASTA.\n')
            sys.stderr.write(traceback.format_exc() + '\n')
            sys.exit(1)
        
        # run rest of codoff (separate function to avoid redundancy between codoff_main_gbk() and codoff_main_coords()
        # Use all genes for simulation pool to match v1.2.1 behavior
        result = _stat_calc_and_simulation(all_cods, cod_freq_dict_focal, cod_freq_dict_background, gene_list, gene_codons, foc_codon_count, all_codon_counts, outfile=outfile, plot_outfile=plot_outfile, verbose=verbose,  num_sims=num_sims)
        return result
    
    
def _stat_calc_and_simulation(all_cods: Set[str], cod_freq_dict_focal: Dict[str, int], cod_freq_dict_background: Dict[str, int], gene_list: List[str], gene_codons: Dict[str, Dict[str, int]], foc_codon_count: int, all_codon_counts: Dict[str, int], outfile: Optional[str] = None, plot_outfile: Optional[str] = None, verbose: bool = True, num_sims: int = 10000) -> Dict[str, Any]:
    """
    Private function that performs the main statistical calculations and Monte Carlo simulations.
    
    SIMULATION LOGIC:
    The simulation uses ALL genes from the genome (both focal and background genes) as the sampling pool.
    For each simulation iteration:
    1. Shuffle the complete list of genes (gene_list)
    2. Select genes sequentially from the shuffled list until we accumulate the same number of codons
       as in the actual focal region (target_codon_count)
    3. Calculate codon frequencies for this simulated "focal" region
    4. Calculate background frequencies as: total_genome_counts - simulated_focal_counts
    5. Compute cosine distance between simulated focal and background frequencies
    6. Count how many simulated distances >= observed distance for empirical P-value
    
    This approach tests whether the observed focal region's codon usage is significantly different
    from what would be expected if we randomly selected the same amount of coding sequence from
    anywhere in the genome.
    """
    # Set random seed for reproducible results
    random.seed(42)
    np.random.seed(42)
    
    cod_order = []
    foc_cod_freqs = []
    bkg_cod_freqs = []
    for cod in sorted(all_cods):
        cod_order.append(cod)
        foc_cod_freqs.append(cod_freq_dict_focal[cod])
        bkg_cod_freqs.append(cod_freq_dict_background[cod])

    foc_cod_freqs = np.array(foc_cod_freqs, dtype=np.float64)
    bkg_cod_freqs = np.array(bkg_cod_freqs, dtype=np.float64)
    
    try:
        assert len(set(foc_cod_freqs)) > 1
        assert len(set(bkg_cod_freqs)) > 1
    except Exception:
        sys.stderr.write('Error: The focal region and/or background genome appear to have no variation in codon usage.\n')
        sys.stderr.write('Focal region codon frequencies: %s\n' % ', '.join([str(x) for x in foc_cod_freqs]))
        sys.stderr.write('Background genome codon frequencies: %s\n' % ', '.join([str(x) for x in bkg_cod_freqs]))
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    # compute stats (align with reference repo behavior)
    rho, cosine_distance = ["NA"]*2
    try:
        # Spearman: return 0.0 if either vector is constant
        if np.std(foc_cod_freqs) == 0 or np.std(bkg_cod_freqs) == 0:
            rho = 0.0
        else:
            rho, _ = stats.spearmanr(foc_cod_freqs, bkg_cod_freqs)
        # Cosine: add small epsilon to avoid division-by-zero artifacts
        eps = 1e-10
        cosine_distance = spatial.distance.cosine(foc_cod_freqs + eps, bkg_cod_freqs + eps)
    except Exception:
        sys.stderr.write('Issues with computing stats!\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    # Prepare data for simulation
    # Convert defaultdict to regular dict for consistency
    gene_codons_dict = {gene: dict(codons) for gene, codons in gene_codons.items()}
    
    # Build codon index map for consistent vector ordering
    cod_to_idx = {cod: i for i, cod in enumerate(cod_order)}
    
    # Precompute per-gene count arrays and totals aligned to cod_order
    gene_counts_arrays = []
    gene_total_codons_list = []
    for gene in gene_list:
        counts = np.zeros(len(cod_order), dtype=np.float64)
        gc = gene_codons_dict[gene]
        for c, cnt in gc.items():
            idx = cod_to_idx.get(c)
            if idx is not None:
                counts[idx] = float(cnt)
        gene_counts_arrays.append(counts)
        gene_total_codons_list.append(int(np.sum(counts)))
    gene_total_codons_arr = np.array(gene_total_codons_list, dtype=np.int64)
    
    # Total codon counts vector for background computation
    # IMPORTANT: This should include ALL genes (focal + background) to match v1.2.1 behavior
    # where simulation background = total_genome_counts - simulated_focal_counts
    total_codon_counts_vec = np.zeros(len(cod_order), dtype=np.float64)
    for i, cod in enumerate(cod_order):
        total_codon_counts_vec[i] = float(all_codon_counts[cod])
    
    # Target number of codons to sample equals focal codon count
    target_codon_count = int(foc_codon_count)
    
    emp_pval = 0
    sim_cosine_distances = []
    
    # Sequential simulation processing
    if verbose:
        sys.stderr.write('Running simulation sequentially...\n')
        util.printProgressBar(0, 10, prefix = 'Progress:', suffix = 'Complete', length = 50)
    
    # Run sequential simulation (vectorized)
    for sim in range(num_sims):
        gene_indices = np.arange(len(gene_counts_arrays))
        np.random.shuffle(gene_indices)
        
        sim_focal_vec = np.zeros(len(cod_order), dtype=np.float64)
        total_collected = 0
        for idx in gene_indices:
            cds_total = gene_total_codons_arr[idx]
            if total_collected + cds_total <= target_codon_count:
                sim_focal_vec += gene_counts_arrays[idx]
                total_collected += cds_total
            else:
                remaining = target_codon_count - total_collected
                if remaining > 0 and cds_total > 0:
                    prop = remaining / float(cds_total)
                    sim_focal_vec += gene_counts_arrays[idx] * prop
                break
        # Compute cosine distance on counts vectors (focal vs background)
        sim_background_vec = total_codon_counts_vec - sim_focal_vec
        # add small epsilon to both to match reference behavior
        eps = 1e-10
        sim_cosine_distance = spatial.distance.cosine(sim_focal_vec + eps, sim_background_vec + eps)
        
        sim_cosine_distances.append(sim_cosine_distance)
        if sim_cosine_distance >= cosine_distance:
            emp_pval += 1
        
        if verbose and num_sims >= 1000 and (sim + 1) % 1000 == 0:
            util.printProgressBar(math.floor((sim + 1)/max(1, num_sims//10)), 10, prefix = 'Progress:', suffix = 'Complete', length = 50)

    emp_pval_freq = (emp_pval+1)/(num_sims+1)

    if outfile == 'stdout':
        sys.stdout.write('Empirical P-value\t%f\n' % emp_pval_freq)
        sys.stdout.write('Cosine Distance\t%f\n' % round(cosine_distance, 3))
        sys.stdout.write('Spearman\'s Rho\t%f\n' % round(rho, 3))
        sys.stdout.write('Codons\t%s\n' % ', '.join(cod_order))
        sys.stdout.write('Focal Region Codon Frequencies\t%s\n' % ', '.join([str(int(x)) for x in foc_cod_freqs]))
        sys.stdout.write('Background Genome Codon Frequencies\t%s\n' % ', '.join([str(int(x)) for x in bkg_cod_freqs]))
    elif outfile is not None:
        try:
            outfile = os.path.abspath(outfile)
            output_handle = open(outfile, 'w')
            output_handle.write('Empirical P-value\t%f\n' % emp_pval_freq)
            output_handle.write('Cosine Distance\t%f\n' % round(cosine_distance, 3))
            output_handle.write('Spearman\'s Rho\t%f\n' % round(rho, 3))
            output_handle.write('Codons\t%s\n' % ', '.join(cod_order))
            output_handle.write('Focal Region Codon Frequencies\t%s\n' % ', '.join([str(int(x)) for x in foc_cod_freqs]))
            output_handle.write('Background Genome Codon Frequencies\t%s\n' % ', '.join([str(int(x)) for x in bkg_cod_freqs]))
            output_handle.close()
        except Exception:
            sys.stderr.write("Issue creating output file!\n")
            sys.stderr.write(traceback.format_exc() + '\n')
            sys.exit(1)    
    
    if plot_outfile is not None:           
        plot_outfile = os.path.abspath(plot_outfile)

        # Create a fresh figure for this plot and close it after saving
        sns.set_theme(style='white')
        plt.figure()
        ax = sns.histplot(sim_cosine_distances, color='grey')
        ax.axvline(cosine_distance, color='#398cbf')
        plt.xlabel('Simulated cosine distances')
        plt.tight_layout()
        plt.savefig(plot_outfile, format='svg')
        plt.close()

    result = {'emp_pval_freq': emp_pval_freq, 
              'cosine_distance': float(cosine_distance), 
              'rho': float(rho), 
              'codon_order': cod_order, 
              'focal_region_codons': [int(x) for x in list(foc_cod_freqs)], 
              'background_genome_codons': [int(x) for x in list(bkg_cod_freqs)]}
    return(result)
