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
from typing import Dict, Any, Optional, List, Set, Tuple

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, AfterPosition
from scipy import stats, spatial

from codoff import util

# Constants
VALID_BASES = {'A', 'C', 'G', 'T'}

def is_partial_feature(feature) -> bool:
    """
    Check if a CDS feature has fuzzy/partial boundaries.
    
    Parameters
    ----------
    feature : Bio.SeqFeature.SeqFeature
        The feature to check for partial boundaries
    
    Returns
    -------
    bool
        True if the feature has partial boundaries (< or > in GenBank notation), 
        False otherwise
    
    Notes
    -----
    Checks both the main location and all parts for compound locations (joins).
    Partial features are typically indicated by BeforePosition or AfterPosition
    in Biopython's location objects.
    """
    # Check the main location
    if isinstance(feature.location.start, (BeforePosition,)) or isinstance(feature.location.end, (AfterPosition,)):
        return True
    
    # Check all parts for compound locations (joins)
    if hasattr(feature.location, 'parts'):
        for part in feature.location.parts:
            if isinstance(part.start, (BeforePosition,)) or isinstance(part.end, (AfterPosition,)):
                return True
    
    return False

def check_data_type(variable: Any, expected_type: type) -> bool:
    """
    Check if a variable matches the expected type.
    
    Parameters
    ----------
    variable : Any
        The variable to check
    expected_type : type
        The expected type (e.g., str, int, list, dict)
    
    Returns
    -------
    bool
        True if the variable is of the expected type, False otherwise
    
    Notes
    -----
    Wrapper around isinstance() for consistent type checking throughout the codebase.
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
        - 'gene_codons': dict mapping locus_tag to codon counts (dict of dicts)
        - 'gene_list': list of all locus_tags
        - 'all_cods': set of all codons found
        - 'all_codon_counts': dict mapping codon to total genome-wide count
        - 'total_cds_length': total CDS length in the genome (int)
        - 'gene_coords': dict mapping locus_tag to (scaffold_name, start, end) tuple
        - 'scaffold_lengths': dict mapping scaffold_name to scaffold length
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
    gene_coords = {}
    scaffold_lengths = {}
    
    try:
        ofgbk = None
        if full_genome_file.endswith('.gz'):
            ofgbk = gzip.open(full_genome_file, 'rt')
        else:
            ofgbk = open(full_genome_file)
        
        for rec in SeqIO.parse(ofgbk, 'genbank'):
            scaffold_name = rec.id
            scaffold_lengths[rec.id] = len(rec.seq)
            for feature in rec.features:
                if not feature.type == 'CDS': 
                    continue
                
                # Skip partial CDS features (those with fuzzy boundaries like < or >)
                if is_partial_feature(feature):
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
                
                all_coords, start, end, direction, is_multi_part = util.parseCDSCoord(feature.location)
                gene_coords[locus_tag] = (scaffold_name, start, end)

                # Use Biopython's robust .extract() method instead of manual slicing
                nucl_seq = str(feature.extract(rec.seq))
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
        'total_cds_length': total_cds_length,
        'gene_coords': gene_coords,
        'scaffold_lengths': scaffold_lengths
    }


def process_bgc_with_cached_data(genome_data: Dict[str, Any], focal_genbank_files: List[str], outfile: Optional[str] = None, plot_outfile: Optional[str] = None, verbose: bool = True, num_sims: int = 10000, seed: Optional[int] = None, max_focal_cds_fraction: float = 0.05) -> Dict[str, Any]:
    """
    Process BGC regions using pre-computed genome data to avoid redundant computation.
    
    Parameters
    ----------
    genome_data : dict
        Pre-computed genome data from extract_genome_codon_data(). Should contain:
        - 'locus_tag_sequences', 'gene_codons', 'gene_list', 'all_cods',
        - 'all_codon_counts', 'total_cds_length', 'gene_coords', 'scaffold_lengths'
    focal_genbank_files : list
        The list of paths to the GenBank files for the focal region.
    outfile : str, optional
        The path to the output file. If not provided, the output will be 
        printed to stdout.
    plot_outfile : str, optional
        The path to the plot output file. If not provided, no plot will 
        be made.
    verbose : bool, optional
        Whether to print progress messages to stderr. Default is True.
    num_sims : int, optional
        Number of simulations to run, by default 10000
    seed : int, optional
        Random seed for reproducible results, by default None
    max_focal_cds_fraction : float, optional
        Maximum allowed fraction of total genome CDS length for focal region, by default 0.05 (5%)
    
    Returns
    -------
    dict
        Result dictionary from _stat_calc_and_simulation() containing:
        - 'empirical_freq', 'cosine_distance', 'rho', 'codon_order',
        - 'focal_region_codons', 'background_genome_codons'
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
    gene_coords = genome_data['gene_coords']
    scaffold_lengths = genome_data['scaffold_lengths']

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
        if size_comparison >= max_focal_cds_fraction:
            sys.stderr.write(f'Error: The focal region aggregate CDS length is >{max_focal_cds_fraction*100:.1f}%% of the total genome aggregate CDS length. This is not appropriate for codoff analysis.\n')
            sys.stderr.write('Focal region(s): %s\n' % str(focal_genbank_files))
            sys.stderr.write('Focal region aggregate CDS length: %d bp\n' % focal_cds_length)
            sys.stderr.write('Total genome aggregate CDS length: %d bp\n' % total_cds_length)
            sys.exit(1)

    except Exception:
        sys.stderr.write('Issues with determining codon usage info for inputs with focal region and genome provided as GenBank.\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    # Determine focal region coordinates for region size
    min_start = float('inf')
    max_end = 0
    for lt in focal_lts:
        if lt in gene_coords:
            _, s, e = gene_coords[lt]
            min_start = min(min_start, s)
            max_end = max(max_end, e)
    region_size = max_end - min_start if max_end > min_start else None

    # Run rest of codoff using pre-computed data
    result = _stat_calc_and_simulation(all_cods, cod_freq_dict_focal, cod_freq_dict_background, gene_list, 
                                       gene_codons, foc_codon_count, all_codon_counts, gene_coords=gene_coords, 
                                       outfile=outfile, plot_outfile=plot_outfile, verbose=verbose,  
                                       num_sims=num_sims, focal_region_size=region_size, 
                                       scaffold_lengths=scaffold_lengths, seed=seed, max_focal_cds_fraction=max_focal_cds_fraction)
    
    # Final writing of results
    if outfile and outfile != 'stdout':
        _write_output_file(result, outfile)
    elif outfile == 'stdout':
        _print_output_stdout(result)

    return result


def codoff_main_gbk(full_genome_file: str, focal_genbank_files: List[str], outfile: Optional[str] = None, plot_outfile: Optional[str] = None, verbose: bool = True, genome_data: Optional[Dict[str, Any]] = None, num_sims: int = 10000, seed: Optional[int] = None, max_focal_cds_fraction: float = 0.05) -> Dict[str, Any]:
    """
    Compare codon usage between a focal genomic region and full genome using GenBank files.
    
    A full genome and a specific region must each be provided in GenBank format, 
    with locus_tags overlapping. Locus tags in the focal region GenBank that are 
    not in the full genome GenBank will be ignored. Uses sequential contiguous-window 
    sampling (v1.2.3+) to perform Monte Carlo simulations.
    
    Parameters
    ----------
    full_genome_file : str
        Path to the full genome file in GenBank format
    focal_genbank_files : List[str]
        List of paths to the GenBank files for the focal region(s)
    outfile : str, optional
        Path to the output file. If not provided, output is printed to stdout
    plot_outfile : str, optional
        Path to the plot output file in SVG format. If not provided, no plot is made
    verbose : bool, optional
        Whether to print progress messages to stderr, by default True
    genome_data : Dict[str, Any], optional
        Pre-computed genome data from extract_genome_codon_data() to avoid 
        redundant computation. If provided, full_genome_file will be ignored.
        Must contain coordinate information (gene_coords, scaffold_lengths)
    num_sims : int, optional
        Number of simulations to run, by default 10000
    seed : int, optional
        Random seed for reproducible results, by default None
    max_focal_cds_fraction : float, optional
        Maximum allowed fraction of total genome CDS length for focal region, by default 0.05 (5%)
    
    Returns
    -------
    Dict[str, Any]
        Results dictionary containing:
        - 'empirical_freq': proportion of simulations as or more extreme than observed
        - 'cosine_distance': cosine distance between focal and background
        - 'rho': Spearman's correlation coefficient
        - 'codon_order': list of codons in order
        - 'focal_region_codons': codon frequencies in focal region
        - 'background_genome_codons': codon frequencies in background
    """
    
    # If genome_data is provided, use it instead of processing the full_genome_file
    if genome_data is not None:
        result = process_bgc_with_cached_data(genome_data, focal_genbank_files, outfile, plot_outfile, verbose, num_sims, seed=seed, max_focal_cds_fraction=max_focal_cds_fraction)
        return result
    
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
    gene_coords = {}
    scaffold_lengths = {}
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
            sys.stderr.write('Unable to find full-genome GenBank file %s.\n' % full_genome_file)
            sys.exit(1)

        ofgbk = None
        if full_genome_file.endswith('.gz'):
            ofgbk = gzip.open(full_genome_file, 'rt')
        else:
            ofgbk = open(full_genome_file)
        for rec in SeqIO.parse(ofgbk, 'genbank'):
            full_sequence = str(rec.seq)
            scaffold_name = rec.id
            scaffold_lengths[rec.id] = len(rec.seq)
            for feature in rec.features:
                if not feature.type == 'CDS': continue
                
                # Skip partial CDS features (those with fuzzy boundaries like < or >)
                if is_partial_feature(feature):
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
                all_coords, start, end, direction, is_multi_part = util.parseCDSCoord(feature.location)
                gene_coords[locus_tag] = (scaffold_name, start, end)
                if end >= len(full_sequence): end = len(full_sequence)
                
                # Use Biopython's robust .extract() method instead of manual slicing
                nucl_seq = str(feature.extract(rec.seq))
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
        if size_comparison >= max_focal_cds_fraction:
            sys.stderr.write(f'Error: The focal region aggregate CDS length is >{max_focal_cds_fraction*100:.1f}%% of the total genome aggregate CDS length. This is not appropriate for codoff analysis.\n')
            sys.stderr.write('Focal region(s): %s\n' % str(focal_genbank_files))
            sys.stderr.write('Focal region aggregate CDS length: %d bp\n' % focal_cds_length)
            sys.stderr.write('Total genome aggregate CDS length: %d bp\n' % total_cds_length)
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

    # Determine region size
    min_start = float('inf')
    max_end = 0
    try:
        for foc_gbk in focal_genbank_files:
            for rec in SeqIO.parse(gzip.open(foc_gbk, 'rt') if foc_gbk.endswith('.gz') else open(foc_gbk), 'genbank'):
                for feature in rec.features:
                    if feature.type == 'CDS':
                        _, start, end, _, _ = util.parseCDSCoord(feature.location)
                        min_start = min(min_start, start)
                        max_end = max(max_end, end)
        region_size = max_end - min_start if max_end > min_start else None
    except Exception as e:
        sys.stderr.write(f"Warning: Could not determine region size: {e}\n")
        region_size = None
 
    # run rest of codoff (separate function to avoid redundancy between codoff_main_gbk() and codoff_main_coords()
    # Use all genes for simulation pool to match v1.2.1 behavior
    result = _stat_calc_and_simulation(
        all_cods, cod_freq_dict_focal, cod_freq_dict_background, gene_list, gene_codons, 
        foc_codon_count, all_codon_counts, gene_coords=gene_coords, 
        outfile=outfile, plot_outfile=plot_outfile, verbose=verbose,  num_sims=num_sims,
        focal_region_size=region_size, scaffold_lengths=scaffold_lengths,
        seed=seed, max_focal_cds_fraction=max_focal_cds_fraction
    )
    
    # Final writing of results
    if outfile and outfile != 'stdout':
        _write_output_file(result, outfile)
    elif outfile == 'stdout':
        _print_output_stdout(result)

    return result


def codoff_main_coords(full_genome_file: str, focal_scaffold: str, focal_start_coord: int, focal_end_coord: int, outfile: Optional[str] = None, plot_outfile: Optional[str] = None, verbose: bool = True, num_sims: int = 10000, seed: Optional[int] = None, max_focal_cds_fraction: float = 0.05) -> Dict[str, Any]:
    """
    Compare codon usage between a focal region and full genome using coordinates.
    
    A full genome file can be provided in either GenBank or FASTA format. If FASTA,
    pyrodigal is used for gene calling (bacteria only). Coordinates provided by users
    determine which CDS features belong to the focal region versus the background genome.
    Uses sequential contiguous-window sampling (v1.2.3+) for Monte Carlo simulations.
    
    Parameters
    ----------
    full_genome_file : str
        Path to the full genome file in GenBank or FASTA format
    focal_scaffold : str
        Scaffold identifier for the focal region
    focal_start_coord : int
        Start coordinate for the focal region (1-based)
    focal_end_coord : int
        End coordinate for the focal region (inclusive)
    outfile : str, optional
        Path to the output file. If not provided, output is printed to stdout
    plot_outfile : str, optional
        Path to the plot output file in SVG format. If not provided, no plot is made
    verbose : bool, optional
        Whether to print progress messages to stderr, by default True
    num_sims : int, optional
        Number of simulations to run, by default 10000
    seed : int, optional
        Random seed for reproducible results, by default None
    max_focal_cds_fraction : float, optional
        Maximum allowed fraction of total genome CDS length for focal region, by default 0.05 (5%)
    
    Returns
    -------
    Dict[str, Any]
        Results dictionary containing:
        - 'empirical_freq': proportion of simulations as or more extreme than observed
        - 'cosine_distance': cosine distance between focal and background
        - 'rho': Spearman's correlation coefficient
        - 'codon_order': list of codons in order
        - 'focal_region_codons': codon frequencies in focal region
        - 'background_genome_codons': codon frequencies in background
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

    # Check assembly quality - N50 must be larger than focal region size
    focal_region_size = focal_end_coord - focal_start_coord
    util.check_assembly_quality(full_genome_file, focal_region_size, verbose=verbose)

    cod_freq_dict_focal = defaultdict(int)
    cod_freq_dict_background = defaultdict(int)
    all_cods = set()
    all_codon_counts = defaultdict(int)
    gene_codons = defaultdict(lambda: defaultdict(int))
    gene_list = []
    foc_codon_count = 0
    gene_coords = {}
    scaffold_lengths = {}

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
                scaffold_name = rec.id
                scaffold_lengths[rec.id] = len(rec.seq)
                for feature in rec.features:
                    if not feature.type == 'CDS': continue
                    
                    # Skip partial CDS features (those with fuzzy boundaries like < or >)
                    if is_partial_feature(feature):
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
                    all_coords, start, end, direction, is_multi_part = util.parseCDSCoord(feature.location)
                    gene_coords[locus_tag] = (scaffold_name, start, end)
                    if end >= len(full_sequence): end = len(full_sequence)
                    
                    # Use Biopython's robust .extract() method instead of manual slicing
                    nucl_seq = str(feature.extract(rec.seq))
                    locus_tag_sequences[locus_tag] = nucl_seq
                    
                    # Correctly identify focal genes based on coordinates
                    if rec.id == focal_scaffold and start >= focal_start_coord and end <= focal_end_coord:
                        focal_lts.add(locus_tag)

                    if len(str(nucl_seq)) % 3 == 0:
                        total_cds_length += len(nucl_seq)
            ofgbk.close()

            # Recalculate focal_cds_length based on identified focal_lts
            for lt in focal_lts:
                if len(locus_tag_sequences.get(lt, '')) % 3 == 0:
                    focal_cds_length += len(locus_tag_sequences[lt])
            
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
            if size_comparison >= max_focal_cds_fraction:
                sys.stderr.write(f'Error: The focal region aggregate CDS length is >{max_focal_cds_fraction*100:.1f}%% of the total genome aggregate CDS length. This is not appropriate for codoff analysis\n')
                sys.stderr.write('Focal region aggregate CDS length: %d bp\n' % focal_cds_length)
                sys.stderr.write('Total genome aggregate CDS length: %d bp\n' % total_cds_length)
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
        
        # run rest of codoff
        result = _stat_calc_and_simulation(
            all_cods, cod_freq_dict_focal, cod_freq_dict_background, gene_list, gene_codons, 
            foc_codon_count, all_codon_counts, gene_coords=gene_coords, 
            outfile=outfile, plot_outfile=plot_outfile, verbose=verbose,  num_sims=num_sims,
            focal_region_size=(focal_end_coord - focal_start_coord),
            scaffold_lengths=scaffold_lengths,
            seed=seed, max_focal_cds_fraction=max_focal_cds_fraction
        )
        
        # Final writing of results
        if outfile and outfile != 'stdout':
            _write_output_file(result, outfile)
        elif outfile == 'stdout':
            _print_output_stdout(result)

        return result

    else:
        try:
            try:
                assert util.confirmFasta(full_genome_file)
            except Exception:
                sys.stderr.write('Input genome does not appear to be either GenBank with CDS features nor a FASTA file.\n')
                sys.stderr.write(traceback.format_exc() + '\n')
                sys.exit(1)            
        
            # parse nucleotide coding sequences from full-genome Genbank
            locus_tag_sequences, focal_lts, gene_coords = util.geneCallUsingPyrodigal(full_genome_file, focal_scaffold, focal_start_coord, focal_end_coord)
            
            # Extract scaffold lengths from FASTA file
            with open(full_genome_file) as off:
                for rec in SeqIO.parse(off, 'fasta'):
                    scaffold_lengths[rec.id] = len(rec.seq)
            
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
            if size_comparison >= max_focal_cds_fraction:
                sys.stderr.write(f'Error: The focal region aggregate CDS length is >{max_focal_cds_fraction*100:.1f}%% of the total genome aggregate CDS length. This is not appropriate for codoff analysis.\n')
                sys.stderr.write('Focal region aggregate CDS length: %d bp\n' % focal_cds_length)
                sys.stderr.write('Total genome aggregate CDS length: %d bp\n' % total_cds_length)
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
        result = _stat_calc_and_simulation(
            all_cods, cod_freq_dict_focal, cod_freq_dict_background, gene_list, gene_codons, 
            foc_codon_count, all_codon_counts, gene_coords=gene_coords, 
            outfile=outfile, plot_outfile=plot_outfile, verbose=verbose,  num_sims=num_sims,
            focal_region_size=(focal_end_coord - focal_start_coord),
            scaffold_lengths=scaffold_lengths,
            seed=seed, max_focal_cds_fraction=max_focal_cds_fraction
        )
        
        # Final writing of results
        if outfile and outfile != 'stdout':
            _write_output_file(result, outfile)
        elif outfile == 'stdout':
            _print_output_stdout(result)

        return result
    
    
def _print_output_stdout(result: Dict[str, Any]) -> None:
    """
    Print analysis results to stdout in tab-separated format.
    
    Parameters
    ----------
    result : Dict[str, Any]
        Results dictionary from _stat_calc_and_simulation()
    """
    percentile = result['empirical_freq'] * 100.0
    sys.stdout.write(f"Discordance Percentile\t{round(percentile, 2)}\n")
    sys.stdout.write(f"Cosine Distance\t{round(result['cosine_distance'], 3)}\n")
    sys.stdout.write(f"Spearman's Rho\t{round(result['rho'], 3)}\n")
    sys.stdout.write(f"Codons\t{', '.join(result['codon_order'])}\n")
    sys.stdout.write(f"Focal Region Codon Frequencies\t{', '.join(map(str, result['focal_region_codons']))}\n")
    sys.stdout.write(f"Background Genome Codon Frequencies\t{', '.join(map(str, result['background_genome_codons']))}\n")


def _write_output_file(result: Dict[str, Any], outfile: str) -> None:
    """
    Write analysis results to a file in tab-separated format.
    
    Parameters
    ----------
    result : Dict[str, Any]
        Results dictionary from _stat_calc_and_simulation()
    outfile : str
        Path to the output file
    """
    try:
        with open(outfile, 'w') as f:
            percentile = result['empirical_freq'] * 100.0
            f.write(f"Discordance Percentile\t{round(percentile, 2)}\n")
            f.write(f"Cosine Distance\t{round(result['cosine_distance'], 3)}\n")
            f.write(f"Spearman's Rho\t{round(result['rho'], 3)}\n")
            f.write(f"Codons\t{', '.join(result['codon_order'])}\n")
            f.write(f"Focal Region Codon Frequencies\t{', '.join(map(str, result['focal_region_codons']))}\n")
            f.write(f"Background Genome Codon Frequencies\t{', '.join(map(str, result['background_genome_codons']))}\n")
    except Exception as e:
        sys.stderr.write(f"Issue creating output file: {e}\n")
        sys.exit(1)


def _stat_calc_and_simulation(all_cods: Set[str], cod_freq_dict_focal: Dict[str, int], cod_freq_dict_background: Dict[str, int], gene_list: List[str], gene_codons: Dict[str, Dict[str, int]], foc_codon_count: int, all_codon_counts: Dict[str, int], outfile: Optional[str] = None, plot_outfile: Optional[str] = None, verbose: bool = True, num_sims: int = 10000, gene_coords: Optional[Dict[str, Tuple[str, int, int]]] = None, focal_region_size: Optional[int] = None, scaffold_lengths: Optional[Dict[str, int]] = None, seed: Optional[int] = None, max_focal_cds_fraction: float = 0.05) -> Dict[str, Any]:
    """
    Private function that performs the main statistical calculations and Monte Carlo simulations.
    
    Parameters
    ----------
    all_cods : Set[str]
        Set of all codons found in the genome
    cod_freq_dict_focal : Dict[str, int]
        Codon frequency counts for the focal region
    cod_freq_dict_background : Dict[str, int]
        Codon frequency counts for the background genome
    gene_list : List[str]
        List of all gene locus tags
    gene_codons : Dict[str, Dict[str, int]]
        Dictionary mapping locus tags to their codon counts
    foc_codon_count : int
        Total number of codons in the focal region
    all_codon_counts : Dict[str, int]
        Genome-wide codon counts
    outfile : str, optional
        Output file path
    plot_outfile : str, optional
        Plot output file path
    verbose : bool, optional
        Whether to print progress messages
    num_sims : int, optional
        Number of simulations to run
    gene_coords : Dict[str, Tuple[str, int, int]], optional
        Dictionary mapping locus tags to (scaffold, start, end)
    focal_region_size : int, optional
        Size of focal region in base pairs
    scaffold_lengths : Dict[str, int], optional
        Dictionary mapping scaffold names to their lengths
    seed : int, optional
        Random seed for reproducibility
    max_focal_cds_fraction : float, optional
        Maximum allowed fraction of total genome CDS length for simulated windows, by default 0.05 (5%)
    
    Returns
    -------
    Dict[str, Any]
        Results dictionary containing empirical_freq, cosine_distance, rho, etc.
    
    Notes
    -----
    SEQUENTIAL CONTIGUOUS-WINDOW SAMPLING:
    - Randomly select genomic windows of the same size as the focal region
    - Extract all genes fully contained within each window
    - Calculate codon frequencies for each simulated window
    - Compare focal region against these random windows
    - Count how many simulated distances >= observed distance
    - Reports: Discordance Percentile
    
    This approach tests whether the observed focal region's codon usage is significantly 
    different from what would be expected for randomly positioned contiguous genomic 
    regions of the same size.
    
    Requires focal_region_size, gene_coords, and scaffold_lengths to be provided. 
    If these are not available, the function will exit with an error.
    """
    # Set random seed for reproducible results
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
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
    gene_codons_dict = {gene: dict(codons) for gene, codons in gene_codons.items()}
    cod_to_idx = {cod: i for i, cod in enumerate(cod_order)}
    
    gene_counts_arrays = []
    valid_gene_list = []
    for gene in gene_list:
        # Skip genes without codon data (e.g., genes filtered out during processing)
        if gene not in gene_codons_dict:
            continue
        counts = np.zeros(len(cod_order), dtype=np.float64)
        gc = gene_codons_dict[gene]
        for c, cnt in gc.items():
            idx = cod_to_idx.get(c)
            if idx is not None:
                counts[idx] = float(cnt)
        gene_counts_arrays.append(counts)
        valid_gene_list.append(gene)
    
    # Use valid_gene_list for simulations
    gene_list = valid_gene_list
    
    total_codon_counts_vec = np.zeros(len(cod_order), dtype=np.float64)
    for i, cod in enumerate(cod_order):
        total_codon_counts_vec[i] = float(all_codon_counts[cod])
    
    # Setup ordered genes for sequential sampling
    sorted_genes = []
    if gene_coords:
        # Only include genes that have coordinate data
        sorted_genes = sorted([g for g in gene_list if g in gene_coords], 
                             key=lambda g: (gene_coords[g][0], gene_coords[g][1]))
    
    # Simulation containers
    sim_cosine_distances: List[float] = []
    extreme_count = 0  # Count of simulations with distance >= observed

    if verbose:
        sys.stderr.write('Running simulation sequentially...\n')
        util.printProgressBar(0, 10, prefix = 'Progress:', suffix = 'Complete', length = 50)

    # Verify we have the required data for contiguous-window sampling
    if focal_region_size is None or not gene_coords or len(sorted_genes) == 0 or not scaffold_lengths:
        sys.stderr.write('Error: Sequential contiguous-window sampling requires genomic coordinates.\n')
        sys.stderr.write('Please ensure your GenBank files contain proper CDS feature locations.\n')
        sys.exit(1)
    
    # Filter scaffolds large enough for the region size
    valid_scaffolds = [(scaf_id, scaf_len) for scaf_id, scaf_len in scaffold_lengths.items() if scaf_len > focal_region_size]
    
    if not valid_scaffolds:
        sys.stderr.write(f'Error: No scaffolds are large enough to contain the focal region size ({focal_region_size} bp).\n')
        sys.exit(1)
    
    # Calculate total genome codons for hyper-dense filtering
    total_genome_codons = np.sum(total_codon_counts_vec)
    
    # Change from 'for' loop to 'while' loop to skip empty windows
    while len(sim_cosine_distances) < num_sims:
        # --- UNIFIED SAMPLING LOGIC ---
        # Pick a random window
        start_scaffold, scaffold_len = random.choice(valid_scaffolds)
        start_coord = random.randint(1, scaffold_len - focal_region_size)
        sim_region_end = start_coord + focal_region_size
        # ---

        sim_focal_vec = np.zeros(len(cod_order), dtype=np.float64)

        # Now, find all genes *fully contained* in this random window
        for gene_locus_tag in sorted_genes:
            scaffold, gene_start, gene_end = gene_coords[gene_locus_tag]

            if scaffold != start_scaffold:
                continue

            # Check if gene is fully contained within the [start_coord, sim_region_end] window
            if gene_start >= start_coord and gene_end <= sim_region_end:
                original_gene_index = gene_list.index(gene_locus_tag)
                sim_focal_vec += gene_counts_arrays[original_gene_index]

            # Optimization: if we're past the window, stop looping over this scaffold
            if gene_start > sim_region_end:
                break

        # Skip empty windows (no genes found) to match observed data filtering
        sim_focal_codons = np.sum(sim_focal_vec)
        if sim_focal_codons == 0:
            continue
        
        # Skip hyper-dense windows to match observed data filtering
        if total_genome_codons > 0:
            size_comparison = sim_focal_codons / total_genome_codons
            if size_comparison >= max_focal_cds_fraction:
                continue

        sim_background_vec = total_codon_counts_vec - sim_focal_vec
        eps = 1e-10
        sim_dist = spatial.distance.cosine(sim_focal_vec + eps, sim_background_vec + eps)
        sim_cosine_distances.append(sim_dist)
        if sim_dist >= cosine_distance:
            extreme_count += 1

        # Update progress bar using current simulation count
        current_sim_count = len(sim_cosine_distances)
        if verbose and num_sims >= 1000 and current_sim_count % 1000 == 0:
            util.printProgressBar(math.floor(current_sim_count / max(1, num_sims // 10)), 10, prefix = 'Progress:', suffix = 'Complete', length = 50)

    # Calculate empirical frequency (with +1 smoothing)
    # This represents the proportion of simulations as or more extreme than observed
    empirical_freq = (extreme_count + 1) / (num_sims + 1)

    result = {
        'empirical_freq': empirical_freq,
        'cosine_distance': float(cosine_distance),
        'rho': float(rho),
        'codon_order': cod_order,
        'focal_region_codons': [int(x) for x in list(foc_cod_freqs)],
        'background_genome_codons': [int(x) for x in list(bkg_cod_freqs)],
    }
    
    if plot_outfile is not None:
        plot_outfile = os.path.abspath(plot_outfile)
        sns.set_theme(style='white')
        plt.figure()
        sns.histplot(sim_cosine_distances, color='grey', alpha=0.5, label='Null distribution')
        plt.axvline(cosine_distance, color='#2d5f8b', label='Observed distance')
        plt.xlabel('Simulated cosine distances')
        plt.legend()
        plt.tight_layout()
        plt.savefig(plot_outfile, format='svg')
        plt.close()

    return result
