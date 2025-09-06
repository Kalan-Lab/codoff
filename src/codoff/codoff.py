"""Codon usage analysis and statistical testing module."""

import gzip
import math
import os
import random
import sys
import traceback
import warnings
from collections import defaultdict, Counter
from multiprocessing import Pool
from operator import itemgetter
from typing import Any, Dict, List, Optional, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from scipy import stats, spatial

from . import util

warnings.filterwarnings("ignore", category=DeprecationWarning)

# Type definitions
ValidBases = Set[str]
CodonList = List[str]
CodonCounts = Dict[str, int]
CodonFrequencies = List[float]
SimulationResults = Tuple[int, List[float]]

# Constants
VALID_BASES: ValidBases = {'A', 'C', 'G', 'T'}

# Global caches for performance optimization
_codon_extraction_cache: Dict[str, CodonList] = {}
_codon_counts_cache: Dict[Tuple[str, Tuple[str, ...]], List[int]] = {}


def _get_cached_codons(sequence: str, locus_tag: str) -> CodonList:
    """Extract codons from sequence with caching.

    Args:
        sequence: DNA sequence string
        locus_tag: unique identifier for the sequence

    Returns:
        List of valid codons (no ambiguous nucleotides)
    """
    # Return cached result if available
    if locus_tag in _codon_extraction_cache:
        return _codon_extraction_cache[locus_tag]

    try:
        # Fast validation and codon extraction
        if len(sequence) % 3 != 0:
            _codon_extraction_cache[locus_tag] = []
            return []

        # Vectorized codon extraction with validation
        codons = []
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i + 3]
            if len(codon) == 3 and all(base in VALID_BASES for base in codon):
                codons.append(codon)

        # Cache the result
        _codon_extraction_cache[locus_tag] = codons
        return codons

    except Exception:
        # Cache empty result for failed extractions
        _codon_extraction_cache[locus_tag] = []
        return []


def _get_cached_codon_counts(sequence: str, locus_tag: str, codon_order: List[str]) -> List[int]:
    """Get codon count array for a sequence.

    Args:
        sequence: DNA sequence string
        locus_tag: unique identifier for the sequence
        codon_order: list of codons in the desired order

    Returns:
        List of codon counts matching codon_order
    """
    # Check if already cached
    cache_key = (locus_tag, tuple(codon_order))
    if cache_key in _codon_counts_cache:
        return _codon_counts_cache[cache_key]

    # Get the codon list (this will populate _codon_extraction_cache if needed)
    codons = _get_cached_codons(sequence, locus_tag)

    if not codons:
        # Cache empty array for failed extractions
        empty_array = [0] * len(codon_order)
        _codon_counts_cache[cache_key] = empty_array
        return empty_array

    # Count codons and create array in codon_order
    codon_counter = Counter(codons)
    codon_counts_array = [codon_counter.get(codon, 0) for codon in codon_order]

    # Cache the result
    _codon_counts_cache[cache_key] = codon_counts_array
    return codon_counts_array


def _clear_caches() -> None:
    """Clear all caches to free memory."""
    global _codon_extraction_cache, _codon_counts_cache
    _codon_extraction_cache.clear()
    _codon_counts_cache.clear()


def _get_cache_stats() -> Dict[str, int]:
    """Get cache statistics for monitoring performance."""
    global _codon_extraction_cache, _codon_counts_cache
    return {
        'codon_extraction_cache_size': len(_codon_extraction_cache),
        'codon_counts_cache_size': len(_codon_counts_cache),
        'total_cached_sequences': len(_codon_extraction_cache)
    }


def _codoff_worker_function(
    cpu_simulations: int,
    region_freqs_list: List[List[float]],
    codon_order: List[str],
    observed_cosine_distance: float,
    all_cds_codon_count_arrays: List[List[int]],
    random_seed: Optional[int] = None
) -> SimulationResults:
    """Codoff simulation worker function for parallel processing.

        cpu_simulations: number of simulations to run on this CPU
        region_freqs_list: region codon frequencies (for target count)
        codon_order: list of codons in order
        observed_cosine_distance: the observed cosine distance to compare against
        all_cds_codon_count_arrays: list of codon count arrays for all CDS features

    Returns:
        tuple of (empirical_count, simulation_distances)
    """
    try:
        import numpy as np
        from scipy import spatial
        import random
    except ImportError:
        return 0, []

    if cpu_simulations <= 0 or not all_cds_codon_count_arrays or not codon_order:
        return 0, []

    # Set random seed for reproducible results
    if random_seed is not None:
        np.random.seed(random_seed)
        random.seed(random_seed)

    emp_count = 0
    sim_cosine_distances: List[float] = []

    # Calculate target codon count from region
    target_codon_count = int(sum(region_freqs_list))

    # Pre-compute total codon counts
    total_codon_counts = np.sum(all_cds_codon_count_arrays, axis=0, dtype=np.float64)

    # Pre-compute individual CDS total codon counts
    cds_total_codons = np.sum(all_cds_codon_count_arrays, axis=1)
    num_cds = len(all_cds_codon_count_arrays)
    cds_indices = np.arange(num_cds)

    for sim_num in range(cpu_simulations):
        try:
            # Shuffle CDS indices
            np.random.shuffle(cds_indices)

            # Initialize simulation focal frequencies
            sim_focal_freqs = np.zeros(len(codon_order), dtype=np.float64)
            total_codons_collected = 0

            # Sample codons from shuffled CDS until target count reached
            for cds_idx in cds_indices:
                cds_codon_count = cds_total_codons[cds_idx]

                if total_codons_collected + cds_codon_count <= target_codon_count:
                    # Take entire CDS
                    sim_focal_freqs += all_cds_codon_count_arrays[cds_idx]
                    total_codons_collected += cds_codon_count
                else:
                    # Partial CDS needed - sample proportionally
                    remaining_codons = target_codon_count - total_codons_collected
                    if remaining_codons > 0 and cds_codon_count > 0:
                        proportion = remaining_codons / cds_codon_count
                        sim_focal_freqs += all_cds_codon_count_arrays[cds_idx] * proportion
                    break

            if np.sum(sim_focal_freqs) == 0:
                sim_cosine_distances.append(1.0)
                continue

            # Calculate background frequencies
            sim_background_freqs = total_codon_counts - sim_focal_freqs

            # Compute cosine distance with proper zero-vector handling
            try:
                # Check for zero vectors to avoid scipy warnings
                if np.sum(sim_focal_freqs) == 0 or np.sum(sim_background_freqs) == 0:
                    sim_distance = 1.0
                else:
                    sim_distance = spatial.distance.cosine(sim_focal_freqs, sim_background_freqs)
                    if not np.isfinite(sim_distance):
                        sim_distance = 1.0
            except Exception:
                sim_distance = 1.0

            sim_cosine_distances.append(sim_distance)

            if sim_distance >= observed_cosine_distance:
                emp_count += 1

        except Exception:
            sim_cosine_distances.append(1.0)
            continue

    return emp_count, sim_cosine_distances


def _run_serial_simulation(
    foc_cod_freqs: List[float],
    all_gene_codon_count_arrays: List[List[int]],
    cod_order: List[str],
    cosine_distance: float,
    foc_codon_count: int,
    gene_list: List[str],
    gene_codons: Dict[str, List[str]],
    verbose: bool
) -> SimulationResults:
    """Run serial simulation for Monte Carlo analysis."""
    emp_pval: int = 0
    if verbose:
        util.print_progress_bar(0, 10, prefix='Progress:', suffix='Complete', length=50)
    sim_cosine_distances: List[float] = []

    # Pre-compute total codon counts for efficiency
    total_codon_counts = np.sum(all_gene_codon_count_arrays, axis=0, dtype=np.float64)

    for sim in range(0, 10000):
        if sim % 1000 == 0 and sim > 0 and verbose:
            util.print_progress_bar(math.floor(sim / 1000) + 1, 10, prefix='Progress:', suffix='Complete', length=50)

        # Shuffle gene indices instead of gene list
        random.shuffle(gene_list)
        limit_hit = False
        cod_count = 0
        foc_codon_counts = np.zeros(len(cod_order), dtype=np.float64)

        for g in gene_list:
            if not limit_hit:
                # Get pre-computed codon count array for this gene
                sequence = ''.join(gene_codons[g])
                gene_codon_count_array = _get_cached_codon_counts(sequence, g, cod_order)
                gene_total_codons = np.sum(gene_codon_count_array)

                if cod_count + gene_total_codons <= foc_codon_count:
                    # Take entire gene
                    foc_codon_counts += np.array(gene_codon_count_array)
                    cod_count += gene_total_codons
                else:
                    # Partial gene needed - sample proportionally
                    remaining_codons = foc_codon_count - cod_count
                    if remaining_codons > 0 and gene_total_codons > 0:
                        proportion = remaining_codons / gene_total_codons
                        foc_codon_counts += np.array(gene_codon_count_array) * proportion
                    break
            else:
                break

        # Calculate background frequencies
        bg_cod_freqs = total_codon_counts - foc_codon_counts

        # Compute cosine distance with proper zero-vector handling
        try:
            # Check for zero vectors to avoid scipy warnings
            if np.sum(foc_codon_counts) == 0 or np.sum(bg_cod_freqs) == 0:
                sim_cosine_distance = 1.0
            else:
                sim_cosine_distance = spatial.distance.cosine(foc_codon_counts, bg_cod_freqs)
                if not np.isfinite(sim_cosine_distance):
                    sim_cosine_distance = 1.0
        except Exception:
            sim_cosine_distance = 1.0

        sim_cosine_distances.append(sim_cosine_distance)
        if sim_cosine_distance >= cosine_distance:
            emp_pval += 1

    return emp_pval, sim_cosine_distances


def check_data_type(variable, expected_type) -> bool:
    """
    Check if the variable is of the expected type.

    Args:
        variable: The variable to check
        expected_type: The expected type of the variable

    Returns:
        True if the variable is of the expected type, False otherwise
    """
    return isinstance(variable, expected_type)


def codoff_main_gbk(
    full_genome_file: str,
    focal_genbank_files: List[str],
    outfile: Optional[str] = None,
    plot_outfile: Optional[str] = None,
    verbose: bool = True,
    threads: int = 1,
    random_seed: Optional[int] = None
) -> Optional[Dict[str, Any]]:
    """
    Analyze codon usage differences between focal regions and background genome using GenBank files.

    A full genome and a specific region must each be provided in GenBank format, with locus_tags
    overlapping. locus_tags in the focal region GenBank that are not in the full genome GenBank
    will be ignored.

    It calls the private function _stat_calc_and_simulation() to perform the main statistical
    calculations and simulations used for inference of the empirical P-value.

    Args:
        full_genome_file: Path to the full genome file in GenBank format
        focal_genbank_files: List of paths to the GenBank files for the focal region
        outfile: Path to the output file. If None, output will be printed to stdout
        plot_outfile: Path to the plot output file. If None, no plot will be made
        verbose: Whether to print progress messages to stderr. Default is True
        threads: Number of threads to use for parallel processing. Default is 1
        random_seed: Random seed for reproducible results. If None, results will be non-deterministic

    Returns:
        Dictionary containing analysis results with keys: 'emp_pval_freq', 'cosine_distance',
        'rho', 'codon_order', 'focal_region_codons', 'background_genome_codons', or None if error
    """

    try:
        assert (check_data_type(full_genome_file, str))
        assert (os.path.isfile(full_genome_file))
        assert (util.check_is_genbank_with_cds(full_genome_file))
    except Exception:
        sys.stderr.write('The full genome file must be a string to a file.\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    try:
        assert (check_data_type(focal_genbank_files, list))
        for reg_gbk in focal_genbank_files:
            assert (check_data_type(reg_gbk, str))
            assert (os.path.isfile(reg_gbk))
            assert (util.check_is_genbank_with_cds(reg_gbk))
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
        assert (check_data_type(verbose, bool))
    except Exception:
        sys.stderr.write('The verbose flag must be a boolean.\n')
        sys.stderr.write('Verbose: %s\n' % str(verbose))
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    # Clear caches at the start of analysis
    _clear_caches()

    cod_freq_dict_focal = defaultdict(int)
    cod_freq_dict_background = defaultdict(int)
    all_cods = set([])
    all_codon_counts = defaultdict(int)
    gene_codons = defaultdict(list)
    gene_list = []
    foc_codon_count = 0
    try:
        # parse GenBank files of focal regions for locus_tags
        focal_lts = set([])
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
                    assert (locus_tag is not None)
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
                assert (locus_tag is not None)
                all_coords, start, end, direction, is_multi_part = util.parse_cds_coord(str(feature.location))
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
                    if locus_tag in focal_lts:
                        focal_cds_length += len(nucl_seq)
                    total_cds_length += len(nucl_seq)
        ofgbk.close()

        if total_cds_length == 0:
            if verbose:
                sys.stderr.write('Error: The genome appears to have no CDS features. Please check the input.\n')
            sys.exit(1)

        if focal_cds_length == 0:
            if verbose:
                sys.stderr.write('Error: The focal region appears to have no CDS features. This might be because locus_tags in the focal region(s) GenBank file do not match locus_tags in the full genome GenBank file. Please check the input.\n')
            sys.exit(1)

        size_comparison = focal_cds_length / total_cds_length
        if size_comparison >= 0.05:
            sys.stderr.write('Error: The size of the focal region is >5%% of the full genome. This is not be appropriate for codoff analysis.\n')
            sys.stderr.write('Focal region(s): %s\n' % str(focal_genbank_files))
            sys.exit(1)

        # get codon frequencies for CDS in BGC and background genome using cached approach
        for locus_tag in locus_tag_sequences:
            gene_list.append(locus_tag)
            seq = locus_tag_sequences[locus_tag]

            # Use cached codon extraction
            codons = _get_cached_codons(str(seq), locus_tag)

            if not codons:
                if verbose:
                    sys.stderr.write("The locus tag %s is ignored because it was not of length 3.\n" % locus_tag)
                continue

            # Store codons for this gene
            gene_codons[locus_tag] = codons

            # Count codons efficiently
            for cod in codons:
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

    # run rest of codoff (separate function to avoid redundancy between codoff_main_gbk() and codoff_main_coords())
    result = _stat_calc_and_simulation(all_cods, cod_freq_dict_focal, cod_freq_dict_background, gene_list, gene_codons, foc_codon_count, all_codon_counts, outfile=outfile, plot_outfile=plot_outfile, verbose=verbose, threads=threads, random_seed=random_seed)
    return result


def codoff_main_coords(
    full_genome_file: str,
    focal_scaffold: str,
    focal_start_coord: int,
    focal_end_coord: int,
    outfile: Optional[str] = None,
    plot_outfile: Optional[str] = None,
    verbose: bool = True,
    threads: int = 1,
    random_seed: Optional[int] = None
) -> Optional[Dict[str, Any]]:
    """
    Analyze codon usage differences between focal regions and background genome using coordinates.

    A full genome file can be provided in either GenBank or FASTA format. If the latter,
    pyrodigal is used for gene calling, so it only works for bacteria.

    Coordinates provided by users for the focal region of interest are used to partition
    which locus tags for CDS features belong to the focal region and which belong to the
    background genome.

    It calls the private function _stat_calc_and_simulation() to perform the main statistical
    calculations and simulations used for inference of the empirical P-value.

    Args:
        full_genome_file: Path to the full genome file in GenBank or FASTA format
        focal_scaffold: Scaffold identifier for the focal region
        focal_start_coord: Start coordinate for the focal region
        focal_end_coord: End coordinate for the focal region
        outfile: Path to the output file. If None, output will be printed to stdout
        plot_outfile: Path to the plot output file. If None, no plot will be made
        verbose: Whether to print progress messages to stderr. Default is True
        threads: Number of threads to use for parallel processing. Default is 1
        random_seed: Random seed for reproducible results. If None, results will be non-deterministic

    Returns:
        Dictionary containing analysis results with keys: 'emp_pval_freq', 'cosine_distance',
        'rho', 'codon_order', 'focal_region_codons', 'background_genome_codons', or None if error
    """

    try:
        assert (check_data_type(full_genome_file, str))
        assert (os.path.isfile(full_genome_file))
    except Exception:
        sys.stderr.write('The full genome file %s must be a string path to an existing file.\n' % full_genome_file)
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    try:
        assert (check_data_type(focal_scaffold, str))
        assert (check_data_type(focal_start_coord, int))
        assert (check_data_type(focal_end_coord, int))
        assert (focal_start_coord < focal_end_coord)
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
        assert (check_data_type(verbose, bool))
    except Exception:
        sys.stderr.write('The verbose flag must be a boolean.\n')
        sys.stderr.write('Verbose: %s\n' % str(verbose))
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    # Clear caches at the start of analysis
    _clear_caches()

    cod_freq_dict_focal = defaultdict(int)
    cod_freq_dict_background = defaultdict(int)
    all_cods = set([])
    all_codon_counts = defaultdict(int)
    gene_codons = defaultdict(list)
    gene_list = []
    foc_codon_count = 0

    if util.check_is_genbank_with_cds(full_genome_file):
        try:
            # parse nucleotide coding sequences from full-genome Genbank
            locus_tag_sequences = {}
            focal_lts = set([])
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
                    assert (locus_tag is not None)
                    all_coords, start, end, direction, is_multi_part = util.parse_cds_coord(str(feature.location))
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
                    if rec.id == focal_scaffold and start >= focal_start_coord and end <= focal_end_coord:
                        focal_lts.add(locus_tag)

                    if len(str(nucl_seq)) % 3 == 0:
                        if locus_tag in focal_lts:
                            focal_cds_length += len(nucl_seq)
                        total_cds_length += len(nucl_seq)
            ofgbk.close()

            if total_cds_length == 0:
                sys.stderr.write('Error: The genome appears to have no CDS features. Please check the input.\n')
                sys.exit(1)

            if focal_cds_length == 0:
                sys.stderr.write('Error: The focal region appears to have no CDS features. Please check the input.\n')
                sys.exit(1)

            size_comparison = focal_cds_length / total_cds_length
            if size_comparison >= 0.05:
                sys.stderr.write('Error: The size of the focal region is >5%% of the full genome. This is not be appropriate for codoff analysis\n')
                sys.stderr.write('Focal region size: %d\n' % focal_cds_length)
                sys.stderr.write('Genome size: %d\n' % total_cds_length)
                sys.exit(1)

            # get codon frequencies for CDS in focal and background genome using cached approach
            for locus_tag in locus_tag_sequences:
                gene_list.append(locus_tag)
                seq = locus_tag_sequences[locus_tag]

                # Use cached codon extraction
                codons = _get_cached_codons(str(seq), locus_tag)

                if not codons:
                    sys.stderr.write("The locus tag %s is ignored because it was not of length 3.\n" % locus_tag)
                    continue

                # Store codons for this gene
                gene_codons[locus_tag] = codons

                # Count codons efficiently
                for cod in codons:
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
                assert (util.confirm_fasta(full_genome_file))
            except Exception:
                sys.stderr.write('Input genome does not appear to be either GenBank with CDS features nor a FASTA file.\n')
                sys.stderr.write(traceback.format_exc() + '\n')
                sys.exit(1)

            # parse nucleotide coding sequences from full-genome Genbank
            locus_tag_sequences, focal_lts = util.gene_call_using_pyrodigal(full_genome_file, focal_scaffold, focal_start_coord, focal_end_coord)

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

            size_comparison = focal_cds_length / total_cds_length
            if size_comparison >= 0.05:
                sys.stderr.write('Error: The size of the focal region is >5%% of the full genome. This is not be appropriate for codoff analysis.\n')
                sys.stderr.write('Focal region size: %s\n' % str(focal_cds_length))
                sys.stderr.write('Genome size: %s\n' % str(total_cds_length))
                sys.exit(1)

            # get codon frequencies for CDS in focal and background genome using cached approach
            for locus_tag in locus_tag_sequences:
                gene_list.append(locus_tag)
                seq = locus_tag_sequences[locus_tag]

                # Use cached codon extraction
                codons = _get_cached_codons(str(seq), locus_tag)

                if not codons:
                    if verbose:
                        sys.stderr.write("The locus tag %s is ignored because it was not of length 3.\n" % locus_tag)
                    continue

                # Store codons for this gene
                gene_codons[locus_tag] = codons

                # Count codons efficiently
                for cod in codons:
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

        # run rest of codoff (separate function to avoid redundancy between codoff_main_gbk() and codoff_main_coords())
        result = _stat_calc_and_simulation(all_cods, cod_freq_dict_focal, cod_freq_dict_background, gene_list, gene_codons, foc_codon_count, all_codon_counts, outfile=outfile, plot_outfile=plot_outfile, verbose=verbose, threads=threads, random_seed=random_seed)
        return result


def _stat_calc_and_simulation(
    all_cods: Set[str],
    cod_freq_dict_focal: CodonCounts,
    cod_freq_dict_background: CodonCounts,
    gene_list: List[str],
    gene_codons: Dict[str, List[str]],
    foc_codon_count: int,
    all_codon_counts: CodonCounts,
    outfile: Optional[str] = None,
    plot_outfile: Optional[str] = None,
    verbose: bool = True,
    threads: int = 1,
    random_seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Perform main statistical calculations and Monte Carlo simulations for codon usage analysis.

    This function calculates cosine distance and Spearman correlation between focal region
    and background genome codon usage, then performs Monte Carlo simulations to determine
    empirical p-values.

    Args:
        all_cods: Set of all unique codons found in the analysis
        cod_freq_dict_focal: Dictionary of codon frequencies in focal region
        cod_freq_dict_background: Dictionary of codon frequencies in background genome
        gene_list: List of gene identifiers
        gene_codons: Dictionary mapping gene IDs to their codon lists
        foc_codon_count: Total number of codons in focal region
        all_codon_counts: Dictionary of total codon counts across all genes
        outfile: Path to output file for results. If None, prints to stdout
        plot_outfile: Path to output file for plots. If None, no plot is generated
        verbose: Whether to print progress messages. Default is True
        threads: Number of threads for parallel processing. Default is 1
        random_seed: Random seed for reproducible results. If None, results will be non-deterministic

    Returns:
        Dictionary containing analysis results with keys: 'emp_pval_freq', 'cosine_distance',
        'rho', 'codon_order', 'focal_region_codons', 'background_genome_codons'
    """
    # Set random seed for reproducible results
    if random_seed is not None:
        np.random.seed(random_seed)
        random.seed(random_seed)

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
        assert (len(set(foc_cod_freqs)) > 1)
        assert (len(set(bkg_cod_freqs)) > 1)
    except Exception:
        sys.stderr.write('Error: The focal region and/or background genome appear to have no variation in codon usage.\n')
        sys.stderr.write('Focal region codon frequencies: %s\n' % ', '.join([str(x) for x in foc_cod_freqs]))
        sys.stderr.write('Background genome codon frequencies: %s\n' % ', '.join([str(x) for x in bkg_cod_freqs]))
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    # compute stats
    rho, cosine_distance = ["NA"] * 2
    try:
        # Check for constant arrays to avoid warnings
        if len(set(foc_cod_freqs)) == 1 or len(set(bkg_cod_freqs)) == 1:
            rho = 0.0  # No correlation for constant arrays
        else:
            rho, _ = stats.spearmanr(foc_cod_freqs, bkg_cod_freqs)

        # Compute cosine distance with proper zero-vector handling
        if np.sum(foc_cod_freqs) == 0 or np.sum(bkg_cod_freqs) == 0:
            cosine_distance = 1.0
        else:
            cosine_distance = spatial.distance.cosine(foc_cod_freqs, bkg_cod_freqs)
            if not np.isfinite(cosine_distance):
                cosine_distance = 1.0
    except Exception:
        sys.stderr.write('Issues with computing stats!\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)

    # Pre-compute codon count arrays for all genes for efficient simulation
    if verbose:
        sys.stderr.write("Pre-computing codon count arrays for efficient simulation...\n")

    all_gene_codon_count_arrays = []
    for g in gene_list:
        # Get codon count array using cached approach
        # We need to reconstruct the sequence from gene_codons for the cache key
        # This is a bit inefficient but maintains compatibility
        sequence = ''.join(gene_codons[g])
        codon_count_array = _get_cached_codon_counts(sequence, g, cod_order)
        if np.sum(codon_count_array) > 0:  # Only include genes with codons
            all_gene_codon_count_arrays.append(codon_count_array)

    if verbose:
        sys.stderr.write("Pre-computed %d gene codon count arrays\n" % len(all_gene_codon_count_arrays))
        cache_stats = _get_cache_stats()
        sys.stderr.write("Cache stats: %d sequences cached, %d count arrays cached\n" %
                         (cache_stats['codon_extraction_cache_size'], cache_stats['codon_counts_cache_size']))

    # Use parallel processing if threads > 1
    if threads > 1 and len(all_gene_codon_count_arrays) > 0:
        if verbose:
            sys.stderr.write("Using %d threads for parallel simulation...\n" % threads)

        # Prepare worker arguments
        num_simulations = 10000
        simulations_per_cpu = max(1, num_simulations // threads)
        remaining_simulations = num_simulations % threads

        worker_args = []
        for cpu_id in range(threads):
            cpu_simulations = simulations_per_cpu + (1 if cpu_id < remaining_simulations else 0)
            if cpu_simulations > 0:
                worker_args.append([
                    cpu_simulations,
                    foc_cod_freqs.tolist(),  # Convert numpy arrays to lists for pickling
                    cod_order,
                    cosine_distance,
                    all_gene_codon_count_arrays,
                    random_seed
                ])

        if worker_args:
            try:
                # Run parallel simulations
                with Pool(processes=threads) as pool:
                    worker_results = pool.starmap(_codoff_worker_function, worker_args)

                # Combine results
                emp_pval: int = 0
                sim_cosine_distances: List[float] = []
                for emp_count, sim_distances in worker_results:
                    emp_pval += emp_count
                    sim_cosine_distances.extend(sim_distances)

            except Exception as e:
                if verbose:
                    sys.stderr.write("Parallel processing failed: %s, falling back to serial processing\n" % str(e))
                # Fall back to serial processing
                emp_pval, sim_cosine_distances = _run_serial_simulation(foc_cod_freqs, all_gene_codon_count_arrays, cod_order, cosine_distance, foc_codon_count, gene_list, gene_codons, verbose)
        else:
            emp_pval, sim_cosine_distances = _run_serial_simulation(foc_cod_freqs, all_gene_codon_count_arrays, cod_order, cosine_distance, foc_codon_count, gene_list, gene_codons, verbose)
    else:
        # Serial processing
        emp_pval, sim_cosine_distances = _run_serial_simulation(foc_cod_freqs, all_gene_codon_count_arrays, cod_order, cosine_distance, foc_codon_count, gene_list, gene_codons, verbose)

    emp_pval_freq = (emp_pval + 1) / 10001

    if outfile == 'stdout':
        sys.stdout.write('Empirical P-value\t % f\n' % emp_pval_freq)
        sys.stdout.write('Cosine Distance\t % f\n' % round(cosine_distance, 3))
        sys.stdout.write('Spearman\'s Rho\t % f\n' % round(rho, 3))
        sys.stdout.write('Codons\t % s\n' % ', '.join(cod_order))
        sys.stdout.write('Focal Region Codon Frequencies\t % s\n' % ', '.join([str(int(x)) for x in foc_cod_freqs]))
        sys.stdout.write('Background Genome Codon Frequencies\t % s\n' % ', '.join([str(int(x)) for x in bkg_cod_freqs]))
    elif outfile is not None:
        try:
            outfile = os.path.abspath(outfile)
            output_handle = open(outfile, 'w')
            output_handle.write('Empirical P-value\t % f\n' % emp_pval_freq)
            output_handle.write('Cosine Distance\t % f\n' % round(cosine_distance, 3))
            output_handle.write('Spearman\'s Rho\t % f\n' % round(rho, 3))
            output_handle.write('Codons\t % s\n' % ', '.join(cod_order))
            output_handle.write('Focal Region Codon Frequencies\t % s\n' % ', '.join([str(int(x)) for x in foc_cod_freqs]))
            output_handle.write('Background Genome Codon Frequencies\t % s\n' % ', '.join([str(int(x)) for x in bkg_cod_freqs]))
            output_handle.close()
        except Exception:
            sys.stderr.write("Issue creating output file!\n")
            sys.stderr.write(traceback.format_exc() + '\n')
            sys.exit(1)

    if plot_outfile is not None:
        plot_outfile = os.path.abspath(plot_outfile)

        sns.set_theme(style='white')
        p = sns.histplot(sim_cosine_distances, color='grey')
        p.axvline(cosine_distance)
        plt.xlabel('Simulated cosine distances')
        plt.savefig(plot_outfile, format='svg')

    result = {'emp_pval_freq': emp_pval_freq,
              'cosine_distance': float(cosine_distance),
              'rho': float(rho),
              'codon_order': cod_order,
              'focal_region_codons': [int(x) for x in list(foc_cod_freqs)],
              'background_genome_codons': [int(x) for x in list(bkg_cod_freqs)]}
    return result
