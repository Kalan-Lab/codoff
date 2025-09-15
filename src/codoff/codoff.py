import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import os
import random
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from scipy import stats, spatial
from codoff import util
from operator import itemgetter
import gzip
import numpy 
import math
import traceback
import seaborn as sns
import matplotlib.pyplot as plt
import multiprocessing as mp
from functools import partial

valid_bases = set(['A', 'C', 'G', 'T'])

def _simulation_batch_worker(shared_data, batch_range):
	"""
	Batch worker function for parallel simulation execution.
	Processes multiple simulations in a single worker to reduce overhead.
	
	Parameters
	----------
	shared_data : dict
		Contains shared data: gene_list, gene_codons, foc_codon_count, bkg_cod_freqs, cod_order, cosine_distance
	batch_range : tuple
		(start_idx, end_idx) for the batch of simulations to process
	
	Returns
	-------
	list
		List of (sim_cosine_distance, exceeds_observed) tuples
	"""
	# Set random seed in each worker process for deterministic behavior
	random.seed(42)
	
	start_idx, end_idx = batch_range
	gene_list = shared_data['gene_list']
	gene_codons = shared_data['gene_codons']
	foc_codon_count = shared_data['foc_codon_count']
	bkg_cod_freqs = shared_data['bkg_cod_freqs']
	cod_order = shared_data['cod_order']
	cosine_distance = shared_data['cosine_distance']
	
	results = []
	
	for sim in range(start_idx, end_idx):
		# Create a local copy of gene_list to avoid shared state issues
		local_gene_list = gene_list.copy()
		random.shuffle(local_gene_list)
		
		limit_hit = False
		cod_count = 0
		foc_codon_counts = defaultdict(int)
		
		for g in local_gene_list:
			if not limit_hit:
				# Calculate how many codons we can take from this gene
				remaining_codons = foc_codon_count - cod_count
				if remaining_codons <= 0:
					limit_hit = True
					break
				
				# Add codons from this gene using proportional sampling
				gene_total_codons = sum(gene_codons[g].values())
				codons_to_sample = min(gene_total_codons, remaining_codons)
				
				# Calculate proportional counts for each codon type
				for c, count in gene_codons[g].items():
					# Calculate proportional count: (original_count / gene_total) * codons_to_sample
					proportional_count = round((count / gene_total_codons) * codons_to_sample)
					foc_codon_counts[c] += proportional_count
					cod_count += proportional_count
				
				remaining_codons -= codons_to_sample
				if remaining_codons <= 0:
					limit_hit = True
			else:
				break

		foc_cod_freqs = []
		for cod in cod_order:
			foc_cod_freqs.append(foc_codon_counts[cod])
		
		# Normalize to frequencies
		total_foc_codons = sum(foc_cod_freqs)
		if total_foc_codons > 0:
			foc_cod_freqs = [x/total_foc_codons for x in foc_cod_freqs]
		else:
			foc_cod_freqs = [0.0] * len(cod_order)
		
		# Calculate cosine distance against the background frequencies
		foc_cod_freqs = numpy.array(foc_cod_freqs, dtype=numpy.float64)
		sim_cosine_distance = spatial.distance.cosine(foc_cod_freqs, bkg_cod_freqs)
		
		exceeds_observed = bool(sim_cosine_distance >= cosine_distance)
		results.append((sim_cosine_distance, exceeds_observed))
	
	return results


def check_data_type(variable, expected_type):
	"""
	Check if the variable is of the expected type.
	
	Parameters
	----------
	variable : any
		The variable to check.
	expected_type : type
		The expected type of the variable.
	
	Returns
	-------
	bool
		True if the variable is of the expected type, False otherwise.
	"""
	return isinstance(variable, expected_type)


def extract_genome_codon_data(full_genome_file, verbose=True):
	"""
	Extract codon usage data from the full genome file once to avoid redundant computation.
	
	Parameters
	----------
	full_genome_file : str
		The path to the full genome file in GenBank format.
	verbose : bool
		Whether to print progress messages to stderr. Default is True.
	
	Returns
	-------
	dict
		Dictionary containing:
		- 'locus_tag_sequences': dict mapping locus_tag to nucleotide sequence
		- 'gene_codons': dict mapping locus_tag to codon counts
		- 'gene_list': list of all locus_tags
		- 'all_cods': set of all codons found
		- 'total_cds_length': total CDS length in the genome
	"""
	try:
		assert(check_data_type(full_genome_file, str))
		assert(os.path.isfile(full_genome_file))
		assert(util.checkIsGenBankWithCDS(full_genome_file))
	except:
		sys.stderr.write('The full genome file must be a string to a file.\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

	locus_tag_sequences = {}
	gene_codons = defaultdict(lambda: defaultdict(int))
	gene_list = []
	all_cods = set([])
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
				except:
					try:
						locus_tag = feature.qualifiers.get('gene')[0]
					except:
						locus_tag = feature.qualifiers.get('protein_id')[0]
				assert(locus_tag != None)
				
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
						if not(len(cod) == 3 and cod[0] in valid_bases and cod[1] in valid_bases and cod[2] in valid_bases): 
							continue
						gene_codons[locus_tag][cod] += 1
						all_cods.add(cod)
		ofgbk.close()
		
		if total_cds_length == 0:
			if verbose:
				sys.stderr.write('Error: The genome appears to have no CDS features. Please check the input.\n')
			sys.exit(1)
			
	except:
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
		'total_cds_length': total_cds_length
	}


def process_bgc_with_cached_data(genome_data, focal_genbank_files, outfile=None, plot_outfile=None, verbose=True, max_jobs=None):
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
	max_jobs : int
		Maximum number of parallel processes to use.
	
	Returns
	-------
	dict
		Result dictionary from _stat_calc_and_simulation()
	"""
	try:
		assert(check_data_type(focal_genbank_files, list))
		for reg_gbk in focal_genbank_files:
			assert(check_data_type(reg_gbk, str))
			assert(os.path.isfile(reg_gbk))
			assert(util.checkIsGenBankWithCDS(reg_gbk))
	except:
		sys.stderr.write('The focal region GenBank files must be a list of string paths to existing files.\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

	if outfile != 'stdout' and outfile != None and os.path.isfile(outfile):
		sys.stderr.write('The outfile must be a string path to a file which does not already exist.\n')
		sys.stderr.write('Outfile: %s\n' % str(outfile))
		sys.exit(1)

	if plot_outfile != 'stdout' and plot_outfile != None and os.path.isfile(plot_outfile):
		sys.stderr.write('The plot outfile must be a string path to a file which does not already exist.\n')
		sys.stderr.write('Plot outfile: %s\n' % str(plot_outfile))
		sys.exit(1)	

	try:
		assert(check_data_type(verbose, bool))
	except:
		sys.stderr.write('The verbose flag must be a boolean.\n')
		sys.stderr.write('Verbose: %s\n' % str(verbose))
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

	# Extract pre-computed data
	locus_tag_sequences = genome_data['locus_tag_sequences']
	gene_codons = genome_data['gene_codons']
	gene_list = genome_data['gene_list']
	all_cods = genome_data['all_cods']
	total_cds_length = genome_data['total_cds_length']

	cod_freq_dict_focal = defaultdict(int)
	cod_freq_dict_background = defaultdict(int)
	all_codon_counts = defaultdict(int)
	foc_codon_count = 0
	focal_cds_length = 0

	try:
		# Parse GenBank files of focal regions for locus_tags
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
					except:
						try:
							locus_tag = feature.qualifiers.get('gene')[0]
						except:
							locus_tag = feature.qualifiers.get('protein_id')[0]
					assert(locus_tag != None)
					focal_lts.add(locus_tag)
			ogbk.close()
		
		# Check for missing focal locus tags and issue warnings
		missing_focal_lts = focal_lts - set(locus_tag_sequences.keys())
		if missing_focal_lts and verbose:
			sys.stderr.write('Warning: The following focal region locus tags were not found in the full genome GenBank file:\n')
			for missing_lt in sorted(missing_focal_lts):
				sys.stderr.write('  %s\n' % missing_lt)
			sys.stderr.write('These locus tags will be ignored in the analysis.\n')
		
		# Calculate focal region statistics using pre-computed data
		for locus_tag in locus_tag_sequences:
			seq = locus_tag_sequences[locus_tag]
			if not len(str(seq)) % 3 == 0:
				if verbose:
					sys.stderr.write("The locus tag %s is ignored because it was not of length 3.\n" % locus_tag)
				continue
			
			if locus_tag in focal_lts:
				focal_cds_length += len(seq)
			
			# Count codons for this gene
			codon_seq = [str(seq)[i:i + 3] for i in range(0, len(str(seq)), 3)]
			for cod in list(codon_seq):
				if not(len(cod) == 3 and cod[0] in valid_bases and cod[1] in valid_bases and cod[2] in valid_bases): 
					continue
				all_codon_counts[cod] += 1
				if locus_tag in focal_lts:
					foc_codon_count += 1
					cod_freq_dict_focal[cod] += 1
				else:
					cod_freq_dict_background[cod] += 1

		if focal_cds_length == 0:
			if verbose:
				sys.stderr.write('Error: The focal region appears to have no CDS features. This might be because locus_tags in the focal region(s) GenBank file do not match locus_tags in the full genome GenBank file. Please check the input.\n')
			sys.exit(1)

		size_comparison = focal_cds_length/total_cds_length
		if size_comparison >= 0.05:
			sys.stderr.write('Error: The size of the focal region is >5%% of the full genome. This is not be appropriate for codoff analysis.\n')
			sys.stderr.write('Focal region(s): %s\n' % str(focal_genbank_files))
			sys.exit(1)

	except:
		sys.stderr.write('Issues with determining codon usage info for inputs with focal region and genome provided as GenBank.\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

	# Run rest of codoff using pre-computed data
	result = _stat_calc_and_simulation(all_cods, cod_freq_dict_focal, cod_freq_dict_background, gene_list, gene_codons, foc_codon_count, all_codon_counts, outfile=outfile, plot_outfile=plot_outfile, verbose=verbose, max_jobs=max_jobs)
	return result


def codoff_main_gbk(full_genome_file, focal_genbank_files, outfile=None, plot_outfile=None, verbose=True, max_jobs=None, genome_data=None):
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
	* max_jobs: int
		Maximum number of parallel processes to use.
	* genome_data: dict, optional
		Pre-computed genome data from extract_genome_codon_data() to avoid 
		redundant computation. If provided, full_genome_file will be ignored.
	"""
	
	# If genome_data is provided, use it instead of processing the full_genome_file
	if genome_data is not None:
		return process_bgc_with_cached_data(genome_data, focal_genbank_files, outfile, plot_outfile, verbose, max_jobs)
	
	try:
		assert(check_data_type(full_genome_file, str))
		assert(os.path.isfile(full_genome_file))
		assert(util.checkIsGenBankWithCDS(full_genome_file))
	except:
		sys.stderr.write('The full genome file must be a string to a file.\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

	try:
		assert(check_data_type(focal_genbank_files, list))
		for reg_gbk in focal_genbank_files:
			assert(check_data_type(reg_gbk, str))
			assert(os.path.isfile(reg_gbk))
			assert(util.checkIsGenBankWithCDS(reg_gbk))
	except:
		sys.stderr.write('The focal region GenBank files must be a list of string paths to existing files.\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)
 
	if outfile != 'stdout' and outfile != None and os.path.isfile(outfile):
		sys.stderr.write('The outfile must be a string path to a file which does not already exist.\n')
		sys.stderr.write('Outfile: %s\n' % str(outfile))
		sys.exit(1)

	if plot_outfile != 'stdout' and plot_outfile != None and os.path.isfile(plot_outfile):
		sys.stderr.write('The plot outfile must be a string path to a file which does not already exist.\n')
		sys.stderr.write('Plot outfile: %s\n' % str(plot_outfile))
		sys.exit(1)	

	try:
		assert(check_data_type(verbose, bool))
	except:
		sys.stderr.write('The verbose flag must be a boolean.\n')
		sys.stderr.write('Verbose: %s\n' % str(verbose))
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

	cod_freq_dict_focal = defaultdict(int)
	cod_freq_dict_background = defaultdict(int)
	all_cods = set([])
	all_codon_counts = defaultdict(int)
	gene_codons = defaultdict(lambda: defaultdict(int))
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
					if not feature.type == 'CDS': continue
					locus_tag = None
					try:
						locus_tag = feature.qualifiers.get('locus_tag')[0]
					except:
						try:
							locus_tag = feature.qualifiers.get('gene')[0]
						except:
							locus_tag = feature.qualifiers.get('protein_id')[0]
					assert(locus_tag != None)
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
				except:
					try:
						locus_tag = feature.qualifiers.get('gene')[0]
					except:
						locus_tag = feature.qualifiers.get('protein_id')[0]
				assert(locus_tag != None)
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

				if len(str(nucl_seq))%3 == 0:
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
			if not len(str(seq))%3 == 0:
				if verbose:
					sys.stderr.write("The locus tag %s is ignored because it was not of length 3.\n" % locus_tag)
				continue
			codon_seq = [str(seq)[i:i + 3] for i in range(0, len(str(seq)), 3)]
			for cod in list(codon_seq):
				if not(len(cod) == 3 and cod[0] in valid_bases and cod[1] in valid_bases and cod[2] in valid_bases): continue
				gene_codons[locus_tag][cod] += 1
				all_codon_counts[cod] += 1
				if locus_tag in focal_lts:
					foc_codon_count += 1
					cod_freq_dict_focal[cod] += 1
				else:
					cod_freq_dict_background[cod] += 1
				all_cods.add(cod)
	except:
		sys.stderr.write('Issues with determinining codon usage info for inputs with focal region and genome provided as GenBank.\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

	# run rest of codoff (separate function to avoid redundancy between codoff_main_gbk() and codoff_main_coords())
	result = _stat_calc_and_simulation(all_cods, cod_freq_dict_focal, cod_freq_dict_background, gene_list, gene_codons, foc_codon_count, all_codon_counts, outfile=outfile, plot_outfile=plot_outfile, verbose=verbose, max_jobs=max_jobs)
	return result


def codoff_main_coords(full_genome_file, focal_scaffold, focal_start_coord, focal_end_coord, outfile=None, plot_outfile=None, verbose=True, max_jobs=None):
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
		assert(check_data_type(full_genome_file, str))
		assert(os.path.isfile(full_genome_file))
	except:
		sys.stderr.write('The full genome file %s must be a string path to an existing file.\n' % full_genome_file)
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

	try:
		assert(check_data_type(focal_scaffold, str))
		assert(check_data_type(focal_start_coord, int))
		assert(check_data_type(focal_end_coord, int))
		assert(focal_start_coord < focal_end_coord)
	except:
		sys.stderr.write('The focal region must be specified via a string for the scaffold and two integers for the start and end coordinates.\n')
		sys.stderr.write('Scaffold: %s\n' % focal_scaffold)
		sys.stderr.write('Start coord: %s\n' % str(focal_start_coord))
		sys.stderr.write('End coord: %s\n' % str(focal_end_coord))
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)
	
	if outfile != 'stdout' and outfile != None and os.path.isfile(outfile):
		sys.stderr.write('The outfile must be a string path to a file which does not already exist.\n')
		sys.stderr.write('Outfile: %s\n' % str(outfile))
		sys.exit(1)

	if plot_outfile != 'stdout' and plot_outfile != None and os.path.isfile(plot_outfile):
		sys.stderr.write('The plot outfile must be a string path to a file which does not already exist.\n')
		sys.stderr.write('Plot outfile: %s\n' % str(plot_outfile))
		sys.exit(1)	

	try:
		assert(check_data_type(verbose, bool))
	except:
		sys.stderr.write('The verbose flag must be a boolean.\n')
		sys.stderr.write('Verbose: %s\n' % str(verbose))
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

	cod_freq_dict_focal = defaultdict(int)
	cod_freq_dict_background = defaultdict(int)
	all_cods = set([])
	all_codon_counts = defaultdict(int)
	gene_codons = defaultdict(lambda: defaultdict(int))
	gene_list = []
	foc_codon_count = 0

	if util.checkIsGenBankWithCDS(full_genome_file):
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
					if not feature.type == 'CDS': continue
					locus_tag = None
					try:
						locus_tag = feature.qualifiers.get('locus_tag')[0]
					except:
						try:
							locus_tag = feature.qualifiers.get('gene')[0]
						except:
							locus_tag = feature.qualifiers.get('protein_id')[0]
					assert(locus_tag != None)
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

					if len(str(nucl_seq))%3 == 0:
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
				if not len(str(seq))%3 == 0:
					sys.stderr.write("The locus tag %s is ignored because it was not of length 3.\n" % locus_tag)
					continue
				codon_seq = [str(seq)[i:i + 3] for i in range(0, len(str(seq)), 3)]
				for cod in list(codon_seq):
					if not(len(cod) == 3 and cod[0] in valid_bases and cod[1] in valid_bases and cod[2] in valid_bases): continue
					gene_codons[locus_tag][cod] += 1
					all_codon_counts[cod] += 1
					if locus_tag in focal_lts:
						foc_codon_count += 1
						cod_freq_dict_focal[cod] += 1
					else:
						cod_freq_dict_background[cod] += 1
					all_cods.add(cod)
		except:
			sys.stderr.write('Issues with determinining codon usage info for inputs with full genome provided as GenBank.\n')
			sys.stderr.write(traceback.format_exc() + '\n')
			sys.exit(1)

	else:
		try:
			try:
				assert(util.confirmFasta(full_genome_file))
			except:
				sys.stderr.write('Input genome does not appear to be either GenBank with CDS features nor a FASTA file.\n')
				sys.stderr.write(traceback.format_exc() + '\n')
				sys.exit(1)			
		
			# parse nucleotide coding sequences from full-genome Genbank
			locus_tag_sequences, focal_lts = util.geneCallUsingPyrodigal(full_genome_file, focal_scaffold, focal_start_coord, focal_end_coord)
			
			total_cds_length = 0
			focal_cds_length = 0
			for locus_tag in locus_tag_sequences:
				seq = locus_tag_sequences[locus_tag]
				if len(str(seq))%3 == 0:
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
				if not len(str(seq))%3 == 0:
					if verbose:
						sys.stderr.write("The locus tag %s is ignored because it was not of length 3.\n" % locus_tag)
					continue
				codon_seq = [str(seq)[i:i + 3] for i in range(0, len(str(seq)), 3)]
				for cod in list(codon_seq):
					if not(len(cod) == 3 and cod[0] in valid_bases and cod[1] in valid_bases and cod[2] in valid_bases): continue
					gene_codons[locus_tag][cod] += 1
					all_codon_counts[cod] += 1
					if locus_tag in focal_lts:
						foc_codon_count += 1
						cod_freq_dict_focal[cod] += 1
					else:
						cod_freq_dict_background[cod] += 1
					all_cods.add(cod)
		except:
			sys.stderr.write('Issues with determinining codon usage info for inputs with full genome provided as FASTA.\n')
			sys.stderr.write(traceback.format_exc() + '\n')
			sys.exit(1)
		
		# run rest of codoff (separate function to avoid redundancy between codoff_main_gbk() and codoff_main_coords())
		result = _stat_calc_and_simulation(all_cods, cod_freq_dict_focal, cod_freq_dict_background, gene_list, gene_codons, foc_codon_count, all_codon_counts, outfile=outfile, plot_outfile=plot_outfile, verbose=verbose, max_jobs=max_jobs)
		return result
	
	
def _stat_calc_and_simulation(all_cods, cod_freq_dict_focal, cod_freq_dict_background, gene_list, gene_codons, foc_codon_count, all_codon_counts, outfile=None, plot_outfile=None, verbose=True, max_jobs=None):
	"""
	private function that performs the main statistical calculations and 
	"""
	# Set random seed for reproducible results
	random.seed(42)
	
	cod_order = []
	foc_cod_freqs = []
	bkg_cod_freqs = []
	for cod in sorted(all_cods):
		cod_order.append(cod)
		foc_cod_freqs.append(cod_freq_dict_focal[cod])
		bkg_cod_freqs.append(cod_freq_dict_background[cod])

	foc_cod_freqs = numpy.array(foc_cod_freqs, dtype=numpy.float64)
	bkg_cod_freqs = numpy.array(bkg_cod_freqs, dtype=numpy.float64)
	
	try:
		assert(len(set(foc_cod_freqs)) > 1)
		assert(len(set(bkg_cod_freqs)) > 1)
	except:
		sys.stderr.write('Error: The focal region and/or background genome appear to have no variation in codon usage.\n')
		sys.stderr.write('Focal region codon frequencies: %s\n' % ', '.join([str(x) for x in foc_cod_freqs]))
		sys.stderr.write('Background genome codon frequencies: %s\n' % ', '.join([str(x) for x in bkg_cod_freqs]))
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

	# compute stats
	rho, cosine_distance = ["NA"]*2
	try:
		rho, _ = stats.spearmanr(foc_cod_freqs, bkg_cod_freqs)
		cosine_distance = spatial.distance.cosine(foc_cod_freqs, bkg_cod_freqs)
	except:
		sys.stderr.write('Issues with computing stats!\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

	# Determine number of processes to use
	if max_jobs is None:
		max_jobs = 1  # Default to sequential processing
	else:
		max_jobs = min(max_jobs, mp.cpu_count())
	
	# For small datasets, force sequential processing to avoid overhead
	# Only use parallel processing for larger datasets where it provides benefit
	if max_jobs > 1 and len(gene_list) < 500:
		if verbose:
			sys.stderr.write(f'Dataset too small for parallel processing (only {len(gene_list)} genes). Using sequential processing.\n')
		max_jobs = 1
	
	# Prepare shared data for parallel simulation (only created once)
	# Convert defaultdict to regular dict for pickling
	gene_codons_dict = {gene: dict(codons) for gene, codons in gene_codons.items()}
	
	# Create shared data dictionary (passed once to all workers)
	# Convert background counts to frequencies for worker functions
	bkg_cod_freqs = []
	for cod in cod_order:
		bkg_cod_freqs.append(cod_freq_dict_background[cod])
	bkg_cod_freqs = numpy.array(bkg_cod_freqs, dtype=numpy.float64)
	
	shared_data = {
		'gene_list': gene_list,
		'gene_codons': gene_codons_dict,
		'foc_codon_count': foc_codon_count,
		'bkg_cod_freqs': bkg_cod_freqs,
		'cod_order': cod_order,
		'cosine_distance': cosine_distance
	}
	
	emp_pval = 0
	sim_cosine_distances = []
	
	if max_jobs == 1:
		# Sequential processing - avoid multiprocessing overhead
		if verbose:
			sys.stderr.write('Running simulation sequentially...\n')
			util.printProgressBar(0, 10, prefix = 'Progress:', suffix = 'Complete', length = 50)
		
		# Run sequential simulation
		for sim in range(10000):
			# Create a local copy of gene_list to avoid shared state issues
			local_gene_list = gene_list.copy()
			random.shuffle(local_gene_list)
			
			limit_hit = False
			cod_count = 0
			foc_codon_counts = defaultdict(int)
			
			for g in local_gene_list:
				if not limit_hit:
					# Calculate how many codons we can take from this gene
					remaining_codons = foc_codon_count - cod_count
					if remaining_codons <= 0:
						limit_hit = True
						break
					
					# Add codons from this gene using proportional sampling
					gene_total_codons = sum(gene_codons_dict[g].values())
					codons_to_sample = min(gene_total_codons, remaining_codons)
					
					# Calculate proportional counts for each codon type
					for c, count in gene_codons_dict[g].items():
						# Calculate proportional count: (original_count / gene_total) * codons_to_sample
						proportional_count = round((count / gene_total_codons) * codons_to_sample)
						foc_codon_counts[c] += proportional_count
						cod_count += proportional_count
					
					remaining_codons -= codons_to_sample
					if remaining_codons <= 0:
						limit_hit = True
				else:
					break

			sim_foc_cod_freqs = []
			for cod in cod_order:
				sim_foc_cod_freqs.append(foc_codon_counts[cod])
			
			# Normalize to frequencies
			total_foc_codons = sum(sim_foc_cod_freqs)
			if total_foc_codons > 0:
				sim_foc_cod_freqs = [x/total_foc_codons for x in sim_foc_cod_freqs]
			else:
				sim_foc_cod_freqs = [0.0] * len(cod_order)
			
			# Calculate cosine distance against the background frequencies
			sim_foc_cod_freqs = numpy.array(sim_foc_cod_freqs, dtype=numpy.float64)
			sim_cosine_distance = spatial.distance.cosine(sim_foc_cod_freqs, bkg_cod_freqs)
			
			sim_cosine_distances.append(sim_cosine_distance)
			if sim_cosine_distance >= cosine_distance:
				emp_pval += 1
			
			if (sim + 1) % 1000 == 0 and verbose:
				util.printProgressBar(math.floor((sim + 1)/1000), 10, prefix = 'Progress:', suffix = 'Complete', length = 50)
	else:
		# Parallel processing with batched approach for better performance
		if verbose:
			sys.stderr.write(f'Running simulation with {max_jobs} parallel processes...\n')
			util.printProgressBar(0, 10, prefix = 'Progress:', suffix = 'Complete', length = 50)
		
		# Batch simulations to reduce multiprocessing overhead
		# Each worker will process multiple simulations
		sims_per_worker = max(1, 10000 // max_jobs)
		batches = []
		for i in range(0, 10000, sims_per_worker):
			batch_end = min(i + sims_per_worker, 10000)
			batches.append((i, batch_end))
		
		# Run parallel simulation with batched approach
		with mp.Pool(processes=max_jobs) as pool:
			# Create a partial function that includes shared_data
			from functools import partial
			worker_func = partial(_simulation_batch_worker, shared_data)
			
			# Process batches in parallel
			results = pool.map(worker_func, batches)
			
			# Collect results from all batches
			completed = 0
			for batch_results in results:
				for sim_cosine_distance, exceeds_observed in batch_results:
					sim_cosine_distances.append(sim_cosine_distance)
					if exceeds_observed:
						emp_pval += 1
					
					completed += 1
					if completed % 1000 == 0 and verbose:
						util.printProgressBar(math.floor(completed/1000), 10, prefix = 'Progress:', suffix = 'Complete', length = 50)

	emp_pval_freq = (emp_pval+1)/10001

	if outfile == 'stdout':
		sys.stdout.write('Empirical P-value\t%f\n' % emp_pval_freq)
		sys.stdout.write('Cosine Distance\t%f\n' % round(cosine_distance, 3))
		sys.stdout.write('Spearman\'s Rho\t%f\n' % round(rho, 3))
		sys.stdout.write('Codons\t%s\n' % ', '.join(cod_order))
		sys.stdout.write('Focal Region Codon Frequencies\t%s\n' % ', '.join([str(int(x)) for x in foc_cod_freqs]))
		sys.stdout.write('Background Genome Codon Frequencies\t%s\n' % ', '.join([str(int(x)) for x in bkg_cod_freqs]))
	elif outfile != None:
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
		except:
			sys.stderr.write("Issue creating output file!\n")
			sys.stderr.write(traceback.format_exc() + '\n')
			sys.exit(1)	
	
	if plot_outfile != None:           
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
	return(result)
