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

valid_bases = set(['A', 'C', 'G', 'T'])
	
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

def codoff_main_gbk(full_genome_file, focal_genbank_files, outfile=None, plot_outfile=None, verbose=True):
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
	"""
	
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
		if foc_gbk.endswith('.gz'):
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
				gene_codons[locus_tag].append(cod)
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
	result = _stat_calc_and_simulation(all_cods, cod_freq_dict_focal, cod_freq_dict_background, gene_list, gene_codons, foc_codon_count, all_codon_counts, outfile=outfile, plot_outfile=plot_outfile, verbose=verbose)
	return result

def codoff_main_coords(full_genome_file, focal_scaffold, focal_start_coord, focal_end_coord, outfile=None, plot_outfile=None, verbose=True):
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
	gene_codons = defaultdict(list)
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
			if foc_gbk.endswith('.gz'):
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
					gene_codons[locus_tag].append(cod)
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
					gene_codons[locus_tag].append(cod)
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
		result = _stat_calc_and_simulation(all_cods, cod_freq_dict_focal, cod_freq_dict_background, gene_list, gene_codons, foc_codon_count, all_codon_counts, outfile=outfile, plot_outfile=plot_outfile, verbose=verbose)
		return result
	
def _stat_calc_and_simulation(all_cods, cod_freq_dict_focal, cod_freq_dict_background, gene_list, gene_codons, foc_codon_count, all_codon_counts, outfile=None, plot_outfile=None, verbose=True):
	"""
	private function that performs the main statistical calculations and 
	"""
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

	emp_pval = 0
	if verbose:
		util.printProgressBar(0, 10, prefix = 'Progress:', suffix = 'Complete', length = 50)
	sim_cosine_distances = []
	for sim in range(0, 10000):
		if sim % 1000 == 0 and sim > 0 and verbose:
			util.printProgressBar(math.floor(sim/1000) + 1, 10, prefix = 'Progress:', suffix = 'Complete', length = 50)
		random.shuffle(gene_list)
		limit_hit = False
		cod_count = 0
		foc_codon_counts = defaultdict(int)
		for g in gene_list:
			if not limit_hit:
				for c in gene_codons[g]:
					foc_codon_counts[c] += 1
					cod_count += 1
					if cod_count >= foc_codon_count:
						limit_hit = True
						break
			else:
				break

		foc_cod_freqs = []
		bg_cod_freqs = []
		for cod in cod_order:
			foc_cod_freqs.append(foc_codon_counts[cod])
			bg_cod_freqs.append(all_codon_counts[cod] - foc_codon_counts[cod])
		foc_cod_freqs = numpy.array(foc_cod_freqs, dtype=numpy.float64)
		bg_cod_freqs = numpy.array(bg_cod_freqs, dtype=numpy.float64)
		sim_cosine_distance = spatial.distance.cosine(foc_cod_freqs, bg_cod_freqs)
		sim_cosine_distances.append(sim_cosine_distance)
		if sim_cosine_distance >= cosine_distance:
			emp_pval += 1

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
