import os
import sys
import traceback
from Bio import SeqIO
from Bio.Seq import Seq
import pyrodigal

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
	"""
	Function from: https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters
    Call in a loop to create terminal progress bar
	@params:
        iteration   - Required  : current iteration (Int)
		total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    
	percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
	filledLength = int(length * iteration // total)
	bar = fill * filledLength + '-' * (length - filledLength)
	sys.stderr.write(f'\r{prefix} |{bar}| {percent}% {suffix}')
    # Print New Line on Complete
	if iteration == total:        
		sys.stderr.write('\n')

def parseCDSCoord(str_gbk_loc):
	try:
			start = None
			end = None
			direction = None
			all_coords = []
			is_multi_part = False
			if not 'join' in str(str_gbk_loc) and not 'order' in str(str_gbk_loc):
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
			return(all_coords, start, end, direction, is_multi_part)
	except Exception as e:
		sys.stderr.write("Issues with GenBank CDS processing.\n")
		raise RuntimeError(traceback.format_exc())

def geneCallUsingPyrodigal(genome_fasta_file, focal_scaffold, focal_start_coord, focal_end_coord):
	locus_tag_sequences = {}
	focal_lts = set([])
	try:
		full_genome_seqs = []
		with open(genome_fasta_file) as off:
			for rec in SeqIO.parse(off, 'fasta'):
				full_genome_seqs.append(str(rec.seq))
		spacer = 'N'*10000
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
					gene_seq = Seq(scaff_seq[start_coord-1:end_coord])
					if strand == -1:
						gene_seq = str(Seq(scaff_seq[start_coord-1:end_coord]).reverse_complement())
					locus_tag_name = 'CDS_' + str(locus_tag_id)
					locus_tag_sequences[locus_tag_name] = gene_seq
					if scaffold == focal_scaffold and start_coord >= focal_start_coord and end_coord <= focal_end_coord:
						focal_lts.add(locus_tag_name)
					locus_tag_id += 1
	except:
		sys.stderr.write("Issues with pyrodigal gene calling.\n")
		raise RuntimeError(traceback.format_exc())
	return([locus_tag_sequences, focal_lts])

def checkIsGenBankWithCDS(gbk_file):
	try:
		gbk_with_cds = False
		with open(gbk_file) as ogf:
			for rec in SeqIO.parse(ogf, 'genbank'):
				for feat in rec.features:
					if feat.type == 'CDS':
						gbk_with_cds = True
						break
		return(gbk_with_cds) 
	except:
		sys.stderr.write("Issues with attemping to parse suspected GenBank file %s" % gbk_file)
		raise RuntimeError(traceback.format_exc())
	
def confirmFasta(fasta_file):
	try:
		fasta_recs = 0
		with open(fasta_file) as off:
			for rec in SeqIO.parse(off, 'fasta'):
				fasta_recs += 1
		return(fasta_recs > 0) 
	except:
		sys.stderr.write("Issues with attemping to parse suspected FASTA file %s" % fasta_file)
		raise RuntimeError(traceback.format_exc())
