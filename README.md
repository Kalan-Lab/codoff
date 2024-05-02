# codoff
codoff: Program to measure the irregularity of the codon usage for a single genomic region (e.g. a BGC, phage, etc.) relative to the full genome. It is primarily designed to work of the output of antiSMASH biosynthetic gene cluster (BGC) predictions - but can be more broadly applied as well. The `-p` argument allows users to infer an empirical p-value as to the significance of observing a codon usage profile for the focal region so discordant with what is observed for the larger genome-wide context.

Inputs are a region-specific (e.g. BGC) Genbank and a full-genome GenBank file. Both types of files are produced by antiSMASH.

codoff is largely adapted from the `compareBGCtoGenomeCodonUsage.py` script in the lsaBGC suite and presented here separately just to make it easy to install via bioconda, therefore, please cite:

```
Evolutionary investigations of the biosynthetic diversity in the skin microbiome using lsaBGC
Rauf Salamzade, J.Z. Alex Cheong, Shelby Sandstrom, Mary Hannah Swaney, Reed M. Stubbendieck,
Nicole Lane Starr, Cameron R. Currie, Anne Marie Singh, and Lindsay R. Kalan
```

If you have suggestions or feedback, please let me know in the GitHub issues.

## Installation

### Bioconda

Coming soon!

### Conda Manual

```
# 1. clone Git repo and change directories into it!
git clone https://github.com/raufs/codoff/
cd codoff/

# 2. create conda environment using yaml file and activate it!
conda env create -f codoff_env.yml -n codoff_env
conda activate codoff_env

# 3. complete python installation with the following command:
pip install -e .
```

### Example of inferring codon usage differences for the *crt* operon in a chromosomal and plasmid context

```
# Download assembly for Staphylococcus warneri st. LK413 / p3-SID855
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/145/625/GCF_025145625.2_ASM2514562v2/GCF_025145625.2_ASM2514562v2_genomic.fna.gz
gunzip GCF_025145625.2_ASM2514562v2_genomic.fna.gz

# Run antiSMASH in environment with antiSMASH available (not a part of codoff conda environment!) or run it on
# their convenient web-server and afterwards download & uncompress the results. Link to antiSMASH webserver:
# https://antismash.secondarymetabolites.org/
antismash --genefinding-tool prodigal --output-dir Sw_LK413 --output-basename LK413 --minimal GCF_025145625.2_ASM2514562v2_genomic.fna

# run codoff for crt (BGC encoding for staphyloxanthin) found in chromosome:
codoff -b Sw_LK413/c00001_NZ_JALX...region001.gbk -g Sw_LK413/LK413.gbk -p -o crt_in_chromosome.txt

# run codoff for crt (BGC encoding for staphyloxanthin) found in plasmid:
codoff -b Sw_LK413/c00002_NZ_JALX...region001.gbk -g Sw_LK413/LK413.gbk -p -o crt_in_plasmid.txt
```

The `crt_in_chromsome.txt` should contain the following results:

```
Cosine_Distance	0.004000
Empirical_Pvalue	0.230000
Spearman_Rho	0.989000
Spearman_Pvalue	0.000000
GCF_Codons	AAA, AAC, AAG, AAT, ACA, ACC, ACG, ACT, AGA, AGC, AGG, AGT, ATA, ATC, ATG, ATT, CAA, CAC, CAG, CAT, CCA, CCC, CCG, CCT, CGA, CGC, CGG, CGT, CTA, CTC, CTG, CTT, GAA, GAC, GAG, GAT, GCA, GCC, GCG, GCT, GGA, GGC, GGG, GGT, GTA, GTC, GTG, GTT, TAA, TAC, TAG, TAT, TCA, TCC, TCG, TCT, TGA, TGC, TGG, TGT, TTA, TTC, TTG, TTT
GCF_Codon_Frequencies	768, 178, 158, 507, 381, 42, 124, 200, 146, 55, 11, 195, 236, 238, 354, 654, 522, 59, 56, 229, 250, 24, 40, 148, 56, 33, 3, 133, 124, 42, 32, 140, 629, 171, 107, 581, 420, 81, 133, 298, 192, 131, 64, 503, 323, 135, 150, 331, 34, 108, 9, 402, 271, 29, 31, 159, 3, 13, 120, 74, 692, 192, 136, 405
Background_Codon_Frequencies	43824, 9285, 8580, 31551, 20161, 2457, 5549, 12786, 9313, 3054, 703, 12462, 13197, 11998, 18709, 35717, 28244, 3360, 3013, 13872, 11636, 935, 1883, 8646, 3218, 1693, 248, 9275, 6537, 1976, 1215, 7724, 40947, 9119, 6862, 33984, 20076, 3810, 5235, 15247, 10203, 6350, 2843, 24522, 17227, 5605, 6841, 17741, 1822, 5458, 343, 21929, 14672, 1550, 1929, 10604, 228, 808, 5314, 3430, 40801, 9128, 6797, 22231
```

The `crt_in_plasmid.txt` should contain the following results:

```
Cosine_Distance	0.021000
Empirical_Pvalue	0.042000
Spearman_Rho	0.967000
Spearman_Pvalue	0.000000
GCF_Codons	AAA, AAC, AAG, AAT, ACA, ACC, ACG, ACT, AGA, AGC, AGG, AGT, ATA, ATC, ATG, ATT, CAA, CAC, CAG, CAT, CCA, CCC, CCG, CCT, CGA, CGC, CGG, CGT, CTA, CTC, CTG, CTT, GAA, GAC, GAG, GAT, GCA, GCC, GCG, GCT, GGA, GGC, GGG, GGT, GTA, GTC, GTG, GTT, TAA, TAC, TAG, TAT, TCA, TCC, TCG, TCT, TGA, TGC, TGG, TGT, TTA, TTC, TTG, TTT
GCF_Codon_Frequencies	425, 74, 102, 337, 138, 16, 47, 93, 74, 32, 17, 125, 172, 60, 167, 240, 213, 32, 34, 133, 80, 14, 17, 80, 25, 12, 1, 53, 48, 16, 18, 56, 240, 67, 75, 279, 130, 23, 51, 122, 94, 66, 14, 172, 126, 57, 68, 144, 11, 58, 5, 223, 103, 17, 14, 102, 6, 8, 60, 27, 274, 56, 75, 197
Background_Codon_Frequencies	44167, 9389, 8636, 31721, 20404, 2483, 5626, 12893, 9385, 3077, 697, 12532, 13261, 12176, 18896, 36131, 28553, 3387, 3035, 13968, 11806, 945, 1906, 8714, 3249, 1714, 250, 9355, 6613, 2002, 1229, 7808, 41336, 9223, 6894, 34286, 20366, 3868, 5317, 15423, 10301, 6415, 2893, 24853, 17424, 5683, 6923, 17928, 1845, 5508, 347, 22108, 14840, 1562, 1946, 10661, 225, 813, 5374, 3477, 41219, 9264, 6858, 22439
```

So we see the cosine distance between the BGC codon usage profile and the background genome is higher in the plasmid instance of crt (which we describe in our manuscript on lsaBGC to correspond to a more recently acquired version of GCF_6) than for the chromosome instance of crt. Further, through calculating an empirical P-value with the -p argument, we see the chances of observing as large a disagreement in the codon distribution of the plasmid BGC to the expected genome-wide distribution is <0.05 which is not the case for the chromosomal BGC.

## Algorithm for computing empirical p-value

To calculate an empirical P-value, we gather codons for each gene across the genome. First the codon frequency distribution of the focal-region/BGC genes is compared to the codon frequency distribution of the background genome (all other genes; genes which have lengths not divisible by 3 are ignored). After, we perform 10,000 simulations where in each simulation we shuffle the full genome-wide set of genes and go through the first N genes until the same number of codons as are present in the focal region/BGC/GCF are observed. The cosine distance between the observed codon frequency is compared to the remainder of the genome-wide codon distribution and checked for whether it is higher than what was actually observed for the BGC; if so, then an empirical P-value counter is appended a count of 1. The final empirical P-value produced is simply this count plus a pseudocount of 1 over 10,001.

![figure](https://github.com/Kalan-Lab/codoff/blob/main/codoff_empirical_pvalue_image.png?raw=true)

## Usage 

Usage is identical to the original `compareBGCtoGenomeCodonUsage.py` script. 

```
usage: codoff [-h] -g FULL_GENOME_GENBANK -b BGC_GENBANKS [BGC_GENBANKS ...] -o OUTPUT [-p]

	Program: codoff
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
		
	This program compares the codon-usage distribution of a BGC (provided as a Genbank) to the codon usage of the 
	background genome (gathered from a full-genome Genbank file - also an output of antiSMASH). It will report
	the cosine distance and Spearman correlation between the two profiles. Only ORFs which are of length divisible 
	by 3 will be considered. It can accept multiple BGC Genbanks in the case that the BGC is fragmented across 
	multiple scaffolds.
	
	Only works for bacterial genomes currently.
	

optional arguments:
  -h, --help            show this help message and exit
  -g FULL_GENOME_GENBANK, --full_genome_genbank FULL_GENOME_GENBANK
                        Path to annotated full-genome Genbank for isolate's genome.
  -b BGC_GENBANKS [BGC_GENBANKS ...], --bgc_genbanks BGC_GENBANKS [BGC_GENBANKS ...]
                        Path to BGC Genbank(s) for isolate. Locus tags must match with tags in full-genome Genbank.
  -o OUTPUT, --output OUTPUT
                        Path to output file.
  -p, --compute_empirical_pval
                        Compute empirical P-value for observed cosine distance being as different than expected to the background codon usage distribution.
```

## License

```
BSD 3-Clause License

Copyright (c) 2024, Kalan-Lab

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```
