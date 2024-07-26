# codoff

**codoff**: A program to measure the irregularity of the codon usage for a focal genomic region (e.g. a BGC, phage, etc.) relative to the full genome. It was primarily designed to work off the output of antiSMASH biosynthetic gene cluster (BGC) predictions - but can be more broadly applied as well, with options allowing users to provide a genome in FASTA format and simply specify the coordinates of the focal region of interest. It computes an empirical P-value as to the significance of observing a codon usage profile for the focal region so discordant with what is observed for the larger genome-wide context.

<img src="https://github.com/Kalan-Lab/codoff/assets/4260723/05b6c7f1-9e8f-4865-b0fd-758c933bab83" width="300"> 


**This is useful because it could indicate that the focal gene cluster has been horizontally transfered. While quick and easy to assess, we encourage users to further investigate if HGT is indeed responsible for any such signal. This is because a discordance between the codon usage for the focal gene cluster and background genome could be due to other reasons, e.g. infrequent expression of the focal gene cluster.**

codoff is largely adapted from the `compareBGCtoGenomeCodonUsage.py` script in the [*lsa*BGC](https://github.com/Kalan-Lab/lsaBGC) suite and presented here separately just to make it easy to install via bioconda, therefore, please cite:

> [Evolutionary investigations of the biosynthetic diversity in the skin microbiome using *lsa*BGC](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000988). *Microbial Genomics 2023.* Rauf Salamzade, J.Z. Alex Cheong, Shelby Sandstrom, Mary Hannah Swaney, Reed M. Stubbendieck, Nicole Lane Starr, Cameron R. Currie, Anne Marie Singh, and Lindsay R. Kalan

If you have suggestions or feedback, please let us know in the GitHub issues.

## Installation

### Bioconda

***Note, (for some setups at least) it is critical to specify the conda-forge channel before the bioconda channel to properly configure priority and lead to a successful installation.***

**Recommended**: For a significantly faster installation process, use mamba in place of conda in the below commands, by installing mamba in your base conda environment.

```bash
conda create -n codoff_env -c conda-forge -c bioconda codoff
conda activate codoff_env
```

### Conda Manual

```bash
# 1. clone Git repo and change directories into it!
git clone https://github.com/raufs/codoff/
cd codoff/

# 2. create conda environment using yaml file and activate it!
conda env create -f codoff_env.yml -n codoff_env
conda activate codoff_env

# 3. complete python installation with the following command:
pip install -e .
```

## Example commands

Uncompress the example inputs (trimmed down to save space) from this git repo for *Staphylococcus warneri* st. 413:

```
# within codoff git repo
tar -zxvf Sw_LK413.tar.gz
```

### Example 1: Providing focal region and genome-wide GenBank files as input (e.g. based on antiSMASH outputs)

(WORKS FOR BOTH EUKARYOTES & BACTERIA) Focal region and full-genome provided as GenBank files with CDS features (compatible with antiSMASH outputs). Multiple focal region GenBank files can be provided, e.g. consider a biosynthetic gene cluster split across multiple scaffolds due to assembly fragmentation. 

```bash
codoff -f Sw_LK413/NZ_JALXLO020000001.1.region001.gbk -g Sw_LK413/LK413.gbk
```

### Example 2: Providing a genome-wide FASTA file and coordinates of a focal region

(WORKS ONLY FOR BACTERIA) Full genome is provided as a FASTA or GenBank file. If CDS features are missing, gene calling is performed using [pyrodigal](https://github.com/althonos/pyrodigal). Afterwards, the focal region is determined through user speciefied coordinates.

```bash
 codoff -s NZ_JALXLO020000001.1 -a 341425 -b 388343 -g Sw_LK413/LK413.fna -p example_plot.svg
```

Here, we also requested the `-p` argument to generate a plot of the simulated distribution of cosine distances for regions of similar size to the focal region and the actual cosine distance for the focal region/cluster (blue vertical line):

![figure](https://raw.githubusercontent.com/Kalan-Lab/codoff/main/codoff_actual_empirical_pvalue_image.svg)

### Detailed example with interpretation

A more detailed example on how homologous instances of the same five-gene operon differ in codon usage between an instance on the plasmid and chromosome can be found on the wiki: [Examples of inferring codon usage for the *crt* operon in a chromosomal and plasmid context](https://github.com/Kalan-Lab/codoff/wiki/Examples-of-inferring-codon-usage-differences-for-the-crt-operon-in-a-chromosomal-and-plasmid-context)

## Algorithm for computing empirical P-value

To calculate an empirical P-value, we gather codons for each gene across the genome. First the codon frequency distribution of the focal-region/BGC genes is compared to the codon frequency distribution of the background genome (all other genes; genes which have lengths not divisible by 3 are ignored). After, we perform 10,000 simulations where in each simulation we shuffle the full genome-wide set of genes and go through the first N genes until the same number of codons as are present in the focal region/BGC/GCF are observed. The cosine distance between the observed codon frequency is compared to the remainder of the genome-wide codon distribution and checked for whether it is higher than what was actually observed for the BGC; if so, then an empirical P-value counter is appended a count of 1. The final empirical P-value produced is simply this count plus a pseudocount of 1 over 10,001.

 <!---![figure](https://github.com/Kalan-Lab/codoff/blob/main/codoff_empirical_pvalue_image.png?raw=true) --->

## Usage 

```
usage: codoff [-h] -g FULL_GENOME [-s SCAFFOLD] [-a START_COORD] [-b END_COORD] [-f FOCAL_GENBANKS [FOCAL_GENBANKS ...]] [-o OUTFILE] [-p PLOT_OUTFILE] [-v]

	Program: codoff
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
		
	This program compares the codon-usage distribution of a focal-region/BGC to the codon usage of 
    the background genome. It will report the cosine distance and Spearman correlation between the 
    two profiles, as well as an empirical P-value for whether the codon usage is siginficantly different
    between the focal region and background genome. Only CDS features which are of length divisible 
    by 3 will be considered. 
                                     
    Two modes of input are supported:
                                     
    1. (WORKS FOR BOTH EUKARYOTES & BACTERIA) Focal region and full-genome provided as GenBank files with 
       CDS features (compatible with antiSMASH outputs). Multiple focal region GenBank files can be provided, 
       e.g. consider a biosynthetic gene cluster split across multiple scaffolds due to assembly fragmentation. 
                                                                         
       Example command: 
                                     
       $ codoff -f Sw_LK413/NZ_JALXLO020000001.1.region001.gbk -g Sw_LK413/LK413.gbk
                                     
    2. (WORKS ONLY FOR BACTERIA) Full genome is provided as a FASTA or GenBank file. If CDS features are missing
       gene calling is performed using pyrodigal. Afterwards, the focal region is determined through user
       speciefied coordinates.
    
       Example command:
                                     
       $ codoff -s NZ_JALXLO020000001.1 -a 341425 -b 388343 -g Sw_LK413/LK413.fna 
 
    

options:
  -h, --help            show this help message and exit
  -g FULL_GENOME, --full_genome FULL_GENOME
                        Path to annotated full-genome in GenBank or FASTA format for isolate's genome.
  -s SCAFFOLD, --scaffold SCAFFOLD
                        Scaffold identifier for focal region.
  -a START_COORD, --start_coord START_COORD
                        Start coordinate for focal region.
  -b END_COORD, --end_coord END_COORD
                        End coordinate for focal region.
  -f FOCAL_GENBANKS [FOCAL_GENBANKS ...], --focal_genbanks FOCAL_GENBANKS [FOCAL_GENBANKS ...]
                        Path to focal region GenBank(s) for isolate. Locus tags must match with tags in full-genome GenBank.
  -o OUTFILE, --outfile OUTFILE
                        Path to output file.
  -p PLOT_OUTFILE, --plot_outfile PLOT_OUTFILE
                        Plot output file name (will be in SVG format). If not provided, no plot will be made.
  -v, --version         Print version and exist
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

## Acknowledgments

codoff utilizes pyrodigal underneath to perform gene calling for bacterial genomes when requested: 

> Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes. *JOSS*, 2022. Martin Larallde.

Earlier work on using codon usage discrepancies to infer genomic islands was previously done by Waack et al. 2006:

> Score-based prediction of genomic islands in prokaryotic genomes using hidden Markov models. _BMC Bioinformatics_, 2006. Stephan Waack, Oliver Keller, Roman Asper, Thomas Brodag, Carsten Damm, Wolfgang Florian Fricke, Katharina Surovcik, Peter Meinicke, and Rainer Merkl.
