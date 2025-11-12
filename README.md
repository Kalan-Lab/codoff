# codoff
[![Bioconda](https://img.shields.io/conda/vn/bioconda/codoff?style=flat-square&maxAge=3600)](https://anaconda.org/bioconda/codoff)
[![PyPI](https://img.shields.io/pypi/v/codoff.svg?style=flat-square&maxAge=3600)](https://pypi.python.org/pypi/codoff)
[![Upload Python Package](https://github.com/Kalan-Lab/codoff/actions/workflows/python-publish.yml/badge.svg)](https://github.com/Kalan-Lab/codoff/actions/workflows/python-publish.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13139677.svg)](https://doi.org/10.5281/zenodo.13139677)

**codoff**: A program to measure the irregularity of the codon usage for a focal genomic region (e.g. a BGC, prophage, etc.) relative to the full genome. It was primarily designed to work off the output of antiSMASH biosynthetic gene cluster (BGC) predictions - but can be more broadly applied as well, with options allowing users to provide a genome in FASTA format and simply specify the coordinates of the focal region of interest. It uses sequential contiguous-window sampling (v1.2.3+) to compare the focal region against randomly positioned genomic windows of the same size, reporting a discordance percentile that indicates how unusual the codon usage is compared to the rest of the genome.

<img src="https://github.com/Kalan-Lab/codoff/assets/4260723/05b6c7f1-9e8f-4865-b0fd-758c933bab83" width="300"> 

> [!IMPORTANT]
> **This is useful because it could indicate that the focal gene cluster has been horizontally transfered. While quick and easy to assess, we encourage users to further investigate if HGT is indeed responsible for any such signal. This is because a discordance between the codon usage for the focal gene cluster and background genome could be due to other reasons, e.g. infrequent/highly-conditional expression of the focal gene cluster.**

codoff is largely adapted from the `compareBGCtoGenomeCodonUsage.py` script in the [*lsa*BGC](https://github.com/Kalan-Lab/lsaBGC) suite and presented here separately just to make it easy to install via Bioconda.

If you have suggestions or feedback, please let us know in the GitHub issues.

## Installation

### Bioconda

***Note, (for some setups at least) it is critical to specify the conda-forge channel before the bioconda channel to properly configure priority and lead to a successful installation.***

**Recommended**: For a significantly faster installation process, use mamba in place of conda in the below commands, by installing mamba in your base conda environment.

```bash
conda create -n codoff_env -c conda-forge -c bioconda codoff
conda activate codoff_env
```

### Pip

```bash
# in some environment with python >=3.10 available
pip install codoff
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

Uncompress the example inputs (trimmed down to save space). The input includes: (1) BGC predictions and (2) the genome of a *Corynebacterium simulans* isolate that we previously reported to feature a non-ribosomal peptide synthase that synthesizes a metallophore in the *lsa*BGC manuscript. Expected results from using v1.2.3 of codoff are also included. 

A demo can be run automatically via the bash script:

```bash
bash run_test_case.sh
```

This should take around half a minute:

```
real   0m34.700s
user   0m32.973s
sys    0m0.558s
```

You can explore the results and also compare the results you produced in the `Csimulans_Data/` folder with the expected results in `Csimulans_Data/expected_results_v1.2.3/` folder.

Let's break down what this test is doing, but first remove the directory: `rm -r Csimulans_Data/` and then uncompress it again: `tar -zxvf Csimulans_Data.tar.gz`.

### First Command: Providing focal region and genome-wide GenBank files as input (e.g. based on antiSMASH outputs)

(WORKS FOR BOTH EUKARYOTES & BACTERIA) Focal region and full-genome provided as GenBank files with CDS features (compatible with antiSMASH outputs). Multiple focal region GenBank files can be provided, e.g. consider a biosynthetic gene cluster split across multiple scaffolds due to assembly fragmentation. 

```bash
codoff -f Csimulans_Data/Coryne_simulans_PES1/NZ_CP014634.1.region001.gbk -g Csimulans_Data/Coryne_simulans_PES1/Coryne_simulans_PES1.gbk
```

> [!NOTE]
> As of v1.2.3, this uses sequential contiguous-window sampling, which provides a more biologically realistic null model by comparing against randomly positioned genomic regions of the same size.

### Second Command: Providing a genome-wide FASTA file and coordinates of a focal region

(WORKS ONLY FOR BACTERIA) Full genome is provided as a FASTA or GenBank file. If CDS features are missing, gene calling is performed using [pyrodigal](https://github.com/althonos/pyrodigal). Afterwards, the focal region is determined through user speciefied coordinates.

```bash
 codoff -s NZ_CP014634.1 -a 1380553 -b 1432929 -g Csimulans_Data/Coryne_simulans_PES1/Coryne_simulans_PES1.gbk -p coordinate_example_plot.svg
```

Here, we also requested the `-p` argument to generate a plot of the simulated distribution of cosine distances for regions of similar size to the focal region and the actual cosine distance for the focal region/cluster (blue vertical line):

![figure](https://raw.githubusercontent.com/Kalan-Lab/codoff/main/images/NZ_CP014634.1.region001.svg)

## Algorithm for computing Discordance Percentile

codoff uses a Monte Carlo simulation approach to calculate a discordance percentile to assess codon usage differences between the focal region(s) of interest and the background genome. The algorithm works as follows:

### 1. Data Preparation
- Extract all CDS features from the genome (genes with lengths divisible by 3)
- Skip CDS features with fuzzy/partial boundaries (indicated by `<` or `>` in GenBank notation, e.g., `<1..500` or `1000..>1500`)
- Calculate codon frequency distributions for:
  - **Focal region**: Codons from genes in the BGC/focal region
  - **Background genome**: All remaining genes in the genome

### 2. Statistical Comparison
- Compute cosine distance between focal and background codon frequency distributions
- Calculate Spearman correlation coefficient between the two distributions

### 3. Monte Carlo Simulation

**As of v1.2.3, codoff uses sequential contiguous-window sampling**, which better preserves the spatial structure of the genome.

#### Sequential Contiguous-Window Sampling

For each of N simulations (default: 10,000, configurable with `--num-sims`):
1. Randomly select a scaffold from those large enough to contain the focal region size
2. Randomly select a starting coordinate on that scaffold
3. Extract all genes fully contained within a contiguous genomic window of the same size as the focal region
4. Calculate codon frequencies for genes in this simulated window
5. Skip windows with no genes or with ≥5% of total genome CDS (to match filtering applied to real data)
6. Calculate background frequencies as: `total_genome_counts - simulated_focal_counts`
7. Compute cosine distance between simulated focal and background frequencies
8. Count how many simulated distances ≥ observed distance

This approach tests whether the observed focal region's codon usage is significantly different from what would be expected for a randomly positioned contiguous genomic region of the same size.

### 4. Result Interpretation

codoff reports a **Discordance Percentile** that indicates where the focal region ranks among all similarly sized genomic windows:

```
Discordance Percentile = (count of simulations with distance ≥ observed distance + 1) / (total simulations + 1) × 100
```

For example, a percentile of 5.0 means the focal region is within the **top 5% most discordant regions** in terms of codon usage compared to the rest of the genome. Lower percentiles indicate more unusual/discordant codon usage, which may be evidence of horizontal gene transfer or other unusual regulatory or evolutionary processes.

> [!NOTE]
> Previously, codoff (<v1.2.2) reported an empirical p-value instead of a dicsocrdance percentile - the two metrics are related in that the discordance percentile is largely just the p-value multiplied by 100. This change was to make the value more easy to interpret and also to simplify the investigation of codoff runs for multiple focal regions (e.g. multiple BGCs across multiple genomes), where expectations of independence across tests might not stand. Note, that in v1.2.3, we also changed the way we perform simulations - so results will also further differ from v1.2.2 based on this.

### 5. Reproducibility

Results can be reproduced exactly using the `--seed/-x` parameter (default: 42). Using the same seed value will produce identical results across multiple runs.

 <!---![figure](https://github.com/Kalan-Lab/codoff/blob/main/codoff_empirical_pvalue_image.png?raw=true) --->

## Usage in Python

Beginning in v1.2.0, codoff can also be used as a function in your Python code. Check out [this wiki page](https://github.com/Kalan-Lab/codoff/wiki/API-and-usage-of-main-functions-in-Python-programs) for info on how to do this and details on the API.

## Commandline usage 

```
usage: codoff [-h] -g FULL_GENOME [-s SCAFFOLD] [-a START_COORD] [-b END_COORD] [-f FOCAL_GENBANKS [FOCAL_GENBANKS ...]] [-o OUTFILE] [-p PLOT_OUTFILE] [-ns NUM_SIMS] [-m MAX_FOCAL_CDS_FRACTION] [-v] [-x SEED]

	Program: codoff
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison

	This program compares the codon-usage distribution of a focal-region/BGC
	to the codon usage of the background genome. It will report the cosine
	distance and Spearman correlation between the two profiles. It uses 
	sequential contiguous-window sampling (v1.2.3+) and reports a discordance
	percentile indicating how unusual the focal region's codon usage is
	compared to similarly sized genomic windows. Only CDS features which
	are of length divisible by 3 will be considered. CDS features with 
	fuzzy/partial boundaries (< or >) are skipped.

	Two modes of input are supported:

	1. (WORKS FOR BOTH EUKARYOTES & BACTERIA) Focal region and full-genome
	   provided as GenBank files with CDS features (compatible with
	   antiSMASH outputs). Multiple focal region GenBank files can be provided,
	   e.g. consider a biosynthetic gene cluster split across multiple
	   scaffolds due to assembly fragmentation.

	   Example command:

	   $ codoff -f Sw_LK413/NZ_JALXLO020000001.1.region001.gbk -g Sw_LK413/LK413.gbk
	   $ codoff -f region.gbk -g genome.gbk --num-sims 5000

	2. (WORKS ONLY FOR BACTERIA) Full genome is provided as a FASTA or GenBank
	   file. If CDS features are missing gene calling is performed using
	   pyrodigal. Afterwards, the focal region is determined through user
	   speciefied coordinates.

	   Example command:

	   $ codoff -s NZ_JALXLO020000001.1 -a 341425 -b 388343 -g Sw_LK413/LK413.fna
	   $ codoff -s scaffold -a 1000 -b 5000 -g genome.fna --num-sims 20000



options:
  -h, --help            show this help message and exit
  -g FULL_GENOME, --full-genome FULL_GENOME
                        Path to a full-genome in GenBank or FASTA format. If GenBank file
                        provided, CDS features are required.
  -s SCAFFOLD, --scaffold SCAFFOLD
                        Scaffold identifier for focal region.
  -a START_COORD, --start-coord START_COORD
                        Start coordinate for focal region.
  -b END_COORD, --end-coord END_COORD
                        End coordinate for focal region.
  -f FOCAL_GENBANKS [FOCAL_GENBANKS ...], --focal-genbanks FOCAL_GENBANKS [FOCAL_GENBANKS ...]
                        Path to focal region GenBank(s) for isolate. Locus tags must match
                        with tags in full-genome GenBank.
  -o OUTFILE, --outfile OUTFILE
                        Path to output file [Default is standard output].
  -p PLOT_OUTFILE, --plot-outfile PLOT_OUTFILE
                        Plot output file name (will be in SVG format). If not provided, no
                        plot will be made.
  -ns NUM_SIMS, --num-sims NUM_SIMS
                        Number of simulations to run [Default: 10000].
  -m MAX_FOCAL_CDS_FRACTION, --max-focal-cds-fraction MAX_FOCAL_CDS_FRACTION
                        Maximum allowed fraction of total genome CDS length for focal region
                        [Default: 0.05].
  -v, --version         Print version and exist
  -x SEED, --seed SEED  Random seed for reproducible results [Default: 42].
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
