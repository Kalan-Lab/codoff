#!/usr/bin/env bash

# Step 0: Uncompress test_case.tar.gz and cd into it.
rm -rf Csimulans_Data/
tar -zxvf Csimulans_Data.tar.gz 
cd Csimulans_Data/

# Step 1: run codoff on the Csimulans_Data directory
codoff -g Coryne_simulans_PES1/Coryne_simulans_PES1.gbk -f Coryne_simulans_PES1/NZ_CP014634.1.region001.gbk -o Coryne_simulans_PES1.gbk.codoff_gbk.txt

# Step 2: run codoff on the Csimulans_Data directory using coordinates
codoff -g Coryne_simulans_PES1/Coryne_simulans_PES1.gbk -s NZ_CP014634.1 -a 1380553 -b 1432929 -o Coryne_simulans_PES1.gbk.codoff_coordinates.txt -p Coryne_simulans_PES1.gbk.codoff.svg

# Step 3: run antismash_codoff on the Csimulans_Data directory
antismash_codoff -a Coryne_simulans_PES1/ -o Coryne_simulans_PES1_antismash_codoff_results/