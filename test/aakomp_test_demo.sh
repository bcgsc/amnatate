#!/bin/bash

set -eu -o pipefail

# Download the H. sapiens genome
echo "Downloading H. sapiens genome..."
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz

# Decompress the genome file
echo "Decompressing genome..."
gunzip -c GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz > GCA_009914755.4_T2T-CHM13v2.0_genomic.fna

# Download and extract primates_odb12.tar.gz
echo "Downloading primates_odb12.tar.gz..."
wget -nc https://zenodo.org/records/15421506/files/primates_odb12.tar.gz

echo "Extracting primates_odb12.tar.gz..."
tar -xzf primates_odb12.tar.gz

# Run aaKomp
echo "Running aaKomp..."
run-aakomp --db-dir ./ --reference primates_odb12/primates_odb12.faa --input GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -o demo

exit 0

