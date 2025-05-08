#!/bin/bash

set -eu -o pipefail

# Download the C. elegans genome
echo "Downloading C. elegans genome..."
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz

# Decompress the genome file
echo "Decompressing genome..."
gunzip -c GCF_000002985.6_WBcel235_genomic.fna.gz > GCF_000002985.6_WBcel235_genomic.fna

# Link the protein reference file
echo "Linking nematoda_odb12.faa..."
ln -sf nematoda_odb12/nematoda_odb12.faa ./

# Run aaKomp
echo "Running aaKomp..."
run-aakomp --db-dir ./ --reference nematoda_odb12.faa --input GCF_000002985.6_WBcel235_genomic.fna -o demo

# Check GFF comparison
echo "Checking for differences from expected output..."
if comm -23 test.gff demo.gff | grep .; then
    echo -e "\nTest failed - output differs from test.gff"
    exit 1
else
    echo -e "\nTest successful - outputs match!"
fi

exit 0

