#!/bin/bash

set -eu -o pipefail

echo "Decompressing genome..."
gunzip -c GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz > GCA_009914755.4_T2T-CHM13v2.0_genomic.fna

echo "Extracting primates_odb12.tar.gz..."
tar -xzf primates_odb12.tar.gz

# Run aaKomp
echo "Running aaKomp..."
run-aakomp --db-dir ./ --reference primates_odb12/primates_odb12.faa --input GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -o demo

exit 0

