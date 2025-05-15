#!/bin/bash

set -eu -o pipefail

echo "Extracting primates_odb12.tar.gz..."
tar -xzf primates_odb12.tar.gz

# Run aaKomp
echo "Running aaKomp..."
run-aakomp --db-dir ./ --reference primates_odb12/primates_odb12.faa --input small_test.fa -t 4  -o demo

exit 0

