#!/bin/bash

set -eu -o pipefail

expected_score="2.00544"

echo "Extracting primates_odb12.tar.gz..."
tar -xzf primates_odb12.tar.gz

# Run aaKomp
echo "Running aaKomp..."
run-aakomp --db-dir ./ --reference primates_odb12/primates_odb12.faa --input small_test.fa -t 4 -o demo

# Compare score
score=$(head -n 1 demo_score.txt)
echo "Observed score: $score"
echo "Expected score: $expected_score"

if (( $(echo "$score == $expected_score" | bc -l) )); then
  echo "Score matches expected value."
else
  echo "Score does not match expected value."
  exit 1
fi

exit 0
