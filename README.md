![GitHub release (latest by date)](https://img.shields.io/github/v/release/bcgsc/aakomp)

# aaKomp

## Description

Assess draft genome completeness using a fast, alignment-free, k-mer hash-based approach (aaKomp). This tool uses amino acid k-mers and a multi-index Bloom filter (miBf) to estimate the completeness of genome assemblies.

## Credits

**Concept**: Johnathan Wong and Rene L. Warren  
**Design and Implementation**: Johnathan Wong

---
## Installing from Conda (Preferred method)
Under construction


## Installing from Source

### Clone from GitHub

```bash
git clone https://github.com/bcgsc/aakomp.git
cd aakomp
meson --prefix /path/to/install build
cd build
ninja install
```

---

## Dependencies

- [GCC 7+](https://gcc.gnu.org/) with [OpenMP](https://www.openmp.org/)
- [Python 3.9+](https://www.python.org/)
- [zlib](https://zlib.net/)
- [meson](https://mesonbuild.com/)
- [ninja](https://github.com/ninja-build/ninja/)
- [tcmalloc](https://google.github.io/tcmalloc/quickstart.html)
- [sdsl-lite](https://github.com/simongog/sdsl-lite)
- [libdivsufsort](https://github.com/y-256/libdivsufsort)
- [btllib](https://github.com/bcgsc/btllib)
- [libsequence](https://github.com/molpopgen/libsequence)
- [gperftools](https://github.com/gperftools/gperftools)
- [boost-cpp](https://www.boost.org/)
- [r-base](https://cran.r-project.org/)
- [r-ggplot2](https://ggplot2.tidyverse.org/)
- [r-dplyr](https://dplyr.tidyverse.org/)
- [r-readr](https://readr.tidyverse.org/)
- [r-cairo](https://www.rdocumentation.org/packages/Cairo/)
- [r-gridextra](https://cran.r-project.org/web/packages/gridExtra/)
- [r-pracma](https://cran.r-project.org/web/packages/pracma/)
- [hmmer=3.1](http://hmmer.org/)
- [pigz](https://zlib.net/pigz/)

### Installing Dependencies with Conda

We recommend creating a fresh conda environment:

```bash
conda create --name aakomp
conda activate aakomp
conda install -c conda-forge -c bioconda --file requirements.txt
```

---

## Running aaKomp

You can run `aaKomp` either directly or using the driver script `run-aakomp`.

### Driver Script: `run-aakomp`

The `run-aakomp` driver automates:

- Downloading BUSCO lineages
- Building a miBf if missing using `make_mibf` with BUSCO lineages or provided references
- Running `aakomp`
- Visualizing with `aakomp_plot.R`

---

## Demo Example

Here are two example usages of `run-aakomp`. In both cases, the `--db-dir` flag controls where the miBf (multi-index Bloom filter) is stored and looked up.

```bash
# Option 1: Run aaKomp using a provided reference file
run-aakomp --db-dir ./ \
  --reference reference.faa \
  --input input.fa \
  -t 4 \
  -o output_ref
  # --visualise optional argument to visualise the cumulative distribution function
```

```bash
# Option 2: Run aaKomp using a lineage name (e.g., "eukaryota")
# The lineage's HMMs will be downloaded and consensus sequences will be extracted to generate a reference
run-aakomp --db-dir ./ \
  --lineage eukaryota \
  --input input.fa \
  -t 4 \
  -o output_eukaryota
```

> **Note:**  
> If the required miBF already exists in the specified --db-dir, it will be reused. Otherwise, run-aakomp will create one using either the provided --reference FASTA or a reference derived from the downloaded lineage.



## Command-line Options

`run-aakomp` options:

| Option                   | Description                                                                 |
|--------------------------|-----------------------------------------------------------------------------|
| `--help-aakomp`          | Show help message for the `aakomp` binary and exit                          |
| `--help-mibf`            | Show help message for the `make_mibf` binary and exit                       |
| `-i`, `--input`          | Input genome file in FASTA format                                           |
| `-o`, `--output`         | Output prefix (default: `_`)                                                |
| `-r`, `--reference`      | Amino acid reference file (e.g., orthologous protein set)                   |
| `-t`, `--threads`        | Number of threads to use (default: 48)                                      |
| `-v`, `--verbose`        | Enable verbose output                                                       |
| `--debug`                | Enable debug mode for internal troubleshooting                              |
| `-H`, `--hash`           | Number of hash functions used in miBF (default: 9)                          |
| `-k`, `--kmer`           | Amino acid k-mer size (default: 9)                                          |
| `-l`, `--lower-bound`    | Minimum occupancy threshold for valid hits (default: 0.7)                   |
| `--rescue-kmer`          | Number of consecutive k-mers to initiate a new seed (default: 4)            |
| `--max-offset`           | Maximum offset allowed when extending a seed during chaining (default: 2)   |
| `--lineage`              | Name of BUSCO lineage to auto-download and use as reference                 |
| `--db-dir`               | Directory for or to store miBf database files (default: `./`)               |
| `--dry-run`              | Print commands that would be executed, but don’t run them                   |
| `--track-time`           | Record and report runtime statistics for each major step                    |
| `--odb-version`          | BUSCO ortholog database version (default: `12`)                             |
| `--list-lineages`        | List all available BUSCO lineages and exit                                  |
| `--visualise`            | Visualise the cumulative distribution function                              |
| `--version`              | Print version of aaKomp                                                     |



## License

aaKomp Copyright (c) 2025  
British Columbia Cancer Agency Branch. All rights reserved.

Licensed under the GNU General Public License v3. See [`LICENSE`](LICENSE) or <http://www.gnu.org/licenses/>.

For commercial licensing inquiries, contact:  
**Patrick Rebstein** – prebstein@bccancer.bc.ca
