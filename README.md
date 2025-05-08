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
- [numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/)

### Installing Dependencies with Conda

We recommend creating a fresh conda environment:

```bash
conda create --name aakomp
conda activate aakomp
conda install -c conda-forge -c bioconda --file requirements.txt
```

---

## Running aaKomp

You can run `aaKomp` either directly or using the wrapper script `run-aakomp`.

### Wrapper Script: `run-aakomp`

The `run-aakomp` wrapper automates:

- Checking for existing miBF
- Building a miBF if missing using `make_mibf`
- Running `aakomp`
- Running post-analysis with `aakomp_score.py`

---

## Demo Example

This demo runs `aaKomp` on the *C. elegans* genome using the `nematoda_odb12` ortholog protein set.

```bash
# Download the C. elegans genome
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz

# Decompress the genome
gunzip -c GCF_000002985.6_WBcel235_genomic.fna.gz > GCF_000002985.6_WBcel235_genomic.fna

# Link the protein reference file to the current directory (assumes you have nematoda_odb12 in the current directory)
ln -sf nematoda_odb12/nematoda_odb12.faa ./

# Run aaKomp through the wrapper
run-aakomp --db-dir ./ \
  --reference nematoda_odb12.faa \
  --input GCF_000002985.6_WBcel235_genomic.fna \
  -o demo
```


## Command-line Options

`run-aakomp` options (partial list):

| Option             | Description                                         |
|--------------------|-----------------------------------------------------|
| `--input` `-i`     | Genome file (FASTA) to assess                      |
| `--reference` `-r` | Protein database (FASTA, amino acid)               |
| `--output` `-o`    | Output prefix                                       |
| `--db-dir`         | Directory to store miBF database                   |
| `--threads` `-t`   | Number of threads (default: 48)                    |
| `--hash` `-H`      | Number of hash functions for miBF (default: 9)    |
| `--kmer` `-k`      | Amino acid k-mer size (default: 9)                 |
| `--lower_bound` `-l` | Minimum occupancy threshold (default: 0.7)        |
| `--rescue_kmer`    | Number of consecutive k-mers to initiate a seed    |
| `--max_offset`     | Max distance to extend seed during chaining        |
| `--track-time`     | Track runtime of each major step                   |
| `--dry-run`        | Print commands only, do not execute                |
| `--verbose` `-v`   | Verbose output                                     |
| `--debug`          | Debug mode for internal troubleshooting            |

---

## License

aaKomp Copyright (c) 2025  
British Columbia Cancer Agency Branch. All rights reserved.

Licensed under the GNU General Public License v3. See [`LICENSE`](LICENSE) or <http://www.gnu.org/licenses/>.

For commercial licensing inquiries, contact:  
**Patrick Rebstein** – prebstein@bccancer.bc.ca
