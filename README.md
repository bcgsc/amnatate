# GoldStandard

## Description
Asses draft genome completeness using a hash based approach.


## Credits
Concept: Johnathan Wong and Rene L. Warren

Design and implementation: Johnathan Wong

### Installing from source code:

#### Github repository main branch
 ```
  git clone https://github.com/bcgsc/goldstandard.git
  cd goldstandard
  meson --prefix /path/to/install build
  cd build
  ninja install
 ```
## Dependencies
 * [GCC 7+](https://gcc.gnu.org/) with [OpenMP](https://www.openmp.org/)
 * [python 3.9+](https://www.python.org/)
 * [zlib](https://zlib.net/)
 * [meson](https://mesonbuild.com/Getting-meson.html)
 * [ninja](https://github.com/ninja-build/ninja/)
 * [tcmalloc](https://google.github.io/tcmalloc/quickstart.html)
 * [sdsl-lite](https://github.com/simongog/sdsl-lite)
 * [libdivsufsort](https://github.com/y-256/libdivsufsort)
 * [btllib](https://github.com/bcgsc/btllib)
 * [libsequence](https://github.com/molpopgen/libsequence)

### Installing Dependencies with Conda

We recommend creating a fresh environment
```
conda create --name goldstandard
```

Installing the dependencies
```
conda install -c conda-forge -c bioconda --file requirements.txt
```

## License
GoldStandard Copyright (c) 2022 British Columbia Cancer Agency Branch. All rights reserved.

GoldStandard is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

For commercial licensing options, please contact Patrick Rebstein (prebstein@bccancer.bc.ca).
