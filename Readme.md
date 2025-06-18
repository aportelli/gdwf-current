# GDWF conserved currents
[![License: GPL v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)

[Grid](https://github.com/paboyle/Grid) implementation of conserved vector, axial, 
and seagull currents for generalized domain-wall fermion actions.

### Building the example
Grid is assumed to be installed in a prefix `<prefix>`, and to have been compiled
with the C & C++ compilers `<grid-cc>` and `<grid-cxx>`.

At the root of the repository, configure the build with
```bash
mkdir build; cd build
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_PREFIX_PATH=<prefix> \
         -DCMAKE_CXX_COMPILER=<grid-cxx> -DCMAKE_C_COMPILER=<grid-cc>
```
Then build with `make`.