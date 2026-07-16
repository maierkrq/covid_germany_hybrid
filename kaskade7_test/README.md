# Kaskade 7 repository

Kaskade 7 is a flexible finite element toolbox implemented in C++17
and based heavily on the Dune libraries.

Current development happens on the master branch. If you're more interested in 
using than in developing Kaskade 7, consider checking out the (more or less)
stable Kaskade7.4 branch.

## Getting started
For new users, follow these steps.

### Getting the code
Clone the Kaskade 7 repository (e.g., using ssh access:
```
git clone git@git.zib.de:numerical-mathematics/computational-anatomy-and-physiology/kaskade7.git
```
The repository contains a master branch (default) that contains current development and release branches which are feature-stable (but are otherwise no full-fledged and stable releases). Depending on how intensely you intend to work with Kaskade 7, choose a suitable branch. E.g., working with version 7.4:
```
git checkout Kaskade7.4
```

### Installing required compilers and libraries
In `InstallDependencies/` there are installer scripts that fetch and install compiler and required libraries in matching versions. This is a semi-automated process - inspect and change the shell scripts as you see fit. Start at `install.sh`, where versions, paths, and URLs can be defined centrally. Start the installation process by
```
cd InstallDependencies; sh install.sh
```
If anything goes wrong, inspect and change the corresponding script. They are intended to be modified as needed, and therefore mostly simple and readable.

> ZIB members: We maintain an installation of required libraries. You can skip the library installation step.

### Building Kaskade 7
Next copy `Makefile.Local-template` to `Makefile.Local` and edit it to reflect the installation paths of Kaskade 7 and the required libraries. After that,
```
make kasklib
```
builds the library, `make test` runs the test suite, and `make tutorial` compies the tutorial and runs the tutorial examples.
