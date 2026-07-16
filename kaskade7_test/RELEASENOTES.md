# Release notes for development version of Kaskade 7.5

## Installation
1. To install required third party libraries, do one of the following:
   - use the installer scripts from InstallDependencies directory
   - rely on the central ZIB installation
2. Edit Makefile.Local to fit your needs
3. Type "make all" - this will build library and documentation, and perform the tests


### Third party libraries
The following libraries are used:

- Dune 2.7 with git version of dune-alugrid, Dune >= 2.7 is required
- boost 1.73
- SuiteSparse 5.3
- MUMPS 5.3.3
- ITSOL 2
- Hypre 2.11.2
- FP from git


### Environment
- A C++17-capable compiler is required.
- GCC 10.2 can be installed by installer script


## Caveats
- Tensor/Tensor3 removed from fem/fixdune.hh, include linalg/tensor.hh directly
  and use Tensor<T,i,j,k> instead of Tensor3<T,i,j,k>.
- FieldMatrix operator * now implemented in Dune, but does not support 1x1 matrix
  times nxm matrix multiplication, interpreting 1x1 as scalar. Instead, there's the
  (experimental) normalForm function, returning a given matrix directly or a scalar
  in case of 1x1 matrix.
- Removed deprecated fem/power.hh. Use utilities/power.hh instead. Removed the deprecated template
  version of power with compile-time fixed integral exponent.
- Removed several unused and undocumented header files.


