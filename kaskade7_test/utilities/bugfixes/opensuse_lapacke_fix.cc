/*
 * This C-source is a fix of a problem which arises from the use of OpenSuse system library liblapacke.so.3
 * When linking the library, the names subsequently defined are referenced, and missing. This is due to
 * fact that these are names if subroutines, which are part of the LAPACK testsuite.
 */
int dlagge_;
int dlatms_;
int zlagsy_;
int dlagsy_;
int clagsy_;
int claghe_;
int zlaghe_;
int slagge_;
int clagge_;
int clatms_;
int slatms_;
int slagsy_;
int zlatms_;
int zlagge_;
