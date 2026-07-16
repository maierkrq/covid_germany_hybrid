/*
 * This C-source is a fix of a problem which arises from the use of Windows 10 Linux subsystem library 
 * liblapacke.so.3
 * When linking the library, the names subsequently defined are referenced, and missing. This is due to
 * fact that these are names if subroutines, which are part of the LAPACK testsuite.
 */
int cgejsv_;
int zggsvd3_;
int cggsvd3_;
int sgghd3_;
int cggev3_;
int zgghd3_;
int cgges3_;
int cgetrf2_;
int zgesvj_;
int dgetrf2_;
int cggsvp3_;
int zggsvp3_;
int sggev3_;
int zggev3_;
int zgesvdx_;
int spotrf2_;
int sgges3_;
int sgetrf2_;
int cgesvdx_;
int zpotrf2_;
int dggsvd3_;
int dgghd3_;
int dggev3_;
int zgejsv_;
int sbdsvdx_;
int dpotrf2_;
int sggsvp3_;
int dgesvdx_;
int dgges3_;
int zgges3_;
int dbdsvdx_;
int zgetrf2_;
int cgesvj_;
int sgesvdx_;
int sggsvd3_;
int dggsvp3_;
int cpotrf2_;
int cgghd3_;