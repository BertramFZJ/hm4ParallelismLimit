int ReadHybMeshNodesBinaryMPI(char *fname, int *NDIST, double **C, FILE *outSTREAM, MPI_Comm comm);

int ReadHybMeshElementsBinaryMPI(char *fname, int *EDIST, int **XE, int **AE, FILE *outSTREAM, MPI_Comm comm);

int ReadHybMeshBoundaryBinaryMPI(char *fname, int *nFACE, char **faceNAME, int *nEB, int **EB, FILE *outSTREAM, MPI_Comm comm);

int ReadHybMeshNameBinaryMPI(char *fname, char **mname, FILE *outSTREAM, MPI_Comm comm);
