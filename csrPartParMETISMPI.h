// опълюъ оюпюккекэмюъ дейнлонгхжхъ цпютю аег беянб бепьхм х пеаеп
int csrPartParMetisV4NoWeightsMPI(int *eDIST, int *csrX, int *csrA, int nPART, int *part, FILE *outSTREAM, MPI_Comm comm);
// опълюъ оюпюккекэмюъ дейнлонгхжхъ цпютю я жекнвхякеммшлх беяюлх бепьхм, пюбмшлх вхякс нохпючыхуяъ мю бепьхмс пеаеп
// беяю пеаеп цпютю ме сярюмюбкхбючряъ
int csrPartParMetisV4NodesWeightsByCsrEdgesMPI(int *eDIST, int *csrX, int *csrA, int nPART, int *part, FILE *outSTREAM, MPI_Comm comm);

// опълюъ онякеднбюрекэмюъ дейнлонгхжхъ цпютю я бнглнфмнярэч сярюмнбйх жекнвхякеммшу беянб бепьхм х пеаеп
// wghtN[Ne] - люяяхб я бевюлх бепьхм цпютю (wghtN == NULL - беянб мер)
// wghtE[csrX[Ne]] - люяяхб я беяюлх пеаеп цпютю (wghtE == NULL - беянб мер)
int csrPartMetisNodesEdgesWeightsSerial(int nE, int *csrX, int *csrA, int *wghtN, int *wghtE, int nPART, int *part, FILE *outSTREAM);
