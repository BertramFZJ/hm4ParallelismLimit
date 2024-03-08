int hm4CreateProcessElementsListMPI(int *eDIST, int *part, int *nePROC, int **ePROC, FILE *outSTREAM, MPI_Comm comm);
int hm4InitElementsLocationMPI(int *eDIST, int *part, int ne, int *elts, int *eltsLOCATION, FILE *outSTREAM, MPI_Comm comm);
int hm4InitElementsLinksMPI(int *eDIST, int *xCSR, int *aCSR, int ne, int *elts, int *neL, int **eltsL, FILE *outSTREAM, MPI_Comm comm);
int hm4InitProcessRecvListMPI(int *eDIST, int *part, int neL, int *eltsL, int *xRECV, int *aRECV, FILE *outSTREAM, MPI_Comm comm);
int hm4InitProcessSendListMPI(int *xRECV, int *aRECV, int **xSEND, int **aSEND, FILE *outSTREAM, MPI_Comm comm);
int hm4CheckSendRecvDataMPI(int *xRECV, int *aRECV, int *xSEND, int *aSEND, FILE *outSTREAM, MPI_Comm comm);
int hm4InitElementsTopologyMPI(int *eDIST, int *xE, int *aE, int neP, int *eP, int **xeP, int **aeP, FILE *outSTREAM, MPI_Comm comm);
int hm4CreateProcessNodesListAndCoordinatesMPI(int nE, int *xE, int *aE, int *vDIST, double *CV, int *nV, int **vLIST, double **vCRD, int nDIM, FILE *outSTREAM, MPI_Comm comm);
