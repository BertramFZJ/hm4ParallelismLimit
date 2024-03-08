int hm4SerialMeshPreprocessorCfdMPI(char *fname, 
                                    int linkLevelsNumber, 
                                    gtypeHm4MeshAttributes *meshAttributes, 
                                    gtypeHm4MeshTopology *gMesh, gtypeHm4MeshTopology *hMesh,
                                    MPI_Comm comm);

int hm4RotateHostMeshCfdMPI(gtypeHm4MeshTopology *hMesh, MPI_Comm comm);

int InterconnectInitializationHostMPI(gtypeHm4MeshTopology *hMesh, MpiHostInterconnectType *hostInterconnect, MPI_Comm comm);

