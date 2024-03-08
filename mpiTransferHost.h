int InitHostsInterconnectRequestsMPI(MpiHostInterconnectType hostInterconnect, MPI_Comm comm);

int HostsInterconnectStartWaitAll(MpiHostInterconnectType hostInterconnect, MPI_Comm comm);

int RebootHostsInterconnectRecvRequestsMPI(MpiHostInterconnectType *hostInterconnect, double *QphRecv, MPI_Comm comm);

int InitHostsInterconnectGradientRequestsMPI(MpiHostInterconnectType *hostInterconnect, double *GradRecv, MPI_Comm comm);

int HostsInterconnectGradientsStartWaitAll(MpiHostInterconnectType hostInterconnect, MPI_Comm comm);
