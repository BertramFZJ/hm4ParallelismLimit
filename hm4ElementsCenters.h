int hm4TetrahedronCenter(double *P1, double *P2, double *P3, double *P4, double *CEL);
int hm4PyramidCenter(double *P0, double *P1, double *P2, double *P3, double *P4, double *CEL);
int hm4PrismCenter(double *P0, double *P1, double *P2, double *P3, double *P4, double *P5, double *CEL);
int hm4HexahedronCenter(double *P0, double *P1, double *P2, double *P3, double *P4, double *P5, double *P6, double *P7, double *CEL);
int hm4InitMeshElementsCenters(int nE, int *xE, int *aE, double *C, double *CEL);
