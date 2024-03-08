double hm4TetrahedronVolume(double *P1, double *P2, double *P3, double *P4);
double hm4PyramidVolume(double *P0, double *P1, double *P2, double *P3, double *P4);
double hm4PrismVolume(double *P0, double *P1, double *P2, double *P3, double *P4, double *P5);
double hm4HexahedronVolume(double *P0, double *P1, double *P2, double *P3, double *P4, double *P5, double *P6, double *P7);
int hm4InitMeshElementsVolumes(int nE, int *xE, int *aE, double *C, double *V);
