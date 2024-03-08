double hm4TetrahedronHeight(double *P1, double *P2, double *P3, double *P4);
double hm4PyramidHeight(double *P0, double *P1, double *P2, double *P3, double *P4);
double hm4PrismHeight(double *P0, double *P1, double *P2, double *P3, double *P4, double *P5);
double hm4HexahedronHeight(double *P0, double *P1, double *P2, double *P3, double *P4, double *P5, double *P6, double *P7);
int hm4InitMeshElementsHeights(int nE, int *xE, int *aE, double *C, double *He);
