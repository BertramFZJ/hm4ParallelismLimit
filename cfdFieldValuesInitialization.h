int hm4ComputePhysForPrimVariables(double *Q, double *Qph);
int hm4ComputePrimForPhysVariables(double *Qph, double *Q);

double hm4CheckCFDDefinedConstants(void);

int hm4SetBoundaryConditionsForBoundaryFaces(char *sfNAMES, int nameL, FILE *stream,
											 int nBCF, gtypeBoundaryCellFace *BCF, 
											 int nBSF, int *BoundaryConditions);

int hm4CalculateResiduals(double DT, double *Qph,
						  int nE, double *cellVOL,
						  int *faceToElementsX, int *faceToElementsA, double *faceCE,
						  FILE *stream);

int hm4CheckMinMaxPhysAndPrimVariables(int nE, double *Q, double *Qph, FILE *stream);

int hm4CfdInitCellsStartValuesForUndisturbedFlow(int nE, double *Q, double *Qph);
