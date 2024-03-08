int mpiCoreCalculateTimeStep(int startI, int finishI, 
							 double *Qph, double *cellH, double *DT, 
							 int maxLevel, int *cellsLevel,
							 MPI_Comm commTASK);

int mpiCoreCheckMinMaxPhysAndPrimVariables(int nE, double *Q, double *Qph, FILE *stream, pm18MpiTopologyType mpiTopology);

int mpiCoreCalculateGaussGradientNoLimiterLevel(int startI, int finishI, int *eLevels, int levelMax,
										        double *Qph, double *CEL, 
								                int *gaussX, int *gaussA, double *kXYZC, double *kXYZS,
										        int *faceToElementX, int *faceToElementA,
										        int nICF, gtypeInteriorCellFace *ICF,
										        int nBCF, gtypeBoundaryCellFace *BCF,
										        double *gaussGRAD,
										        int omp_nested);
int mpiCoreCalculateGaussGradientNoLimiter(int startI, int finishI,
										   double *Qph, double *CEL, 
								           int *gaussX, int *gaussA, double *kXYZC, double *kXYZS,
										   int *faceToElementX, int *faceToElementA,
										   int nICF, gtypeInteriorCellFace *ICF,
										   int nBCF, gtypeBoundaryCellFace *BCF,
										   double *gaussGRAD,
										   int omp_nested);

int mpiCoreCalculateEulerInteriorFacesFluxesAc1(int startI, int finishI, gtypeInteriorCellFace *ICF, double *faceCE,
											    double *Qph, int omp_nested);
int mpiCoreCalculateEulerInteriorFacesFluxesAc2(int startI, int finishI, gtypeInteriorCellFace *ICF, double *faceCE,
											    double *Qph, double *dQph, double *cellCEL,
											    int omp_nested);

int mpiCoreCalculateEulerBoundaryFacesFluxes(int nICF, int startI, int finishI, gtypeBoundaryCellFace *BCF, 
											 double *faceCE, double *Qph, int omp_nested);

int mpiCoreUpdateCellsValuesForElementsLevel(int startI, int finishI, double DT, double alfa, int UpDateQ,
										     int maxLevel, int *cellsLevel,
										     double *Q, double *Qph, double *cellVOL,
										     int *faceToElementsX, int *faceToElementsA,
										     double *faceCE,
										     int omp_nested);
int mpiCoreUpdateCellsValuesForElements(int startI, int finishI, double DT, double alfa, int UpDateQ,
										double *Q, double *Qph, double *cellVOL,
										int *faceToElementsX, int *faceToElementsA,
										double *faceCE,
										int omp_nested);

#if 0
int mpiCoreCalculateResidualsAMR(double DT, double *Qph,
							     int nE, double *cellVOL,
							     int *faceToElementsX, int *faceToElementsA, double *faceCE,
							     FILE *stream,
								 amrTypeMPITopology mpiTopology);

int mpiCoreCalculateForcesAMR(int nBCF, gtypeBoundaryCellFace *BCF,
						      double *Qph, double *cellCEL,
						      int nSurfaceTag, int *surfaceTags, int useVisc,
							  double MT, char *fName,
							  amrTypeMPITopology mpiTopology);

int mpiCoreCheckMinMaxPhysAndPrimVariablesAMR(int nE, double *Q, double *Qph, FILE *stream, amrTypeMPITopology mpiTopology);
#endif
