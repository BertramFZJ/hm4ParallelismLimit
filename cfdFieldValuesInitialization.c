#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "constantsCFD.h"
#include "cfdHeader.h"

#include "cfdFieldValuesInitialization.h"

int hm4ComputePhysForPrimVariables(double *Q, double *Qph)
{
	double rho1 = 1.0 / Q[0];

	Qph[0] = Q[0];
	Qph[1] = Q[1] * rho1;
	Qph[2] = Q[2] * rho1;
	Qph[3] = Q[3] * rho1;
	Qph[4] = _GAMMA_1_ * (Q[4] - 0.5 * (Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) * rho1);

	return 0;
}

int hm4ComputePrimForPhysVariables(double *Qph, double *Q)
{
	Q[0] = Qph[0];
	Q[1] = Qph[0] * Qph[1];
	Q[2] = Qph[0] * Qph[2];
	Q[3] = Qph[0] * Qph[3];
	Q[4] = _1_DIV_GAMMA_1_ * Qph[4] + 0.5 * Qph[0] * (Qph[1]*Qph[1] + Qph[2]*Qph[2] + Qph[3]*Qph[3]);

	return 0;
}

double hm4CheckCFDDefinedConstants(void)
{
	double diffMAX = 0.0;
	double diffCUR;

	diffCUR = (_GAMMA_ - 1.0) - _GAMMA_1_; diffCUR = fabs(diffCUR);
	if(diffMAX < diffCUR) diffMAX = diffCUR;

	diffCUR = _GAMMA_ / (_GAMMA_ - 1.0) - _GAMMA_DIV_GAMMA_1_; diffCUR = fabs(diffCUR);
	if(diffMAX < diffCUR) diffMAX = diffCUR;

	diffCUR = 1.0 / (_GAMMA_ - 1.0) - _1_DIV_GAMMA_1_; diffCUR = fabs(diffCUR);
	if(diffMAX < diffCUR) diffMAX = diffCUR;

	diffCUR = 1.0 / (_GAMMA_ * _MACH_NUMBER_ * _MACH_NUMBER_) - _1_DIV_GAMMA_DIV_MACH2_; diffCUR = fabs(diffCUR);
	if(diffMAX < diffCUR) diffMAX = diffCUR;

	diffCUR = 1.0 / _RE_ - _MU_; diffCUR = fabs(diffCUR);
	if(diffMAX < diffCUR) diffMAX = diffCUR;

	diffCUR = _GAMMA_ / (_PR_ * (_GAMMA_ - 1.0)) - _GAMMA_DIV_PR_DIV_GAMMA_1_; diffCUR = fabs(diffCUR);
	if(diffMAX < diffCUR) diffMAX = diffCUR;

	diffCUR = _2_DIV_RE_DIV_PR_ - 2.0 / _RE_ / _PR_;
	diffCUR = fabs(diffCUR);
	if(diffMAX < diffCUR) diffMAX = diffCUR;

	return diffMAX;
}

int hm4CfdInitCellsStartValuesForUndisturbedFlow(int nE, double *Q, double *Qph) 
{
	int i, j;

	double QphInf[5], Qinf[5];

	QphInf[0] = 1.0;
	QphInf[1] = 1.0;
	QphInf[2] = 0.0;
	QphInf[3] = 0.0;
	QphInf[4] = _1_DIV_GAMMA_DIV_MACH2_;

	hm4ComputePrimForPhysVariables(QphInf, Qinf);

	for(i=0; i<nE; i++)
		for(j=0; j<5; j++)
		{
			Qph[i*5 + j] = QphInf[j];
			Q[i*5 + j]   = Qinf[j]  ;
		} // for i	

	return 0;
}

int hm4SetBoundaryConditionsForBoundaryFaces(char *sfNAMES, int nameL, FILE *stream,
											 int nBCF, gtypeBoundaryCellFace *BCF, 
											 int nBSF, int *BoundaryConditions)
{
	int i;

	for(i=0; i<nBCF; i++)
	{
		if(BCF[i].surfaceIndex >= nBSF) { fprintf(stderr, "hm4SetBoundaryConditionsForBoundaryFaces --> *** ERROR *** WRONG SURFACE ID = %d FOR BOUNDARY FACE %d OF ELELEMNT %d\n", BCF[i].surfaceIndex, i, BCF[i].eL); exit(0); }
		BCF[i].boundaryConditionType = BoundaryConditions[BCF[i].surfaceIndex];	
	} // for i

	if(stream != NULL)
	{
		fprintf(stream, "BOUNDARY CONDITIONS\n");

		for(i=0; i<nBSF; i++)
		{
			if(BoundaryConditions[i] == 0) fprintf(stream, "INLET FREE STREAM ");
			if(BoundaryConditions[i] == 1) fprintf(stream, "       REFLECTION ");
			if(BoundaryConditions[i] == 2) fprintf(stream, "           NOSLIP ");
			if(BoundaryConditions[i] == 3) fprintf(stream, "      OUTLET CHAR ");
			if(BoundaryConditions[i] == 4) fprintf(stream, "      OUTLET COPY ");			
			
			if(sfNAMES != NULL) fprintf(stream, "%s\n", sfNAMES + i*nameL);
			               else fprintf(stream, "NONAME %5d\n", i);
		} // for i

		fflush(stream);
	}	

	return 0;
}

#if 0
int hm4CalculateResiduals(double DT, double *Qph,
						  int nE, double *cellVOL,
						  int *faceToElementsX, int *faceToElementsA, double *faceCE,
						  FILE *stream)
{
	double residuals[5]  = {0.0, 0.0, 0.0, 0.0, 0.0};
		
	int i;
	
	for(i=0; i<nE; i++)
	{
		double CE[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
		double kSUM = DT / cellVOL[i];

		int j, n;

		for(j=faceToElementsX[i]; j<faceToElementsX[i + 1]; j++)
		{
			if(faceToElementsA[j] > 0)
			{
				int idFace = faceToElementsA[j] - 1;
				for(n=0; n<5; n++) CE[n] -= faceCE[idFace*5 + n];
			}
			else
			{
				int idFace = - faceToElementsA[j] - 1;
				for(n=0; n<5; n++) CE[n] += faceCE[idFace*5 + n];
			}
		} // for j
				
		for(j=0; j<5; j++) residuals[j] += (CE[j] * kSUM)*(CE[j] * kSUM) * cellVOL[i];
	} // for i
	
	for(i=0; i<5; i++) residuals[i] = sqrt(residuals[i]) / DT;
	fprintf(stream, "RESIDUALS: RO = %15g E = %15g\n", residuals[0], residuals[4]);
	fprintf(stream, "\n"); fflush(stream);
	
	return 0;
}
#endif
