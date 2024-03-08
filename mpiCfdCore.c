#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "pm18MpiHeader.h"

// #include <cfdTYPES.h>
// #include <mpiTopologyHeader.h>

#include "constantsCFD.h"
#include "cfdHeader.h"
#include "cfdFieldValuesInitialization.h"

#include "solverEuler3D.h"

#include "mpiCfdCore.h"

int mpiCoreCalculateTimeStep(int startI, int finishI, 
							 double *Qph, double *cellH, double *DT, 
							 int maxLevel, int *cellsLevel,
							 MPI_Comm commTASK)
{
	double DTD;

	int i;

	for(i=startI; i<finishI; i++) if(cellsLevel[i] <= maxLevel)
	{
		double *dPTR = Qph + i*5;
		double u, c, edt;
		
		u = sqrt(dPTR[1]*dPTR[1] + dPTR[2]*dPTR[2] + dPTR[3]*dPTR[3]);
		c = sqrt(_GAMMA_ * dPTR[4] / dPTR[0]);
		
#if 1 // ÝÉËÅÐ
		edt = cellH[i] / (u + c);
#else // ÍÀÂÜÅ-ÑÒÎÊÑ
		edt =  cellH[i] * cellH[i] 
			/ (cellH[i] * (u + c) + _2_DIV_RE_DIV_PR_);
#endif

		if( (i == 0) || (edt < DTD) ) DTD = edt;		
	} // for i

	DTD *= _CFL_;
	MPI_Allreduce(&DTD, DT, 1, MPI_DOUBLE, MPI_MIN, commTASK);		

	return 0;
}

int mpiCoreCheckMinMaxPhysAndPrimVariables(int nE, double *Q, double *Qph, FILE *stream, pm18MpiTopologyType mpiTopology)
{
	double fminQG[5], fmaxQG[5];
	double fminQphG[5], fmaxQphG[5];

	double fminQGM[5], fmaxQGM[5];
	double fminQphGM[5], fmaxQphGM[5];
	
	int i, j;
	
	for(i=0; i<5; i++)
		fminQG[i]   = fmaxQG[i]   = Q[i],
		fminQphG[i] = fmaxQphG[i] = Qph[i];

	for(i=0; i<nE; i++)
		for(j=0; j<5; j++)
		{
			if(fminQG[j] > Q[i*5 + j]) fminQG[j] = Q[i*5 + j];
			if(fmaxQG[j] < Q[i*5 + j]) fmaxQG[j] = Q[i*5 + j];

			if(fminQphG[j] > Qph[i*5 + j]) fminQphG[j] = Qph[i*5 + j];
			if(fmaxQphG[j] < Qph[i*5 + j]) fmaxQphG[j] = Qph[i*5 + j];
		} // for i

	MPI_Allreduce(fminQG, fminQGM, 5, MPI_DOUBLE, MPI_MIN, mpiTopology.commTASK);
	MPI_Allreduce(fmaxQG, fmaxQGM, 5, MPI_DOUBLE, MPI_MAX, mpiTopology.commTASK);
	MPI_Allreduce(fminQphG, fminQphGM, 5, MPI_DOUBLE, MPI_MIN, mpiTopology.commTASK);
	MPI_Allreduce(fmaxQphG, fmaxQphGM, 5, MPI_DOUBLE, MPI_MAX, mpiTopology.commTASK);

	if(mpiTopology.rankTASK == 0)
	{
		fprintf(stream, "PHYS AND PRIM VARIABLES VALUES\n");
		fprintf(stream, "  Q   MIN: %12g %12g %12g %12g %12g\n", fminQGM[0]  , fminQGM[1]  , fminQGM[2]  , fminQGM[3]  , fminQGM[4]  );
		fprintf(stream, "  Q   MAX: %12g %12g %12g %12g %12g\n", fmaxQGM[0]  , fmaxQGM[1]  , fmaxQGM[2]  , fmaxQGM[3]  , fmaxQGM[4]  );
		fprintf(stream, "  QPH MIN: %12g %12g %12g %12g %12g\n", fminQphGM[0], fminQphGM[1], fminQphGM[2], fminQphGM[3], fminQphGM[4]);
		fprintf(stream, "  QPH MAX: %12g %12g %12g %12g %12g\n", fmaxQphGM[0], fmaxQphGM[1], fmaxQphGM[2], fmaxQphGM[3], fmaxQphGM[4]);
		fprintf(stream, "\n"); fflush(stream);
	}

	MPI_Barrier(mpiTopology.commTASK);
	
	return 0;
}

int mpiCoreCalculateGaussGradientNoLimiterLevel(int startI, int finishI, int *eLevels, int levelMax,
										        double *Qph, double *CEL, 
								                int *gaussX, int *gaussA, double *kXYZC, double *kXYZS,
										        int *faceToElementX, int *faceToElementA,
										        int nICF, gtypeInteriorCellFace *ICF,
										        int nBCF, gtypeBoundaryCellFace *BCF,
										        double *gaussGRAD,
										        int omp_nested)
{
	int i;

	for(i=startI; i<finishI; i++) if(eLevels[i] <= levelMax)
	{
		int nST = gaussX[i + 1] - gaussX[i];
		double *kC = kXYZC + i*6;
		double *kS = kXYZS + gaussX[i]*3;
		double *pQ = Qph + i*5;
		int *indexS = gaussA + gaussX[i];

		double grad[5][3];
		
		int j, k, n;

		// ÑÎÁÑÒÂÅÍÍÛÉ ÂÊËÀÄ		
		for(j=0; j<3; j++)
		{
			grad[0][j] = pQ[0] * (kC[j] + kC[j + 3]);
			grad[1][j] = pQ[1] * (kC[j]            );
			grad[2][j] = pQ[2] * (kC[j]            );
			grad[3][j] = pQ[3] * (kC[j]            );
			grad[4][j] = pQ[4] * (kC[j] + kC[j + 3]);
		} // for j
		// ÑÎÁÑÒÂÅÍÍÛÉ ÂÊËÀÄ

		// ÂÊËÀÄ ØÀÁËÎÍÀ
		for(n=0; n<nST; n++)
		{
			pQ = Qph + indexS[n]*5;
			for(k=0; k<5; k++) for(j=0; j<3; j++) grad[k][j] += pQ[k] * kS[n*3 + j];			
		} // for n
		// ÂÊËÀÄ ØÀÁËÎÍÀ		

		// ÊÎÏÈÐÎÂÀÍÈÅ ÐÅÇÓËÜÒÀÒÀ
		for(k=0; k<5; k++) for(j=0; j<3; j++) gaussGRAD[i*15 + k*3 + j] = grad[k][j];		
	} // for i
	
	return 0;
}

int mpiCoreCalculateGaussGradientNoLimiter(int startI, int finishI,
										   double *Qph, double *CEL, 
								           int *gaussX, int *gaussA, double *kXYZC, double *kXYZS,
										   int *faceToElementX, int *faceToElementA,
										   int nICF, gtypeInteriorCellFace *ICF,
										   int nBCF, gtypeBoundaryCellFace *BCF,
										   double *gaussGRAD,
										   int omp_nested)
{
	int i;

	for(i=startI; i<finishI; i++)
	{
		int nST = gaussX[i + 1] - gaussX[i];
		double *kC = kXYZC + i*6;
		double *kS = kXYZS + gaussX[i]*3;
		double *pQ = Qph + i*5;
		int *indexS = gaussA + gaussX[i];

		double grad[5][3];
		
		int j, k, n;

		// ÑÎÁÑÒÂÅÍÍÛÉ ÂÊËÀÄ		
		for(j=0; j<3; j++)
		{
			grad[0][j] = pQ[0] * (kC[j] + kC[j + 3]);
			grad[1][j] = pQ[1] * (kC[j]            );
			grad[2][j] = pQ[2] * (kC[j]            );
			grad[3][j] = pQ[3] * (kC[j]            );
			grad[4][j] = pQ[4] * (kC[j] + kC[j + 3]);
		} // for j
		// ÑÎÁÑÒÂÅÍÍÛÉ ÂÊËÀÄ

		// ÂÊËÀÄ ØÀÁËÎÍÀ
		for(n=0; n<nST; n++)
		{
			pQ = Qph + indexS[n]*5;
			for(k=0; k<5; k++) for(j=0; j<3; j++) grad[k][j] += pQ[k] * kS[n*3 + j];			
		} // for n
		// ÂÊËÀÄ ØÀÁËÎÍÀ		

		// ÊÎÏÈÐÎÂÀÍÈÅ ÐÅÇÓËÜÒÀÒÀ
		for(k=0; k<5; k++) for(j=0; j<3; j++) gaussGRAD[i*15 + k*3 + j] = grad[k][j];		
	} // for i
	
	return 0;
}

int mpiCoreCalculateEulerInteriorFacesFluxesAc1(int startI, int finishI, gtypeInteriorCellFace *ICF, double *faceCE,
											    double *Qph, int omp_nested)
{
	int i;

	for(i=startI; i<finishI; i++)
	{
		double fluxEuler[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
		
		double *QphL = Qph + ICF[i].eL * 5;
		double *QphR = Qph + ICF[i].eR * 5;
		
		int j;

		EULERSOLVER3D(QphR, QphL,
					  fluxEuler, 
					  ICF[i].n[0], ICF[i].n[1], ICF[i].n[2],
					  ICF[i].s[0], ICF[i].s[1], ICF[i].s[2],
					  ICF[i].t[0], ICF[i].t[1], ICF[i].t[2]);
		
		for(j=0; j<5; j++) faceCE[i*5 + j] = fluxEuler[j] * ICF[i].fS;
	} // for i
	// ÏÅÐÂÛÉ ÏÎÐßÄÎÊ

	return 0;
}

int mpiCoreCalculateEulerInteriorFacesFluxesAc2(int startI, int finishI, gtypeInteriorCellFace *ICF, double *faceCE,
											    double *Qph, double *dQph, double *cellCEL,
											    int omp_nested)
{
	int i;

	for(i=startI; i<finishI; i++)
	{
		double fluxEuler[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
		double QPHL[5], QPHR[5];

		double *QphL = Qph + ICF[i].eL * 5;
		double *dQphL = dQph + ICF[i].eL * 15;
		
		double *QphR = Qph + ICF[i].eR * 5;
		double *dQphR = dQph + ICF[i].eR * 15;
		
		double dxL = ICF[i].fMC[0] - cellCEL[ICF[i].eL*3 + 0];
		double dyL = ICF[i].fMC[1] - cellCEL[ICF[i].eL*3 + 1];
		double dzL = ICF[i].fMC[2] - cellCEL[ICF[i].eL*3 + 2];
		double dxR = ICF[i].fMC[0] - cellCEL[ICF[i].eR*3 + 0];
		double dyR = ICF[i].fMC[1] - cellCEL[ICF[i].eR*3 + 1];
		double dzR = ICF[i].fMC[2] - cellCEL[ICF[i].eR*3 + 2];

		int j;


		for(j=0; j<5; j++) QPHL[j] = QphL[j] + 1.0 * (dQphL[j*3 + 0] * dxL + dQphL[j*3 + 1] * dyL + dQphL[j*3 + 2] * dzL);
		for(j=0; j<5; j++) QPHR[j] = QphR[j] + 1.0 * (dQphR[j*3 + 0] * dxR + dQphR[j*3 + 1] * dyR + dQphR[j*3 + 2] * dzR);

		EULERSOLVER3D(QPHR, QPHL,
					  fluxEuler, 
					  ICF[i].n[0], ICF[i].n[1], ICF[i].n[2],
					  ICF[i].s[0], ICF[i].s[1], ICF[i].s[2],
					  ICF[i].t[0], ICF[i].t[1], ICF[i].t[2]);
		
		for(j=0; j<5; j++) faceCE[i*5 + j] = fluxEuler[j] * ICF[i].fS;
	} // for i

	return 0;
}

int mpiCoreCalculateEulerBoundaryFacesFluxes(int nICF, int startI, int finishI, gtypeBoundaryCellFace *BCF, 
											 double *faceCE, double *Qph, int omp_nested)
{
	int i;

	for(i=startI; i<finishI; i++)
	{
		// ÂÕÎÄÍÀß ÃÐÀÍÈÖÀ
		if(BCF[i].boundaryConditionType == 0)
		{
			double qphINF[5] = {1.0, 1.0, 0.0, 0.0, _1_DIV_GAMMA_DIV_MACH2_};
			double fluxEuler[5];

			int j;

			EULERSOLVER3D(qphINF, &Qph[BCF[i].eL * 5],
						  fluxEuler, 
				          BCF[i].n[0], BCF[i].n[1], BCF[i].n[2],
						  BCF[i].s[0], BCF[i].s[1], BCF[i].s[2],
						  BCF[i].t[0], BCF[i].t[1], BCF[i].t[2]);
			for(j=0; j<5; j++) faceCE[(i + nICF)*5 + j] = fluxEuler[j] * BCF[i].fS;			
		}
		// ÂÕÎÄÍÀß ÃÐÀÍÈÖÀ

		// ÓÑËÎÂÈÅ ÎÒÐÀÆÅÍÈß
		if(BCF[i].boundaryConditionType == 1)
		{
			double fluxEuler[5] = {0.0, 
				                   Qph[BCF[i].eL*5 + 4] * BCF[i].n[0], 
								   Qph[BCF[i].eL*5 + 4] * BCF[i].n[1], 
								   Qph[BCF[i].eL*5 + 4] * BCF[i].n[2],
								   0.0};

			int j;

			for(j=0; j<5; j++) faceCE[(i + nICF)*5 + j] = fluxEuler[j] * BCF[i].fS;			
		}
		// ÓÑËÎÂÈÅ ÎÒÐÀÆÅÍÈß

		// ÓÑËÎÂÈÅ ÏÐÈËÈÏÀÍÈß
		if(BCF[i].boundaryConditionType == 2)
		{
			fprintf(stderr, "ERROR NOSLIP FOR EULRER\n"); fflush(stderr);
			exit(0);
		}
		// ÓÑËÎÂÈÅ ÏÐÈËÈÏÀÍÈß

		// ÑÂÎÁÎÄÍÀß ÃÐÀÍÈÖÀ ÑÎ ÑÍÎÑÎÌ
		if(BCF[i].boundaryConditionType == 3)
		{
			double fluxEuler[5];

			int j;

			EULERSOLVER3D(Qph + BCF[i].eL*5, Qph + BCF[i].eL*5, 
						  fluxEuler, 
				          BCF[i].n[0], BCF[i].n[1], BCF[i].n[2],
						  BCF[i].s[0], BCF[i].s[1], BCF[i].s[2],
						  BCF[i].t[0], BCF[i].t[1], BCF[i].t[2]);
					
			for(j=0; j<5; j++) faceCE[(i + nICF)*5 + j] = fluxEuler[j] * BCF[i].fS;			
		}
		// ÑÂÎÁÎÄÍÀß ÃÐÀÍÈÖÀ ÑÎ ÑÍÎÑÎÌ

		if(BCF[i].boundaryConditionType == 4)
		{
			fprintf(stderr, "ERROR BC 4 FOR EULRER\n"); fflush(stderr);
			exit(0);
		}
	} // for i

	return 0;
}

int mpiCoreUpdateCellsValuesForElementsLevel(int startI, int finishI, double DT, double alfa, int UpDateQ,
										     int maxLevel, int *cellsLevel,
										     double *Q, double *Qph, double *cellVOL,
										     int *faceToElementsX, int *faceToElementsA,
										     double *faceCE,
										     int omp_nested)
{
	double alfaDT = alfa * DT;

	int i;

	for(i=startI; i<finishI; i++) if(cellsLevel[i] <= maxLevel)
	{
		double CE[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
		double QPLUS[5];
		double kSUM = alfaDT / cellVOL[i];
		
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

		for(n=0; n<5; n++) QPLUS[n] = Q[i*5 + n] + CE[n] * kSUM;
		
		if(QPLUS[0] < 0.0) { fprintf(stderr, "Update: ERROR --> ELEMENT %d : DENSITY < 0.0\n", i);  exit(0); }
		if(QPLUS[4] < 0.0) { fprintf(stderr, "Update: ERROR --> ELEMENT %d :  ENERGY < 0.0\n", i);  exit(0); }

		hm4ComputePhysForPrimVariables(QPLUS, Qph + i*5);				
		if(Qph[i*5 + 4] <= 0.0) { fprintf(stderr, "Update: ERROR --> ELEMENT %d : PRESSURE < 0.0\n", i);  exit(0); }	

		if(UpDateQ == 1) for(n=0; n<5; n++) Q[i*5 + n] = QPLUS[n];		
	} // for i

	return 0;
}

int mpiCoreUpdateCellsValuesForElements(int startI, int finishI, double DT, double alfa, int UpDateQ,
										double *Q, double *Qph, double *cellVOL,
										int *faceToElementsX, int *faceToElementsA,
										double *faceCE,
										int omp_nested)
{
	double alfaDT = alfa * DT;

	int i;

	for(i=startI; i<finishI; i++)
	{
		double CE[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
		double QPLUS[5];
		double kSUM = alfaDT / cellVOL[i];
		
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

		for(n=0; n<5; n++) QPLUS[n] = Q[i*5 + n] + CE[n] * kSUM;
		
		if(QPLUS[0] < 0.0) { fprintf(stderr, "Update: ERROR --> ELEMENT %d : DENSITY < 0.0\n", i);  exit(0); }
		if(QPLUS[4] < 0.0) { fprintf(stderr, "Update: ERROR --> ELEMENT %d :  ENERGY < 0.0\n", i);  exit(0); }

		hm4ComputePhysForPrimVariables(QPLUS, Qph + i*5);				
		if(Qph[i*5 + 4] <= 0.0) { fprintf(stderr, "Update: ERROR --> ELEMENT %d : PRESSURE < 0.0\n", i);  exit(0); }	

		if(UpDateQ == 1) for(n=0; n<5; n++) Q[i*5 + n] = QPLUS[n];		
	} // for i

	return 0;
}

#if 0
int mpiCoreCalculateResidualsAMR(double DT, double *Qph,
							     int nE, double *cellVOL,
							     int *faceToElementsX, int *faceToElementsA, double *faceCE,
							     FILE *stream,
								 amrTypeMPITopology mpiTopology)
{
	double  residuals[5]  = {0.0, 0.0, 0.0, 0.0, 0.0};
	double residualsG[5]  = {0.0, 0.0, 0.0, 0.0, 0.0};
		
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
	
	MPI_Allreduce(residuals, residualsG, 5, MPI_DOUBLE, MPI_SUM, mpiTopology.commTASK);
	for(i=0; i<5; i++) residualsG[i] = sqrt(residualsG[i]) / DT;

	if(mpiTopology.runkTASK == 0)
	{
		fprintf(stream, "RESIDUALS: RO = %15g E = %15g\n", residualsG[0], residualsG[4]);
		fprintf(stream, "\n"); fflush(stream);
	}

	MPI_Barrier(mpiTopology.commTASK);
	
	return 0;
}

int mpiCoreCalculateForcesAMR(int nBCF, gtypeBoundaryCellFace *BCF,
						      double *Qph, double *cellCEL,
						      int nSurfaceTag, int *surfaceTags, int useVisc,
							  double MT, char *fName,
							  amrTypeMPITopology mpiTopology)
{
	double   ForcesPV[6]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double ForcesGLOB[6]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	int i, j;

	for(i=0; i<nBCF; i++)
	{
		for(j=0; j<nSurfaceTag; j++) if(BCF[i].surfaceIndex == surfaceTags[j]) break;

		if(j < nSurfaceTag)
		{
			int eL = BCF[i].eL;
			double *ptrQPH = Qph + eL*5;
			double *n = BCF[i].n;
			double *s = BCF[i].s;
			double *t = BCF[i].t;
			double  S = BCF[i].fS;

			ForcesPV[0] += ptrQPH[4] * S * n[0]; 
			ForcesPV[1] += ptrQPH[4] * S * n[1];
			ForcesPV[2] += ptrQPH[4] * S * n[2];				             

			if(useVisc == 1)
			{
				double detX, detY, detZ, detG;
				double txx, txy, txz, tyy, tyz, tzz;
				double dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz;
				double u, v, w, divU;
				double dx, dy, dz, nX, nY, nZ, sX, sY, sZ, tX, tY, tZ;					

				double SVISC = 120.0 / 293.15;
				double G = _GAMMA_ * _MACH_NUMBER_ * _MACH_NUMBER_ * ptrQPH[4] / ptrQPH[0];
				double kMU = pow(G, 1.5) * (1.0 + SVISC) / (G + SVISC);		
				
				u = ptrQPH[1]; 
				v = ptrQPH[2]; 
				w = ptrQPH[3];
				
				dx = BCF[i].fMC[0] - cellCEL[eL*3 + 0];
				dy = BCF[i].fMC[1] - cellCEL[eL*3 + 1];
				dz = BCF[i].fMC[2] - cellCEL[eL*3 + 2];
				
				nX = n[0]; nY = n[1]; nZ = n[2];
				sX = s[0]; sY = s[1]; sZ = s[2];
				tX = t[0]; tY = t[1]; tZ = t[2];

#define _det3_(a11,a12,a13,a21,a22,a23,a31,a32,a33) (a11 * (a22*a33 - a23*a32) - a12 * (a21*a33 - a23*a31) + a13 * (a21*a32 - a22*a31))
				detG = _det3_(dx,dy,dz, sX,sY,sZ, tX,tY,tZ); detG = 1.0 / detG;

				detX = _det3_(-u,dy,dz, 0.0,sY,sZ, 0.0,tY,tZ); dux = detX * detG;
				detY = _det3_(dx,-u,dz, sX,0.0,sZ, tX,0.0,tZ); duy = detY * detG;
				detZ = _det3_(dx,dy,-u, sX,sY,0.0, tX,tY,0.0); duz = detZ * detG;

				detX = _det3_(-v,dy,dz, 0.0,sY,sZ, 0.0,tY,tZ); dvx = detX * detG;
				detY = _det3_(dx,-v,dz, sX,0.0,sZ, tX,0.0,tZ); dvy = detY * detG;
				detZ = _det3_(dx,dy,-v, sX,sY,0.0, tX,tY,0.0); dvz = detZ * detG;

				detX = _det3_(-w,dy,dz, 0.0,sY,sZ, 0.0,tY,tZ); dwx = detX * detG;
				detY = _det3_(dx,-w,dz, sX,0.0,sZ, tX,0.0,tZ); dwy = detY * detG;
				detZ = _det3_(dx,dy,-w, sX,sY,0.0, tX,tY,0.0); dwz = detZ * detG;
#undef _det3_

				divU = dux + dvy + dwz;

				txx = _MU_ * (2.0*dux - 2.0*divU/3.0);
				tyy = _MU_ * (2.0*dvy - 2.0*divU/3.0);
				tzz = _MU_ * (2.0*dwz - 2.0*divU/3.0);

				txy = _MU_ * (duy + dvx);
				txz = _MU_ * (duz + dwx);
				tyz = _MU_ * (dvz + dwy);

				ForcesPV[3] += - kMU * (nX * txx + nY * txy + nZ * txz) * S;
				ForcesPV[4] += - kMU * (nX * txy + nY * tyy + nZ * tyz) * S;
				ForcesPV[5] += - kMU * (nX * txz + nY * tyz + nZ * tzz) * S;
			} // if useVisc
		} // if
	} // for i

	ForcesPV[3] += ForcesPV[0]; ForcesPV[4] += ForcesPV[1]; ForcesPV[5] += ForcesPV[2];

	MPI_Allreduce(ForcesPV, ForcesGLOB, 6, MPI_DOUBLE, MPI_SUM, mpiTopology.commTASK);

	if(mpiTopology.runkTASK == 0)
	{
		FILE *stream = fopen(fName, "a");
		
		if(stream == NULL) exit(0);
		fprintf(stream, "TIME %8.4lf FORCES %12.6E %12.6E %12.6E %12.6E %12.6E %12.6E\n",
			             MT,
						 ForcesGLOB[0], ForcesGLOB[3], ForcesGLOB[1], ForcesGLOB[2], ForcesGLOB[4], ForcesGLOB[5]);
		fclose(stream);
	}

	MPI_Barrier(mpiTopology.commTASK);
	
	return 0;
}
#endif
