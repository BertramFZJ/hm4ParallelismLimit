#define _OMP_NESTED_ 1

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>
#include <omp.h>

#include "cfdHeader.h"
#include "pm18MpiHeader.h"
#include "mpiTransferHost.h"
#include "mpiCfdCore.h"
#include "cfdFieldValuesInitialization.h"

#include "mpiCfdCoreMain.h"

double cfdCoreModelMainLevelAc1(int nITERATION, gtypeHm4CfdHostCore *cfdCore, pm18MpiTopologyType mpiTopology)
{
	double timerStart, timerFinish, taskTime, totalTime;

	int iLevel = -1;
	int ITERATION = 0;

	MPI_Barrier(mpiTopology.commTASK);
	if(mpiTopology.rankTASK == 0)
	{ printf("\n"); printf("============== CFD CORE ==============\n"); printf("\n"); fflush(stdout); }
	MPI_Barrier(mpiTopology.commTASK);

	timerStart = MPI_Wtime();
	for(ITERATION=0; ITERATION<nITERATION; ITERATION++)
	{
		double DT;

		// ����� �������
		if(iLevel == -1)
		{
			int ii;		
			
			// ���������� ������ �������� ������
			for(ii=0; ii<mpiTopology.hostInterconnect.sendBufferTotalLength; ii++)
			{ 
				double *dPtr = cfdCore->cellQph + mpiTopology.hostInterconnect.SendALocalIndex[ii]*5;
				int j;

				for(j=0; j<5; j++) 
					mpiTopology.hostInterconnect.sendValuesBuffer[ii*5 + j] = dPtr[j]; 
			} // for ii
			
			// ����������� ����� �������
			HostsInterconnectStartWaitAll(mpiTopology.hostInterconnect, mpiTopology.commTASK);

			iLevel = cfdCore->levelMax - 1;

			for(ii=cfdCore->nEinterior; ii<cfdCore->nE; ii++)
				if(cfdCore->cellLevels[ii] <= iLevel)
					hm4ComputePrimForPhysVariables(cfdCore->cellQph + ii*5, cfdCore->cellQ + ii*5);
		}
		// ����� �������

		// ��� �� �������
		DT = 0.0000025;
		// ��� �� �������

		// ���������� ������
		mpiCoreCalculateEulerInteriorFacesFluxesAc1(0, cfdCore->icfOffset[iLevel], cfdCore->ICF, cfdCore->faceCE,
			                                        cfdCore->cellQph, 1);
		// ���������� ������

		// ��������� ������
		mpiCoreCalculateEulerBoundaryFacesFluxes(cfdCore->nICF, 0, cfdCore->bcfOffset[iLevel], cfdCore->BCF,
												 cfdCore->faceCE, cfdCore->cellQph, 1);											 
		// ��������� ������

		// ����������
		mpiCoreUpdateCellsValuesForElementsLevel(0, cfdCore->nE, DT, 1.0, 1,
			                                     iLevel, cfdCore->cellLevels,
											     cfdCore->cellQ, cfdCore->cellQph, cfdCore->cellVOL,
											     cfdCore->faceToElementsX, cfdCore->faceToElementsA, cfdCore->faceCE,
											     1);
		// ����������

		iLevel -= 1;
	} // for ITERATION

	MPI_Barrier(mpiTopology.commTASK);
	timerFinish = MPI_Wtime();
	taskTime = timerFinish - timerStart;

	mpiCoreCheckMinMaxPhysAndPrimVariables(cfdCore->nEinterior, cfdCore->cellQ, cfdCore->cellQph, stdout, mpiTopology);

	MPI_Allreduce(&taskTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, mpiTopology.commTASK);
	MPI_Barrier(mpiTopology.commTASK);
	if(mpiTopology.rankTASK == 0) printf("TOTAL TIME %.4lf\n", totalTime); fflush(stdout);
	MPI_Barrier(mpiTopology.commTASK);

	return totalTime;
}

double cfdCoreModelMainAc1(int nITERATION, gtypeHm4CfdHostCore *cfdCore, pm18MpiTopologyType mpiTopology)
{
	double timerStart, timerFinish, taskTime, totalTime;
	double localStart, localFinish;
	double localTransferTime = 0.0;
	double localComputeTime = 0.0;
	double localComputeSpeed = 0.0;
	double transferTime[2], computeTime[2], computeSpeed[2];
	int nInteriorFace[2];

	int ITERATION = 0;

	MPI_Barrier(mpiTopology.commTASK);
	if(mpiTopology.rankTASK == 0)
	{ printf("\n"); printf("============== CFD CORE ==============\n"); printf("\n"); fflush(stdout); }
	MPI_Barrier(mpiTopology.commTASK);

	timerStart = MPI_Wtime();
	for(ITERATION=0; ITERATION<nITERATION; ITERATION++)
	{
		double DT;

		// ����� �������
		localStart = MPI_Wtime();
		{
			int ii;		
			
			// ���������� ������ �������� ������
			for(ii=0; ii<mpiTopology.hostInterconnect.sendBufferTotalLength; ii++)
			{ 
				double *dPtr = cfdCore->cellQph + mpiTopology.hostInterconnect.SendALocalIndex[ii]*5;
				int j;

				for(j=0; j<5; j++) 
					mpiTopology.hostInterconnect.sendValuesBuffer[ii*5 + j] = dPtr[j]; 
			} // for ii
			
			// ����������� ����� �������
			HostsInterconnectStartWaitAll(mpiTopology.hostInterconnect, mpiTopology.commTASK);			
		}
		localFinish = MPI_Wtime(); localTransferTime += localFinish - localStart;
		// ����� �������

		// ��� �� �������
		DT = 0.0000025;
		// ��� �� �������

		localStart = MPI_Wtime();
		// ���������� ������
		mpiCoreCalculateEulerInteriorFacesFluxesAc1(0, cfdCore->icfOffset[0], cfdCore->ICF, cfdCore->faceCE,
			                                        cfdCore->cellQph, 1);
		// ���������� ������

		// ��������� ������
		mpiCoreCalculateEulerBoundaryFacesFluxes(cfdCore->nICF, 0, cfdCore->bcfOffset[0], cfdCore->BCF,
												 cfdCore->faceCE, cfdCore->cellQph, 1);											 
		// ��������� ������

		// ����������
		mpiCoreUpdateCellsValuesForElements(0, cfdCore->nEinterior, DT, 1.0, 1,
			                                cfdCore->cellQ, cfdCore->cellQph, cfdCore->cellVOL,
											cfdCore->faceToElementsX, cfdCore->faceToElementsA, cfdCore->faceCE,
											1);
		// ����������
		localFinish = MPI_Wtime(); localComputeTime += localFinish - localStart;

	} // for ITERATION

	MPI_Barrier(mpiTopology.commTASK);
	timerFinish = MPI_Wtime();
	taskTime = timerFinish - timerStart;

	mpiCoreCheckMinMaxPhysAndPrimVariables(cfdCore->nEinterior, cfdCore->cellQ, cfdCore->cellQph, stdout, mpiTopology);

	localComputeSpeed = ((double)ITERATION * (double)cfdCore->nICF) / localComputeTime;

	MPI_Allreduce(&taskTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, mpiTopology.commTASK);
	MPI_Allreduce(&localComputeTime, &computeTime[0], 1, MPI_DOUBLE, MPI_MIN, mpiTopology.commTASK);
	MPI_Allreduce(&localComputeTime, &computeTime[1], 1, MPI_DOUBLE, MPI_MAX, mpiTopology.commTASK);
	MPI_Allreduce(&localTransferTime, &transferTime[0], 1, MPI_DOUBLE, MPI_MIN, mpiTopology.commTASK);
	MPI_Allreduce(&localTransferTime, &transferTime[1], 1, MPI_DOUBLE, MPI_MAX, mpiTopology.commTASK);
	MPI_Allreduce(&cfdCore->nICF, &nInteriorFace[0], 1, MPI_INT, MPI_MIN, mpiTopology.commTASK);
	MPI_Allreduce(&cfdCore->nICF, &nInteriorFace[1], 1, MPI_INT, MPI_MAX, mpiTopology.commTASK);
	MPI_Allreduce(&localComputeSpeed, &computeSpeed[0], 1, MPI_DOUBLE, MPI_MIN, mpiTopology.commTASK);
	MPI_Allreduce(&localComputeSpeed, &computeSpeed[1], 1, MPI_DOUBLE, MPI_MAX, mpiTopology.commTASK);

	MPI_Barrier(mpiTopology.commTASK);
	if(mpiTopology.rankTASK == 0)
	{
		printf("   TOTAL TIME: %12.4lf\n", totalTime);
		printf(" COMPUTE TIME: %12.4lf <> %12.4lf ==> %10.4lf\n", computeTime[0], computeTime[1], computeTime[1] / computeTime[0]);
		printf("COMPUTE SPEED: %12.4E <> %12.4E ==> %10.4lf\n", computeSpeed[0], computeSpeed[1], computeSpeed[1] / computeSpeed[0]);
		printf("TRANSFER TIME: %12.4lf <> %12.4lf\n", transferTime[0], transferTime[1]);
		printf("   ICF NUMBER: %12d <> %12d ==> %10.4lf\n", nInteriorFace[0], nInteriorFace[1], (double)nInteriorFace[1] / (double)nInteriorFace[0]);
		fflush(stdout);
	}
	MPI_Barrier(mpiTopology.commTASK);

	return totalTime;
}

double cfdCoreModelMainLevelAc2(int nITERATION, gtypeHm4CfdHostCore *cfdCore, pm18MpiTopologyType mpiTopology)
{
	double timerStart, timerFinish, taskTime, totalTime;

	int iLevel = -1;
	int ITERATION = 0;

	MPI_Barrier(mpiTopology.commTASK);
	if(mpiTopology.rankTASK == 0)
	{ printf("\n"); printf("============== CFD CORE ==============\n"); printf("\n"); fflush(stdout); }
	MPI_Barrier(mpiTopology.commTASK);

	timerStart = MPI_Wtime();
	for(ITERATION=0; ITERATION<nITERATION; ITERATION++)
	{
		double DT;

#if 0
		MPI_Barrier(mpiTopology.commTASK);
		if(mpiTopology.rankTASK == 0) printf("ITERATION %d\n", ITERATION); fflush(stdout);
		MPI_Barrier(mpiTopology.commTASK);
#endif

		// ����� �������
		if(iLevel < 0)
		{
			int ii;		
			
			// ���������� ������ �������� ������
			for(ii=0; ii<mpiTopology.hostInterconnect.sendBufferTotalLength; ii++)
			{ 
				double *dPtr = cfdCore->cellQph + mpiTopology.hostInterconnect.SendALocalIndex[ii]*5;
				int j;

				for(j=0; j<5; j++) 
					mpiTopology.hostInterconnect.sendValuesBuffer[ii*5 + j] = dPtr[j]; 
			} // for ii
			
			// ����������� ����� �������
			HostsInterconnectStartWaitAll(mpiTopology.hostInterconnect, mpiTopology.commTASK);

			iLevel = cfdCore->levelMax - 2;

			for(ii=cfdCore->nEinterior; ii<cfdCore->nE; ii++)
				if(cfdCore->cellLevels[ii] <= iLevel)
					hm4ComputePrimForPhysVariables(cfdCore->cellQph + ii*5, cfdCore->cellQ + ii*5);
		}
		// ����� �������

		// ��� �� �������
		DT = 0.00001;
		// ��� �� �������

		// ���������
		mpiCoreCalculateGaussGradientNoLimiterLevel(0, cfdCore->nE, cfdCore->cellLevels, iLevel + 1,
			                                        cfdCore->cellQph, cfdCore->cellCEL, 
											        cfdCore->gaussX, cfdCore->gaussA, cfdCore->kXYZC, cfdCore->kXYZS,
											        cfdCore->faceToElementsX, cfdCore->faceToElementsA,
											        cfdCore->nICF, cfdCore->ICF,
											        cfdCore->nBCF, cfdCore->BCF,
											        cfdCore->celldQph, 
											        1);										   
		// ���������

		// ���������� ������
		mpiCoreCalculateEulerInteriorFacesFluxesAc2(0, cfdCore->icfOffset[iLevel], cfdCore->ICF, cfdCore->faceCE,
		                                            cfdCore->cellQph, cfdCore->celldQph, cfdCore->cellCEL, 1);
		// ���������� ������

		// ��������� ������
		mpiCoreCalculateEulerBoundaryFacesFluxes(cfdCore->nICF, 0, cfdCore->bcfOffset[iLevel], cfdCore->BCF,
												 cfdCore->faceCE, cfdCore->cellQph, 1);											 
		// ��������� ������

		// ����������
		mpiCoreUpdateCellsValuesForElementsLevel(0, cfdCore->nE, DT, 1.0, 1,
			                                     iLevel, cfdCore->cellLevels,
											     cfdCore->cellQ, cfdCore->cellQph, cfdCore->cellVOL,
											     cfdCore->faceToElementsX, cfdCore->faceToElementsA, cfdCore->faceCE,
											     1);
		// ����������

		iLevel -= 2;
	} // for ITERATION

	MPI_Barrier(mpiTopology.commTASK);
	timerFinish = MPI_Wtime();
	taskTime = timerFinish - timerStart;

	mpiCoreCheckMinMaxPhysAndPrimVariables(cfdCore->nEinterior, cfdCore->cellQ, cfdCore->cellQph, stdout, mpiTopology);

	MPI_Allreduce(&taskTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, mpiTopology.commTASK);
	MPI_Barrier(mpiTopology.commTASK);
	if(mpiTopology.rankTASK == 0) printf("TOTAL TIME %.4lf\n", totalTime); fflush(stdout);
	MPI_Barrier(mpiTopology.commTASK);

	return totalTime;
}

double cfdCoreModelMainAc2(int nITERATION, gtypeHm4CfdHostCore *cfdCore, pm18MpiTopologyType mpiTopology)
{
	double timerStart, timerFinish, taskTime, totalTime;

	int ITERATION = 0;

	MPI_Barrier(mpiTopology.commTASK);
	if(mpiTopology.rankTASK == 0)
	{ printf("\n"); printf("============== CFD CORE ==============\n"); printf("\n"); fflush(stdout); }
	MPI_Barrier(mpiTopology.commTASK);

	timerStart = MPI_Wtime();
	for(ITERATION=0; ITERATION<nITERATION; ITERATION++)
	{
		double DT;

#if 0
		MPI_Barrier(mpiTopology.commTASK);
		if(mpiTopology.rankTASK == 0) printf("ITERATION %d\n", ITERATION); fflush(stdout);
		MPI_Barrier(mpiTopology.commTASK);
#endif

		// ����� �������
		{
			int ii;		
			
			// ���������� ������ �������� ������
			for(ii=0; ii<mpiTopology.hostInterconnect.sendBufferTotalLength; ii++)
			{ 
				double *dPtr = cfdCore->cellQph + mpiTopology.hostInterconnect.SendALocalIndex[ii]*5;
				int j;

				for(j=0; j<5; j++) 
					mpiTopology.hostInterconnect.sendValuesBuffer[ii*5 + j] = dPtr[j]; 
			} // for ii
			
			// ����������� ����� �������
			HostsInterconnectStartWaitAll(mpiTopology.hostInterconnect, mpiTopology.commTASK);			
		}
		// ����� �������

		// ��� �� �������
		DT = 0.00001;
		// ��� �� �������

		// ���������
		mpiCoreCalculateGaussGradientNoLimiter(0, cfdCore->nE,
			                                   cfdCore->cellQph, cfdCore->cellCEL, 
											   cfdCore->gaussX, cfdCore->gaussA, cfdCore->kXYZC, cfdCore->kXYZS,
											   cfdCore->faceToElementsX, cfdCore->faceToElementsA,
											   cfdCore->nICF, cfdCore->ICF,
											   cfdCore->nBCF, cfdCore->BCF,
											   cfdCore->celldQph, 
											   1);										   
		// ���������

		// ���������� ������
		mpiCoreCalculateEulerInteriorFacesFluxesAc2(0, cfdCore->icfOffset[0], cfdCore->ICF, cfdCore->faceCE,
		                                            cfdCore->cellQph, cfdCore->celldQph, cfdCore->cellCEL, 1);
		// ���������� ������

		// ��������� ������
		mpiCoreCalculateEulerBoundaryFacesFluxes(cfdCore->nICF, 0, cfdCore->bcfOffset[0], cfdCore->BCF,
												 cfdCore->faceCE, cfdCore->cellQph, 1);											 
		// ��������� ������

		// ����������
		mpiCoreUpdateCellsValuesForElements(0, cfdCore->nEinterior, DT, 1.0, 1,
			                                cfdCore->cellQ, cfdCore->cellQph, cfdCore->cellVOL,
											cfdCore->faceToElementsX, cfdCore->faceToElementsA, cfdCore->faceCE,
											1);
		// ����������		
	} // for ITERATION

	MPI_Barrier(mpiTopology.commTASK);
	timerFinish = MPI_Wtime();
	taskTime = timerFinish - timerStart;

	mpiCoreCheckMinMaxPhysAndPrimVariables(cfdCore->nEinterior, cfdCore->cellQ, cfdCore->cellQph, stdout, mpiTopology);

	MPI_Allreduce(&taskTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, mpiTopology.commTASK);
	MPI_Barrier(mpiTopology.commTASK);
	if(mpiTopology.rankTASK == 0) printf("TOTAL TIME %.4lf\n", totalTime); fflush(stdout);
	MPI_Barrier(mpiTopology.commTASK);

	return totalTime;
}

double cfdCoreModelMainTransfer2Ac2(int nITERATION, gtypeHm4CfdHostCore *cfdCore, pm18MpiTopologyType mpiTopology)
{
	double timerStart, timerFinish, taskTime, totalTime;

	int ITERATION = 0;

	MPI_Barrier(mpiTopology.commTASK);
	if(mpiTopology.rankTASK == 0)
	{ printf("\n"); printf("============== CFD CORE ==============\n"); printf("\n"); fflush(stdout); }
	MPI_Barrier(mpiTopology.commTASK);

	timerStart = MPI_Wtime();
	for(ITERATION=0; ITERATION<nITERATION; ITERATION++)
	{
		double DT;

#if 0
		MPI_Barrier(mpiTopology.commTASK);
		if(mpiTopology.rankTASK == 0) printf("ITERATION %d\n", ITERATION); fflush(stdout);
		MPI_Barrier(mpiTopology.commTASK);
#endif

		// ����� ���������� ���������� ����������
#if 1
		{
			int ii;		
			
			// ���������� ������ �������� ������
#if 1
			for(ii=0; ii<mpiTopology.hostInterconnect.sendBufferTotalLength; ii++)
			{ 
				double *dPtr = cfdCore->cellQph + mpiTopology.hostInterconnect.SendALocalIndex[ii]*5;
				int j;

				for(j=0; j<5; j++) 
					mpiTopology.hostInterconnect.sendValuesBuffer[ii*5 + j] = dPtr[j]; 
			} // for ii
#endif
			
			// ����������� ����� �������
			HostsInterconnectStartWaitAll(mpiTopology.hostInterconnect, mpiTopology.commTASK);			
		}
#endif
		// ����� ���������� ���������� ����������

		// ��� �� �������
		   DT = 0.00001;
		// DT = 0.00000001;
		// ��� �� �������

		// ���������
		mpiCoreCalculateGaussGradientNoLimiter(0, cfdCore->nEinterior,
			                                   cfdCore->cellQph, cfdCore->cellCEL, 
											   cfdCore->gaussX, cfdCore->gaussA, cfdCore->kXYZC, cfdCore->kXYZS,
											   cfdCore->faceToElementsX, cfdCore->faceToElementsA,
											   cfdCore->nICF, cfdCore->ICF,
											   cfdCore->nBCF, cfdCore->BCF,
											   cfdCore->celldQph, 
											   1);										   
		// ���������

		// ����� ���������� ���������� ���������� ����������
#if 1
		{
			int ii;		
			
			// ���������� ������ �������� ������
#if 1
			for(ii=0; ii<mpiTopology.hostInterconnect.sendBufferTotalLength; ii++)
			{ 
				double *dPtr = cfdCore->celldQph + mpiTopology.hostInterconnect.SendALocalIndex[ii]*15;
				int j;

				for(j=0; j<15; j++) 
					mpiTopology.hostInterconnect.sendValuesBufferGrad[ii*15 + j] = dPtr[j]; 
			} // for ii		
#endif

			// ����������� ����� �������
			HostsInterconnectGradientsStartWaitAll(mpiTopology.hostInterconnect, mpiTopology.commTASK);			
		}
#endif
		// ����� ���������� ���������� ���������� ����������

		// ���������� ������
		mpiCoreCalculateEulerInteriorFacesFluxesAc2(0, cfdCore->icfOffset[0], cfdCore->ICF, cfdCore->faceCE,
		                                            cfdCore->cellQph, cfdCore->celldQph, cfdCore->cellCEL, 1);
		// ���������� ������

		// ��������� ������
		mpiCoreCalculateEulerBoundaryFacesFluxes(cfdCore->nICF, 0, cfdCore->bcfOffset[0], cfdCore->BCF,
												 cfdCore->faceCE, cfdCore->cellQph, 1);											 
		// ��������� ������

		// ����������
		mpiCoreUpdateCellsValuesForElements(0, cfdCore->nEinterior, DT, 1.0, 1,
			                                cfdCore->cellQ, cfdCore->cellQph, cfdCore->cellVOL,
											cfdCore->faceToElementsX, cfdCore->faceToElementsA, cfdCore->faceCE,
											1);
		// ����������		
	} // for ITERATION

	MPI_Barrier(mpiTopology.commTASK);
	timerFinish = MPI_Wtime();
	taskTime = timerFinish - timerStart;

	mpiCoreCheckMinMaxPhysAndPrimVariables(cfdCore->nEinterior, cfdCore->cellQ, cfdCore->cellQph, stdout, mpiTopology);

	MPI_Allreduce(&taskTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, mpiTopology.commTASK);
	MPI_Barrier(mpiTopology.commTASK);
	if(mpiTopology.rankTASK == 0) printf("TOTAL TIME %.4lf\n", totalTime); fflush(stdout);
	MPI_Barrier(mpiTopology.commTASK);

	return totalTime;
}
