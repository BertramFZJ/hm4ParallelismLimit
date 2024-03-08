#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "typesMeshHeader.h"
#include "cfdHeader.h"

#include "hm4ElementsVolumes.h"
#include "hm4ElementsHeights.h"
#include "hm4ElementsCenters.h"
#include "hm4ElementsFaces.h"
#include "hm4DualGraphCSR.h"
#include "sortIcfList.h"

#include "cfdGeometryInitialization.h"

int cfdGeometryInitializationHost(gtypeHm4MeshTopology *hMesh, gtypeHm4CfdHostCore *cfdCore, MPI_Comm comm)
{
	int *dualX, *dualA;
	int rankMPI, sizeMPI;

	int i, j, k;

	MPI_Comm_rank(comm, &rankMPI);
	MPI_Comm_size(comm, &sizeMPI);
	MPI_Barrier(comm);

	cfdCore->nE = hMesh->meshElementsNumber;
		
	cfdCore->cellVOL = (double *)malloc(cfdCore->nE * sizeof(double));
	if(cfdCore->cellVOL == NULL) exit(0);
	hm4InitMeshElementsVolumes(hMesh->meshElementsNumber, hMesh->meshEX, hMesh->meshEA, hMesh->meshNodesCoordinates, cfdCore->cellVOL);

	MPI_Barrier(comm);
	if(rankMPI == 0) { printf("GEOMETRY >> CELL VOLUMES\n"); fflush(stdout); }
	MPI_Barrier(comm);

	cfdCore->cellH = (double *)malloc(hMesh->meshElementsNumber * sizeof(double));
	if(cfdCore->cellH == NULL) exit(0);
	hm4InitMeshElementsHeights(hMesh->meshElementsNumber, hMesh->meshEX, hMesh->meshEA, hMesh->meshNodesCoordinates, cfdCore->cellH);
	MPI_Barrier(comm);
	if(rankMPI == 0) { printf("GEOMETRY >> CELL HEIGHTS\n"); fflush(stdout); }
	MPI_Barrier(comm);	

	cfdCore->cellCEL = (double *)malloc(hMesh->meshElementsNumber * 3 * sizeof(double));
	if(cfdCore->cellCEL == NULL) exit(0);
	hm4InitMeshElementsCenters(hMesh->meshElementsNumber, hMesh->meshEX, hMesh->meshEA, hMesh->meshNodesCoordinates, cfdCore->cellCEL);
	MPI_Barrier(comm);
	if(rankMPI == 0) { printf("GEOMETRY >> CELL CENTERS\n"); fflush(stdout); }
	MPI_Barrier(comm);

	hm4MeshDualGraphCSRSerial(hMesh->meshElementsNumber, hMesh->meshEX, hMesh->meshEA, &dualX, &dualA);
	MPI_Barrier(comm);
	if(rankMPI == 0) { printf("GEOMETRY >> DUAL GRAPH\n"); fflush(stdout); }
	MPI_Barrier(comm);

	// спнбмх ябъгмнярх ъвеей он дсюкэмнлс цпютс
	cfdCore->levelMax = 0;
	cfdCore->cellLevels = (int *)malloc(cfdCore->nE * sizeof(int));
	if(cfdCore->cellLevels == NULL) exit(0);
	for(j=0, i=0; i<cfdCore->nE; i++)
		if(hMesh->partitionUP[i] == rankMPI) { cfdCore->cellLevels[i] = 0; j += 1; }
		else cfdCore->cellLevels[i] = -1;
	while(j != cfdCore->nE)
	{
		int nSet = 0;

		cfdCore->levelMax += 1;

		for(i=0; i<cfdCore->nE; i++) if(cfdCore->cellLevels[i] == cfdCore->levelMax - 1)
			for(k=dualX[i]; k<dualX[i+1]; k++)
			{
				if(cfdCore->cellLevels[dualA[k]] == -1)
				{
					cfdCore->cellLevels[dualA[k]] = cfdCore->levelMax;
					j += 1;
					nSet += 1;
				} // if
			} // for k

		if( (j != cfdCore->nE) && (nSet == 0) )
		{ fprintf(stderr, "[%d]: Error cell level initializaetion\n", rankMPI); fflush(stderr); exit(0); }
	} // while

	MPI_Barrier(comm);
	{
		int minLEVEL, maxLEVEL;

		MPI_Allreduce(&cfdCore->levelMax, &minLEVEL, 1, MPI_INT, MPI_MIN, comm);
		MPI_Allreduce(&cfdCore->levelMax, &maxLEVEL, 1, MPI_INT, MPI_MAX, comm);
						
		if(rankMPI == 0)
		{
			fprintf(stdout, "\n");
			fprintf(stdout, "CELL LEVEL\n");
			fprintf(stdout, "MINIMUM CELL LEVEL: %8d\n", minLEVEL);
			fprintf(stdout, "MAXIMUM CELL LEVEL: %8d\n", maxLEVEL);
			fprintf(stdout, "\n"); fflush(stdout);
		}

		MPI_Barrier(comm);
		if(minLEVEL != maxLEVEL) exit(0);		
	}
	MPI_Barrier(comm);
	// спнбмх ябъгмнярх ъвеей он дсюкэмнлс цпютс

	// люяяхб бмсрпеммху цпюмеи ъвеей
	cfdCore->nICF = dualX[cfdCore->nE] / 2;
	cfdCore->ICF = (gtypeInteriorCellFace *)malloc(cfdCore->nICF * sizeof(gtypeInteriorCellFace)); if(cfdCore->ICF == NULL) exit(0);
	// тнплхпнбюмхе яохяйю цпюмеи
	for(k=0, i=0; i<cfdCore->nE; i++)
	{
		for(j=dualX[i]; j<dualX[i+1]; j++)
		{
			int eL = i;
			int eR = dualA[j];

			if(eL < eR)
			{
				cfdCore->ICF[k].eL = eL;
				cfdCore->ICF[k].eR = eR;
				k++;
			} // if eL < eR
		} // for j
	} // for i

	// янпрхпнбйю цпюмеи он спнбмъл ъвеей
	wsortIcfListLevel(cfdCore->ICF, 0, cfdCore->nICF - 1, cfdCore->cellLevels);
	for(i=0; i<cfdCore->nICF; i++)
	{
		int dlevel = cfdCore->cellLevels[cfdCore->ICF[i].eL] - cfdCore->cellLevels[cfdCore->ICF[i].eR];
		if(dlevel < 0) dlevel = - dlevel;
		if(dlevel > 1) { fprintf(stderr, "[%d]: ERROR ICF dLevel > 1\n", rankMPI); exit(0); }
	} // for i

	// хмхжхюкхгюжхъ ценлерпхх цпюмеи
	for(i=0; i<cfdCore->nICF; i++)
	{
		int idL = cfdCore->ICF[i].eL;
		int eL = hMesh->meshEX[idL + 1] - hMesh->meshEX[idL];
		int *vL = &hMesh->meshEA[hMesh->meshEX[idL]];

		int idR = cfdCore->ICF[i].eR;
		int eR = hMesh->meshEX[idR + 1] - hMesh->meshEX[idR];
		int *vR = &hMesh->meshEA[hMesh->meshEX[idR]];
				
		int nF, vF[4];

		hm4ElementsSameFace(eL, vL, eR, vR, &nF, vF);

		// рпесцнкэмюъ цпюмэ
		if(nF == 3)
		{
			hm4TriangleCenter(hMesh->meshNodesCoordinates + vF[0]*3, hMesh->meshNodesCoordinates + vF[1]*3, 
				              hMesh->meshNodesCoordinates + vF[2]*3, cfdCore->ICF[i].fMC);
			hm4TriangleFace(hMesh->meshNodesCoordinates + vF[0]*3, hMesh->meshNodesCoordinates + vF[1]*3, 
				            hMesh->meshNodesCoordinates + vF[2]*3, cfdCore->ICF[i].n);
		}

		// вершпеусцнкэмюъ цпюмэ
		if(nF == 4)
		{
			hm4QuadCenter(hMesh->meshNodesCoordinates + vF[0]*3, hMesh->meshNodesCoordinates + vF[1]*3, 
				          hMesh->meshNodesCoordinates + vF[2]*3, hMesh->meshNodesCoordinates + vF[3]*3,
						  cfdCore->ICF[i].fMC);
			hm4QuadFace(hMesh->meshNodesCoordinates + vF[0]*3, hMesh->meshNodesCoordinates + vF[1]*3, 
				        hMesh->meshNodesCoordinates + vF[2]*3, hMesh->meshNodesCoordinates + vF[3]*3,
						cfdCore->ICF[i].n);
		}

		// пюгбнпнр бейрнпю
		if(cfdCore->ICF[i].eL > cfdCore->ICF[i].eR)
		{ fprintf(stderr, "[%d]: ERROR ICF eL < eR\n", rankMPI); exit(0); }

		// бшвхякемхе окныюдх цпюмх
		cfdCore->ICF[i].fS = sqrt(cfdCore->ICF[i].n[0] * cfdCore->ICF[i].n[0] 
		                        + cfdCore->ICF[i].n[1] * cfdCore->ICF[i].n[1] 
								+ cfdCore->ICF[i].n[2] * cfdCore->ICF[i].n[2]);

		// мнплхпнбюмхе бейрнпю окныюдх цпюмх
		for(j=0; j<3; j++) cfdCore->ICF[i].n[j] /= cfdCore->ICF[i].fS;

		// хмхжхюкхгюжхъ рпнийх бейрнпнб опюбни яхярелш йннпдхмюр мю нямнбе бейрнпю мнплюкх й цпюмх
		hm4FaceCoordinateSystem(cfdCore->ICF[i].n, cfdCore->ICF[i].s, cfdCore->ICF[i].t);
	} // for i

	MPI_Barrier(comm);
	if(rankMPI == 0) { printf("GEOMETRY >> INTERIOR FACES\n"); fflush(stdout); }
	MPI_Barrier(comm);
	// люяяхб бмсрпеммху цпюмеи ъвеей

	// люяяхб цпюмхвмшу цпюмеи
	cfdCore->nBCF = hMesh->meshBoundaryFacesNumber;

	if(cfdCore->nBCF > 0)
	{
		cfdCore->BCF = (gtypeBoundaryCellFace *)malloc(cfdCore->nBCF * sizeof(gtypeBoundaryCellFace)); if(cfdCore->BCF == NULL) exit(0);
			
		// хмхжхюкхгюжхъ яохяйю цпюмеи
		for(i=0; i<hMesh->meshBoundaryFacesNumber; i++)
		{
			int nF, vF[4];

			cfdCore->BCF[i].eL           = hMesh->meshBoundaryFacesTopology[i*3 + 0];
			cfdCore->BCF[i].faceIndex    = hMesh->meshBoundaryFacesTopology[i*3 + 1];
			cfdCore->BCF[i].surfaceIndex = hMesh->meshBoundaryFacesTopology[i*3 + 2];
			cfdCore->BCF[i].boundaryConditionType = -1;

			hm4ElementFaceNodes(hMesh->meshEX[cfdCore->BCF[i].eL + 1] - hMesh->meshEX[cfdCore->BCF[i].eL],
				                hMesh->meshEA + hMesh->meshEX[cfdCore->BCF[i].eL],
				                cfdCore->BCF[i].faceIndex, &nF, vF);

			// рпесцнкэмюъ цпюмэ
			if(nF == 3)
			{
				hm4TriangleCenter(hMesh->meshNodesCoordinates + vF[0]*3, hMesh->meshNodesCoordinates + vF[1]*3, 
					              hMesh->meshNodesCoordinates + vF[2]*3, cfdCore->BCF[i].fMC);
				  hm4TriangleFace(hMesh->meshNodesCoordinates + vF[0]*3, hMesh->meshNodesCoordinates + vF[1]*3, 
					              hMesh->meshNodesCoordinates + vF[2]*3, cfdCore->BCF[i].n  );
			}

			// вершпеусцнкэмюъ цпюмэ
			if(nF == 4)
			{
				hm4QuadCenter(hMesh->meshNodesCoordinates + vF[0]*3, hMesh->meshNodesCoordinates + vF[1]*3, 
					          hMesh->meshNodesCoordinates + vF[2]*3, hMesh->meshNodesCoordinates + vF[3]*3, cfdCore->BCF[i].fMC);
				  hm4QuadFace(hMesh->meshNodesCoordinates + vF[0]*3, hMesh->meshNodesCoordinates + vF[1]*3, 
					          hMesh->meshNodesCoordinates + vF[2]*3, hMesh->meshNodesCoordinates + vF[3]*3, cfdCore->BCF[i].n  );
			}

			// бшвхякемхе окныюдх цпюмх
			cfdCore->BCF[i].fS = sqrt(cfdCore->BCF[i].n[0] * cfdCore->BCF[i].n[0] 
				                    + cfdCore->BCF[i].n[1] * cfdCore->BCF[i].n[1] 
									+ cfdCore->BCF[i].n[2] * cfdCore->BCF[i].n[2]);

			// мнплхпнбюмхе бейрнпю окныюдх цпюмх
			for(j=0; j<3; j++) cfdCore->BCF[i].n[j] /= cfdCore->BCF[i].fS;

			// хмхжхюкхгюжхъ рпнийх бейрнпнб опюбни яхярелш йннпдхмюр мю нямнбе бейрнпю мнплюкх й цпюмх
			hm4FaceCoordinateSystem(cfdCore->BCF[i].n, cfdCore->BCF[i].s, cfdCore->BCF[i].t);
		} // for i
		// хмхжхюкхгюжхъ яохяйю цпюмеи

		// янпрхпнбйю яохяйю цпюмеи
		wsortBcfListLevel(cfdCore->BCF, 0, cfdCore->nBCF - 1, cfdCore->cellLevels);
	} // if
	else
		cfdCore->BCF = NULL;
	
	MPI_Barrier(comm);
	if(rankMPI == 0) { printf("GEOMETRY >> BOUNDARY FACES\n"); fflush(stdout); }
	MPI_Barrier(comm);
	// люяяхб цпюмхвмшу цпюмеи

	// йюпрю опхмюдкефмнярх цпюмеи ъвеийюл
	cfdCore->faceToElementsX = (int *)malloc((cfdCore->nE + 1) * sizeof(int));
	if(cfdCore->faceToElementsX == NULL) { fprintf(stderr, "[%d]: ERROR FTE <000>\n", rankMPI); exit(0); }
	for(i=0; i<=cfdCore->nE; i++) cfdCore->faceToElementsX[i] = 0;

	for(i=0; i<cfdCore->nICF; i++)
	{
		cfdCore->faceToElementsX[cfdCore->ICF[i].eL + 1] += 1;
		cfdCore->faceToElementsX[cfdCore->ICF[i].eR + 1] += 1;
	} // for i
	for(i=0; i<cfdCore->nBCF; i++)
	cfdCore->faceToElementsX[cfdCore->BCF[i].eL + 1] += 1;		
	for(i=2; i<=cfdCore->nE; i++) cfdCore->faceToElementsX[i] += cfdCore->faceToElementsX[i - 1];

	if(cfdCore->faceToElementsX[cfdCore->nE] != 2 * cfdCore->nICF + cfdCore->nBCF)
	{ fprintf(stderr, "[%d]: ERROR FTE <001>\n", rankMPI); exit(0); }

	cfdCore->faceToElementsA = (int *)malloc(cfdCore->faceToElementsX[cfdCore->nE] * sizeof(int));
	if(cfdCore->faceToElementsA == NULL) { fprintf(stderr, "[%d]: ERROR FTE <002>\n", rankMPI); exit(0); }
	for(i=0; i<cfdCore->faceToElementsX[cfdCore->nE]; i++) cfdCore->faceToElementsA[i] = 0;

	for(i=0; i<cfdCore->nICF; i++)
	{
		cfdCore->faceToElementsA[cfdCore->faceToElementsX[cfdCore->ICF[i].eL]] =    i + 1 ; cfdCore->faceToElementsX[cfdCore->ICF[i].eL] += 1;
		cfdCore->faceToElementsA[cfdCore->faceToElementsX[cfdCore->ICF[i].eR]] = - (i + 1); cfdCore->faceToElementsX[cfdCore->ICF[i].eR] += 1;
	} // for i
	for(i=0; i<cfdCore->nBCF; i++)
	{
		cfdCore->faceToElementsA[cfdCore->faceToElementsX[cfdCore->BCF[i].eL]] = cfdCore->nICF + i + 1 ; 
		cfdCore->faceToElementsX[cfdCore->BCF[i].eL] += 1;
	} // for i
		
	if(cfdCore->faceToElementsX[cfdCore->nE] != cfdCore->faceToElementsX[cfdCore->nE - 1]) { fprintf(stderr, "[%d]: ERROR FTE <003>\n", rankMPI); exit(0); }
	for(i=0; i<cfdCore->faceToElementsX[cfdCore->nE]; i++) if(cfdCore->faceToElementsA[i] == 0) { fprintf(stderr, "[%d]: ERROR FTE <004> %d\n", i, rankMPI); exit(0); }

	for(i=cfdCore->nE - 1; i>0; i--) cfdCore->faceToElementsX[i] = cfdCore->faceToElementsX[i - 1];
	cfdCore->faceToElementsX[0] = 0;

	MPI_Barrier(comm);
	if(rankMPI == 0) { printf("GEOMETRY >> FACE TO CELL MAP\n"); fflush(stdout); }
	MPI_Barrier(comm);
	// йюпрю опхмюдкефмнярх цпюмеи ъвеийюл

	// опнбепйю ценлерпхвеяйни йнмяепбюрхбмнярх он цпюмъл
	{
		double dGLOB[3] = {0.0, 0.0, 0.0};
		double dMAX[3]  = {0.0, 0.0, 0.0};

		for(i=0; i<cfdCore->nE; i++) if(cfdCore->cellLevels[i] < cfdCore->levelMax)
		{
			double dCUR[3] = {0.0, 0.0, 0.0};

			for(j=cfdCore->faceToElementsX[i]; j<cfdCore->faceToElementsX[i + 1]; j++)
			{
				int idFACE = cfdCore->faceToElementsA[j];

				if(idFACE < 0)
				{
					int index = - idFACE - 1;
					for(k=0; k<3; k++) dCUR[k] -= cfdCore->ICF[index].n[k] * cfdCore->ICF[index].fS;					
				}

				if(idFACE > 0)
				{
					int index = idFACE - 1;
					if(index >= cfdCore->nICF) for(k=0; k<3; k++) dCUR[k] += cfdCore->BCF[index - cfdCore->nICF].n[k] * cfdCore->BCF[index - cfdCore->nICF].fS;
					else for(k=0; k<3; k++) dCUR[k] += cfdCore->ICF[index].n[k] * cfdCore->ICF[index].fS;
				}
			} // for j
			
			for(k=0; k<3; k++) dCUR[k] = fabs(dCUR[k]);
			for(k=0; k<3; k++) if(dMAX[k] < dCUR[k]) dMAX[k] = dCUR[k];
		} // for i

		MPI_Allreduce(&dMAX[0], &dGLOB[0], 1, MPI_DOUBLE, MPI_MAX, comm);
		MPI_Allreduce(&dMAX[1], &dGLOB[1], 1, MPI_DOUBLE, MPI_MAX, comm);
		MPI_Allreduce(&dMAX[2], &dGLOB[2], 1, MPI_DOUBLE, MPI_MAX, comm);
		
		MPI_Barrier(comm);
		if(rankMPI == 0)
		{
			printf("\n");
			printf("CELLS ICF & BCF GEOMETRY TEST\n");
			printf("      X: %20e\n", dGLOB[0]);
			printf("      Y: %20e\n", dGLOB[1]);
			printf("      Z: %20e\n", dGLOB[2]);
			fflush(stdout);
		}
		MPI_Barrier(comm);
	}
	// опнбепйю ценлерпхвеяйни йнмяепбюрхбмнярх он цпюмъл

	// хмхжхюкхгюжхъ йнщттхжхемрнб онкхмнлхюкэмни пейнмярпсйжхх
	if(1)
	{
		int iCell, iFace;
		
		cfdCore->gaussX = (int *)malloc((cfdCore->nE + 1) * sizeof(int));
		if(cfdCore->gaussX == NULL) { fprintf(stderr, "[%d]: ERROR GAUSS <000>\n", rankMPI); exit(0); }
		for(i=0; i<=cfdCore->nE; i++) cfdCore->gaussX[i] = 0;

		for(i=0; i<cfdCore->nE; i++)	
			cfdCore->gaussX[i + 1] = cfdCore->faceToElementsX[i + 1] - cfdCore->faceToElementsX[i];
		for(i=0; i<cfdCore->nBCF; i++) cfdCore->gaussX[cfdCore->BCF[i].eL + 1] -= 1;
		for(i=2; i<=cfdCore->nE; i++) cfdCore->gaussX[i] += cfdCore->gaussX[i - 1];
		if(cfdCore->gaussX[cfdCore->nE] != cfdCore->nICF * 2) { fprintf(stderr, "[%d]: ERROR GAUSS <001>\n", rankMPI); exit(0); }

		cfdCore->gaussA = (int *)malloc(cfdCore->gaussX[cfdCore->nE] * sizeof(int));
		if(cfdCore->gaussA == NULL) { fprintf(stderr, "[%d]: ERROR GAUSS <002>\n", rankMPI); exit(0); }
		for(i=0; i<cfdCore->gaussX[cfdCore->nE]; i++) cfdCore->gaussA[i] = -1;

		cfdCore->kXYZC = (double *)malloc(cfdCore->nE * 6 * sizeof(double));
		if(cfdCore->kXYZC == NULL) { fprintf(stderr, "[%d]: ERROR GAUSS <003>\n", rankMPI); exit(0); }
		for(i=0; i<cfdCore->nE * 6; i++) cfdCore->kXYZC[i] = 0.0;
		cfdCore->kXYZS = (double *)malloc(cfdCore->gaussX[cfdCore->nE] * 3 * sizeof(double));
		if(cfdCore->kXYZS == NULL) { fprintf(stderr, "[%d]: ERROR GAUSS <004>\n", rankMPI); exit(0); }
		for(i=0; i<cfdCore->gaussX[cfdCore->nE] * 3; i++) cfdCore->kXYZS[i] = 0.0;

		for(iCell=0; iCell<cfdCore->nE; iCell++)
		{
			int offset = cfdCore->gaussX[iCell];
			int nST = 0;			

			for(iFace=cfdCore->faceToElementsX[iCell]; iFace<cfdCore->faceToElementsX[iCell + 1]; iFace++)
			{
				int idFACE = cfdCore->faceToElementsA[iFace];

				if(idFACE < 0)
				{
					int index = - idFACE - 1;
					cfdCore->gaussA[offset + nST] = cfdCore->ICF[index].eL;					
					for(k=0; k<3; k++) cfdCore->kXYZS[(offset + nST)*3 + k]  = - 0.5 * (cfdCore->ICF[index].n[k] * cfdCore->ICF[index].fS);
					for(k=0; k<3; k++) cfdCore->kXYZC[        iCell *6 + k] -=   0.5 * (cfdCore->ICF[index].n[k] * cfdCore->ICF[index].fS);
					
					nST += 1;
				}

				if(idFACE > 0)
				{
					int index = idFACE - 1;

					if(index < cfdCore->nICF)
					{
						cfdCore->gaussA[offset + nST] = cfdCore->ICF[index].eR;					
						for(k=0; k<3; k++) cfdCore->kXYZS[(offset + nST)*3 + k]  =   0.5 * (cfdCore->ICF[index].n[k] * cfdCore->ICF[index].fS);
						for(k=0; k<3; k++) cfdCore->kXYZC[        iCell *6 + k] +=   0.5 * (cfdCore->ICF[index].n[k] * cfdCore->ICF[index].fS);

						nST += 1;
					}
					else
					{
						int indexB = index - cfdCore->nICF;
						for(k=0; k<3; k++) cfdCore->kXYZC[iCell *6 + 0 + k] += 1.0 * (cfdCore->BCF[indexB].n[k] * cfdCore->BCF[indexB].fS);
					}
				}
			} // for iFace
			
			if(nST != cfdCore->gaussX[iCell + 1] - cfdCore->gaussX[iCell])
			{ fprintf(stderr, "[%d]: ERROR GAUSS <005>\n", rankMPI); exit(0); }

			for(i=cfdCore->gaussX[iCell]*3; i<cfdCore->gaussX[iCell + 1]*3; i++) 
				cfdCore->kXYZS[i] /= cfdCore->cellVOL[iCell];
			for(i=0; i<6; i++) cfdCore->kXYZC[iCell*6 + i] /= cfdCore->cellVOL[iCell];				
		} // for iCell

		MPI_Barrier(comm);
		if(rankMPI == 0) { printf("GEOMETRY >> SET AC2 GEOMETRY STRUCTURES\n"); fflush(stdout); }
		MPI_Barrier(comm);

		// опнбепйю ценлерпхвеяйни йнмяепбюрхбмнярэ пейнмярпсйжхнммшу ьюакнмнб
		{
			double geomGLOB[3] = {0.0, 0.0, 0.0};
			double geomMAX[3]  = {0.0, 0.0, 0.0};

			for(i=0; i<cfdCore->nE; i++) if(cfdCore->cellLevels[i] < cfdCore->levelMax)
			{
				double geomLOCAL[3] = {0.0, 0.0, 0.0};

				for(j=0; j<3; j++) geomLOCAL[j] = cfdCore->kXYZC[i*6 + j] + cfdCore->kXYZC[i*6 + j + 3];

				for(k=cfdCore->gaussX[i]; k<cfdCore->gaussX[i + 1]; k++)
					for(j=0; j<3; j++) geomLOCAL[j] += cfdCore->kXYZS[k*3 + j];

				for(j=0; j<3; j++) geomLOCAL[j] = fabs(geomLOCAL[j]);
				for(j=0; j<3; j++) if(geomMAX[j] < geomLOCAL[j]) geomMAX[j] = geomLOCAL[j];
			} // for i

			MPI_Allreduce(&geomMAX[0], &geomGLOB[0], 1, MPI_DOUBLE, MPI_MAX, comm);
			MPI_Allreduce(&geomMAX[1], &geomGLOB[1], 1, MPI_DOUBLE, MPI_MAX, comm);
			MPI_Allreduce(&geomMAX[2], &geomGLOB[2], 1, MPI_DOUBLE, MPI_MAX, comm);
		
			MPI_Barrier(comm);
			if(rankMPI == 0)
			{
				printf("CELLS GAUSS RECONSTRUCTION GEOMETRY TEST\n");
				printf("      X: %20e\n", geomGLOB[0]);
				printf("      Y: %20e\n", geomGLOB[1]);
				printf("      Z: %20e\n", geomGLOB[2]);
				printf("\n");
				fflush(stdout);
			}
			MPI_Barrier(comm);			
		}
		// опнбепйю ценлерпхвеяйни йнмяепбюрхбмнярэ пейнмярпсйжхнммшу ьюакнмнб
	}
	// хмхжхюкхгюжхъ йнщттхжхемрнб онкхмнлхюкэмни пейнмярпсйжхх

	// ондявер вхякю ярпнцн бмсрпеммху ъвеей
	{
		int minN, maxN;

		cfdCore->nEinterior = 0;
		for(i=0; i<cfdCore->nE; i++)
			if(cfdCore->cellLevels[i] != 0) break;
		cfdCore->nEinterior = i;

		MPI_Allreduce(&cfdCore->nEinterior, &minN, 1, MPI_INT, MPI_MIN, comm);
		MPI_Allreduce(&cfdCore->nEinterior, &maxN, 1, MPI_INT, MPI_MAX, comm);

		MPI_Barrier(comm);
		if(rankMPI == 0)
		{
			printf("INTERIOR CELLS NUMBER\n");
			printf(" MIN: %d\n", minN);
			printf(" MAX: %d\n", maxN);			
			printf("\n");
			fflush(stdout);
		}
		MPI_Barrier(comm);
	}
	// ондявер вхякю ярпнцн бмсрпеммху ъвеей

	// пюяопедекемхе цпюмеи он рхоюл
	cfdCore->icfOffset = (int *)malloc(cfdCore->levelMax * sizeof(int));
	if(cfdCore->icfOffset == NULL) exit(0);
	for(i=0; i<cfdCore->levelMax; i++) cfdCore->icfOffset[i] = -1;
	for(i=0; i<cfdCore->levelMax; i++)
	{
		for(j=0; j<cfdCore->nICF; j++)
			if( (cfdCore->cellLevels[cfdCore->ICF[j].eL] > i) && (cfdCore->cellLevels[cfdCore->ICF[j].eR] > i) )
				break;
		cfdCore->icfOffset[i] = j;
	} // for i
	for(i=0; i<cfdCore->levelMax; i++) if(cfdCore->icfOffset[i] == -1)
	{ fprintf(stderr, "[%d]: ERROR OFFSET ICF\n", rankMPI); fflush(stderr); exit(0); }

	cfdCore->bcfOffset = (int *)malloc(cfdCore->levelMax * sizeof(int));
	if(cfdCore->bcfOffset == NULL) exit(0);
	for(i=0; i<cfdCore->levelMax; i++) cfdCore->bcfOffset[i] = -1;
	for(i=0; i<cfdCore->levelMax; i++)
	{
		for(j=0; j<cfdCore->nBCF; j++)
			if(cfdCore->cellLevels[cfdCore->BCF[j].eL] > i)
				break;
		cfdCore->bcfOffset[i] = j;
	} // for i
	for(i=0; i<cfdCore->levelMax; i++) if(cfdCore->bcfOffset[i] == -1)
	{ fprintf(stderr, "[%d]: ERROR OFFSET BCF\n", rankMPI); fflush(stderr); exit(0); }
	// пюяопедекемхе цпюмеи он рхоюл

	return 0;
}
