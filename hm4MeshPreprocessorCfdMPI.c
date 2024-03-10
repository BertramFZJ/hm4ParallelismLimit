#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "typesMeshHeader.h"
#include "pm18MpiHeader.h"
#include "mpiTransferHost.h"

#include "iohm3DMPI.h"
#include "csrHmeshDualMPI.h"
#include "csrPartParMETISMPI.h"
#include "hm4PreprocessorMPI.h"
#include "processIntegerLIST.h"
#include "hm4DualGraphCSR.h"
#include "csrProcessCuthillMcKee.h"

#include "hm4MeshPreprocessorCfdMPI.h"

int hm4SerialMeshPreprocessorCfdMPI(char *fname, 
                                    int linkLevelsNumber, 
                                    gtypeHm4MeshAttributes *meshAttributes, 
                                    gtypeHm4MeshTopology *gMesh, gtypeHm4MeshTopology *hMesh,
                                    MPI_Comm comm)
{
	// ﾈﾄﾅﾍﾒﾈﾔﾈﾊﾀﾒﾎﾐ MPI-ﾏﾐﾎﾖﾅﾑﾑﾀ ﾈ ﾐﾀﾇﾌﾅﾐ ﾃﾐﾓﾏﾏﾛ
	int rankMPI, sizeMPI;
	// ﾈﾄﾅﾍﾒﾈﾔﾈﾊﾀﾒﾎﾐ MPI-ﾏﾐﾎﾖﾅﾑﾑﾀ ﾈ ﾐﾀﾇﾌﾅﾐ ﾃﾐﾓﾏﾏﾛ

	char *meshName; // ﾇﾀﾃﾎﾋﾎﾂﾎﾊ ﾑﾅﾒﾊﾈ

	// ﾐﾀﾂﾍﾎﾌﾅﾐﾍﾎﾅ ﾐﾀﾑﾏﾐﾅﾄﾅﾋﾅﾍﾈﾅ ﾒﾎﾏﾎﾋﾎﾃﾈﾈ ﾑﾅﾒﾎﾗﾍﾛﾕ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾃﾋﾎﾁﾀﾋﾜﾍﾎﾉ ﾑﾅﾒﾊﾈ
	int *eDIST, *xE, *aE;
	int *vDIST; double *CV;
	// ﾐﾀﾂﾍﾎﾌﾅﾐﾍﾎﾅ ﾐﾀﾑﾏﾐﾅﾄﾅﾋﾅﾍﾈﾅ ﾒﾎﾏﾎﾋﾎﾃﾈﾈ ﾑﾅﾒﾎﾗﾍﾛﾕ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾃﾋﾎﾁﾀﾋﾜﾍﾎﾉ ﾑﾅﾒﾊﾈ

	// ﾎﾏﾈﾑﾀﾍﾈﾅ ﾏﾎﾑﾒﾀﾍﾎﾂﾊﾈ ﾃﾐﾀﾍﾈﾗﾍﾛﾕ ﾓﾑﾋﾎﾂﾈﾉ ﾈ ﾒﾎﾏﾎﾋﾎﾃﾈﾈ ﾑﾅﾒﾊﾈ ﾍﾀ ﾏﾎﾂﾅﾐﾕﾍﾎﾑﾒﾈ ﾐﾀﾑﾗﾅﾒﾍﾎﾉ ﾎﾁﾋﾀﾑﾒﾈ
	int nFACE;
	char *faceNAME;
	int nEB, *EB;
	// ﾎﾏﾈﾑﾀﾍﾈﾅ ﾏﾎﾑﾒﾀﾍﾎﾂﾊﾈ ﾃﾐﾀﾍﾈﾗﾍﾛﾕ ﾓﾑﾋﾎﾂﾈﾉ ﾈ ﾒﾎﾏﾎﾋﾎﾃﾈﾈ ﾑﾅﾒﾊﾈ ﾍﾀ ﾏﾎﾂﾅﾐﾕﾍﾎﾑﾒﾈ ﾐﾀﾑﾗﾅﾒﾍﾎﾉ ﾎﾁﾋﾀﾑﾒﾈ

	// ﾐﾀﾑﾏﾐﾅﾄﾅﾋﾅﾍﾍﾛﾉ ﾄﾓﾀﾋﾜﾍﾛﾉ ﾃﾐﾀﾔ ﾃﾋﾎﾁﾀﾋﾜﾍﾎﾉ ﾑﾅﾒﾊﾈ
	int *xCSR, *aCSR;
	// ﾐﾀﾑﾏﾐﾅﾄﾅﾋﾅﾍﾍﾛﾉ ﾄﾓﾀﾋﾜﾍﾛﾉ ﾃﾐﾀﾔ ﾃﾋﾎﾁﾀﾋﾜﾍﾎﾉ ﾑﾅﾒﾊﾈ

	// ﾐﾀﾑﾏﾐﾅﾄﾅﾋﾅﾍﾍﾀﾟ ﾄﾅﾊﾎﾌﾏﾎﾇﾈﾖﾈﾟ ﾃﾋﾎﾁﾀﾋﾜﾍﾎﾉ ﾑﾅﾒﾊﾈ
	int nPART, *part;
	// ﾐﾀﾑﾏﾐﾅﾄﾅﾋﾅﾍﾍﾀﾟ ﾄﾅﾊﾎﾌﾏﾎﾇﾈﾖﾈﾟ ﾃﾋﾎﾁﾀﾋﾜﾍﾎﾉ ﾑﾅﾒﾊﾈ

	// ﾗﾈﾑﾋﾎ ﾈ ﾑﾏﾈﾑﾎﾊ ﾈﾍﾄﾅﾊﾑﾎﾂ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ, ﾂﾎﾘﾅﾄﾘﾈﾕ ﾂ HOST ﾄﾎﾌﾅﾍ
	int nePROC, *ePROC;
	// ﾗﾈﾑﾋﾎ ﾈ ﾑﾏﾈﾑﾎﾊ ﾈﾍﾄﾅﾊﾑﾎﾂ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ, ﾂﾎﾘﾅﾄﾘﾈﾕ ﾂ HOST ﾄﾎﾌﾅﾍ

	// ﾗﾈﾑﾋﾎ ﾈ ﾑﾏﾈﾑﾎﾊ ﾃﾋﾎﾁﾀﾋﾜﾍﾛﾕ ﾈﾍﾄﾅﾊﾑﾎﾂ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ, ﾏﾎﾏﾀﾂﾘﾈﾕ ﾂ ﾁﾓﾔﾅﾐﾍﾓﾞ ﾇﾎﾍﾓ HOST-ﾄﾎﾌﾅﾍﾀ
	int neLINK, *eLINK;
	// ﾗﾈﾑﾋﾎ ﾈ ﾑﾏﾈﾑﾎﾊ ﾃﾋﾎﾁﾀﾋﾜﾍﾛﾕ ﾈﾍﾄﾅﾊﾑﾎﾂ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ, ﾏﾎﾏﾀﾂﾘﾈﾕ ﾂ ﾁﾓﾔﾅﾐﾍﾓﾞ ﾇﾎﾍﾓ HOST-ﾄﾎﾌﾅﾍﾀ

    // ﾊﾀﾐﾒﾀ ﾏﾐﾈﾅﾌﾀ ﾄﾀﾍﾍﾛﾕ
	int *xRECV, *aRECV;
	// ﾊﾀﾐﾒﾀ ﾏﾐﾈﾅﾌﾀ ﾄﾀﾍﾍﾛﾕ

    // ﾒﾎﾏﾎﾋﾎﾃﾈﾟ ﾑﾅﾒﾎﾗﾍﾛﾕ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾈ ﾊﾎﾎﾐﾄﾈﾍﾀﾒﾛ ﾓﾇﾋﾎﾂ HOST ﾄﾎﾌﾅﾍﾀ
	int *xePROC, *aePROC;
	int nvPROC, *vPROC; double *cvPROC;
	// ﾒﾎﾏﾎﾋﾎﾃﾈﾟ ﾑﾅﾒﾎﾗﾍﾛﾕ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾈ ﾊﾎﾎﾐﾄﾈﾍﾀﾒﾛ ﾓﾇﾋﾎﾂ HOST ﾄﾎﾌﾅﾍﾀ

	MPI_Comm_rank(comm, &rankMPI);
	MPI_Comm_size(comm, &sizeMPI);
	  MPI_Barrier(comm);	

	ReadHybMeshNameBinaryMPI(fname, &meshName, stdout, comm); fflush(stdout);
	
	eDIST = (int *)malloc( (sizeMPI + 1) * sizeof(int));
	if(eDIST == NULL) { fprintf(stderr, "HM4GHP [%4d]: ERROR --> Can't allocate memory < %d * INT>\n", rankMPI, sizeMPI); MPI_Abort(comm, 0); }
	ReadHybMeshElementsBinaryMPI(fname, eDIST, &xE, &aE, stdout, comm);

	vDIST = (int *)malloc( (sizeMPI + 1) * sizeof(int));
	if(vDIST == NULL) { fprintf(stderr, "HM4GHP [%4d]: ERROR --> Can't allocate memory < %d * INT>\n", rankMPI, sizeMPI); MPI_Abort(comm, 0); }
	ReadHybMeshNodesBinaryMPI(fname, vDIST, &CV, stdout, comm);

	ReadHybMeshBoundaryBinaryMPI(fname, &nFACE, &faceNAME, &nEB, &EB, stdout, comm);

	fflush(stdout); fflush(stderr);
	MPI_Barrier(comm);
	if(rankMPI == 0)
	{ printf("\n"); printf("READ DISTRIBUTED TOPOLOGY FROM SERIAL FILE\n"); fflush(stdout); }
	MPI_Barrier(comm);

// ﾂﾛﾂﾎﾄ ﾒﾀﾁﾋﾈﾖ ﾐﾀﾑﾏﾐﾅﾄﾅﾋﾅﾍﾈﾟ ﾑﾅﾒﾎﾗﾍﾛﾕ ﾓﾇﾋﾎﾂ ﾈ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾍﾀ ﾊﾎﾍﾑﾎﾋﾜ
#if 1
	if(rankMPI == 0) 
	{ 
		int pflag, i;
			
		pflag = 0; fprintf(stdout, "eDIST = { %d ", eDIST[1] - eDIST[0]);
		for(i=1; i<sizeMPI - 1; i++)
		{
			if(eDIST[i+1] - eDIST[i] != eDIST[i] - eDIST[i-1]) { fprintf(stdout, ", %d ", eDIST[i+1] - eDIST[i]); pflag = 0; }
			else if(pflag == 0) { fprintf(stdout, ", ... "); pflag = 1; }
		} // for i
		fprintf(stdout, ", %d} eTOTAL = %d\n", eDIST[i+1] - eDIST[i], eDIST[sizeMPI]); fflush(stdout);

		pflag = 0; fprintf(stdout, "vDIST = { %d ", vDIST[1] - vDIST[0]);
		for(i=1; i<sizeMPI - 1; i++)
		{
			if(vDIST[i+1] - vDIST[i] != vDIST[i] - vDIST[i-1]) { fprintf(stdout, ", %d ", vDIST[i+1] - vDIST[i]); pflag = 0; }
			else if(pflag == 0) { fprintf(stdout, ", ... "); pflag = 1; }
		} // for i
		fprintf(stdout, ", %d} vTOTAL = %d\n", vDIST[i+1] - vDIST[i], vDIST[sizeMPI]); fflush(stdout);			
	}
	fflush(stdout); MPI_Barrier(comm);
#endif
	// ﾂﾛﾂﾎﾄ ﾒﾀﾁﾋﾈﾖ ﾐﾀﾑﾏﾐﾅﾄﾅﾋﾅﾍﾈﾟ ﾑﾅﾒﾎﾗﾍﾛﾕ ﾓﾇﾋﾎﾂ ﾈ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾍﾀ ﾊﾎﾍﾑﾎﾋﾜ

	// ﾈﾍﾈﾖﾈﾀﾋﾈﾇﾀﾖﾈﾟ ﾐﾀﾑﾏﾐﾅﾄﾅﾋﾅﾍﾍﾎﾃﾎ ﾄﾓﾀﾋﾜﾍﾎﾃﾎ ﾑﾅﾒﾎﾗﾍﾎﾃﾎ ﾃﾐﾀﾔﾀ
	{
		double startTIME, finishTIME, procedureTIME;

		MPI_Barrier(comm);
		startTIME = MPI_Wtime();
		if(rankMPI == 0) fprintf(stdout, "\n");

	    csrHmeshDualMPI(eDIST, xE, aE, &xCSR, &aCSR, NULL, comm);

		MPI_Barrier(comm);
		finishTIME = MPI_Wtime();
		procedureTIME = finishTIME - startTIME;

		{
			int totalEdges;
			double maxTIME;

			MPI_Allreduce(&procedureTIME, &maxTIME, 1, MPI_DOUBLE, MPI_MAX, comm);
			MPI_Allreduce(&xCSR[eDIST[rankMPI + 1] - eDIST[rankMPI]], &totalEdges, 1, MPI_INT, MPI_SUM, comm);
			totalEdges /= 2;

			if(rankMPI == 0) fprintf(stdout, "MESH DUAL CSR GRAPH edgesTOTAL = %d TIME = %.3lf\n\n", totalEdges, maxTIME);
		}

		fflush(stdout); MPI_Barrier(comm);
	}
	// ﾈﾍﾈﾖﾈﾀﾋﾈﾇﾀﾖﾈﾟ ﾐﾀﾑﾏﾐﾅﾄﾅﾋﾅﾍﾍﾎﾃﾎ ﾄﾓﾀﾋﾜﾍﾎﾃﾎ ﾑﾅﾒﾎﾗﾍﾎﾃﾎ ﾃﾐﾀﾔﾀ

	// ﾏﾐﾟﾌﾀﾟ ﾄﾅﾊﾎﾌﾏﾎﾇﾈﾖﾈﾟ ﾐﾀﾑﾏﾐﾅﾄﾅﾋﾅﾍﾍﾎﾃﾎ ﾄﾓﾀﾋﾜﾍﾎﾃﾎ ﾑﾅﾒﾎﾗﾍﾎﾃﾎ ﾃﾐﾀﾔﾀ
	{
		int nE = eDIST[rankMPI + 1] - eDIST[rankMPI];
		int i;
		
		nPART = sizeMPI;

		part = (int *)malloc( nE * sizeof(int));
		if(part == NULL) { fprintf(stderr, "HM4GHP [%4d]: ERROR --> Can't allocate memory < %d * INT>\n", rankMPI, nE); MPI_Abort(comm, 0); }
		for(i=0; i<nE; i++) part[i] = -1;		

		// ﾄﾅﾊﾎﾌﾏﾎﾇﾈﾖﾈﾟ ﾃﾐﾀﾔﾀ ﾁﾅﾇ ﾂﾅﾑﾎﾂ
		if(0)
		csrPartParMetisV4NoWeightsMPI(eDIST, xCSR, aCSR, nPART, part, stderr, comm);
		// ﾄﾅﾊﾎﾌﾏﾎﾇﾈﾖﾈﾟ ﾃﾐﾀﾔﾀ ﾑ ﾂﾅﾑﾀﾌﾈ ﾂﾅﾐﾘﾈﾍ, ﾐﾀﾂﾍﾛﾌﾈ ﾗﾈﾑﾋﾓ ﾎﾏﾈﾐﾀﾞﾙﾈﾕﾑﾟ ﾍﾀ ﾂﾅﾐﾘﾈﾍﾓ ﾐﾅﾁﾅﾐ ﾃﾐﾀﾔﾀ
		if(1)
		{
			// CASA
			fprintf(stderr, "CASA: --- WeightDecomposition ---\n");
			csrPartParMetisV4NodesWeightsByCsrEdgesMPI(eDIST, xCSR, aCSR, nPART, part, stderr, comm);
			fprintf(stderr, "CASA: *** WeightDecomposition ***\n");
			// CASA
		}		


		fflush(stdout); fflush(stderr); MPI_Barrier(comm);
	}
	// ﾏﾐﾟﾌﾀﾟ ﾄﾅﾊﾎﾌﾏﾎﾇﾈﾖﾈﾟ ﾐﾀﾑﾏﾐﾅﾄﾅﾋﾅﾍﾍﾎﾃﾎ ﾄﾓﾀﾋﾜﾍﾎﾃﾎ ﾑﾅﾒﾎﾗﾍﾎﾃﾎ ﾃﾐﾀﾔﾀ

	// ﾑﾎﾑﾒﾀﾂﾋﾅﾍﾈﾅ ﾑﾏﾈﾑﾊﾎﾂ ﾈﾍﾄﾅﾊﾑﾎﾂ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾑﾎﾁﾑﾒﾂﾅﾍﾍﾎﾉ ﾏﾎﾄﾎﾁﾋﾀﾑﾒﾈ ﾏﾐﾎﾖﾅﾑﾑﾀ
	{
		double startTIME, finishTIME;

		startTIME = MPI_Wtime();
		hm4CreateProcessElementsListMPI(eDIST, part, &nePROC, &ePROC, NULL, comm);
		MPI_Barrier(comm);
		finishTIME = MPI_Wtime();
		if(rankMPI == 0) { fprintf(stdout, "hm4CreateProcessElementsListMPI TIME = %.2lf\n", finishTIME - startTIME); fflush(stdout); }
		fflush(stdout); MPI_Barrier(comm);
	}
	// ﾑﾎﾑﾒﾀﾂﾋﾅﾍﾈﾅ ﾑﾏﾈﾑﾊﾎﾂ ﾈﾍﾄﾅﾊﾑﾎﾂ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾑﾎﾁﾑﾒﾂﾅﾍﾍﾎﾉ ﾏﾎﾄﾎﾁﾋﾀﾑﾒﾈ ﾏﾐﾎﾖﾅﾑﾑﾀ

	// ﾈﾍﾈﾖﾈﾀﾋﾈﾇﾀﾖﾈﾟ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ-ﾑﾂﾟﾇﾅﾉ ﾄﾎ ﾇﾀﾄﾀﾍﾍﾎﾃﾎ ﾓﾐﾎﾂﾍﾟ
	{
		double startTIME, finishTIME;

		int iL;
		
		// ﾏﾅﾐﾂﾛﾉ ﾓﾐﾎﾂﾅﾍﾜ ﾑﾂﾟﾇﾍﾎﾑﾒﾈ
		// ﾊﾀﾊ ﾂﾑﾅ ﾝﾒﾎ ﾐﾀﾁﾎﾒﾀﾅﾒ ﾑ ﾂﾛﾐﾎﾆﾄﾅﾍﾍﾛﾌﾈ ﾌﾀﾑﾑﾈﾂﾀﾌﾈ (ﾍﾓﾋﾅﾂﾀﾟ ﾄﾋﾈﾍﾀ ﾈ ﾒ.ﾄ.)
		MPI_Barrier(comm);
		startTIME = MPI_Wtime();
		hm4InitElementsLinksMPI(eDIST, xCSR, aCSR, nePROC, ePROC, &neLINK, &eLINK, NULL, comm);
		MPI_Barrier(comm);
		finishTIME = MPI_Wtime();
		if(rankMPI == 0) fprintf(stderr, "hm4InitElementsLinksMPI TIME = %.2lf\n", finishTIME - startTIME); fflush(stderr);
		
		// ﾏﾎﾑﾋﾅﾄﾓﾞﾙﾈﾅ ﾓﾐﾎﾂﾍﾈ ﾑﾂﾟﾇﾍﾎﾑﾒﾈ
		for(iL=1; iL<linkLevelsNumber; iL++)
		{
			int neL, *eL;

			MPI_Barrier(comm);
			startTIME = MPI_Wtime();
			hm4InitElementsLinksMPI(eDIST, xCSR, aCSR, neLINK, eLINK, &neL, &eL, NULL, comm);
			MPI_Barrier(comm);
			finishTIME = MPI_Wtime();
			if(rankMPI == 0) { fprintf(stdout, "hm4InitElementsLinksMPI TIME = %.2lf\n", finishTIME - startTIME); fflush(stdout); }

			DeleteListFromOrderedList(nePROC, ePROC, neL, eL, &neL);
			DeleteListFromOrderedList(neLINK, eLINK, neL, eL, &neL);

			if(neL != 0)
			{
				{
					int *irPTR = (int *)realloc(eLINK, (neLINK + neL) * sizeof(int));
					if(irPTR == NULL) { fprintf(stderr, "HM4GHP --> %4d: *** ERROR *** Can't reallocate memory <%d * INT>\n", rankMPI, neL + neLINK); MPI_Abort(comm, 0); }
					if(irPTR != eLINK) eLINK = irPTR;
				}

				{
					int i;
					for(i=0; i<neL; i++) eLINK[neLINK + i] = eL[i];
					neLINK += neL;

					free(eL); neL = 0; eL = NULL;
				}				
			} // if
		} // for iL
		// ﾏﾎﾑﾋﾅﾄﾓﾞﾙﾈﾅ ﾓﾐﾎﾂﾍﾈ ﾑﾂﾟﾇﾍﾎﾑﾒﾈ

		// ﾎﾁﾚﾅﾄﾈﾍﾅﾍﾈﾅ ﾑﾏﾈﾑﾊﾎﾂ
		{
			int i;

			int *irPTR = (int *)realloc(ePROC, (neLINK + nePROC) * sizeof(int));
			if(irPTR == NULL) { fprintf(stderr, "HM4GHP --> %4d: *** ERROR *** Can't reallocate memory <%d * INT>\n", rankMPI, nePROC + neLINK); MPI_Abort(comm, 0); }
			if(irPTR != eLINK) ePROC = irPTR;

			for(i=0; i<neLINK; i++) ePROC[nePROC + i] = eLINK[i];
			free(eLINK); eLINK = ePROC + nePROC;
		}
		// ﾎﾁﾚﾅﾄﾈﾍﾅﾍﾈﾅ ﾑﾏﾈﾑﾊﾎﾂ

		// ﾑﾒﾀﾒﾈﾑﾒﾈﾊﾀ
		{
			int minLINKS, maxLINKS;
			
			MPI_Allreduce(&neLINK, &minLINKS, 1, MPI_INT, MPI_MIN, comm);
			MPI_Allreduce(&neLINK, &maxLINKS, 1, MPI_INT, MPI_MAX, comm);
						
			if(rankMPI == 0) 
			{
				fprintf(stdout, "\n");
				fprintf(stdout, "   NUMBER OF LINKS LEVELS: %d\n", linkLevelsNumber);
				fprintf(stdout, "MINIMUM LINKS LIST LENTGH: %8d\n", minLINKS);
				fprintf(stdout, "MAXIMUM LINKS LIST LENTGH: %8d\n", maxLINKS);
				fprintf(stdout, "\n"); fflush(stdout);
			}
			MPI_Barrier(comm);			
		}
		// ﾑﾒﾀﾒﾈﾑﾒﾈﾊﾀ
	}
	// ﾈﾍﾈﾖﾈﾀﾋﾈﾇﾀﾖﾈﾟ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ-ﾑﾂﾟﾇﾅﾉ ﾄﾎ ﾇﾀﾄﾀﾍﾍﾎﾃﾎ ﾓﾐﾎﾂﾍﾟ

    // ﾈﾍﾈﾖﾈﾀﾋﾈﾇﾀﾖﾈﾟ ﾏﾐﾈﾍﾀﾄﾋﾅﾆﾍﾎﾑﾒﾈ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾁﾓﾔﾅﾐﾍﾎﾉ ﾇﾎﾍﾛ HOST ﾄﾎﾌﾅﾍﾀ ﾗﾅﾐﾅﾇ ﾑﾒﾐﾓﾊﾒﾓﾐﾓ ﾏﾐﾈﾅﾌﾀ ﾄﾀﾍﾍﾛﾕ
	{
		double startTIME, finishTIME;

		int i;

        xRECV = (int *)malloc( (sizeMPI + 1) * sizeof(int));
        if(xRECV == NULL) { fprintf(stderr, "HM4GHP [%4d]: ERROR --> Can't allocate memory <%d * INT>\n", rankMPI, sizeMPI + 1); MPI_Abort(comm, 0); }
		aRECV = (int *)malloc( neLINK * sizeof(int));
        if(aRECV == NULL) { fprintf(stderr, "HM4GHP [%4d]: ERROR --> Can't allocate memory <%d * INT>\n", rankMPI, neLINK); MPI_Abort(comm, 0); }

		MPI_Barrier(comm); startTIME = MPI_Wtime();
		hm4InitProcessRecvListMPI(eDIST, part, neLINK, eLINK, xRECV, aRECV, NULL, comm);
		finishTIME = MPI_Wtime(); MPI_Barrier(comm);
        if(rankMPI == 0) { fprintf(stdout, "hm4InitProcessRecvListMPI TIME = %.2lf\n", finishTIME - startTIME); fflush(stdout); }

		// ﾒﾐﾀﾍﾑﾔﾎﾐﾌﾀﾖﾈﾟ ﾑﾏﾈﾑﾊﾀ ﾑﾂﾟﾇﾅﾉ ﾂ ﾑﾎﾎﾒﾂﾅﾒﾑﾒﾂﾈﾈ ﾑ ﾌﾀﾑﾑﾈﾂﾎﾌ ﾏﾐﾈﾅﾌﾀ ﾄﾀﾍﾍﾛﾕ
		for(i=0; i<neLINK; i++) ePROC[nePROC + i] = aRECV[i];
		MPI_Barrier(comm);		
	}
	// ﾈﾍﾈﾖﾈﾀﾋﾈﾇﾀﾖﾈﾟ ﾏﾐﾈﾍﾀﾄﾋﾅﾆﾍﾎﾑﾒﾈ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾁﾓﾔﾅﾐﾍﾎﾉ ﾇﾎﾍﾛ HOST ﾄﾎﾌﾅﾍﾀ ﾗﾅﾐﾅﾇ ﾑﾒﾐﾓﾊﾒﾓﾐﾓ ﾏﾐﾈﾅﾌﾀ ﾄﾀﾍﾍﾛﾕ

    // ﾔﾎﾐﾌﾈﾐﾎﾂﾀﾍﾈﾅ ﾌﾀﾑﾑﾈﾂﾎﾂ ﾑ ﾎﾏﾈﾑﾀﾍﾈﾅﾌ ﾒﾎﾏﾎﾋﾎﾃﾈﾈ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾈ ﾊﾎﾎﾐﾄﾈﾍﾀﾒﾀﾌﾈ ﾓﾇﾋﾎﾂ HOST ﾄﾎﾌﾅﾍﾀ
	{
		double startTIME, finishTIME;

		MPI_Barrier(comm); startTIME = MPI_Wtime();
		hm4InitElementsTopologyMPI(eDIST, xE, aE, nePROC + neLINK, ePROC, &xePROC, &aePROC, NULL, comm);
		MPI_Barrier(comm); finishTIME = MPI_Wtime();
        if(rankMPI == 0) { fprintf(stdout, "hm4InitElementsTopologyMPI TIME = %.2lf\n", finishTIME - startTIME); fflush(stdout); }

		MPI_Barrier(comm); startTIME = MPI_Wtime();
		hm4CreateProcessNodesListAndCoordinatesMPI(nePROC + neLINK, xePROC, aePROC, vDIST, CV, &nvPROC, &vPROC, &cvPROC, 3, NULL, comm);
		MPI_Barrier(comm); finishTIME = MPI_Wtime();
        if(rankMPI == 0) { fprintf(stdout, "hm4CreateProcessNodesListAndCoordinatesMPI TIME = %.2lf\n", finishTIME - startTIME); fflush(stdout); }

		// ﾑﾒﾀﾒﾈﾑﾒﾈﾊﾀ ﾏﾎ ﾑﾓﾌﾌﾀﾐﾍﾎﾌﾓ ﾗﾈﾑﾋﾓ ﾓﾇﾋﾎﾂ
		{
			int minNV, maxNV;
			
			MPI_Allreduce(&nvPROC, &minNV, 1, MPI_INT, MPI_MIN, comm);
			MPI_Allreduce(&nvPROC, &maxNV, 1, MPI_INT, MPI_MAX, comm);
						
            if(rankMPI == 0) 
			{
				fprintf(stdout, "\n");
				fprintf(stdout, "MINIMUM HOST DOMEN NODES NUMBER: %8d\n", minNV);
				fprintf(stdout, "MAXIMUM HOST DOMEN NODES NUMBER: %8d\n", maxNV);
				fprintf(stdout, "\n"); fflush(stdout);
			}
			MPI_Barrier(comm);
		}
		// ﾑﾒﾀﾒﾈﾑﾒﾈﾊﾀ ﾏﾎ ﾑﾓﾌﾌﾀﾐﾍﾎﾌﾓ ﾗﾈﾑﾋﾓ ﾓﾇﾋﾎﾂ

		// ﾏﾅﾐﾅﾕﾎﾄ ﾊ ﾋﾎﾊﾀﾋﾜﾍﾎﾉ ﾈﾍﾄﾅﾊﾑﾀﾖﾈﾈ ﾓﾇﾋﾎﾂ ﾂ ﾌﾀﾑﾑﾈﾂﾅ ﾑﾎ ﾑﾏﾈﾑﾊﾀﾌﾈ ﾈﾍﾄﾅﾊﾑﾎﾂ ﾂﾅﾐﾘﾈﾍ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ
		{
			int i;

			// ﾏﾐﾎﾂﾅﾐﾊﾀ ﾓﾏﾎﾐﾟﾄﾎﾗﾅﾍﾍﾎﾑﾒﾈ
			for(i=1; i<nvPROC; i++)
			if(vPROC[i] <= vPROC[i - 1]) 
            { fprintf(stderr, "HM4GHP [%4d]: ERROR --> 2378\n", rankMPI); MPI_Abort(comm, 0); }

			for(i=0; i<xePROC[nePROC + neLINK]; i++)
			{
				int pos = DichotomySearchElementInOrderIntegerListUnit(vPROC, nvPROC, aePROC[i]);
				
				if(pos == -1)
                { fprintf(stderr, "HM4GHP [%4d]: ERROR --> 3455\n", rankMPI); MPI_Abort(comm, 0); }

				aePROC[i] = pos;
			} // for i

			MPI_Barrier(comm);
            if(rankMPI == 0) { fprintf(stdout, "SET LOCAL NODES INDEXES\n"); fflush(stdout); }
			MPI_Barrier(comm);
		}
		// ﾏﾅﾐﾅﾕﾎﾄ ﾊ ﾋﾎﾊﾀﾋﾜﾍﾎﾉ ﾈﾍﾄﾅﾊﾑﾀﾖﾈﾈ ﾓﾇﾋﾎﾂ ﾂ ﾌﾀﾑﾑﾈﾂﾅ ﾑﾎ ﾑﾏﾈﾑﾊﾀﾌﾈ ﾈﾍﾄﾅﾊﾑﾎﾂ ﾂﾅﾐﾘﾈﾍ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ
	}
	// ﾔﾎﾐﾌﾈﾐﾎﾂﾀﾍﾈﾅ ﾌﾀﾑﾑﾈﾂﾎﾂ ﾑ ﾎﾏﾈﾑﾀﾍﾈﾅﾌ ﾒﾎﾏﾎﾋﾎﾃﾈﾈ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾈ ﾊﾎﾎﾐﾄﾈﾍﾀﾒﾀﾌﾈ ﾓﾇﾋﾎﾂ HOST ﾄﾎﾌﾅﾍﾀ

    // ﾇﾀﾏﾎﾋﾍﾅﾍﾈﾅ ﾏﾎﾋﾅﾉ ﾑﾒﾐﾓﾊﾒﾓﾐ ﾑ ﾎﾏﾈﾑﾀﾍﾈﾅﾌ ﾒﾎﾏﾎﾋﾎﾃﾈﾈ ﾃﾋﾎﾁﾀﾋﾜﾍﾎﾉ ﾑﾅﾒﾊﾈ ﾈ ﾅﾅ ﾀﾒﾐﾈﾁﾓﾒﾎﾂ
    meshAttributes->meshName = meshName;
    meshAttributes->meshNodesNumber = vDIST[sizeMPI];
    meshAttributes->meshElementsNumber = eDIST[sizeMPI];
    meshAttributes->meshBoundaryFacesNumber = nEB;
    meshAttributes->boundarySurfaceNames = faceNAME;
    meshAttributes->meshBoundarySurfacesNumber = nFACE;
	
    gMesh->meshNodesNumber    = vDIST[sizeMPI];
	gMesh->meshNodesCoordinates = NULL;
    gMesh->meshElementsNumber = eDIST[sizeMPI];	
	gMesh->meshEX = gMesh->meshEA = NULL;
	gMesh->partitionUP = gMesh->partitionDOWN = NULL;
	gMesh->indexesUP = gMesh->indexesDOWN = NULL;
	gMesh->meshBoundaryFacesNumber = nEB; gMesh->meshBoundaryFacesTopology = NULL;	
	// ﾇﾀﾏﾎﾋﾍﾅﾍﾈﾅ ﾏﾎﾋﾅﾉ ﾑﾒﾐﾓﾊﾒﾓﾐ ﾑ ﾎﾏﾈﾑﾀﾍﾈﾅﾌ ﾒﾎﾏﾎﾋﾎﾃﾈﾈ ﾃﾋﾎﾁﾀﾋﾜﾍﾎﾉ ﾑﾅﾒﾊﾈ ﾈ ﾅﾅ ﾀﾒﾐﾈﾁﾓﾒﾎﾂ

    // ﾇﾀﾏﾎﾋﾍﾅﾍﾈﾅ ﾏﾎﾋﾅﾉ ﾑﾒﾐﾓﾊﾒﾓﾐﾛ ﾑ ﾎﾏﾈﾑﾀﾍﾈﾅﾌ ﾑﾅﾒﾊﾈ HOST ﾄﾎﾌﾅﾍﾀ
	{
		int i, k;		

		hMesh->meshNodesNumber = nvPROC;
		hMesh->meshNodesCoordinates = cvPROC;

		hMesh->meshElementsNumber = nePROC + neLINK;
		hMesh->meshEX = xePROC; hMesh->meshEA = aePROC;

		hMesh->partitionDOWN = NULL; hMesh->indexesDOWN = NULL;

		hMesh->indexesUP = ePROC;
		hMesh->partitionUP = (int *)malloc(hMesh->meshElementsNumber * sizeof(int)); if(hMesh->partitionUP == NULL) exit(0);
        for(i=0; i<nePROC; i++) hMesh->partitionUP[i] = rankMPI;
        for(k=0; k<sizeMPI; k++)
			for(i=xRECV[k]; i<xRECV[k + 1]; i++) hMesh->partitionUP[i + nePROC] = k;
        		
		{
			double startTIME, finishTIME;

			int *indexesUP = (int *)malloc(hMesh->meshElementsNumber * 2 * sizeof(int));			
			if(indexesUP == NULL) exit(0);

			startTIME = MPI_Wtime();

			for(i=0; i<hMesh->meshElementsNumber; i++)
			indexesUP[i*2    ] = hMesh->indexesUP[i],
			indexesUP[i*2 + 1] =                  i ;			
			wsortIntegerListSizeUnit(indexesUP, 0, hMesh->meshElementsNumber - 1, 2);

			hMesh->meshBoundaryFacesNumber = 0;
			for(i=0; i<nEB; i++)
			{
				k = DichotomySearchElementInOrderIntegerListSizeUnit(indexesUP, hMesh->meshElementsNumber, EB[i*3], 2);
				if(k != -1)	k = indexesUP[k*2 + 1];			

				if(k != -1)
				{
					if(hMesh->indexesUP[k] != EB[i*3]) 
					{ fprintf(stderr, "EEEEEEE\n"); fflush(stderr); exit(0); }

					EB[hMesh->meshBoundaryFacesNumber*3 + 0] = k;
					EB[hMesh->meshBoundaryFacesNumber*3 + 1] = EB[i*3 + 1];
					EB[hMesh->meshBoundaryFacesNumber*3 + 2] = EB[i*3 + 2];
					hMesh->meshBoundaryFacesNumber += 1;
				}				
			} // for i

			free(indexesUP);

			finishTIME = MPI_Wtime();
			if(0) { printf("BOUNDARY SEARCH TIME = %g\n", finishTIME - startTIME); fflush(stdout); }
			MPI_Barrier(comm);

			if(hMesh->meshBoundaryFacesNumber == 0)
			{ free(EB); }
			else
			{
				if(hMesh->meshBoundaryFacesNumber != nEB)
				{
					int *irPTR = (int *)realloc(EB, hMesh->meshBoundaryFacesNumber * 3 * sizeof(int));
                    if(irPTR == NULL) { fprintf(stderr, "HM4GHP --> %4d: *** ERROR *** Can't reallocate memory <%d * INT>\n", rankMPI, hMesh->meshBoundaryFacesNumber); MPI_Abort(comm, 0); }
					if(irPTR != EB) EB = irPTR;
				}

				hMesh->meshBoundaryFacesTopology = EB;
			}
		}

		MPI_Barrier(comm);
        if(rankMPI == 0) 
		{ fprintf(stdout, "EXTRACT HOST MESH TOPOLOGY\n"); fflush(stdout); }
		MPI_Barrier(comm);
	}
	// ﾇﾀﾏﾎﾋﾍﾅﾍﾈﾅ ﾏﾎﾋﾅﾉ ﾑﾒﾐﾓﾊﾒﾓﾐﾛ ﾑ ﾎﾏﾈﾑﾀﾍﾈﾅﾌ ﾑﾅﾒﾊﾈ HOST ﾄﾎﾌﾅﾍﾀ

    // ﾎﾑﾂﾎﾁﾎﾆﾄﾅﾍﾈﾅ ﾂﾑﾏﾎﾌﾎﾃﾀﾒﾅﾋﾜﾍﾎﾉ ﾏﾀﾌﾟﾒﾈ
	{
		free(eDIST); free(xE); free(aE);
		free(vDIST); free(CV);
		free(xCSR); free(aCSR);
		free(part);

		MPI_Barrier(comm);
        if(rankMPI == 0) 
		{ fprintf(stdout, "FREE GLOBAL MESH MEMORY\n"); fflush(stdout); }
		MPI_Barrier(comm);
	}
	// ﾎﾑﾂﾎﾁﾎﾆﾄﾅﾍﾈﾅ ﾂﾑﾏﾎﾌﾎﾃﾀﾒﾅﾋﾜﾍﾎﾉ ﾏﾀﾌﾟﾒﾈ

	// ﾂﾛﾂﾎﾄ ﾏﾀﾐﾀﾌﾅﾒﾐﾎﾂ ﾃﾋﾎﾁﾀﾋﾜﾍﾎﾉ ﾑﾅﾒﾊﾈ
	MPI_Barrier(comm);
    if(rankMPI == 0)
	{
		int i;

		printf("\n");
		printf("******************************** MESH ********************************\n");
		printf("\n");
		printf("                         FILE: %s\n", fname);
        printf("                    MESH NAME: %s\n", meshAttributes->meshName);
		printf("            MESH NODES NUMBER: %15d\n", meshAttributes->meshNodesNumber);
		printf("         MESH ELEMENTS NUMBER: %15d\n", meshAttributes->meshElementsNumber);
		printf("\n");
		printf("   MESH BOUNDARY FACES NUMBER: %15d\n", meshAttributes->meshBoundaryFacesNumber);
		printf("MESH BOUNDARY SURFACES NUMBER: %15d\n", meshAttributes->meshBoundarySurfacesNumber);
		for(i=0; i<meshAttributes->meshBoundarySurfacesNumber; i++) 
		printf("                   SURFACE %2d: %s\n", i, meshAttributes->boundarySurfaceNames + i*81);
		printf("\n");
		printf("**********************************************************************\n");
		printf("\n"); fflush(stdout);
	}
	MPI_Barrier(comm);
	// ﾂﾛﾂﾎﾄ ﾏﾀﾐﾀﾌﾅﾒﾐﾎﾂ ﾃﾋﾎﾁﾀﾋﾜﾍﾎﾉ ﾑﾅﾒﾊﾈ
  
	// ﾑﾒﾀﾒﾈﾑﾒﾈﾊﾀ ﾏﾎ ﾑﾅﾒﾊﾀﾌ HOST ﾄﾎﾌﾅﾍﾎﾂ
    if(sizeMPI > 1)
	{
        int ehAV = gMesh->meshElementsNumber / sizeMPI;
		int ehMIN, ehMAX;

		MPI_Barrier(comm);
		MPI_Allreduce(&hMesh->meshElementsNumber, &ehMIN, 1, MPI_INT, MPI_MIN, comm);
		MPI_Allreduce(&hMesh->meshElementsNumber, &ehMAX, 1, MPI_INT, MPI_MAX, comm);
        if(rankMPI == 0)
		{
			printf("HOSTS MESH PARAMETERS\n");
			printf("eMIN: %15d (%8.2lf %%)\n", ehMIN, 100.0 * (double)ehMIN/(double)ehAV - 100.0);
			printf("eMAX: %15d (%8.2lf %%)\n", ehMAX, 100.0 * (double)ehMAX/(double)ehAV - 100.0);
			printf("\n"); fflush(stdout);			
		}
		MPI_Barrier(comm);
	}
	// ﾑﾒﾀﾒﾈﾑﾒﾈﾊﾀ ﾏﾎ ﾑﾅﾒﾊﾀﾌ HOST ﾄﾎﾌﾅﾍﾎﾂ

	MPI_Barrier(comm);
    if(rankMPI == 0)
	{ printf("FINISH hm4SerialMeshPreprocessorCfdMPI\n"); fflush(stdout); }
	MPI_Barrier(comm);

	return 0;
}

int hm4RotateHostMeshCfdMPI(gtypeHm4MeshTopology *hMesh, MPI_Comm comm)
{
    int rankMPI, sizeMPI;
    int nePROC = 0;

    int *dualX, *dualA;
    int firstNodeIndex = -1;

	int *baseIndexFromRotateIndex = NULL;
	int *rotateIndexFromBaseIndex = NULL;

	int *meshEX, *meshEA;
	int *meshBoundaryFacesTopology;
	int *partitionUP;
	int *indexesUP;

    int i;

    MPI_Comm_rank(comm, &rankMPI);
	MPI_Comm_size(comm, &sizeMPI);
	MPI_Barrier(comm);

	// ﾏﾎﾄﾑﾗﾅﾒ ﾗﾈﾑﾋﾀ ﾑﾎﾁﾑﾒﾂﾅﾍﾍﾛﾕ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾐﾀﾑﾗﾅﾒﾍﾎﾉ ﾎﾁﾋﾀﾑﾒﾈ ﾏﾐﾎﾖﾅﾑﾑﾀ
    for(i=0; i<hMesh->meshElementsNumber; i++)
		if(hMesh->partitionUP[i] != rankMPI)
        {
			nePROC = i;
            break;
        }

	// ﾃﾅﾍﾅﾐﾀﾖﾈﾟ ﾄﾓﾀﾋﾜﾍﾎﾃﾎ ﾃﾐﾀﾔﾀ ﾑﾎﾁﾑﾒﾂﾅﾍﾍﾛﾕ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾐﾀﾑﾗﾅﾒﾍﾎﾉ ﾎﾁﾋﾀﾑﾒﾈ
	hm4MeshDualGraphCSRSerial(nePROC, hMesh->meshEX, hMesh->meshEA, &dualX, &dualA);
	MPI_Barrier(comm);
	if(rankMPI == 0) { printf("Create local csr graphs\n"); fflush(stdout); }
	MPI_Barrier(comm);

	// ﾐﾅﾈﾍﾄﾅﾊﾑﾀﾖﾈﾟ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾀﾋﾃﾎﾐﾈﾒﾌﾎﾌ ﾊﾀﾒﾕﾈﾋﾋﾀ-ﾌﾀﾊﾊﾈ
#if 0
	// ﾎﾏﾐﾅﾄﾅﾋﾅﾍﾈﾅ ﾈﾍﾄﾅﾊﾑﾀ ﾏﾑﾅﾂﾄﾎﾏﾅﾐﾈﾔﾅﾐﾈﾉﾍﾎﾉ ﾂﾅﾐﾘﾈﾍﾛ ﾃﾐﾀﾔﾀ ﾏﾎﾄﾎﾁﾋﾀﾑﾒﾈ ﾑﾎﾁﾑﾒﾂﾅﾍﾍﾛﾕ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ
	firstNodeIndex = PseudoPeripheralCsrGraphVertex(nePROC, dualX, dualA);
	MPI_Barrier(comm);
    if(rankMPI == 0) { printf("Set pseudo peripheral csr graph vertex\n"); fflush(stdout); }
    MPI_Barrier(comm);

	// ﾐﾅﾈﾍﾄﾅﾊﾑﾀﾖﾈﾟ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾀﾋﾃﾎﾐﾈﾒﾌﾎﾌ ﾊﾀﾒﾕﾈﾋﾋﾀ-ﾌﾀﾊﾊﾈ
	baseIndexFromRotateIndex = (int *)malloc(hMesh->meshElementsNumber * sizeof(int));
    if(baseIndexFromRotateIndex == NULL) exit(0);
    for(i=0; i<hMesh->meshElementsNumber; i++) baseIndexFromRotateIndex[i] = i;
	csrGraphRenameIndexCuthillMcKee(nePROC, firstNodeIndex, dualX, dualA, baseIndexFromRotateIndex);
    MPI_Barrier(comm);
    if(rankMPI == 0) { printf("Set Cuthill-McKee indexes\n"); fflush(stdout); }
    MPI_Barrier(comm);
	// ﾏﾅﾐﾅﾕﾎﾄ ﾊ ﾎﾁﾐﾀﾒﾍﾎﾉ ﾈﾍﾄﾅﾊﾑﾀﾖﾈﾈ ﾊﾀﾒﾕﾈﾋﾋﾀ-ﾌﾀﾊﾊﾈ
	{
		int *index = (int *)malloc(nePROC * sizeof(int));
		if(index == NULL) exit(0);

		for(i=0; i<nePROC; i++) index[i] = baseIndexFromRotateIndex[nePROC - 1 - i];
		for(i=0; i<nePROC; i++) baseIndexFromRotateIndex[i] = index[i];

		free(index);

		MPI_Barrier(comm);
		if(rankMPI == 0) { printf("Set reverse indexes\n"); fflush(stdout); }
		MPI_Barrier(comm);
	} // if
#endif
	// ﾐﾅﾈﾍﾄﾅﾊﾑﾀﾖﾈﾟ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾀﾋﾃﾎﾐﾈﾒﾌﾎﾌ ﾊﾀﾒﾕﾈﾋﾋﾀ-ﾌﾀﾊﾊﾈ

	// ﾐﾅﾈﾍﾄﾅﾊﾑﾀﾖﾈﾟ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾀﾋﾃﾎﾐﾈﾒﾌﾎﾌ ﾊﾀﾒﾕﾈﾋﾋﾀ-ﾌﾀﾊﾊﾈ ﾂ ﾍﾅﾑﾂﾟﾇﾍﾎﾉ ﾎﾁﾋﾀﾑﾒﾈ
	baseIndexFromRotateIndex = (int *)malloc(hMesh->meshElementsNumber * sizeof(int));
	if(baseIndexFromRotateIndex == NULL) exit(0);
	for(i=0; i<hMesh->meshElementsNumber; i++) baseIndexFromRotateIndex[i] = i;
	csrGraphRenameIndexCuthillMcKeeMultiCluster(nePROC, dualX, dualA, baseIndexFromRotateIndex, 1);
	MPI_Barrier(comm);
	if(rankMPI == 0) { printf("Set reverse Cuthill-McKee indexing\n"); fflush(stdout); }
	MPI_Barrier(comm);
	// ﾐﾅﾈﾍﾄﾅﾊﾑﾀﾖﾈﾟ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾀﾋﾃﾎﾐﾈﾒﾌﾎﾌ ﾊﾀﾒﾕﾈﾋﾋﾀ-ﾌﾀﾊﾊﾈ ﾂ ﾍﾅﾑﾂﾟﾇﾍﾎﾉ ﾎﾁﾋﾀﾑﾒﾈ

	// ﾌﾀﾒﾐﾈﾖﾀ ﾏﾅﾐﾅﾕﾎﾄﾀ ﾎﾒ ﾈﾑﾕﾎﾄﾍﾎﾉ ﾈﾍﾄﾅﾊﾑﾀﾖﾈﾈ ﾊ ﾎﾏﾒﾈﾌﾈﾇﾈﾐﾎﾂﾀﾍﾍﾎﾉ
	rotateIndexFromBaseIndex = (int *)malloc(hMesh->meshElementsNumber * sizeof(int));
    for(i=0; i<hMesh->meshElementsNumber; i++)
		rotateIndexFromBaseIndex[baseIndexFromRotateIndex[i]] = i;

	free(dualX); free(dualA);

	MPI_Barrier(comm);
	if(rankMPI == 0) { printf("Set indexes down\n"); fflush(stdout); }
	MPI_Barrier(comm);

	// ﾒﾐﾀﾍﾑﾔﾎﾐﾌﾀﾖﾈﾟ ﾑﾅﾒﾎﾗﾍﾎﾉ ﾒﾎﾏﾎﾋﾎﾃﾈﾈ ﾂ ﾑﾎﾎﾒﾂﾅﾒﾑﾒﾂﾈﾈ ﾑ ﾍﾎﾂﾎﾉ ﾈﾍﾄﾅﾊﾑﾀﾖﾈﾅﾉ ﾟﾗﾅﾅﾊ
	meshEX = (int *)malloc((hMesh->meshElementsNumber + 1) * sizeof(int)); if(meshEX == NULL) exit(0);
	meshEA = (int *)malloc(hMesh->meshEX[hMesh->meshElementsNumber] * sizeof(int)); if(meshEA == NULL) exit(0);

	meshEX[0] = 0;
	for(i=0; i<hMesh->meshElementsNumber; i++)
	{
		int idBase = baseIndexFromRotateIndex[i];
		int nVtx = hMesh->meshEX[idBase + 1] - hMesh->meshEX[idBase];
		int offsetBase = hMesh->meshEX[idBase];

		int j;

		meshEX[i + 1] = meshEX[i];
		for(j=0; j<nVtx; j++)
		{
			meshEA[meshEX[i + 1]] = hMesh->meshEA[offsetBase + j];
			meshEX[i + 1] += 1;
		} // for j
	} // for i

   if(meshEX[hMesh->meshElementsNumber] != hMesh->meshEX[hMesh->meshElementsNumber]) 
   { fprintf(stderr, "ROTATE EX ERROR\n"); exit(0); }

   memcpy(hMesh->meshEX, meshEX, (hMesh->meshElementsNumber + 1)   * sizeof(int));
   memcpy(hMesh->meshEA, meshEA, meshEX[hMesh->meshElementsNumber] * sizeof(int));
   free(meshEX); free(meshEA);

   meshBoundaryFacesTopology = (int *)malloc(hMesh->meshBoundaryFacesNumber * 3 * sizeof(int));
   if(meshBoundaryFacesTopology == NULL) exit(0);
   memcpy((void *)meshBoundaryFacesTopology, (void *)hMesh->meshBoundaryFacesTopology, hMesh->meshBoundaryFacesNumber * 3 * sizeof(int));
   for(i=0; i<hMesh->meshBoundaryFacesNumber; i++) 
	   meshBoundaryFacesTopology[i*3] = rotateIndexFromBaseIndex[meshBoundaryFacesTopology[i*3]];

   memcpy(hMesh->meshBoundaryFacesTopology, meshBoundaryFacesTopology, hMesh->meshBoundaryFacesNumber * 3 * sizeof(int));
   free(meshBoundaryFacesTopology);

   MPI_Barrier(comm);
   if(rankMPI == 0) { printf("Rotate topology\n"); fflush(stdout); }
   MPI_Barrier(comm);
   // ﾒﾐﾀﾍﾑﾔﾎﾐﾌﾀﾖﾈﾟ ﾑﾅﾒﾎﾗﾍﾎﾉ ﾒﾎﾏﾎﾋﾎﾃﾈﾈ ﾂ ﾑﾎﾎﾒﾂﾅﾒﾑﾒﾂﾈﾈ ﾑ ﾍﾎﾂﾎﾉ ﾈﾍﾄﾅﾊﾑﾀﾖﾈﾅﾉ ﾟﾗﾅﾅﾊ   

   partitionUP = (int *)malloc(hMesh->meshElementsNumber * sizeof(int)); if(partitionUP == NULL) exit(0);
   for(i=0; i<hMesh->meshElementsNumber; i++) partitionUP[i] = -1;
   for(i=0; i<hMesh->meshElementsNumber; i++) partitionUP[i] = hMesh->partitionUP[baseIndexFromRotateIndex[i]];
   for(i=0; i<hMesh->meshElementsNumber; i++) if(partitionUP[i] == -1) exit(0);

   memcpy(hMesh->partitionUP, partitionUP, hMesh->meshElementsNumber * sizeof(int));
   free(partitionUP);

   MPI_Barrier(comm);
   if(rankMPI == 0) { printf("Set partition up\n"); fflush(stdout); }
   MPI_Barrier(comm);

   indexesUP = (int *)malloc(hMesh->meshElementsNumber * sizeof(int)); if(indexesUP == NULL) exit(0);
   for(i=0; i<hMesh->meshElementsNumber; i++) indexesUP[i] = -1;
   for(i=0; i<hMesh->meshElementsNumber; i++) indexesUP[i] = hMesh->indexesUP[baseIndexFromRotateIndex[i]];
   for(i=0; i<hMesh->meshElementsNumber; i++) if(indexesUP[i] == -1) exit(0);

   memcpy(hMesh->indexesUP, indexesUP, hMesh->meshElementsNumber * sizeof(int));
   free(indexesUP);

   MPI_Barrier(comm);
   if(rankMPI == 0) { printf("Set indexes up\n"); fflush(stdout); }
   MPI_Barrier(comm);

   free(rotateIndexFromBaseIndex); free(baseIndexFromRotateIndex);
   MPI_Barrier(comm);

   return 0;
}

int InterconnectInitializationHostMPI(gtypeHm4MeshTopology *hMesh, MpiHostInterconnectType *hostInterconnect, MPI_Comm comm)
{
	int i, j;

	MPI_Comm_size(comm, &hostInterconnect->hostsNumber);
	MPI_Comm_rank(comm, &hostInterconnect->hostID);

	// ﾈﾍﾈﾖﾈﾀﾋﾈﾇﾀﾖﾈﾟ ﾑﾒﾐﾓﾊﾒﾓﾐ ﾏﾐﾈﾅﾌﾀ ﾄﾀﾍﾍﾛﾕ ﾂ ﾃﾋﾎﾁﾀﾋﾜﾍﾎﾉ ﾈﾍﾄﾅﾊﾑﾀﾖﾈﾈ
	// ﾏﾐﾅﾄﾏﾎﾋﾀﾃﾀﾅﾒﾑﾟ ﾏﾐﾅﾄﾂﾀﾐﾈﾒﾅﾋﾜﾍﾎﾅ ﾓﾏﾎﾐﾟﾄﾎﾗﾈﾂﾀﾍﾈﾅ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ ﾏﾐﾈ ﾂﾛﾄﾅﾋﾅﾍﾈﾈ HOST-ﾄﾎﾌﾅﾍﾀ
	// ﾇﾀﾏﾐﾀﾘﾈﾂﾀﾅﾌﾛﾅ ﾎﾒ ﾄﾐﾓﾃﾈﾕ HOST ﾇﾍﾀﾗﾅﾍﾈﾟ ﾄﾎﾋﾆﾍﾛ ﾁﾛﾒﾜ ﾑﾍﾅﾑﾅﾍﾛ ﾂ ﾊﾎﾍﾅﾖ ﾋﾎﾊﾀﾋﾜﾍﾎﾉ ﾈﾍﾄﾅﾊﾑﾀﾖﾈﾈ
	// ﾑ ﾓﾏﾎﾐﾟﾄﾎﾗﾈﾂﾀﾍﾈﾅﾌ ﾏﾎ ﾂﾎﾇﾐﾀﾑﾒﾀﾍﾈﾞ ﾈﾄﾅﾍﾒﾈﾔﾈﾊﾀﾒﾎﾐﾎﾂ HOST'ﾎﾂ, ﾀ ﾂﾍﾓﾒﾐﾈ ﾋﾎﾊﾀﾋﾜﾍﾛﾕ ﾑﾏﾈﾑﾊﾎﾂ HOST'ﾎﾂ
	// ﾏﾎ ﾂﾎﾇﾐﾀﾑﾒﾀﾍﾈﾞ ﾃﾋﾎﾁﾀﾋﾜﾍﾛﾕ ﾈﾍﾄﾅﾊﾑﾎﾂ ﾝﾋﾅﾌﾅﾍﾒﾎﾂ
	hostInterconnect->recvBufferTotalLength = 0;
	for(i=0; i<hMesh->meshElementsNumber; i++) if(hMesh->partitionUP[i] != hostInterconnect->hostID) break;
	for(   ; i<hMesh->meshElementsNumber; i++)
	if(hMesh->partitionUP[i] != hostInterconnect->hostID) hostInterconnect->recvBufferTotalLength += 1;
													 else break;
	if(i != hMesh->meshElementsNumber)
	{ fprintf(stderr, "[%4d] : ERROR (InterconnectInitializationHostMPI) --> Elements not ordered 1\n", hostInterconnect->hostID); exit(0); }
	for(i=hMesh->meshElementsNumber - hostInterconnect->recvBufferTotalLength + 1; i<hMesh->meshElementsNumber; i++)
	if(hMesh->partitionUP[i] < hMesh->partitionUP[i-1])
	{ fprintf(stderr, "[%4d] : ERROR (InterconnectInitializationHostMPI) --> Elements not ordered 2\n", hostInterconnect->hostID); exit(0); }

	hostInterconnect->RecvX = (int *)malloc((hostInterconnect->hostsNumber + 1) * sizeof(int)); if(hostInterconnect->RecvX == NULL) exit(0);
	for(i=0; i<=hostInterconnect->hostsNumber; i++) hostInterconnect->RecvX[i] = 0;
	for(i=hMesh->meshElementsNumber - hostInterconnect->recvBufferTotalLength; i<hMesh->meshElementsNumber; i++)
	hostInterconnect->RecvX[hMesh->partitionUP[i] + 1] += 1;
	for(i=2; i<=hostInterconnect->hostsNumber; i++) hostInterconnect->RecvX[i] += hostInterconnect->RecvX[i-1];
	if(hostInterconnect->RecvX[hostInterconnect->hostsNumber] != hostInterconnect->recvBufferTotalLength)
	{ fprintf(stderr, "[%4d] : ERROR (InterconnectInitializationHostMPI) --> RecvX 1\n", hostInterconnect->hostID); exit(0); }
	if(hostInterconnect->RecvX[hostInterconnect->hostID + 1] - hostInterconnect->RecvX[hostInterconnect->hostID] != 0)
	{ fprintf(stderr, "[%4d] : ERROR (InterconnectInitializationHostMPI) --> RecvX 2\n", hostInterconnect->hostID); exit(0); }

	if(hostInterconnect->recvBufferTotalLength > 0)
	{ hostInterconnect->RecvAGlobalIndex = (int *)malloc(hostInterconnect->recvBufferTotalLength * sizeof(int)); if(hostInterconnect->RecvAGlobalIndex == NULL) exit(0); }
	else hostInterconnect->RecvAGlobalIndex = NULL;
	for(j=0, i=hMesh->meshElementsNumber - hostInterconnect->recvBufferTotalLength; i<hMesh->meshElementsNumber; i++)
	{ hostInterconnect->RecvAGlobalIndex[j] = hMesh->indexesUP[i]; j += 1; }

	// ﾈﾍﾈﾖﾈﾀﾋﾈﾇﾀﾖﾈﾟ ﾑﾒﾐﾓﾊﾒﾓﾐ ﾏﾅﾐﾅﾄﾀﾗﾈ ﾄﾀﾍﾍﾛﾕ ﾂ ﾃﾋﾎﾁﾀﾋﾜﾍﾎﾉ ﾈﾍﾄﾅﾊﾑﾀﾖﾈﾈ
	{
		int sendN = 50000;

		hostInterconnect->SendX = (int *)malloc((hostInterconnect->hostsNumber + 1) * sizeof(int)); if(hostInterconnect->SendX == NULL) exit(0);
		hostInterconnect->SendX[0] = 0; for(i=1; i<=hostInterconnect->hostsNumber; i++) hostInterconnect->SendX[i] = -1;
		hostInterconnect->SendAGlobalIndex = (int *)malloc(sendN * sizeof(int)); if(hostInterconnect->SendAGlobalIndex == NULL) exit(0);

		for(i=0; i<hostInterconnect->hostsNumber; i++)
		{
			if(hostInterconnect->hostID == i)
			{
				hostInterconnect->SendX[i + 1] = hostInterconnect->SendX[i];
			
				for(j=hostInterconnect->hostID + 1; j<hostInterconnect->hostsNumber; j++)
				{
					int nSEND;
					int nRECV = hostInterconnect->RecvX[j + 1] - hostInterconnect->RecvX[j];
					MPI_Status status;
				
					MPI_Send(&nRECV, 1, MPI_INT, j, 1184, comm);
					if(nRECV > 0) MPI_Send(hostInterconnect->RecvAGlobalIndex + hostInterconnect->RecvX[j], nRECV, MPI_INT, j, 1283, comm);

					MPI_Recv(&nSEND, 1, MPI_INT, j, 1382, comm, &status);
					if(nSEND > 0)
					{
						hostInterconnect->SendX[j + 1] = hostInterconnect->SendX[j] + nSEND;
						if(hostInterconnect->SendX[j + 1] > sendN)
						{
							int *irPTR = (int *)realloc(hostInterconnect->SendAGlobalIndex, (hostInterconnect->SendX[j + 1] + 1000) * sizeof(int)); if(irPTR == NULL) exit(0);
							if(irPTR != hostInterconnect->SendAGlobalIndex) hostInterconnect->SendAGlobalIndex = irPTR;
							sendN = hostInterconnect->SendX[j + 1] + 1000;
						} // realloc
						MPI_Recv(hostInterconnect->SendAGlobalIndex + hostInterconnect->SendX[j], nSEND, MPI_INT, j, 1481, comm, &status);
					}
					else hostInterconnect->SendX[j + 1] = hostInterconnect->SendX[j];				
				} // for j			
			} // rank == i

			if(hostInterconnect->hostID > i)
			{
				int nSEND;
				int nRECV = hostInterconnect->RecvX[i+1] - hostInterconnect->RecvX[i];
				MPI_Status status;
			
				MPI_Recv(&nSEND, 1, MPI_INT, i, 1184, comm, &status);
				if(nSEND > 0)
				{
					hostInterconnect->SendX[i + 1] = hostInterconnect->SendX[i] + nSEND;
					if(hostInterconnect->SendX[i + 1] > sendN)
					{
						int *irPTR = (int *)realloc(hostInterconnect->SendAGlobalIndex, (hostInterconnect->SendX[i + 1] + 1000) * sizeof(int)); if(irPTR == NULL) exit(0);
						if(irPTR != hostInterconnect->SendAGlobalIndex) hostInterconnect->SendAGlobalIndex = irPTR; sendN = hostInterconnect->SendX[i + 1] + 1000;
					} // realloc
					MPI_Recv(hostInterconnect->SendAGlobalIndex + hostInterconnect->SendX[i], nSEND, MPI_INT, i, 1283, comm, &status);
				}
				else hostInterconnect->SendX[i + 1] = hostInterconnect->SendX[i];

				MPI_Send(&nRECV, 1, MPI_INT, i, 1382, comm);
				if(nRECV > 0) MPI_Send(hostInterconnect->RecvAGlobalIndex + hostInterconnect->RecvX[i], nRECV, MPI_INT, i, 1481, comm);
			} // rank > i
		} // for i

		if( (hostInterconnect->SendX[hostInterconnect->hostsNumber] > 0) && (hostInterconnect->SendX[hostInterconnect->hostsNumber] > sendN) )
		{
			int *irPTR = (int *)realloc(hostInterconnect->SendAGlobalIndex, hostInterconnect->SendX[hostInterconnect->hostsNumber] * sizeof(int)); if(irPTR == NULL) exit(0);
			if(irPTR != hostInterconnect->SendAGlobalIndex) hostInterconnect->SendAGlobalIndex = irPTR;
		} // realloc

		if(hostInterconnect->SendX[hostInterconnect->hostsNumber] == 0) { free(hostInterconnect->SendAGlobalIndex); hostInterconnect->SendAGlobalIndex = NULL; }

		hostInterconnect->sendBufferTotalLength = hostInterconnect->SendX[hostInterconnect->hostsNumber];
	}

	if(hostInterconnect->sendBufferTotalLength > 0)
	{
		int neProc = hMesh->meshElementsNumber - hostInterconnect->recvBufferTotalLength;

		int *globalVSlocal = (int *)malloc(neProc * 2 * sizeof(int));
		if(globalVSlocal == NULL) { fprintf(stderr, "[%4d] : ERROR (InterconnectInitializationHostMPI) --> globalVSlocal 1\n", hostInterconnect->hostID); exit(0); }
		for(i=0; i<neProc; i++)
		{
			globalVSlocal[i*2 + 0] = hMesh->indexesUP[i];
			globalVSlocal[i*2 + 1] = i;
		} // for i
		wsortIntegerListSizeUnit(globalVSlocal, 0, neProc - 1, 2);
		for(i=1; i<neProc; i++) if(globalVSlocal[(i-1)*2] >= globalVSlocal[i*2]) { fprintf(stderr, "[%4d] : ERROR (InterconnectInitializationHostMPI) --> globalVSlocal 2\n", hostInterconnect->hostID); exit(0); }

		hostInterconnect->SendALocalIndex = (int *)malloc(hostInterconnect->sendBufferTotalLength * sizeof(int));
		if(hostInterconnect->SendALocalIndex == NULL) { fprintf(stderr, "[%4d] : ERROR (InterconnectInitializationHostMPI) --> SendALocalIndex 1\n", hostInterconnect->hostID); exit(0); }

		for(i=0; i<hostInterconnect->sendBufferTotalLength; i++)
		{
			j = DichotomySearchElementInOrderIntegerListSizeUnit(globalVSlocal, neProc, hostInterconnect->SendAGlobalIndex[i], 2);
			if(j == -1) { fprintf(stderr, "[%4d] : ERROR (InterconnectInitializationHostMPI) --> Search Local Index 1\n", hostInterconnect->hostID); exit(0); }

			hostInterconnect->SendALocalIndex[i] = globalVSlocal[j*2 + 1];
		} // for i

		free(globalVSlocal);
	}
	else
		hostInterconnect->SendALocalIndex = NULL;

	// ﾏﾐﾎﾂﾅﾐﾊﾀ ﾊﾎﾐﾐﾅﾊﾒﾍﾎﾑﾒﾈ ﾏﾅﾐﾅﾄﾀﾂﾀﾅﾌﾛﾕ ﾄﾀﾍﾍﾛﾕ
	
	// ﾈﾍﾈﾖﾈﾀﾋﾈﾇﾀﾖﾈﾟ ﾖﾅﾋﾎﾗﾈﾑﾋﾅﾍﾍﾛﾕ ﾑﾒﾐﾓﾊﾒﾓﾐ ﾄﾋﾟ ﾏﾅﾐﾅﾄﾀﾗﾈ ﾄﾀﾍﾍﾛﾕ
	hostInterconnect->sendMessagesNumber = 0;
	for(i=0; i<hostInterconnect->hostsNumber; i++) if(hostInterconnect->SendX[i + 1] - hostInterconnect->SendX[i] > 0) hostInterconnect->sendMessagesNumber += 1;
	if(hostInterconnect->hostsNumber > 0)
	{
		hostInterconnect->sendMessageHostsIDs = (int *)malloc(hostInterconnect->sendMessagesNumber * sizeof(int)); if(hostInterconnect->sendMessageHostsIDs == NULL) exit(0);
		hostInterconnect->sendMessageLengths  = (int *)malloc(hostInterconnect->sendMessagesNumber * sizeof(int)); if(hostInterconnect->sendMessageLengths  == NULL) exit(0);
	}
	else
	{
		hostInterconnect->sendMessageHostsIDs = NULL;
		hostInterconnect->sendMessageLengths  = NULL;
	}
	for(j=0, i=0; i<hostInterconnect->hostsNumber; i++) if(hostInterconnect->SendX[i + 1] - hostInterconnect->SendX[i] > 0) 
	{
		hostInterconnect->sendMessageHostsIDs[j] = i;
		hostInterconnect->sendMessageLengths[j] = hostInterconnect->SendX[i + 1] - hostInterconnect->SendX[i];
		j += 1;
	} // for i

	// ﾈﾍﾈﾖﾈﾀﾋﾈﾇﾀﾖﾈﾟ ﾑﾒﾐﾓﾊﾒﾓﾐ ﾄﾋﾟ ﾏﾅﾐﾅﾄﾀﾗﾈ ﾄﾀﾍﾍﾛﾕ ﾂﾅﾙﾅﾑﾒﾂﾅﾍﾍﾎﾃﾎ ﾒﾈﾏﾀ
	if(hostInterconnect->sendMessagesNumber > 0)
	{
		hostInterconnect->sendBufferPointers = (void **)malloc(hostInterconnect->sendMessagesNumber * sizeof(void *)); if(hostInterconnect->sendBufferPointers == NULL) exit(0);
		hostInterconnect->sendValuesBuffer = (double *)malloc(hostInterconnect->sendBufferTotalLength * 5 * sizeof(double)); if(hostInterconnect->sendValuesBuffer == NULL) exit(0);
		for(j=0, i=0; i<hostInterconnect->sendMessagesNumber; i++)
		{
			hostInterconnect->sendBufferPointers[i] = (void *)(&hostInterconnect->sendValuesBuffer[j * 5]);
			j += hostInterconnect->sendMessageLengths[i];
		} // for i
	}
	else
	{
		hostInterconnect->sendBufferPointers = NULL;
		hostInterconnect->sendValuesBuffer = NULL;
	}

	// ﾈﾍﾈﾖﾈﾀﾋﾈﾇﾀﾖﾈﾟ ﾖﾅﾋﾎﾗﾈﾑﾋﾅﾍﾍﾛﾕ ﾑﾒﾐﾓﾊﾒﾓﾐ ﾄﾋﾟ ﾏﾐﾈﾅﾌﾀ ﾄﾀﾍﾍﾛﾕ
	hostInterconnect->recvMessagesNumber = 0;
	for(i=0; i<hostInterconnect->hostsNumber; i++) if(hostInterconnect->RecvX[i + 1] - hostInterconnect->RecvX[i] > 0) hostInterconnect->recvMessagesNumber += 1;
	if(hostInterconnect->recvMessagesNumber > 0)
	{
		hostInterconnect->recvMessageHostsIDs = (int *)malloc(hostInterconnect->recvMessagesNumber * sizeof(int)); if(hostInterconnect->recvMessageHostsIDs == NULL) exit(0);
		hostInterconnect->recvMessageLengths  = (int *)malloc(hostInterconnect->recvMessagesNumber * sizeof(int)); if(hostInterconnect->recvMessageLengths  == NULL) exit(0);
		for(j=0, i=0; i<hostInterconnect->hostsNumber; i++) if(hostInterconnect->RecvX[i + 1] - hostInterconnect->RecvX[i] > 0) 
		{
			hostInterconnect->recvMessageHostsIDs[j] = i;
			hostInterconnect->recvMessageLengths[j] = hostInterconnect->RecvX[i + 1] - hostInterconnect->RecvX[i];
			j += 1;
		} // for i
	}
	else
	{
		hostInterconnect->recvMessageHostsIDs = NULL;
		hostInterconnect->recvMessageLengths = NULL;
	}

	// ﾈﾍﾈﾖﾈﾀﾋﾈﾇﾀﾖﾈﾟ ﾑﾒﾐﾓﾊﾒﾓﾐ ﾄﾋﾟ ﾏﾐﾈﾅﾌﾀ ﾄﾀﾍﾍﾛﾕ ﾂﾅﾙﾅﾑﾒﾂﾅﾍﾍﾎﾃﾎ ﾒﾈﾏﾀ
	if(hostInterconnect->recvMessagesNumber > 0)
	{
		hostInterconnect->recvBufferPointers = (void **)malloc(hostInterconnect->recvMessagesNumber * sizeof(void *)); if(hostInterconnect->recvBufferPointers == NULL) exit(0);
		hostInterconnect->recvValuesBuffer = (double *)malloc(hostInterconnect->recvBufferTotalLength * 5 * sizeof(double)); if(hostInterconnect->recvValuesBuffer == NULL) exit(0);
		for(j=0, i=0; i<hostInterconnect->recvMessagesNumber; i++)
		{
			hostInterconnect->recvBufferPointers[i] = (void *)(&hostInterconnect->recvValuesBuffer[j * 5]);
			j += hostInterconnect->recvMessageLengths[i];
		} // for i
	}
	else
	{
		hostInterconnect->recvBufferPointers = NULL;
		hostInterconnect->recvValuesBuffer = NULL;
	}

	// ﾈﾍﾈﾖﾈﾀﾋﾈﾇﾀﾖﾈﾟ ﾈﾄﾅﾍﾒﾈﾔﾈﾊﾀﾒﾎﾐﾎﾂ ﾃﾐﾓﾏﾏﾎﾂﾛﾕ ﾎﾁﾌﾅﾍﾎﾂ ﾄﾀﾍﾍﾛﾌﾈ
	InitHostsInterconnectRequestsMPI(hostInterconnect[0], comm);
	MPI_Barrier(MPI_COMM_WORLD);

	// ﾒﾅﾑﾒﾎﾂﾛﾉ ﾎﾁﾌﾅﾍ ﾄﾀﾍﾍﾛﾌﾈ
	if(hostInterconnect->hostsNumber > 1)
	{
		int diffGLOB;
		int diffABS = - 1000;

		for(i=0; i<hostInterconnect->sendBufferTotalLength; i++)
		{
			int baseVAL = hMesh->indexesUP[hostInterconnect->SendALocalIndex[i]];
			int *iPtr = (int *)(hostInterconnect->sendValuesBuffer + i*5);
			iPtr[0] = baseVAL; iPtr[1] = 0;
			iPtr[2] = baseVAL; iPtr[3] = 1;
			iPtr[4] = baseVAL; iPtr[5] = 2;
			iPtr[6] = baseVAL; iPtr[7] = 3;
			iPtr[8] = baseVAL; iPtr[9] = 4;			
		} // for i

		for(i=0; i<hostInterconnect->recvBufferTotalLength * 5; i++) 	hostInterconnect->recvValuesBuffer[i] = -1.235e-8;

		HostsInterconnectStartWaitAll(hostInterconnect[0], comm);
		MPI_Barrier(comm);
		
		for(i=0; i<hostInterconnect->recvBufferTotalLength; i++)
		{
			int baseVAL = hostInterconnect->RecvAGlobalIndex[i];
			int  QEX[10] = {baseVAL, 0, baseVAL, 1, baseVAL, 2, baseVAL, 3, baseVAL, 4};
			int *iPtr = (int *)(hostInterconnect->recvValuesBuffer + i*5);
			for(j=0; j<10; j++)
			{
				int diff = iPtr[j] - QEX[j];
				if(diff < 0) diff = - diff;

				if(diffABS < diff) diffABS = diff;
			} // for j
		} // for i

		MPI_Barrier(comm);
		MPI_Allreduce(&diffABS, &diffGLOB, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		if(hostInterconnect->hostID == 0)
		{ printf("MPI HOST DEBUG TRINSFER ERROR CODE: %d\n", diffGLOB); fflush(stdout); }
		MPI_Barrier(comm);
	}
	// ﾒﾅﾑﾒﾎﾂﾛﾉ ﾎﾁﾌﾅﾍ ﾄﾀﾍﾍﾛﾌﾈ

	MPI_Barrier(comm);
			
	return 0;
}
