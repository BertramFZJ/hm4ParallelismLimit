#include <stdio.h>
#include <stdlib.h>

#include <metis.h>

#include "csrProcessCuthillMcKee.h"

#define _DEBUG_PRINTF_ 1
int csrGraphRenameIndexCuthillMcKeeMultiLevel(int n, int subdomenSize, int *x, int *a, int *index, int reverseFlag)
{
	int *part     = (int *)malloc(n * sizeof(int));
	int *localID  = (int *)malloc(n * sizeof(int));
	int *globalID = (int *)malloc(n * sizeof(int));
	
	int npart = 0;
	int nset = 0;

	int iPART;

	if(part     == NULL) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiLevel: ERROR #1\n"); exit(0); }
	if(localID  == NULL) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiLevel: ERROR #2\n"); exit(0); }
	if(globalID == NULL) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiLevel: ERROR #3\n"); exit(0); }

	npart = n / subdomenSize; if(npart < 2) npart = 2;
	if(_DEBUG_PRINTF_) printf("SUBDOMENS NUMBER: %d\n", npart);

	// РАЗБИЕНИЕ НА ДОМЕНЫ
	if(npart > 1)
	{
		int wgtflag = 0; // флаг наличия весов ребер и вершин графа
		int *vwgt = NULL, *adjwgt = NULL; // массивы весов вершин графа и ребер соответственно
		int numflag = 0; // флаг стиля нумерации узлов				
		int options[5] = {0, 0, 0, 0, 0}; // параметры алгоритма разбиения
		int edgecut = 0; // число разрезанных в результате разбиения ребер сетки		
		int Nedge = x[n] / 2;
		int nPart = npart;

		if(_DEBUG_PRINTF_) { printf(" Start METIS_PartGraphKway .....\n"); fflush(stdout); }
		METIS_PartGraphKway(&n, x, a, vwgt, adjwgt, &wgtflag, &numflag, &npart, options, &edgecut, part); 
		if(_DEBUG_PRINTF_) printf("Finish METIS_PartGraphKway      \n");
		if(_DEBUG_PRINTF_) printf("edgecut = %d (%.2lf %%)\n", edgecut, 100.0 * (double)edgecut/(double)Nedge);		
		if(_DEBUG_PRINTF_) { printf("\n"); fflush(stdout); }
	}
	// РАЗБИЕНИЕ НА ДОМЕНЫ

	for(iPART=0; iPART<npart; iPART++)
	{
		int i, j;

		int *xx, *aa, *localIndex;
		int localN = 0;
		int firstID;
			
		for(i=0; i<n; i++) localID[i] = globalID[i] = -1;
		for(j=0, i=0; i<n; i++) if(part[i] == iPART)
		{
			 localID[i]      = localN; 
			globalID[localN] = i;
			localN += 1; 

			j += x[i+1] - x[i];
		} // if in for

		if(_DEBUG_PRINTF_) printf("SUB %4d: NE = %8d\n", iPART, localN);

		        xx = (int *)malloc((localN + 1) * sizeof(int)); if(xx         == NULL) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiLevel: ERROR #4.1\n"); exit(0); }
		        aa = (int *)malloc(j            * sizeof(int)); if(aa         == NULL) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiLevel: ERROR #4.2\n"); exit(0); }
		localIndex = (int *)malloc(localN       * sizeof(int)); if(localIndex == NULL) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiLevel: ERROR #4.3\n"); exit(0); }

		xx[0] = 0;
		for(i=0; i<n; i++) if(part[i] == iPART)
		{
			xx[localID[i] + 1] = xx[localID[i]];

			for(j=x[i]; j<x[i+1]; j++) if(part[a[j]] == iPART)
			{ 
				aa[xx[localID[i] + 1]] = localID[a[j]]; 
				xx[localID[i] + 1] += 1; 
			}

			if(xx[localID[i] + 1] == xx[localID[i]]) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiLevel: ERROR #5 (NO CONVEX DOMEN)\n"); exit(0); }
		} // if in for

		firstID = PseudoPeripheralCsrGraphVertex(localN, xx, aa);
		if(_DEBUG_PRINTF_) printf("FIRST ID = %8d LOCAL AND %8d GLOBAL\n", firstID, globalID[firstID]);

		csrGraphRenameIndexCuthillMcKee(localN, firstID, xx, aa, localIndex);

		if(reverseFlag == 1)
		{
			int *indexR = (int *)malloc(localN * sizeof(int));
			if(indexR == NULL) exit(0);

			for(i=0; i<localN; i++)
				indexR[i] = localIndex[localN - 1 - i];
			
			for(i=0; i<localN; i++) localIndex[i] = indexR[i];

			free(indexR);
		} // if reverseFlag

		for(i=0; i<localN; i++)
		{
			index[nset] = globalID[localIndex[i]];
			nset += 1;
		} // for i

		free(xx); free(aa); free(localIndex);
	} // for iPART

	if(nset != n) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiLevel: ERROR #6\n"); exit(0); }

	free(part); free(localID); free(globalID);

	return 0;
}
#undef _DEBUG_PRINTF_

#define _DEBUG_PRINTF_ 0
int csrGraphRenameIndexCuthillMcKeeMultiCluster(int n, int *x, int *a, int *index, int reverseFlag)
{
	int *cluster  = (int *)malloc(n * sizeof(int));
	int *localID  = (int *)malloc(n * sizeof(int));
	int *globalID = (int *)malloc(n * sizeof(int));
	
	int ncluster = 0;
	int nset = 0;

	int iCLUSTER;

	if(cluster  == NULL) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiCluster: ERROR #1\n"); exit(0); }
	if(localID  == NULL) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiCluster: ERROR #2\n"); exit(0); }
	if(globalID == NULL) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiCluster: ERROR #3\n"); exit(0); }

	// ВЫДЕЛЕНИЕ КЛАСТЕРОВ
	{
		int *queue = localID;
		int i, j;
			
		for(i=0; i<n; i++) cluster[i] = -1;

		while(nset != n)
		{
			int *queue = localID;
			int lqueue = 0;

			for(i=0; i<n; i++) if(cluster[i] == -1) break;			
			if(i == n) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiCluster: ERROR #M\n"); exit(0); }

			queue[0] = i; cluster[i] = ncluster; lqueue = 1; nset += 1;

			for(i=0; i<lqueue; i++)
			{
				for(j=x[queue[i]]; j<x[queue[i] + 1]; j++)
					if(cluster[a[j]] == -1)
					{
						cluster[a[j]] = ncluster; nset += 1;
						queue[lqueue] = a[j]; lqueue += 1;
					} // if
			} // for i

			ncluster += 1;
		} // while

		for(i=0; i<n; i++) queue[i] = -1;
		for(i=0; i<n; i++) if(cluster[i] == -1) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiCluster: ERROR #N\n"); exit(0); }

		nset = 0;
	}
	// ВЫДЕЛЕНИЕ КЛАСТЕРОВ

	if(_DEBUG_PRINTF_) printf("CLUSTER NUMBER: %d\n", ncluster);	

	for(iCLUSTER=0; iCLUSTER<ncluster; iCLUSTER++)
	{
		int i, j;

		int *xx, *aa, *localIndex;
		int localN = 0;
		int firstID;
			
		for(i=0; i<n; i++) localID[i] = globalID[i] = -1;
		for(j=0, i=0; i<n; i++) if(cluster[i] == iCLUSTER)
		{
			 localID[i]      = localN; 
			globalID[localN] = i;
			localN += 1; 

			j += x[i+1] - x[i];
		} // if in for

		if(_DEBUG_PRINTF_) printf("SUB %4d: NE = %8d\n", iCLUSTER, localN);

		if(localN > 1)
		{
			        xx = (int *)malloc((localN + 1) * sizeof(int)); if(xx         == NULL) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiCluster: ERROR #4.1\n"); exit(0); }
		    	    aa = (int *)malloc(j            * sizeof(int)); if(aa         == NULL) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiCluster: ERROR #4.2\n"); exit(0); }
			localIndex = (int *)malloc(localN       * sizeof(int)); if(localIndex == NULL) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiCluster: ERROR #4.3\n"); exit(0); }

			xx[0] = 0;
			for(i=0; i<n; i++) if(cluster[i] == iCLUSTER)
			{
				xx[localID[i] + 1] = xx[localID[i]];

				for(j=x[i]; j<x[i+1]; j++) if(cluster[a[j]] == iCLUSTER)
				{ 
					aa[xx[localID[i] + 1]] = localID[a[j]]; 
					xx[localID[i] + 1] += 1; 
				}

				if(xx[localID[i] + 1] == xx[localID[i]]) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiCluster: ERROR #5 (NO CONVEX DOMEN)\n"); exit(0); }
			} // if in for

			firstID = PseudoPeripheralCsrGraphVertex(localN, xx, aa);
			if(_DEBUG_PRINTF_) printf("FIRST ID = %8d LOCAL AND %8d GLOBAL\n", firstID, globalID[firstID]);

			csrGraphRenameIndexCuthillMcKee(localN, firstID, xx, aa, localIndex);

			if(reverseFlag == 1)
			{
				int *indexR = (int *)malloc(localN * sizeof(int));
				if(indexR == NULL) exit(0);

				for(i=0; i<localN; i++)
					indexR[i] = localIndex[localN - 1 - i];
			
				for(i=0; i<localN; i++) localIndex[i] = indexR[i];

				free(indexR);
			} // if reverseFlag

			for(i=0; i<localN; i++)
			{
				index[nset] = globalID[localIndex[i]];
				nset += 1;
			} // for i

			free(xx); free(aa); free(localIndex);
		} // if localN > 1
		else
		{
			if(_DEBUG_PRINTF_) printf("SINGLE CELL SUBDOMAIN\n");
			index[nset] = globalID[0];
			nset += 1;
		}
	} // for iPART

	if(nset != n) { fprintf(stderr, "csrGraphRenameIndexCuthillMcKeeMultiCluster: ERROR #6\n"); exit(0); }

	free(cluster); free(localID); free(globalID);

	return 0;
}
#undef _DEBUG_PRINTF_


/*
if(1)
	{
		int *newEID = (int *)malloc(Ne * sizeof(int));

		int i, j, k;

		if(newEID == NULL) exit(0);

		for(k=0; k<_CLUSTER_NUMBER_; k++)
		{
			int *xx, *aa, *pp;
			int Nec = 0;
			
			for(i=0; i<Ne; i++) newEID[i] = -1;

			j = 0;
			for(i=0; i<Ne; i++) if( (part[i] == k) && (eClusterType[i] == -1) )
			{
				newEID[i] = Nec; Nec += 1; 
				j += xadjE[i+1] - xadjE[i];
			} // if in for

			printf("k = %d Nec = %d L = %d\n", k, Nec, j);

			xx = (int *)malloc((Nec + 1) * sizeof(int)); if(xx == NULL) exit(0);
			aa = (int *)malloc(j * sizeof(int)); if(aa == NULL) exit(0);
			pp = (int *)malloc(j * sizeof(int)); if(pp == NULL) exit(0);

			xx[0] = 0;
			for(i=0; i<Ne; i++) if( (part[i] == k) && (eClusterType[i] == -1) )
			{
				xx[newEID[i]+1] = xx[newEID[i]];
				for(j=xadjE[i]; j<xadjE[i+1]; j++) if( (eClusterType[adjncyE[j]] == -1) && (part[adjncyE[j]] == k) )
				{ aa[xx[newEID[i]+1]] = newEID[adjncyE[j]]; xx[newEID[i]+1] += 1; }

				if(xx[newEID[i]+1] == xx[newEID[i]]) { fprintf(stderr, "ERROR 111\n"); exit(0); }
			} // if in for

			printf("k = %d Nec = %d L = %d\n", k, Nec, xx[Nec]);

			{
				int wgtflag = 0; // флаг наличия весов ребер и вершин графа
				int *vwgt = NULL, *adjwgt = NULL; // массивы весов вершин графа и ребер соответственно
				int numflag = 0; // флаг стиля нумерации узлов				
				int options[5] = {0, 0, 0, 0, 0}; // параметры алгоритма разбиения
				int edgecut = 0; // число разрезанных в результате разбиения ребер сетки		
				int Nedge;
				int nnn = _GPU_NUMBER_;

				printf(" Start METIS_PartGraphKway .....\n");
				METIS_PartGraphKway(&Nec, xx, aa, vwgt, adjwgt, &wgtflag, &numflag, &nnn, options, &edgecut, pp); 
				Nedge = xadjE[Nec] / 2;
				printf("Finish METIS_PartGraphKway      \n");
				printf("edgecut = %d (%.2lf %%)\n", edgecut, 100.0 * (double)edgecut/(double)Nedge);		
				printf("\n");
			}

			for(i=0; i<Ne; i++) if(newEID[i] >= 0) partM[i] = nparts + k*_GPU_NUMBER_ + pp[newEID[i]];
			// for(i=0; i<Ne; i++) if(newEID[i] >= 0) partM[i] = k*_GPU_NUMBER_ + pp[newEID[i]];

			free(xx); free(aa); free(pp);
		} // for k

		free(newEID);

		for(i=0; i<Ne; i++) if(eClusterType[i] != -1) partM[i] = part[i];

		if(1) VisualHybMeshPart2D(Np, CV, Ne, XE, AE, nparts*(_GPU_NUMBER_ + 1), partM, "./output/MeshPartGPU42L2.plt");
		// if(1) VisualHybMeshPart2D(Np, CV, Ne, XE, AE, nparts*_GPU_NUMBER_, partM, "./output/MeshPartCPU42.plt");
		if(0) VisualCSRGraphPart2D(Ne, CEL, xadjE, adjncyE, nparts*(_GPU_NUMBER_ + 1), partM, "./output/GraphPartGPU.plt");
	}
	*/
