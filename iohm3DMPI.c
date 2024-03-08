#define _DEBUG_OUTPUT_ 0

#define _CRT_SECURE_NO_WARNINGS
#define HMSH_STR_LENGTH 81

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "iohm3DMPI.h"

int ReadHybMeshNodesBinaryMPI(char *fname, int *NDIST, double **C, FILE *outSTREAM, MPI_Comm comm)
{
	int rank, size;

	  MPI_Barrier(comm);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	  MPI_Barrier(comm);

	if( (outSTREAM != NULL) && (rank == 0) )
	{ fprintf(outSTREAM, "ReadHybMeshNodesBinaryMPI --> MAIN: START PROCEDURE IN %d MPI PROCESSES TOPOLOGY\n", size); }

	if(rank == 0)
	{
		FILE *file = NULL;
		int offsetREAD;
		int get;
		int i, j, k;

		file = fopen(fname, "rb+");
		if(file == NULL) { fprintf(stderr, "ReadHybMeshNodesBinaryMPI --> MAIN: *** ERROR *** Can't open data file \"%s\"\n", fname); MPI_Abort(comm, 0); }
		else if(outSTREAM != NULL) { fprintf(outSTREAM, "ReadHybMeshNodesBinaryMPI --> MAIN: Open data file \"%s\"\n", fname); }

		// Õ¿«¬¿Õ»≈ —≈“ »
		{
			char meshNAME[HMSH_STR_LENGTH];
			for(i=0; i<HMSH_STR_LENGTH; i++) meshNAME[i] = '\0';
			get = (int)fread(meshNAME, sizeof(char), HMSH_STR_LENGTH, file);
			if(get != HMSH_STR_LENGTH) { fprintf(stderr, "ReadHybMeshNodesBinaryMPI --> MAIN: *** ERROR *** Can't read data from file \"%s\"\n", fname); MPI_Abort(comm, 0); }
			if(outSTREAM != NULL) { fprintf(outSTREAM, "ReadHybMeshNodesBinaryMPI --> MAIN: MESH NAME = \"%s\"\n", meshNAME); }
		}
		// Õ¿«¬¿Õ»≈ —≈“ »

		// —◊»“€¬¿Õ»≈ ’¿–¿ “≈–»—“»◊≈— »’ œ¿–¿Ã≈“–Œ¬ » »Õ»÷»¿À»«¿÷»ﬂ –¿—œ–≈ƒ≈À≈Õ»ﬂ ”«ÀŒ¬
		{
			int SP[10];
			
			get = (int)fread(SP, sizeof(int), 10, file);
			if(get != 10) { fprintf(stderr, "ReadHybMeshNodesBinaryMPI --> MAIN: *** ERROR *** Can't read data from file \"%s\"\n", fname); MPI_Abort(comm, 0); }

			j = SP[0] / size; k = SP[0] % size;
			NDIST[0] = 0;
			for(i=0; i<size; i++)
			{
				NDIST[i + 1] = NDIST[i] + j;
				if(i < k) NDIST[i + 1] += 1;
			} // for i

			if(size > 1) MPI_Bcast(NDIST, size + 1, MPI_INT, 0, comm);			

			offsetREAD = HMSH_STR_LENGTH*sizeof(char) + 10*sizeof(int) + SP[6]*HMSH_STR_LENGTH*sizeof(char);
		}
		// —◊»“€¬¿Õ»≈ ’¿–¿ “≈–»—“»◊≈— »’ œ¿–¿Ã≈“–Œ¬ » »Õ»÷»¿À»«¿÷»ﬂ –¿—œ–≈ƒ≈À≈Õ»ﬂ ”«ÀŒ¬

		// —◊»“€¬¿Õ»≈  ŒŒ–ƒ»Õ¿“ —Œ¡—“¬≈ÕÕ€’ ”«ÀŒ¬
		{
			C[0] = (double *)malloc(NDIST[1] * 3 * sizeof(double));
			if(C[0] == NULL) { fprintf(stderr, "ReadHybMeshNodesBinaryMPI --> MAIN: *** ERROR *** Can't allocate memory <%d * DOUBLE>\n", 3 * NDIST[1]); MPI_Abort(comm, 0); }
			fseek(file, offsetREAD, SEEK_SET);
			get = (int)fread(C[0], sizeof(double), NDIST[1] * 3, file);
			if(get != NDIST[1] * 3) { fprintf(stderr, "ReadHybMeshNodesBinaryMPI --> MAIN: *** ERROR *** Can't read data from file \"%s\"\n", fname); MPI_Abort(comm, 0); }			

#if(_DEBUG_OUTPUT_ == 1)
			if(outSTREAM != NULL) { fprintf(outSTREAM, "ReadHybMeshNodesBinaryMPI --> MAIN: READ COORDINATES DATA NV = %8d\n", NDIST[1]); }
#endif
		}
		// —◊»“€¬¿Õ»≈  ŒŒ–ƒ»Õ¿“ —Œ¡—“¬≈ÕÕ€’ ”«ÀŒ¬

		// –¿——€À ¿  ŒŒ–ƒ»Õ¿“ ”«ÀŒ¬ Œ—“¿À‹Õ€Ã œ–Œ÷≈——¿Ã
		if(size > 1)
		{
			double *cBUF = (double *)malloc(NDIST[1] * 3 * sizeof(double));
			if(cBUF == NULL) { fprintf(stderr, "ReadHybMeshNodesBinaryMPI --> MAIN: *** ERROR *** Can't allocate memory <%d * DOUBLE>\n", 3 * NDIST[1]); MPI_Abort(comm, 0); }
			

			for(i=1; i<size; i++)
			{
				int n = NDIST[i + 1] - NDIST[i];
				
				fseek(file, offsetREAD + NDIST[i] * 3 * sizeof(double), SEEK_SET);
				get = (int)fread(cBUF, sizeof(double), n * 3, file);
				if(get != n * 3) { fprintf(stderr, "ReadHybMeshNodesBinaryMPI --> MAIN: *** ERROR *** Can't read data from file \"%s\"\n", fname); MPI_Abort(comm, 0); }
				MPI_Send(cBUF, n * 3, MPI_DOUBLE, i, 1000 + i, comm);
								
#if(_DEBUG_OUTPUT_ == 1)
				if(outSTREAM != NULL) { fprintf(outSTREAM, "ReadHybMeshNodesBinaryMPI --> MAIN: SEND COORDINATES DATA NV = %8d TO %4d\n", n, i); }
#endif
			} // for i

			free(cBUF);
		}
		// –¿——€À ¿  ŒŒ–ƒ»Õ¿“ ”«ÀŒ¬ Œ—“¿À‹Õ€Ã œ–Œ÷≈——¿Ã

		fclose(file);
		if(outSTREAM != NULL) { fprintf(outSTREAM, "ReadHybMeshNodesBinaryMPI --> MAIN: Close data file \"%s\"\n", fname); }
	} // rank == 0
	
	if(rank != 0)
	{
		MPI_Status status;
		
		MPI_Bcast(NDIST, size + 1, MPI_INT, 0, comm);

		C[0] = (double *)malloc((NDIST[rank + 1] - NDIST[rank]) * 3 * sizeof(double));
		if(C[0] == NULL) { fprintf(stderr, "ReadHybMeshNodesBinaryMPI --> %4d: *** ERROR *** Can't allocate memory <%d * DOUBLE>\n", rank, (NDIST[rank + 1] - NDIST[rank]) * 3); MPI_Abort(comm, 0); }
		MPI_Recv(C[0], (NDIST[rank + 1] - NDIST[rank]) * 3, MPI_DOUBLE, 0, 1000 + rank, comm, &status);		
	} // rank != 0

	 MPI_Barrier(comm);

	return 0;
}

int ReadHybMeshElementsBinaryMPI(char *fname, int *EDIST, int **XE, int **AE, FILE *outSTREAM, MPI_Comm comm)
{
	int rank, size;

	  MPI_Barrier(comm);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	  MPI_Barrier(comm);

	if( (outSTREAM != NULL) && (rank == 0) )
	{ fprintf(outSTREAM, "ReadHybMeshElementsBinaryMPI --> MAIN: START PROCEDURE IN %d MPI PROCESSES TOPOLOGY\n", size); }

	if(rank == 0)
	{
		FILE *file = NULL;
		int offsetREAD;
		int get;
		int i, j, k;

		file = fopen(fname, "rb+");
		if(file == NULL) { fprintf(stderr, "ReadHybMeshElementsBinaryMPI --> MAIN: *** ERROR *** Can't open data file \"%s\"\n", fname); MPI_Abort(comm, 0); }
		else if(outSTREAM != NULL) { fprintf(outSTREAM, "ReadHybMeshElementsBinaryMPI --> MAIN: Open data file \"%s\"\n", fname); }

		// Õ¿«¬¿Õ»≈ —≈“ »
		{
			char meshNAME[HMSH_STR_LENGTH];
			for(i=0; i<HMSH_STR_LENGTH; i++) meshNAME[i] = '\0';
			get = (int)fread(meshNAME, sizeof(char), HMSH_STR_LENGTH, file);
			if(get != HMSH_STR_LENGTH) { fprintf(stderr, "ReadHybMeshElementsBinaryMPI --> MAIN: *** ERROR *** Can't read data from file \"%s\"\n", fname); MPI_Abort(comm, 0); }
			if(outSTREAM != NULL) { fprintf(outSTREAM, "ReadHybMeshElementsBinaryMPI --> MAIN: MESH NAME = \"%s\"\n", meshNAME); }
		}
		// Õ¿«¬¿Õ»≈ —≈“ »

		// —◊»“€¬¿Õ»≈ ’¿–¿ “≈–»—“»◊≈— »’ œ¿–¿Ã≈“–Œ¬ » »Õ»÷»¿À»«¿÷»ﬂ –¿—œ–≈ƒ≈À≈Õ»ﬂ ›À≈Ã≈Õ“Œ¬
		{
			int SP[10];
			
			get = (int)fread(SP, sizeof(int), 10, file);
			if(get != 10) { fprintf(stderr, "ReadHybMeshElementsBinaryMPI --> MAIN: *** ERROR *** Can't read data from file \"%s\"\n", fname); MPI_Abort(comm, 0); }

			j = SP[1] / size; k = SP[1] % size;
			EDIST[0] = 0;
			for(i=0; i<size; i++)
			{
				EDIST[i + 1] = EDIST[i] + j;
				if(i < k) EDIST[i + 1] += 1;
			} // for i

			if(size > 1) MPI_Bcast(EDIST, size + 1, MPI_INT, 0, comm);			

			offsetREAD = HMSH_STR_LENGTH*sizeof(char) + 10*sizeof(int) + SP[6]*HMSH_STR_LENGTH*sizeof(char) + SP[0]*3*sizeof(double);
		}
		// —◊»“€¬¿Õ»≈ ’¿–¿ “≈–»—“»◊≈— »’ œ¿–¿Ã≈“–Œ¬ » »Õ»÷»¿À»«¿÷»ﬂ –¿—œ–≈ƒ≈À≈Õ»ﬂ ›À≈Ã≈Õ“Œ¬

		// —◊»“€¬¿Õ»≈ “ŒœŒÀŒ√»» —Œ¡—“¬≈ÕÕ€’ ›À≈Ã≈Õ“Œ¬
		{
			XE[0] = (int *)malloc((EDIST[1] + 1) * sizeof(int));
			if(XE[0] == NULL) { fprintf(stderr, "ReadHybMeshElementsBinaryMPI --> MAIN: *** ERROR *** Can't allocate memory <%d * INT>\n", EDIST[1]); MPI_Abort(comm, 0); }
			fseek(file, offsetREAD, SEEK_SET);
			get = (int)fread(XE[0], sizeof(int), EDIST[1] + 1, file);
			if(get != EDIST[1] + 1) { fprintf(stderr, "ReadHybMeshElementsBinaryMPI --> MAIN: *** ERROR *** Can't read data from file \"%s\"\n", fname); MPI_Abort(comm, 0); }

			AE[0] = (int *)malloc(XE[0][EDIST[1]] * sizeof(int));
			if(AE[0] == NULL) { fprintf(stderr, "ReadHybMeshElementsBinaryMPI --> MAIN: *** ERROR *** Can't allocate memory <%d * INT>\n", XE[0][EDIST[1]]); MPI_Abort(comm, 0); }
			fseek(file, offsetREAD + (EDIST[size] + 1) * sizeof(int), SEEK_SET);
			get = (int)fread(AE[0], sizeof(int), XE[0][EDIST[1]], file);
			if(get != XE[0][EDIST[1]]) { fprintf(stderr, "ReadHybMeshElementsBinaryMPI --> MAIN: *** ERROR *** Can't read data from file \"%s\"\n", fname); MPI_Abort(comm, 0); }

#if(_DEBUG_OUTPUT_ == 1)
			if(outSTREAM != NULL) { fprintf(outSTREAM, "ReadHybMeshElementsBinaryMPI --> MAIN: READ DATA NE = %8d LA = %8d\n", EDIST[1], XE[0][EDIST[1]]); }
#endif
		}
		// —◊»“€¬¿Õ»≈ “ŒœŒÀŒ√»» —Œ¡—“¬≈ÕÕ€’ ›À≈Ã≈Õ“Œ¬

		// –¿——€À ¿ “ŒœŒÀŒ√»» ›À≈Ã≈Õ“Œ¬ Œ—“¿À‹Õ€Ã œ–Œ÷≈——¿Ã
		if(size > 1)
		{
			int *xBUF = (int *)malloc((EDIST[1] + 1) * sizeof(int));
			int *aBUF = (int *)malloc( EDIST[1] * 8  * sizeof(int));

			if(xBUF == NULL) { fprintf(stderr, "ReadHybMeshElementsBinaryMPI --> MAIN: *** ERROR *** Can't allocate memory <%d * INT>\n", EDIST[1] + 1); MPI_Abort(comm, 0); }
			if(aBUF == NULL) { fprintf(stderr, "ReadHybMeshElementsBinaryMPI --> MAIN: *** ERROR *** Can't allocate memory <%d * INT>\n", EDIST[1] * 8); MPI_Abort(comm, 0); }

			for(i=1; i<size; i++)
			{
				int xn = EDIST[i + 1] - EDIST[i];
				int an;

				fseek(file, offsetREAD + EDIST[i]*sizeof(int), SEEK_SET);
				get = (int)fread(xBUF, sizeof(int), xn + 1, file);
				if(get != xn + 1) { fprintf(stderr, "ReadHybMeshElementsBinaryMPI --> MAIN: *** ERROR *** Can't read data from file \"%s\"\n", fname); MPI_Abort(comm, 0); }
				MPI_Send(xBUF, xn + 1, MPI_INT, i, 1000 + i, comm);
								
				an = xBUF[xn] - xBUF[0];
				fseek(file, offsetREAD + (EDIST[size] + 1 + xBUF[0]) * sizeof(int), SEEK_SET);
				get = (int)fread(aBUF, sizeof(int), an, file);
				if(get != an) { fprintf(stderr, "ReadHybMeshElementsBinaryMPI --> MAIN: *** ERROR *** Can't read data from file \"%s\"\n", fname); MPI_Abort(comm, 0); }					
				MPI_Send(aBUF, an, MPI_INT, i, 5000 + i, comm);

#if(_DEBUG_OUTPUT_ == 1)
				if(outSTREAM != NULL) { fprintf(outSTREAM, "ReadHybMeshElementsBinaryMPI --> MAIN: SEND DATA NE = %8d LA = %8d TO %4d\n", xn, an, i); }
#endif
			} // for i

			free(xBUF); free(aBUF);
		}
		// –¿——€À ¿ “ŒœŒÀŒ√»» ›À≈Ã≈Õ“Œ¬ Œ—“¿À‹Õ€Ã œ–Œ÷≈——¿Ã

		fclose(file);
		if(outSTREAM != NULL) { fprintf(outSTREAM, "ReadHybMeshElementsBinaryMPI --> MAIN: Close data file \"%s\"\n", fname); }
	} // rank == 0
	
	if(rank != 0)
	{
		MPI_Status status;
		int xn, an;
		int i;

		MPI_Bcast(EDIST, size + 1, MPI_INT, 0, comm);

		xn = EDIST[rank + 1] - EDIST[rank];
		XE[0] = (int *)malloc((xn + 1) * sizeof(int));
		if(XE[0] == NULL) { fprintf(stderr, "ReadHybMeshElementsBinaryMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, xn + 1); MPI_Abort(comm, 0); }
		MPI_Recv(XE[0], xn + 1, MPI_INT, 0, 1000 + rank, comm, &status);

		an = XE[0][xn] - XE[0][0];
		AE[0] = (int *)malloc(an * sizeof(int));
		if(AE[0] == NULL) { fprintf(stderr, "ReadHybMeshElementsBinaryMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, an); MPI_Abort(comm, 0); }
		MPI_Recv(AE[0], an, MPI_INT, 0, 5000 + rank, comm, &status);

		for(i=xn; i>=0; i--) XE[0][i] -= XE[0][0];
	} // rank != 0

	 MPI_Barrier(comm);

	return 0;
}

int ReadHybMeshBoundaryBinaryMPI(char *fname, int *nFACE, char **faceNAME, int *nEB, int **EB, FILE *outSTREAM, MPI_Comm comm)
{
	int rank, size;

	  MPI_Barrier(comm);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	  MPI_Barrier(comm);

	if( (outSTREAM != NULL) && (rank == 0) )
	{ fprintf(outSTREAM, "ReadHybMeshBoundaryBinaryMPI --> MAIN: START PROCEDURE IN %d MPI PROCESSES TOPOLOGY\n", size); }

	if(rank == 0)
	{
		FILE *file = NULL;
		int get;
		int i;

		file = fopen(fname, "rb+");
		if(file == NULL) { fprintf(stderr, "ReadHybMeshBoundaryBinaryMPI --> MAIN: *** ERROR *** Can't open data file \"%s\"\n", fname); MPI_Abort(comm, 0); }
		else if(outSTREAM != NULL) { fprintf(outSTREAM, "ReadHybMeshBoundaryBinaryMPI --> MAIN: Open data file \"%s\"\n", fname); }

		// Õ¿«¬¿Õ»≈ —≈“ »
		{
			char meshNAME[HMSH_STR_LENGTH];
			for(i=0; i<HMSH_STR_LENGTH; i++) meshNAME[i] = '\0';
			get = (int)fread(meshNAME, sizeof(char), HMSH_STR_LENGTH, file);
			if(get != HMSH_STR_LENGTH) { fprintf(stderr, "ReadHybMeshBoundaryBinaryMPI --> MAIN: *** ERROR *** Can't read data from file \"%s\"\n", fname); MPI_Abort(comm, 0); }
			if(outSTREAM != NULL) { fprintf(outSTREAM, "ReadHybMeshBoundaryBinaryMPI --> MAIN: MESH NAME = \"%s\"\n", meshNAME); }
		}
		// Õ¿«¬¿Õ»≈ —≈“ »

		// —◊»“€¬¿Õ»≈ ’¿–¿ “≈–»—“»◊≈— »’ œ¿–¿Ã≈“–Œ¬
		{
			int SP[10];
			
			get = (int)fread(SP, sizeof(int), 10, file);
			if(get != 10) { fprintf(stderr, "ReadHybMeshBoundaryBinaryMPI --> MAIN: *** ERROR *** Can't read data from file \"%s\"\n", fname); MPI_Abort(comm, 0); }
						
			nFACE[0] = SP[6];			
			  nEB[0] = SP[7];			
		}
		// —◊»“€¬¿Õ»≈ ’¿–¿ “≈–»—“»◊≈— »’ œ¿–¿Ã≈“–Œ¬

		// —◊»“€¬¿Õ»≈ Õ¿«¬¿Õ»… √–¿Õ»◊Õ€’ œŒ¬≈–’ÕŒ—“≈…
		{
			faceNAME[0] = (char *)malloc(HMSH_STR_LENGTH * nFACE[0] * sizeof(char));
			if(faceNAME[0] == NULL) { fprintf(stderr, "ReadHybMeshBoundaryBinaryMPI --> MAIN: *** ERROR *** Can't allocate memory <%d * CHAR>\n", HMSH_STR_LENGTH * nFACE[0]); MPI_Abort(comm, 0); }

			get = (int)fread(faceNAME[0], sizeof(char), HMSH_STR_LENGTH * nFACE[0], file);
			if(get != HMSH_STR_LENGTH * nFACE[0]) { fprintf(stderr, "ReadHybMeshBoundaryBinaryMPI --> MAIN: *** ERROR *** Can't read data from file \"%s\"\n", fname); MPI_Abort(comm, 0); }
		}
		// —◊»“€¬¿Õ»≈ Õ¿«¬¿Õ»… √–¿Õ»◊Õ€’ œŒ¬≈–’ÕŒ—“≈…

		// —◊»“€¬¿Õ»≈ Œœ»—¿Õ»… œŒ¬≈–’ÕŒ—“Õ€’ √–¿Õ≈… —≈“Œ◊Õ€’ ›À≈Ã≈Õ“Œ¬
		{
			EB[0] = (int *)malloc(3 * nEB[0] * sizeof(int));
			if(EB[0] == NULL) { fprintf(stderr, "ReadHybMeshBoundaryBinaryMPI --> MAIN: *** ERROR *** Can't allocate memory <%d * INT>\n", 3 * nEB[0]); MPI_Abort(comm, 0); }

			fseek(file, - nEB[0] * 3 * sizeof(int), SEEK_END);
			get = (int)fread(EB[0], sizeof(int), nEB[0] * 3, file);
			if(get != nEB[0] * 3) { fprintf(stderr, "ReadHybMeshBoundaryBinaryMPI --> MAIN: *** ERROR *** Can't read data from file \"%s\"\n", fname); MPI_Abort(comm, 0); }
		}
		// —◊»“€¬¿Õ»≈ Œœ»—¿Õ»… œŒ¬≈–’ÕŒ—“Õ€’ √–¿Õ≈… —≈“Œ◊Õ€’ ›À≈Ã≈Õ“Œ¬		

		// –¿——€À ¿ ƒ¿ÕÕ€’ Œ—“¿À‹Õ€Ã œ–Œ÷≈——¿Ã
		if(size > 1)
		{
			MPI_Bcast(      nFACE,                          1,  MPI_INT, 0, comm);
			MPI_Bcast(        nEB,                          1,  MPI_INT, 0, comm);
			MPI_Bcast(faceNAME[0], HMSH_STR_LENGTH * nFACE[0], MPI_CHAR, 0, comm);
			MPI_Bcast(      EB[0],                 nEB[0] * 3,  MPI_INT, 0, comm);
		}
		// –¿——€À ¿ ƒ¿ÕÕ€’ Œ—“¿À‹Õ€Ã œ–Œ÷≈——¿Ã

		fclose(file);
		if(outSTREAM != NULL) { fprintf(outSTREAM, "ReadHybMeshBoundaryBinaryMPI --> MAIN: Close data file \"%s\"\n", fname); }
	} // rank == 0
	
	if(rank != 0)
	{
		MPI_Bcast(      nFACE,                          1,  MPI_INT, 0, comm);
		MPI_Bcast(        nEB,                          1,  MPI_INT, 0, comm);

		faceNAME[0] = (char *)malloc(HMSH_STR_LENGTH * nFACE[0] * sizeof(char));
		if(faceNAME[0] == NULL) { fprintf(stderr, "ReadHybMeshBoundaryBinaryMPI --> [%d]: *** ERROR *** Can't allocate memory <%d * CHAR>\n", rank, HMSH_STR_LENGTH * nFACE[0]); MPI_Abort(comm, 0); }
		EB[0] = (int *)malloc(3 * nEB[0] * sizeof(int));
		if(EB[0] == NULL) { fprintf(stderr, "ReadHybMeshBoundaryBinaryMPI --> [%d]: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, 3 * nEB[0]); MPI_Abort(comm, 0); }

		MPI_Bcast(faceNAME[0], HMSH_STR_LENGTH * nFACE[0], MPI_CHAR, 0, comm);
		MPI_Bcast(      EB[0],                 nEB[0] * 3,  MPI_INT, 0, comm);
	} // rank != 0

	 MPI_Barrier(comm);

	return 0;
}

int ReadHybMeshNameBinaryMPI(char *fname, char **mname, FILE *outSTREAM, MPI_Comm comm)
{
	int rank, size;

	  MPI_Barrier(comm);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	  MPI_Barrier(comm);	

	if( (outSTREAM != NULL) && (rank == 0) )
	{ fprintf(outSTREAM, "ReadHybMeshNameBinaryMPI --> MAIN: START PROCEDURE IN %d MPI PROCESSES TOPOLOGY\n", size); }

	mname[0] = NULL; mname[0] = (char *)malloc(HMSH_STR_LENGTH * sizeof(char));
	if(mname[0] == NULL) { fprintf(stderr, "ReadHybMeshNameBinaryMPI --> MAIN: *** ERROR *** Can't allocate memory <%d * CHAR>\n", HMSH_STR_LENGTH); MPI_Abort(comm, 0); }
		
	if(rank == 0)
	{
		FILE *file = NULL;
		int get;
		
		file = fopen(fname, "rb+");
		if(file == NULL) { fprintf(stderr, "ReadHybMeshNameBinaryMPI --> MAIN: *** ERROR *** Can't open data file \"%s\"\n", fname); MPI_Abort(comm, 0); }
		else if(outSTREAM != NULL) { fprintf(outSTREAM, "ReadHybMeshNameBinaryMPI --> MAIN: Open data file \"%s\"\n", fname); }
		
		// Õ¿«¬¿Õ»≈ —≈“ »
		{
			int i; for(i=0; i<HMSH_STR_LENGTH; i++) mname[0][i] = '\0';
			get = (int)fread(mname[0], sizeof(char), HMSH_STR_LENGTH, file);
			if(get != HMSH_STR_LENGTH) { fprintf(stderr, "ReadHybMeshNameBinaryMPI --> MAIN: *** ERROR *** Can't read data from file \"%s\"\n", fname); MPI_Abort(comm, 0); }
			if(outSTREAM != NULL) { fprintf(outSTREAM, "ReadHybMeshNameBinaryMPI --> MAIN: MESH NAME = \"%s\"\n", mname[0]); }
		}
		// Õ¿«¬¿Õ»≈ —≈“ »

		fclose(file);
		if(outSTREAM != NULL) { fprintf(outSTREAM, "ReadHybMeshNameBinaryMPI --> MAIN: Close data file \"%s\"\n", fname); }
	}

	MPI_Barrier(comm);
	MPI_Bcast(mname[0], HMSH_STR_LENGTH, MPI_CHAR, 0, comm);	

	return 0;
}

