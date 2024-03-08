#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include <parmetis.h>
#include <metis.h>

#include "csrPartParMETISMPI.h"

int csrPartParMetisV4NoWeightsMPI(int *eDIST, int *csrX, int *csrA, int nPART, int *part, FILE *outSTREAM, MPI_Comm comm)
{
	int *vwgt, *adjwgt;
	// опциональные массивы весов вершин и ребер сетки. При этом каждая вершина и 
	// каждое ребро может иметь несколько весов.

	int wgtflag;
	// параметр, определяющий налицие весов у узлов и ребер сетки.
	// Может принимать 4 фиксированных значения:
	// 0 - нет весов (vwgt = adjwgt = NULL)
	// 1 - веса только у ребер сетки (vwgt = NULL)
	// 2 - веса только у узлов сетки (adjwgt = NULL)
	// 3 - веса у узлов и ребер сетки одновременно

	int numflag;
	// параметр, определяющий стиль нумерации узлов сетки
	// 0 - нумерация в стиле C (идентификаторы элементов начинаются с 0)
	// 1 - нумерация в стиле Fortran (идентификаторы элементов начинаются с 1)

	int ncon;
	// параметр, определяющий число весов, которое имеют вершины сетки	

	float *tpwgts;
	// массив размера ncon x nparts, в котором задается часть веса вершин, которая
	// будет относиться к каждой получаемой в результате разбиения подобласти сетки.
	// Если целью ставится получить подобласти одинакового размера с одинаковыми
	// весами, то значение каждого элемента массива устанавливается равным
	// 1/nparts. Сумма весов подобластей по каждому из ncon параметров должна
	// быть равна 1.

	float *ubvec;
	// массив из ncon элементов, в котором задается "дисбаланс и толерантность"
	// по каждому весу. Рекомендуемое разработчиками библиотеки подпрограмм
	// значение элементов массива равняется 1.05

	int *options;
	// массив целых чисел (integer) для задания дополнительных параметров функций
	// options[0] = 0 - значения параметров по умолчанию
	// options[0] = 1 - задание собственных параметров
	// задаваемые пользователем параметры различаются для различных функций
	// бибилиотеки, но можно выделить общие для функций, выполняющих разбиения
	// options[1] - задает уровень вывода "на консоль" в процессе выполнения
	// подпрограмм. Значение по умолчанию - 0, 1 - позволяет выводить информацию
	// о времени работы подпрограмм. Остальную дополнительную информацию можно
	// посмотреть или настроить в файле defs.h
	// options[2] - число, которое инициализирует некий генератор случайных
	// чисел. Значение по умолчанию равно 15.

	int edgecut;
	// в данную переменную записывается число разрезанных в результате разбиения
	// ребер исходной сетки. Судя по всему на всех процессорах будет записано
	// одно суммарное значение.

	int rank, size;
	double startTIME, finishTIME;
	int i;

	  MPI_Barrier(comm);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	  MPI_Barrier(comm);

	// веса узлов и ребер сетки
	vwgt = adjwgt = NULL;
	wgtflag = 0;

	// стиль нумерации узлов сетки
	numflag = 0;

	// параметры процедуры разбиения
	ncon = 1;
		
	tpwgts = (float *)malloc(ncon * nPART * sizeof(float));
	if(!tpwgts) { fprintf(stderr, "csrPartParMetisV4NoWeightsMPI --> %4d: *** ERROR *** Can't allocate memory <%d * float>\n", rank, ncon * nPART); MPI_Abort(comm, 0); }
	for(i=0; i<(ncon * nPART); i++) tpwgts[i] = (float)1.0 / ((float)nPART);

	ubvec = (float *)malloc(ncon * sizeof(float)); 
	if(!ubvec) { fprintf(stderr, "csrPartParMetisV4NoWeightsMPI --> %4d: *** ERROR *** Can't allocate memory <%d * float>\n", rank, ncon); MPI_Abort(comm, 0); }
	
	// for(i=0; i<ncon; i++) ubvec[i] = (float)1.05;
	   for(i=0; i<ncon; i++) ubvec[i] = (float)1.01;

	options = (int *)malloc(3 * sizeof(int));
	if(!options) { fprintf(stderr, "csrPartParMetisV4NoWeightsMPI --> %4d: *** ERROR *** Can't allocate memory <%d * int>\n", rank, 3); MPI_Abort(comm, 0); }
	options[0] = 1; options[1] = 1; options[2] = 15;

	MPI_Barrier(MPI_COMM_WORLD);
	if( (outSTREAM != NULL) && (rank == 0) )
	{ 
		fprintf(outSTREAM, "csrPartParMetisV4NoWeightsMPI --> MAIN:  START PARTITION PROCEDURE\n");
		fprintf(outSTREAM, "csrPartParMetisV4NoWeightsMPI --> MAIN:  NUMBER OF PARTS = %d\n", nPART);
		fflush(outSTREAM); 
	}
	startTIME = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);

	ParMETIS_V3_PartKway(eDIST, csrX, csrA, vwgt, adjwgt,
		                 &wgtflag, &numflag, &ncon, &nPART, tpwgts, ubvec,
				         options, &edgecut, part, &comm);

	MPI_Barrier(MPI_COMM_WORLD);
	finishTIME = MPI_Wtime();
	if( (outSTREAM != NULL) && (rank == 0) )
	{ fprintf(outSTREAM, "csrPartParMetisV4NoWeightsMPI --> MAIN: FINISH PARTITION PROCEDURE\n"); fflush(outSTREAM); }		
	MPI_Barrier(MPI_COMM_WORLD);

	if( (outSTREAM != NULL) && (rank == 0) )
	{ fprintf(outSTREAM, "csrPartParMetisV4NoWeightsMPI --> TIME: %4.2lf sec\n", finishTIME - startTIME); fflush(outSTREAM); }

	free(ubvec); free(tpwgts); free(options);

	return 0;
}

int csrPartParMetisV4NodesWeightsByCsrEdgesMPI(int *eDIST, int *csrX, int *csrA, int nPART, int *part, FILE *outSTREAM, MPI_Comm comm)
{
	int *vwgt, *adjwgt;
	// опциональные массивы весов вершин и ребер сетки. При этом каждая вершина и 
	// каждое ребро может иметь несколько весов.

	int wgtflag;
	// параметр, определяющий налицие весов у узлов и ребер сетки.
	// Может принимать 4 фиксированных значения:
	// 0 - нет весов (vwgt = adjwgt = NULL)
	// 1 - веса только у ребер сетки (vwgt = NULL)
	// 2 - веса только у узлов сетки (adjwgt = NULL)
	// 3 - веса у узлов и ребер сетки одновременно

	int numflag;
	// параметр, определяющий стиль нумерации узлов сетки
	// 0 - нумерация в стиле C (идентификаторы элементов начинаются с 0)
	// 1 - нумерация в стиле Fortran (идентификаторы элементов начинаются с 1)

	int ncon;
	// параметр, определяющий число весов, которое имеют вершины сетки	

	float *tpwgts;
	// массив размера ncon x nparts, в котором задается часть веса вершин, которая
	// будет относиться к каждой получаемой в результате разбиения подобласти сетки.
	// Если целью ставится получить подобласти одинакового размера с одинаковыми
	// весами, то значение каждого элемента массива устанавливается равным
	// 1/nparts. Сумма весов подобластей по каждому из ncon параметров должна
	// быть равна 1.

	float *ubvec;
	// массив из ncon элементов, в котором задается "дисбаланс и толерантность"
	// по каждому весу. Рекомендуемое разработчиками библиотеки подпрограмм
	// значение элементов массива равняется 1.05

	int *options;
	// массив целых чисел (integer) для задания дополнительных параметров функций
	// options[0] = 0 - значения параметров по умолчанию
	// options[0] = 1 - задание собственных параметров
	// задаваемые пользователем параметры различаются для различных функций
	// бибилиотеки, но можно выделить общие для функций, выполняющих разбиения
	// options[1] - задает уровень вывода "на консоль" в процессе выполнения
	// подпрограмм. Значение по умолчанию - 0, 1 - позволяет выводить информацию
	// о времени работы подпрограмм. Остальную дополнительную информацию можно
	// посмотреть или настроить в файле defs.h
	// options[2] - число, которое инициализирует некий генератор случайных
	// чисел. Значение по умолчанию равно 15.

	int edgecut;
	// в данную переменную записывается число разрезанных в результате разбиения
	// ребер исходной сетки. Судя по всему на всех процессорах будет записано
	// одно суммарное значение.

	int rank, size;
	double startTIME, finishTIME;
	int i;

	  MPI_Barrier(comm);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	  MPI_Barrier(comm);

	// веса узлов и ребер сетки
	wgtflag = 2;
	vwgt = (int *)malloc((eDIST[rank + 1] - eDIST[rank]) * sizeof(int));
	if(!vwgt) { fprintf(stderr, "csrPartParMetisV4NodesWeightsByCsrEdgesMPI --> %4d: *** ERROR *** Can't allocate memory <%d * int>\n", rank, eDIST[rank + 1] - eDIST[rank]); MPI_Abort(comm, 0); }
	for(i=0; i<eDIST[rank + 1] - eDIST[rank]; i++) vwgt[i] = csrX[i + 1] - csrX[i];
	adjwgt = NULL;	

	// стиль нумерации узлов сетки
	numflag = 0;

	// параметры процедуры разбиения
	ncon = 1;
		
	tpwgts = (float *)malloc(ncon * nPART * sizeof(float));
	if(!tpwgts) { fprintf(stderr, "csrPartParMetisV4NodesWeightsByCsrEdgesMPI --> %4d: *** ERROR *** Can't allocate memory <%d * float>\n", rank, ncon * nPART); MPI_Abort(comm, 0); }
	for(i=0; i<(ncon * nPART); i++) tpwgts[i] = (float)1.0 / ((float)nPART);

	ubvec = (float *)malloc(ncon * sizeof(float)); 
	if(!ubvec) { fprintf(stderr, "csrPartParMetisV4NodesWeightsByCsrEdgesMPI --> %4d: *** ERROR *** Can't allocate memory <%d * float>\n", rank, ncon); MPI_Abort(comm, 0); }
	
	// for(i=0; i<ncon; i++) ubvec[i] = (float)1.05;
	   for(i=0; i<ncon; i++) ubvec[i] = (float)1.01;

	options = (int *)malloc(3 * sizeof(int));
	if(!options) { fprintf(stderr, "csrPartParMetisV4NodesWeightsByCsrEdgesMPI --> %4d: *** ERROR *** Can't allocate memory <%d * int>\n", rank, 3); MPI_Abort(comm, 0); }
	options[0] = 1; options[1] = 1; options[2] = 15;

	MPI_Barrier(MPI_COMM_WORLD);
	if( (outSTREAM != NULL) && (rank == 0) )
	{
		fprintf(outSTREAM, "csrPartParMetisV4NodesWeightsByCsrEdgesMPI --> MAIN:  START PARTITION PROCEDURE\n"); 
		fprintf(outSTREAM, "csrPartParMetisV4NodesWeightsByCsrEdgesMPI --> MAIN:  NUMBER OF PARTS = %d\n", nPART); 
	}
	startTIME = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);

	ParMETIS_V3_PartKway(eDIST, csrX, csrA, vwgt, adjwgt,
		                 &wgtflag, &numflag, &ncon, &nPART, tpwgts, ubvec,
				         options, &edgecut, part, &comm);

	MPI_Barrier(MPI_COMM_WORLD);
	finishTIME = MPI_Wtime();
	if( (outSTREAM != NULL) && (rank == 0) )
	{ fprintf(outSTREAM, "csrPartParMetisV4NodesWeightsByCsrEdgesMPI --> MAIN: FINISH PARTITION PROCEDURE\n"); }		
	MPI_Barrier(MPI_COMM_WORLD);

	if( (outSTREAM != NULL) && (rank == 0) )
	{ 
		fprintf(outSTREAM, "csrPartParMetisV4NodesWeightsByCsrEdgesMPI --> TIME: %4.2lf sec\n", finishTIME - startTIME); 
	}

	free(vwgt); free(ubvec); free(tpwgts); free(options);

	return 0;
}

int csrPartMetisNodesEdgesWeightsSerial(int nE, int *csrX, int *csrA, int *wghtN, int *wghtE, int nPART, int *part, FILE *outSTREAM)
{
	int wgtflag = -1;
	// параметр, определяющий налицие весов у узлов и ребер сетки.
	// Может принимать 4 фиксированных значения:
	// 0 - нет весов (vwgt = adjwgt = NULL)
	// 1 - веса только у ребер сетки (vwgt = NULL)
	// 2 - веса только у узлов сетки (adjwgt = NULL)
	// 3 - веса у узлов и ребер сетки одновременно

	int *vwgt = NULL, *adjwgt = NULL; // массивы весов вершин графа и ребер соответственно
	int numflag = 0; // флаг стиля нумерации узлов				
	int options[5] = {0, 0, 0, 0, 0}; // параметры алгоритма разбиения
	int edgecut = 0; // число разрезанных в результате разбиения ребер сетки
	int nparts = nPART;
	
	if( (wgtflag == -1) && (wghtN == NULL) && (wghtE == NULL) ) wgtflag = 0;
	if( (wgtflag == -1) && (wghtN == NULL) && (wghtE != NULL) ) wgtflag = 1;
	if( (wgtflag == -1) && (wghtN != NULL) && (wghtE == NULL) ) wgtflag = 2;
	if( (wgtflag == -1) && (wghtN != NULL) && (wghtE != NULL) ) wgtflag = 3;

	  vwgt = wghtN;
	adjwgt = wghtE;
		
	if(outSTREAM != NULL) 
	{
		fprintf(outSTREAM, "\nStart METIS_PartGraphKway .....\n");
		if( (wgtflag == 2) || (wgtflag == 3) ) fprintf(outSTREAM, "VERTEX WEIGHTS: YES\n");		                                       
		                                  else fprintf(outSTREAM, "VERTEX WEIGHTS:  NO\n");
		if( (wgtflag == 1) || (wgtflag == 3) ) fprintf(outSTREAM, " EDGES WEIGHTS: YES\n");
		                                  else fprintf(outSTREAM, " EDGES WEIGHTS:  NO\n");
		fprintf(outSTREAM, "NUMBER OF DOMAINS: %8d\n", nparts);
	}
	
	METIS_PartGraphKway(&nE, csrX, csrA, vwgt, adjwgt, &wgtflag, &numflag, &nparts, options, &edgecut, part); 
		
	if(outSTREAM != NULL) 
	{
		int *wghtPART, wghtMAX, wghtTOTAL;
		int wghtEdgesTotal, wghtEdgesCut;
		int i;

		if( (wgtflag == 2) || (wgtflag == 3) )
		{
			wghtPART = (int *)malloc(nparts * sizeof(int));
			if(wghtPART == NULL) { fprintf(stderr, "csrPartMetisNodesEdgesWeightsSerial --> *** ERROR *** Can't allocate memory <%d * int>\n", nparts); exit(0); }
			wghtTOTAL = 0; for(i=0; i<nPART; i++) wghtPART[i] = 0;
			for(i=0; i<nE; i++) { wghtTOTAL += vwgt[i]; wghtPART[part[i]] += vwgt[i]; }
		} // if
		else
		{
			wghtPART = (int *)malloc(nparts * sizeof(int));
			if(wghtPART == NULL) { fprintf(stderr, "csrPartMetisNodesEdgesWeightsSerial --> *** ERROR *** Can't allocate memory <%d * int>\n", nparts); exit(0); }
			wghtTOTAL = 0; for(i=0; i<nPART; i++) wghtPART[i] = 0;
			for(i=0; i<nE; i++) { wghtTOTAL += 1; wghtPART[part[i]] += 1; }
		} // else

		wghtMAX = wghtPART[0];
		for(i=1; i<nparts; i++) if(wghtMAX < wghtPART[i]) wghtMAX = wghtPART[i];		
		fprintf(outSTREAM, "wMAX = %d wTOTAL = %d\n", wghtMAX, wghtTOTAL);
		fprintf(outSTREAM, "DISBALANCE: %6.3lf\n", (double)wghtMAX / ((double)wghtTOTAL / (double)nparts) );
		free(wghtPART);
		
		wghtEdgesTotal = wghtEdgesCut = 0;
		if( (wgtflag == 1) || (wgtflag == 3) )
		{
			for(i=0; i<csrX[nE]; i++) wghtEdgesTotal += wghtE[i]; wghtEdgesTotal /= 2;
			for(i=0; i<nE; i++)
			{
				int j; for(j=csrX[i]; j<csrX[i+1]; j++) if(part[i] != part[csrA[j]]) wghtEdgesCut += adjwgt[j];		
			} // for i
			wghtEdgesCut /= 2;
		} // if
		else
		{
			wghtEdgesTotal = csrX[nE] / 2;
			for(i=0; i<nE; i++)
			{
				int j; for(j=csrX[i]; j<csrX[i+1]; j++) if(part[i] != part[csrA[j]]) wghtEdgesCut += 1;		
			} // for i
			wghtEdgesCut /= 2;			
		} // else

		fprintf(outSTREAM, "          edgecut: %8d (%.2lf %%)\n", wghtEdgesCut, 100.0 * (double)wghtEdgesCut/(double)wghtEdgesTotal);
		fprintf(outSTREAM, "      REF edgecut: %8d           \n", edgecut);		
	} // if

	return 0;
}
