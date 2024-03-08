#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "processIntegerLIST.h"
#include "hm4PreprocessorMPI.h"

int hm4CreateProcessElementsListMPI(int *eDIST, int *part, int *nePROC, int **ePROC, FILE *outSTREAM, MPI_Comm comm)
{
	int rank, size;

	int ne, *elts, eltsSIZE;
	int *buf, bufSIZE;

	int i, j;
	
	  MPI_Barrier(comm);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	  MPI_Barrier(comm);

	if( (outSTREAM != NULL) && (rank == 0) )
	{ fprintf(outSTREAM, "hm4CreateProcessElementsListMPI --> MAIN: START PROCEDURE IN %d MPI PROCESSES TOPOLOGY\n", size); }

	// буфер для приема частей распределенного массива "раскраски" от других MPI-процессов
	// размер буфера - максимальная длина из частей распределенного массива
	bufSIZE = eDIST[1] - eDIST[0];
	for(i=1; i<size; i++) if(bufSIZE < eDIST[i+1] - eDIST[i]) bufSIZE = eDIST[i+1] - eDIST[i];
	buf = (int *)malloc(bufSIZE * sizeof(int));
	if(buf == NULL) { fprintf(stderr, "hm4CreateProcessElementsListMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, bufSIZE); MPI_Abort(comm, 0); }

	// буфер для формирования списка индексов собственных элементов части сетки MPI-процесса
	// начальный размер буфера - удвоенное среднее арифметическое число сеточных элементов на процесс
	eltsSIZE = 2 * eDIST[size] / size;
	elts = (int *)malloc(eltsSIZE * sizeof(int));
	if(elts == NULL) { fprintf(stderr, "hm4CreateProcessElementsListMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, eltsSIZE); MPI_Abort(comm, 0); }
	ne = 0;

	MPI_Barrier(comm);
	
	// цикл по числу параллельных процессов
	for(i=0; i<size; i++)
	{
		int *iPTR; // указатель на обрабатываемую часть массива "раскраски"

		// рассылка обрабатываемой части массива
		// инициализация указателя на обрабатываемую часть массива
		if(i == rank)
		{ iPTR = part; MPI_Bcast(part, eDIST[i+1] - eDIST[i], MPI_INT, i, comm); }
		else
		{ iPTR =  buf; MPI_Bcast( buf, eDIST[i+1] - eDIST[i], MPI_INT, i, comm); }

		// цикл по элементам обрабатываемой части
		for(j=0; j<eDIST[i+1] - eDIST[i]; j++)
		// добавление подходящего элемента в список
		if(iPTR[j] == rank)
		{
			// дополнительное выделение памяти
			if(ne == eltsSIZE)
			{
				int *irPTR = (int *)realloc(elts, (eltsSIZE + 1000) * sizeof(int));
				if(irPTR == NULL) { fprintf(stderr, "hm4CreateProcessElementsListMPI --> %4d: *** ERROR *** Can't reallocate memory <%d * INT>\n", rank, eltsSIZE + 1000); MPI_Abort(comm, 0); }
				if(irPTR != elts) elts = irPTR; eltsSIZE += 1000;
			}
			// добавление элемента
			elts[ne] = eDIST[i] + j;
			ne += 1;		
		} // for j
	} // for i

	// выравнивание распределенной под массив памяти к фактической длине массива
	if(ne != eltsSIZE)
	{
		int *irPTR = (int *)realloc(elts, ne * sizeof(int));
		if(irPTR == NULL) { fprintf(stderr, "hm4CreateProcessElementsListMPI --> %4d: *** ERROR *** Can't reallocate memory <%d * INT>\n", rank, ne); MPI_Abort(comm, 0); }
		if(irPTR != elts) elts = irPTR;
	}

	// освобождение выделенный под буфер приема данных памяти
	free(buf);

	// инициализация выходных параметров подпрограммы
	nePROC[0] = ne; ePROC[0] = elts;
	
	MPI_Barrier(comm);

	// вывод статистической информации относительно баланса распределения сеточных элементов
	if(rank == 0)
	{
		int nMIN, nMAX, nTOTAL;
		int *ePDIST = (int *)malloc(size * sizeof(int));
		if(ePDIST == NULL) { fprintf(stderr, "hm4CreateProcessElementsListMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, size); MPI_Abort(comm, 0); }

		ePDIST[0] = ne;
		for(i=1; i<size; i++)
		{
			MPI_Status status;
			MPI_Recv(&ePDIST[i], 1, MPI_INT, i, 752, comm, &status);
		} // for i

		nMIN = nMAX = nTOTAL = ePDIST[0];
		for(i=1; i<size; i++)
		{
			if(nMIN > ePDIST[i]) nMIN    = ePDIST[i];
			if(nMAX < ePDIST[i]) nMAX    = ePDIST[i];
			                     nTOTAL += ePDIST[i];
		} // for i

		if(nTOTAL != eDIST[size]) { fprintf(stderr, "hm4CreateProcessElementsListMPI --> %4d: *** ERROR *** nTOTAL <%d> != eDIST[size] <%d>\n", rank, nTOTAL, eDIST[size]); MPI_Abort(comm, 0); }
		if(outSTREAM != NULL)
		{
			fprintf(outSTREAM, "\n");
			fprintf(outSTREAM, "ELEMENTS DISTRIBUTION FOR %d SUBAREAS\n", size);
			fprintf(outSTREAM, "nMIN = %8d nMAX = %d\n", nMIN, nMAX);
			for(i=0; i<size; i++)
			{
				fprintf(outSTREAM, "%8d ", ePDIST[i]);
				if((i + 1) % 10 == 0) fprintf(outSTREAM, "\n");
			} // for i
			if(i % 10 != 0) fprintf(outSTREAM, "\n");
		} // if

		free(ePDIST);
	}
	else
	{
		MPI_Send(&ne, 1, MPI_INT, 0, 752, comm);
	}

	MPI_Barrier(comm);

	return 0;
}

int hm4InitElementsLocationMPI(int *eDIST, int *part, int ne, int *elts, int *eltsLOCATION, FILE *outSTREAM, MPI_Comm comm)
{
	int rank, size;

	int *buf, bufSIZE;

	int i, j;
	
	  MPI_Barrier(comm);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	  MPI_Barrier(comm);

	if( (outSTREAM != NULL) && (rank == 0) )
	{ fprintf(outSTREAM, "hm4InitElementsLocationMPI --> MAIN: START PROCEDURE IN %d MPI PROCESSES TOPOLOGY\n", size); }

	bufSIZE = eDIST[1] - eDIST[0];
	for(i=1; i<size; i++) if(bufSIZE < eDIST[i+1] - eDIST[i]) bufSIZE = eDIST[i+1] - eDIST[i];
	buf = (int *)malloc(bufSIZE * sizeof(int));
	if(buf == NULL) { fprintf(stderr, "hm4InitElementsLocationMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, bufSIZE); MPI_Abort(comm, 0); }

	for(i=0; i<ne; i++) eltsLOCATION[i] = -1;

	MPI_Barrier(comm);
	
	for(i=0; i<size; i++)
	{
		int *iPTR;

		if(i == rank)
		{ iPTR = part; MPI_Bcast(part, eDIST[i+1] - eDIST[i], MPI_INT, i, comm); }
		else
		{ iPTR =  buf; MPI_Bcast( buf, eDIST[i+1] - eDIST[i], MPI_INT, i, comm); }

		for(j=0; j<ne; j++) if( (elts[j] >= eDIST[i]) && (elts[j] < eDIST[i + 1]) ) eltsLOCATION[j] = iPTR[elts[j] - eDIST[i]];		
	} // for i

	free(buf);

	MPI_Barrier(comm);
	
	for(j=0, i=0; i<ne; i++) if(eltsLOCATION[i] == -1) j += 1;
	if(j != 0) { fprintf(stderr, "hm4InitElementsLocationMPI --> %4d: *** ERROR *** Not all (%d) elements of eltsLOCAL array inited\n", rank, j); MPI_Abort(comm, 0); }
		
	MPI_Barrier(comm);	

	return 0;
}

// ПРОВЕРИТЬ: ДОПУСКАЕТСЯ ЛИ АЛГОРИТМОМ СИТУАЦИЯ, КОГДА ИНИЦИАЛИЗИРУЕМЫЙ МАССИВ eltsL БУДЕТ ПУСТЫМ (neL = 0) ????????
// АЛГОРИТМ: КАЖДЫЙ РАССЫЛАЕТ КАЖДОМУ СВОЙ СПИСОК, СОБИРАЕТ РЕЗУЛЬТАТ И УДАЛЯЕТ ПОВТОРЕНИЯ
int hm4InitElementsLinksMPI(int *eDIST, int *xCSR, int *aCSR, int ne, int *elts, int *neL, int **eltsL, FILE *outSTREAM, MPI_Comm comm)
{
	int rank, size;

	int *requestPROC; // массив с длинами запросов (длины списков элементов) от каждого процесса
	int *requestELTS; // массив для приема списка запрашиваемых элементов
	int replyN, *replyELTS; // массив для формирования списка связей для запрашиваемых элементов
	int eltsLLengh; // фактическая длина выделенного под запись результата массива

	int i;

	  MPI_Barrier(comm);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	  MPI_Barrier(comm);

	if( (outSTREAM != NULL) && (rank == 0) )
	{ fprintf(outSTREAM, "hm4InitElementsLinksMPI --> MAIN: START PROCEDURE IN %d MPI PROCESSES TOPOLOGY\n", size); fflush(outSTREAM); }

	// ИНИЦИАЛИЗАЦИЯ ДЛИН ЗАПРОСОВ
	{
		requestPROC = (int *)malloc(size * sizeof(int));
		if(requestPROC == NULL) { fprintf(stderr, "hm4InitElementsLinksMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, size); MPI_Abort(comm, 0); }

		if(rank == 0)
		{
			requestPROC[0] = ne;
			for(i=1; i<size; i++)
			{
				MPI_Status status;
				MPI_Recv(&requestPROC[i], 1, MPI_INT, i, 1049, comm, &status);
			} // for i
			MPI_Bcast(requestPROC, size, MPI_INT, 0, comm);
		}
		else
		{
			MPI_Send(&ne, 1, MPI_INT, 0, 1049, comm);
			MPI_Bcast(requestPROC, size, MPI_INT, 0, comm);
		}
	}
	// ИНИЦИАЛИЗАЦИЯ ДЛИН ЗАПРОСОВ

	// ВЫДЕЛЕНИЕ ВСПОМОГАТЕЛЬНОЙ ПАМЯТИ
	{
		int requestMAX = requestPROC[0];
		
		for(i=1; i<size; i++) if(requestMAX < requestPROC[i]) requestMAX = requestPROC[i];
		if( (outSTREAM != NULL) && (rank == 0) ){ fprintf(outSTREAM, "hm4InitElementsLinksMPI --> MAIN: MAX REQUEST SIZE = %d ELEMENTS\n", requestMAX); }

		requestELTS = (int *)malloc(requestMAX * sizeof(int));
		if(requestPROC == NULL) { fprintf(stderr, "hm4InitElementsLinksMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, requestMAX); MPI_Abort(comm, 0); }

		replyELTS = (int *)malloc(xCSR[eDIST[rank + 1] - eDIST[rank]] * sizeof(int));
		if(replyELTS == NULL) { fprintf(stderr, "hm4InitElementsLinksMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, xCSR[eDIST[rank + 1] - eDIST[rank]]); MPI_Abort(comm, 0); }
	}
	// ВЫДЕЛЕНИЕ ВСПОМОГАТЕЛЬНОЙ ПАМЯТИ

	// ОСНОВНОЙ ЦИКЛ
	for(i=0; i<size; i++)
	{
		int j, k;

		if(rank == i)
		{
			MPI_Bcast(elts, requestPROC[i], MPI_INT, i, comm);

			replyN = 0;
			for(j=0; j<requestPROC[i]; j++)
				if( (elts[j] >= eDIST[rank]) && (elts[j] < eDIST[rank + 1]) )
				{
					int offset = elts[j] - eDIST[rank];
					for(k=xCSR[offset]; k<xCSR[offset + 1]; k++) 
					{
						replyELTS[replyN] = aCSR[k];						
						replyN += 1;
					}
				} // if

			eltsLLengh = requestPROC[i] * 6;
			eltsL[0] = (int *)malloc(eltsLLengh * sizeof(int));
			if(eltsL[0] == NULL) { fprintf(stderr, "hm4InitElementsLinksMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, eltsLLengh); MPI_Abort(comm, 0); }
			neL[0] = 0;

			if(replyN > 0)
			{
				SortAndClearIntegerListUnit(replyN, replyELTS, &replyN);
				DeleteListFromOrderedList(requestPROC[i], elts, replyN, replyELTS, &replyN);
			}

			if(replyN != 0)
			{
				if(neL[0] + replyN > eltsLLengh)
				{
					int *irPTR = (int *)realloc(eltsL[0], (neL[0] + replyN + 10000) * sizeof(int));
					if(irPTR == NULL) { fprintf(stderr, "hm4InitElementsLinksMPI --> %4d: *** ERROR *** Can't reallocate memory <%d * INT>\n", rank, neL[0] + replyN + 10000); MPI_Abort(comm, 0); }
					if(irPTR != eltsL[0]) eltsL[0] = irPTR; eltsLLengh = neL[0] + replyN + 10000;
				} // if

				for(j=0; j<replyN; j++) eltsL[0][j] = replyELTS[j];
				neL[0] = replyN;				
			}

			for(j=0; j<size; j++) if(j != rank)
			{
				int nRECV;
				MPI_Status status;
				
				MPI_Recv(&nRECV, 1, MPI_INT, j, 2048, comm, &status);
				if(nRECV > 0)
				{
					if(neL[0] + nRECV > eltsLLengh)
					{
						int *irPTR = (int *)realloc(eltsL[0], (neL[0] + nRECV + 10000) * sizeof(int));
						if(irPTR == NULL) { fprintf(stderr, "hm4InitElementsLinksMPI --> %4d: *** ERROR *** Can't reallocate memory <%d * INT>\n", rank, neL[0] + nRECV + 10000); MPI_Abort(comm, 0); }
						if(irPTR != eltsL[0]) eltsL[0] = irPTR; eltsLLengh = neL[0] + nRECV + 10000;
					} // if
					
					MPI_Recv(eltsL[0] + neL[0], nRECV, MPI_INT, j, 3047, comm, &status);
					neL[0] += nRECV;
				} // if
			} // for j

			if(neL[0] > 0)
			{
				SortAndClearIntegerListUnit(neL[0], eltsL[0], neL);
				if(neL[0] != eltsLLengh)
				{
					int *irPTR = (int *)realloc(eltsL[0], neL[0] * sizeof(int));
					if(irPTR == NULL) { fprintf(stderr, "hm4InitElementsLinksMPI --> %4d: *** ERROR *** Can't reallocate memory <%d * INT>\n", rank, neL[0]); MPI_Abort(comm, 0); }
					if(irPTR != eltsL[0]) eltsL[0] = irPTR;
				} // if
			}
			else
			{
				free(eltsL[0]); neL[0] = 0;
			}
		}
		else
		{
			MPI_Bcast(requestELTS, requestPROC[i], MPI_INT, i, comm);

			replyN = 0;
			for(j=0; j<requestPROC[i]; j++)
				if( (requestELTS[j] >= eDIST[rank]) && (requestELTS[j] < eDIST[rank + 1]) )
				{
					int offset = requestELTS[j] - eDIST[rank];
					for(k=xCSR[offset]; k<xCSR[offset + 1]; k++)
					{
						replyELTS[replyN] = aCSR[k];						
						replyN += 1;
					}
				} // if

			if(replyN > 0)
			{
				SortAndClearIntegerListUnit(replyN, replyELTS, &replyN);
				DeleteListFromOrderedList(requestPROC[i], requestELTS, replyN, replyELTS, &replyN);
			}

			MPI_Send(&replyN, 1, MPI_INT, i, 2048, comm);
			if(replyN > 0) MPI_Send(replyELTS, replyN, MPI_INT, i, 3047, comm);
		}
	} // for i
	// ОСНОВНОЙ ЦИКЛ

	// ОСВОБОЖДЕНИЕ ПАМЯТИ
	{
		free(requestPROC);
		free(requestELTS);
		free(replyELTS);
	}
	// ОСВОБОЖДЕНИЕ ПАМЯТИ
	
	MPI_Barrier(comm);

	return 0;
}

int hm4InitProcessRecvListMPI(int *eDIST, int *part, int neL, int *eltsL, int *xRECV, int *aRECV, FILE *outSTREAM, MPI_Comm comm)
{
	int rank, size;
	int *eltsLOCATION;

	int i;

	  MPI_Barrier(comm);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	  MPI_Barrier(comm);

	if( (outSTREAM != NULL) && (rank == 0) )
	{ fprintf(outSTREAM, "hm4InitProcessRecvListMPI --> MAIN: START PROCEDURE IN %d MPI PROCESSES TOPOLOGY\n", size); }

	eltsLOCATION = (int *)malloc(neL * sizeof(int));
	if(eltsLOCATION == NULL) { fprintf(stderr, "hm4InitProcessRecvListMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, neL); MPI_Abort(comm, 0); }
	hm4InitElementsLocationMPI(eDIST, part, neL, eltsL, eltsLOCATION, NULL, comm);
	for(i=0; i<neL; i++) if( (eltsLOCATION[i] < 0) || (eltsLOCATION[i] == rank) || (eltsLOCATION[i] >= size) )
	{ fprintf(stderr, "hm4InitProcessRecvListMPI --> %4d: *** ERROR *** (eltsLOCATION[i] < 0) || (eltsLOCATION[i] == rank) || (eltsLOCATION[i] >= size)\n", 
	                   rank); MPI_Abort(comm, 0); }

	for(i=0; i<=size; i++) xRECV[i] =  0;
	for(i=0; i<  neL; i++) aRECV[i] = -1;
	for(i=0; i<neL; i++) xRECV[eltsLOCATION[i] + 1] += 1;
	if(xRECV[rank + 1] != 0) { fprintf(stderr, "hm4InitProcessRecvListMPI --> %4d: *** ERROR *** xRECV[rank + 1] != 0 ( %d != 0 )\n", rank, xRECV[rank + 1]); MPI_Abort(comm, 0); }
	for(i=2; i<=size; i++) xRECV[i] += xRECV[i-1];
	if(xRECV[size] != neL) { fprintf(stderr, "hm4InitProcessRecvListMPI --> %4d: *** ERROR *** xRECV[size] != neL ( %d != %d )\n", rank, xRECV[size], neL); MPI_Abort(comm, 0); }
	for(i=0; i<neL; i++) aRECV[xRECV[eltsLOCATION[i]]] = eltsL[i], xRECV[eltsLOCATION[i]] += 1;
	if(xRECV[size] != xRECV[size - 1]) { fprintf(stderr, "hm4InitProcessRecvListMPI --> %4d: *** ERROR *** xRECV[size] != xRECV[size - 1] ( %d != %d )\n", rank, xRECV[size], xRECV[size - 1]); MPI_Abort(comm, 0); }
	for(i=size-1; i>0; i--) xRECV[i] = xRECV[i-1]; xRECV[0] = 0;
	for(i=0; i<size; i++) if( (i != rank) && (xRECV[i+1] - xRECV[i] > 0) ) wsortIntegerListUnit(aRECV, xRECV[i], xRECV[i+1] - 1);
	
	free(eltsLOCATION);
	
	MPI_Barrier(comm);

	return 0;
}

int hm4InitProcessSendListMPI(int *xRECV, int *aRECV, int **xSEND, int **aSEND, FILE *outSTREAM, MPI_Comm comm)
{
	int rank, size;
	int sendN = 50000;
		
	int i, j;

	  MPI_Barrier(comm);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	  MPI_Barrier(comm);

	if( (outSTREAM != NULL) && (rank == 0) )
	{ fprintf(outSTREAM, "hm4InitProcessSendListMPI --> MAIN: START PROCEDURE IN %d MPI PROCESSES TOPOLOGY\n", size); }

	xSEND[0] = (int *)malloc((size + 1) * sizeof(int));
	if(xSEND[0] == NULL) { fprintf(stderr, "hm4InitProcessSendListMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, size + 1); MPI_Abort(comm, 0); }
	xSEND[0][0] = 0; for(i=1; i<=size; i++) xSEND[0][i] = -1;
	aSEND[0] = (int *)malloc(sendN * sizeof(int));
	if(aSEND[0] == NULL) { fprintf(stderr, "hm4InitProcessSendListMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, sendN); MPI_Abort(comm, 0); }

	for(i=0; i<size; i++)
	{
		if(rank == i)
		{
			xSEND[0][i + 1] = xSEND[0][i];
			
			for(j=rank+1; j<size; j++)
			{
				int nSEND;
				int nRECV = xRECV[j + 1] - xRECV[j];
				MPI_Status status;
				
				MPI_Send(&nRECV, 1, MPI_INT, j, 1184, comm);
				if(nRECV > 0) MPI_Send(aRECV + xRECV[j], nRECV, MPI_INT, j, 1283, comm);

				MPI_Recv(&nSEND, 1, MPI_INT, j, 1382, comm, &status);
				if(nSEND > 0)
				{
					xSEND[0][j + 1] = xSEND[0][j] + nSEND;
					if(xSEND[0][j + 1] > sendN)
					{
						int *irPTR = (int *)realloc(aSEND[0], (xSEND[0][j + 1] + 1000) * sizeof(int));
						if(irPTR == NULL) { fprintf(stderr, "hm4InitProcessSendListMPI --> %4d: *** ERROR *** Can't reallocate memory <%d * INT>\n", rank, xSEND[0][j + 1] + 1000); MPI_Abort(comm, 0); }
						if(irPTR != aSEND[0]) aSEND[0] = irPTR; sendN = xSEND[0][j + 1] + 1000;
					} // realloc
					MPI_Recv(aSEND[0] + xSEND[0][j], nSEND, MPI_INT, j, 1481, comm, &status);
				}
				else xSEND[0][j + 1] = xSEND[0][j];				
			} // for j			
		} // rank == i

		if(rank > i)
		{
			int nSEND;
			int nRECV = xRECV[i+1] - xRECV[i];
			MPI_Status status;
			
			MPI_Recv(&nSEND, 1, MPI_INT, i, 1184, comm, &status);			
			if(nSEND > 0)
			{
				xSEND[0][i + 1] = xSEND[0][i] + nSEND;
				if(xSEND[0][i + 1] > sendN)
				{
					int *irPTR = (int *)realloc(aSEND[0], (xSEND[0][i + 1] + 1000) * sizeof(int));
					if(irPTR == NULL) { fprintf(stderr, "hm4InitProcessSendListMPI --> %4d: *** ERROR *** Can't reallocate memory <%d * INT>\n", rank, xSEND[0][i + 1] + 1000); MPI_Abort(comm, 0); }
					if(irPTR != aSEND[0]) aSEND[0] = irPTR; sendN = xSEND[0][i + 1] + 1000;
				} // realloc
				MPI_Recv(aSEND[0] + xSEND[0][i], nSEND, MPI_INT, i, 1283, comm, &status);
			}
			else xSEND[0][i + 1] = xSEND[0][i];

			MPI_Send(&nRECV, 1, MPI_INT, i, 1382, comm);
			if(nRECV > 0) MPI_Send(aRECV + xRECV[i], nRECV, MPI_INT, i, 1481, comm);
		} // rank > i
	} // for i

	if(xSEND[0][size] > sendN)
	{
		int *irPTR = (int *)realloc(aSEND[0], xSEND[0][size] * sizeof(int));
		if(irPTR == NULL) { fprintf(stderr, "hm4InitProcessSendListMPI --> %4d: *** ERROR *** Can't reallocate memory <%d * INT>\n", rank, xSEND[0][size]); MPI_Abort(comm, 0); }
		if(irPTR != aSEND[0]) aSEND[0] = irPTR;
	} // realloc

	return 0;
}

int hm4CheckSendRecvDataMPI(int *xRECV, int *aRECV, int *xSEND, int *aSEND, FILE *outSTREAM, MPI_Comm comm)
{
	int rank, size;

	int i, j;

	  MPI_Barrier(comm);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	  MPI_Barrier(comm);

	if( (outSTREAM != NULL) && (rank == 0) )
	{ fprintf(outSTREAM, "hm4CheckSendRecvDataMPI --> MAIN: START PROCEDURE IN %d MPI PROCESSES TOPOLOGY\n", size); }

	if(rank == 0)
	{
		int *xRALL, *xSALL;
		int *aRALL, *aSALL;
		int *offsetR, *offsetS;

		xRALL = (int *)malloc(size * (size + 1) * sizeof(int));
		if(xRALL == NULL) { fprintf(stderr, "hm4CheckSendRecvDataMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, size * (size + 1)); MPI_Abort(comm, 0); }
		xSALL = (int *)malloc(size * (size + 1) * sizeof(int));
		if(xSALL == NULL) { fprintf(stderr, "hm4CheckSendRecvDataMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, size * (size + 1)); MPI_Abort(comm, 0); }

		offsetR = (int *)malloc((size + 1) * sizeof(int));
		if(offsetR == NULL) { fprintf(stderr, "hm4CheckSendRecvDataMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, size + 1); MPI_Abort(comm, 0); }
		offsetS = (int *)malloc((size + 1) * sizeof(int));
		if(offsetS == NULL) { fprintf(stderr, "hm4CheckSendRecvDataMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, size + 1); MPI_Abort(comm, 0); }

		for(i=0; i<=size; i++) xRALL[i] = xRECV[i], xSALL[i] = xSEND[i];
		offsetR[0] = 0; offsetR[1] = xRECV[size];
		offsetS[0] = 0; offsetS[1] = xSEND[size];

		for(i=1; i<size; i++)
		{
			MPI_Status status;
			
			MPI_Recv(xRALL + i*(size + 1), size + 1, MPI_INT, i, 4028 + i, comm, &status);
			MPI_Recv(xSALL + i*(size + 1), size + 1, MPI_INT, i, 5027 + i, comm, &status);

			offsetR[i + 1] = offsetR[i] + xRALL[i*(size + 1) + size];
			offsetS[i + 1] = offsetS[i] + xSALL[i*(size + 1) + size];			
		} // for i

		aRALL = (int *)malloc(offsetR[size] * sizeof(int));
		if(aRALL == NULL) { fprintf(stderr, "hm4CheckSendRecvDataMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, offsetR[size]); MPI_Abort(comm, 0); }
		for(i=0; i<xRALL[size]; i++) aRALL[i] = aRECV[i];

		aSALL = (int *)malloc(offsetS[size] * sizeof(int));
		if(aSALL == NULL) { fprintf(stderr, "hm4CheckSendRecvDataMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, offsetS[size]); MPI_Abort(comm, 0); }
		for(i=0; i<xSALL[size]; i++) aSALL[i] = aSEND[i];

		for(i=1; i<size; i++)
		{
			MPI_Status status;
			
			MPI_Recv(aRALL + offsetR[i], offsetR[i + 1] - offsetR[i], MPI_INT, i, 6026 + i, comm, &status);
			MPI_Recv(aSALL + offsetS[i], offsetS[i + 1] - offsetS[i], MPI_INT, i, 7025 + i, comm, &status);			
		} // for i

		for(i=0; i<size-1; i++)
		{
			int *xRMAIN = xRALL + i*(size + 1);
			int *aRMAIN = aRALL + offsetR[i];
			int *xSMAIN = xSALL + i*(size + 1);
			int *aSMAIN = aSALL + offsetS[i];

			for(j=i+1; j<size; j++)
			{
				int *xRLINK = xRALL + j*(size + 1);
				int *aRLINK = aRALL + offsetR[j];
				int *xSLINK = xSALL + j*(size + 1);
				int *aSLINK = aSALL + offsetS[j];

				int nMainFromLink = xRMAIN[j + 1] - xRMAIN[j];
				int   nMainToLink = xSMAIN[j + 1] - xSMAIN[j];
				int nLinkFromMain = xRLINK[i + 1] - xRLINK[i];
				int   nLinkToMain = xSLINK[i + 1] - xSLINK[i];

				int k;

				if( (nMainFromLink != nLinkToMain)  || (nMainToLink != nLinkFromMain) )				
				{ 
					fprintf(stderr, "hm4CheckSendRecvDataMPI --> MAIN: *** ERROR *** %04d <--> %04d\n", i, j);
					fprintf(stderr, "hm4CheckSendRecvDataMPI --> MAIN: nMainFromLink (%d) != nLinkToMain (%d)\n", nMainFromLink, nLinkToMain);
					fprintf(stderr, "hm4CheckSendRecvDataMPI --> MAIN: nMainToLink (%d) != nLinkFromMain (%d)\n", nMainToLink, nLinkFromMain);
					MPI_Abort(comm, 0);					
				}

				for(k=0; k<nMainFromLink; k++) if(aRMAIN[xRMAIN[j] + k] != aSLINK[xSLINK[i] + k])
				{
					fprintf(stderr, "hm4CheckSendRecvDataMPI --> MAIN: *** ERROR *** %04d <--> %04d\n", i, j);
					fprintf(stderr, "hm4CheckSendRecvDataMPI --> MAIN: RECV %04d LIST[%d] = %d\n", i, k, aRMAIN[xRMAIN[j] + k]);
					fprintf(stderr, "hm4CheckSendRecvDataMPI --> MAIN: SEND %04d LIST[%d] = %d\n", j, k, aSLINK[xSLINK[i] + k]);
					MPI_Abort(comm, 0);
				} // if in for k

				for(k=0; k<nMainToLink; k++) if(aSMAIN[xSMAIN[j] + k] != aRLINK[xRLINK[i] + k])
				{
					fprintf(stderr, "hm4CheckSendRecvDataMPI --> MAIN: *** ERROR *** %04d <--> %04d\n", i, j);
					fprintf(stderr, "hm4CheckSendRecvDataMPI --> MAIN: SEND %04d LIST[%d] = %d\n", i, k, aSMAIN[xSMAIN[j] + k]);
					fprintf(stderr, "hm4CheckSendRecvDataMPI --> MAIN: RECV %04d LIST[%d] = %d\n", j, k, aRLINK[xRLINK[i] + k]);
					MPI_Abort(comm, 0);
				} // if in for k
			} // for j			
		} // for i

		  free(xRALL);   free(xSALL);
		  free(aRALL);   free(aSALL);
		free(offsetR); free(offsetS);		
	}
	else
	{
		MPI_Send(xRECV,    size + 1, MPI_INT, 0, 4028 + rank, comm);
		MPI_Send(xSEND,    size + 1, MPI_INT, 0, 5027 + rank, comm);
		
		MPI_Send(aRECV, xRECV[size], MPI_INT, 0, 6026 + rank, comm);
		MPI_Send(aSEND, xSEND[size], MPI_INT, 0, 7025 + rank, comm);
	}

	MPI_Barrier(comm);

	return 0;
}

// ПОРЯДОК ПЕРЕЧИСЛЕНИЯ ИНДЕКСОВ В МАССИВЕ eP ПРОИЗВОЛЬНЫЙ
// ПОДПРОГРАММА МОЖЕТ РАБОТАТЬ С ЭЛЕМЕНТАМИ С ПРОИЗВОЛЬНЫМ ЧИСЛОМ ВЕРШИН
// АЛГОРИТМ: КАЖДЫЙ РАССЫЛАЕТ КАЖДОМУ СВОЮ ЧАСТЬ ТОПОЛОГИЧЕСКОГО МАССИВА ДЛЯ ЛОКАЛЬНОЙ ОБРАБОТКИ
// ИЗНАЧАЛЬНО РЕЗУЛЬТИРУЮЩИЙ МАССИВ ВЫДЕЛЯЕТСЯ ПОД МАКСИМАЛЬНЫЙ РАЗМЕР СРЕДИ ВСЕХ СЕТОЧНЫХ ЭЛЕМЕНТОВ
// В ФИНАЛЕ ПРОЦЕДУРЫ ДЛИНА МАССИВА КОРРЕКТИРУЕТСЯ ДО АКТУАЛЬНОГО ЗНАЧЕНИЯ
int hm4InitElementsTopologyMPI(int *eDIST, int *xE, int *aE, int neP, int *eP, int **xeP, int **aeP, FILE *outSTREAM, MPI_Comm comm)
{
	int rank, size;

	int *lengthMSG;
	int *bufXE, *bufAE;

	int i, j, k;
	
	  MPI_Barrier(comm);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	  MPI_Barrier(comm);

	if( (outSTREAM != NULL) && (rank == 0) )
	{ fprintf(outSTREAM, "hm4InitElementsTopologyMPI --> MAIN: START PROCEDURE IN %d MPI PROCESSES TOPOLOGY\n", size); }

	lengthMSG = (int *)malloc((size + 1) * sizeof(int));
	if(lengthMSG == NULL) { fprintf(stderr, "hm4InitElementsTopologyMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, size + 1); MPI_Abort(comm, 0); }

	{
		int lengthTOTAL = xE[eDIST[rank + 1] - eDIST[rank]];
		int eltMAX = xE[1] - xE[0];
		int buf2[2];

		for(i=1; i<eDIST[rank + 1] - eDIST[rank]; i++)
		{
			int n = xE[i + 1] - xE[i];
			if(eltMAX < n) eltMAX = n;
		} // for i

		if(rank == 0)
		{
			lengthMSG[0] = lengthTOTAL; lengthMSG[size] = eltMAX;

			for(i=1; i<size; i++)
			{
				MPI_Status status;
				MPI_Recv(buf2, 2, MPI_INT, i, 1076 + i, comm, &status);
				lengthMSG[i] = buf2[0];
				if(lengthMSG[size] < buf2[1]) lengthMSG[size] = buf2[1];
			} // for i

			MPI_Bcast(lengthMSG, size + 1, MPI_INT, 0, comm);
		}
		else
		{
			buf2[0] = lengthTOTAL; buf2[1] = eltMAX;
			MPI_Send(buf2, 2, MPI_INT, 0, 1076 + rank, comm);
			MPI_Bcast(lengthMSG, size + 1, MPI_INT, 0, comm);
		}
	}

	{
		int nMAX = eDIST[1] - eDIST[0];
		int lMAX = lengthMSG[0];

		for(i=1; i<size; i++)
		{
			int n = eDIST[i + 1] - eDIST[i];
			if(nMAX < n) nMAX = n;
			if(lMAX < lengthMSG[i]) lMAX = lengthMSG[i];
		} // for i

		bufXE = (int *)malloc((nMAX + 1) * sizeof(int));
		if(bufXE == NULL) { fprintf(stderr, "hm4InitElementsTopologyMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, nMAX + 1); MPI_Abort(comm, 0); }
		bufAE = (int *)malloc(lMAX * sizeof(int));
		if(bufAE == NULL) { fprintf(stderr, "hm4InitElementsTopologyMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, lMAX); MPI_Abort(comm, 0); }

		xeP[0] = (int *)malloc((neP + 1) * sizeof(int));
		if(xeP[0] == NULL) { fprintf(stderr, "hm4InitElementsTopologyMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, neP + 1); MPI_Abort(comm, 0); }
		for(i=0; i<=neP; i++) xeP[0][i] = 0;
		aeP[0] = (int *)malloc(neP * lengthMSG[size] * sizeof(int));
		if(aeP[0] == NULL) { fprintf(stderr, "hm4InitElementsTopologyMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, neP * lengthMSG[size]); MPI_Abort(comm, 0); }
		for(i=0; i<neP * lengthMSG[size]; i++) aeP[0][i] = -1;
	}

	for(i=0; i<size; i++)
	{
		int *ixPTR, *iaPTR;
		int n = eDIST[i + 1] - eDIST[i];

		if(rank == i)
		{
			ixPTR = xE; iaPTR = aE;
			MPI_Bcast(ixPTR, n + 1, MPI_INT, i, comm);
			MPI_Bcast(iaPTR, lengthMSG[i], MPI_INT, i, comm);
		}
		else
		{
			ixPTR = bufXE; iaPTR = bufAE;
			MPI_Bcast(ixPTR, n + 1, MPI_INT, i, comm);
			MPI_Bcast(iaPTR, lengthMSG[i], MPI_INT, i, comm);
		}

		for(j=0; j<neP; j++) if( (eP[j] >= eDIST[i]) && (eP[j] < eDIST[i + 1]) )
		{
			int offset = eP[j] - eDIST[i];
			xeP[0][j + 1] = ixPTR[offset + 1] - ixPTR[offset];
			for(k=0; k<xeP[0][j + 1]; k++)
			{
				aeP[0][lengthMSG[size]*j + k] = iaPTR[ixPTR[offset] + k];
			} // for k 
		} // if in for j
	} // for i

	for(i=1; i<=neP; i++) if(xeP[0][i] == 0)
	{ fprintf(stderr, "hm4InitElementsTopologyMPI --> %4d: *** ERROR *** Element %d global id = %d topology not inited\n", rank, i-1, eP[i-1]); MPI_Abort(comm, 0); }

	for(i=2; i<=neP; i++) xeP[0][i] += xeP[0][i - 1];
	for(j=0, i=0; i<neP*lengthMSG[size]; i++) if(aeP[0][i] != -1) aeP[0][j] = aeP[0][i], j += 1;
	if(j != xeP[0][neP]) { fprintf(stderr, "hm4InitElementsTopologyMPI --> %4d: *** ERROR *** xeP[neP] != real aeP length (%d <--> %d)\n", rank, xeP[0][neP], j); MPI_Abort(comm, 0); }

	if(xeP[0][neP] != neP * lengthMSG[size])
	{
		int *irPTR = (int *)realloc(aeP[0], xeP[0][neP] * sizeof(int));
		if(irPTR == NULL) { fprintf(stderr, "hm4InitElementsTopologyMPI --> %4d: *** ERROR *** Can't reallocate memory <%d * INT>\n", rank, xeP[0][neP]); MPI_Abort(comm, 0); }
		if(irPTR != aeP[0]) aeP[0] = irPTR;
	}

	free(lengthMSG);
	free(bufXE); free(bufAE);
	
	MPI_Barrier(comm);

	return 0;
}

// ПОРЯДОК ПЕРЕЧИСЛЕНИЯ ИНДЕКСОВ В МАССИВЕ vLIST ПО ВОЗРАСТАНИЮ
// ПОРЯДОК ПЕРЕЧИСЛЕНИЯ КООРДИНАТ В МАССИВЕ vCRD СООТВЕТСТВУЕТ ПОРЯДКУ ПЕРЕЧИСЛЕНИЯ ИНДЕКСОВ В МАССИВЕ vLIST
// АЛГОРИТМ: ЛОКАЛЬНАЯ ИНИЦИАЛИЗАЦИЯ УПОРЯДОЧЕННЫХ ПО ВОЗРАСТАНИБ ИНДЕКСОВ СПИСКОВ СЕТОЧНЫХ УЗЛОВ
// ДАЛЕЕ КАЖДЫЙ РАССЫЛАЕТ КАЖДОМУ СВОЮ ЧАСТЬ МАССИВА КООРДИНАТ ЛОКАЛЬНОЙ ОБРАБОТКИ
int hm4CreateProcessNodesListAndCoordinatesMPI(int nE, int *xE, int *aE, int *vDIST, double *CV, int *nV, int **vLIST, double **vCRD, int nDIM, FILE *outSTREAM, MPI_Comm comm)
{
	int rank, size;

	double *bufC;

	int i, j, k;
	
	  MPI_Barrier(comm);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	  MPI_Barrier(comm);

	if( (outSTREAM != NULL) && (rank == 0) )
	{ fprintf(outSTREAM, "hm4CreateProcessNodesListAndCoordinatesMPI --> MAIN: START PROCEDURE IN %d MPI PROCESSES TOPOLOGY\n", size); }

	vLIST[0] = (int *)malloc(xE[nE] * sizeof(int));
	if(vLIST[0] == NULL) { fprintf(stderr, "hm4CreateProcessNodesListAndCoordinatesMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, xE[nE]); MPI_Abort(comm, 0); }
	for(i=0; i<xE[nE]; i++) vLIST[0][i] = aE[i];
	SortAndClearIntegerListUnit(xE[nE], vLIST[0], nV);
	if(nV[0] != xE[nE])
	{
		int *iPTR = (int *)realloc(vLIST[0], nV[0] * sizeof(int));
		if(iPTR == NULL) { fprintf(stderr, "hm4CreateProcessNodesListAndCoordinatesMPI --> %4d: *** ERROR *** Can't reallocate memory <%d * INT>\n", rank, nV[0]); MPI_Abort(comm, 0); }
		if(iPTR != vLIST[0]) vLIST[0] = iPTR;
	}
	vCRD[0] = (double *)malloc(nV[0] * nDIM * sizeof(double));
	if(vCRD[0] == NULL) { fprintf(stderr, "hm4CreateProcessNodesListAndCoordinatesMPI --> %4d: *** ERROR *** Can't allocate memory <%d * DOUBLE>\n", rank, nV[0] * nDIM); MPI_Abort(comm, 0); }

	{
		int nMAX = vDIST[1] - vDIST[0];
		for(i=1; i<size; i++)
		{
			int n = vDIST[i + 1] - vDIST[i];
			if(nMAX < n) nMAX = n;
		} // for i
		bufC = (double *)malloc(nMAX * nDIM * sizeof(double));
		if(bufC == NULL) { fprintf(stderr, "hm4CreateProcessNodesListAndCoordinatesMPI --> %4d: *** ERROR *** Can't allocate memory <%d * DOUBLE>\n", rank, nMAX * nDIM); MPI_Abort(comm, 0); }
	}

	for(i=0; i<size; i++)
	{
		double *dPTR;
		int n = vDIST[i + 1] - vDIST[i];

		if(rank == i)
		{
			dPTR = CV;
			MPI_Bcast(dPTR, n * nDIM, MPI_DOUBLE, i, comm);			
		}
		else
		{
			dPTR = bufC;
			MPI_Bcast(bufC, n * nDIM, MPI_DOUBLE, i, comm);			
		}

		for(j=0; j<nV[0]; j++) if( (vLIST[0][j] >= vDIST[i]) && (vLIST[0][j] < vDIST[i + 1]) )
		{
			int offset = vLIST[0][j] - vDIST[i];
			for(k=0; k<nDIM; k++) vCRD[0][j*nDIM + k] = dPTR[offset*nDIM + k];			
		} // if in for j
	} // for i

	free(bufC);

	MPI_Barrier(comm);
	
	return 0;
}
