#define _DEBUG_OUTPUT_ 1

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "csrHmeshDualMPI.h"

// ВНУТРЕННИЕ ЛОКАЛЬНЫЕ ТИПЫ ДАННЫХ
typedef struct
{
	int idE[2];
	int idV[4];
} t_Face;
// ВНУТРЕННИЕ ЛОКАЛЬНЫЕ ТИПЫ ДАННЫХ

static int eFaceStructureInit(int eTYPE, int *ePTR, int idELT, int idFACE, t_Face *face);
static int PuzSort(int *P, int N);
static int compareFACES(t_Face face1, t_Face face2);
static void wsortFACES(t_Face *A, int l, int r);
static int FindFaceInOrderList(t_Face face, t_Face *FACES, int n);

int csrHmeshDualMPI(int *eDIST, int *xE, int *aE, int **x, int **a, FILE *outSTREAM, MPI_Comm comm)
{
	int rank, size;

	int nIF, nPF, nNF;
	t_Face *initFACE, *posFACE, *negFACE;

	int i, j, k, m;

	  MPI_Barrier(comm);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	  MPI_Barrier(comm);

	if( (outSTREAM != NULL) && (rank == 0) )
	{ fprintf(outSTREAM, "csrHmeshDualMPI --> MAIN: START PROCEDURE IN %d MPI PROCESSES TOPOLOGY\n", size); }

	// ИНИЦИАЛИЗАЦИЯ ОБЩЕГО СПИСКА ГРАНЕЙ
	// Инициализируется массив структур с параметрами всех граней всех элементов из списка данного параллельного процесса
	{
		int nE = eDIST[rank + 1] - eDIST[rank];

		// подсчет общего числа граней
		nIF = 0;
		for(i=0; i<nE; i++)
		{
			int eTYPE = xE[i + 1] - xE[i];
			if(eTYPE == 4) nIF += 4; if(eTYPE == 5) nIF += 5; if(eTYPE == 6) nIF += 5; if(eTYPE == 8) nIF += 6;
		} // for i

		// выделение памяти под инициализируемый массив
		initFACE = (t_Face *)malloc(nIF * sizeof(t_Face));
		if(initFACE == NULL) { fprintf(stderr, "csrHmeshDualMPI --> %4d: *** ERROR *** Can't allocate memory <%d * t_Face>\n", rank, nIF); MPI_Abort(comm, 0); }

		// непосредственно инициализация массива
		// цикл по сеточным элементам
		for(i=0, k=0; i<nE; i++)
		{
			int eTYPE = xE[i + 1] - xE[i];
			int *ePTR = aE + xE[i];
			int enFACE;

			if(eTYPE == 4) enFACE = 4; if(eTYPE == 5) enFACE = 5; if(eTYPE == 6) enFACE = 5; if(eTYPE == 8) enFACE = 6;
		
			// цикл по граням сеточного элемента
			for(j=0; j<enFACE; j++)
			{
				// подпрограмма инициализирует структуру с параметрами ориентированной грани
				// после определения ориентации (позитивная/негативная грань) индексы вершин грани упорядочиваются по возрастанию
				// индекс четвертой вершины треугольной грани приравнивается -1
				// идентификатор элемента - локальный idELT -> [0; nE)
				eFaceStructureInit(eTYPE, ePTR, i, j, initFACE + k);				
				k += 1;
			} // for j
		} // for i

		if(k != nIF) { fprintf(stderr, "csrHmeshDualMPI --> %4d: *** ERROR *** nIF (%d) != local nIF (%d)\n", rank, nIF, k); MPI_Abort(comm, 0); }
	}
	// ИНИЦИАЛИЗАЦИЯ ОБЩЕГО СПИСКА ГРАНЕЙ

	// СОРТИРОВКА ОБЩЕГО СПИСКА ГРАНЕЙ
	// в данном случае грани сортируются по возрастанию идентификаторов списков принадлежащих грани вершин
	// индексы сеточных элементов, которым принадлежит грань, и ориентация грани не учитывается
	wsortFACES(initFACE, 0, nIF - 1);	
	// СОРТИРОВКА ОБЩЕГО СПИСКА ГРАНЕЙ

	// ДЕЛЕНИЕ ОБЩЕГО СПИСКА НА ИНИЦИАЛИЗИРОВАННЫЕ, ПОЗИТИВНЫЕ И НЕГАТИВНЫЕ ГРАНИ
	{
		int nIFDUPL = nIF;
		nIF = nPF = nNF = 0;

		for(i=0; i<nIFDUPL; )
		{
			if(i != nIFDUPL - 1)
			{
				// списки индексов вершин двух граней эквивалентны - инициализированная грань
				// грань с индексом i полностью инициализируется (инициализируются идексы двух сеточных элементов)
				// грань с индексом i+1 обнуляется
				if(0 == compareFACES(initFACE[i], initFACE[i + 1]))
				{
					if(    ( (initFACE[i].idE[0] >= 0) && (initFACE[i + 1].idE[0] >= 0) )
						||
						   ( (initFACE[i].idE[1] >= 0) && (initFACE[i + 1].idE[1] >= 0) )
					  )
					{ fprintf(stderr, "csrHmeshDualMPI --> %4d: *** ERROR *** initFACE 2 equal oriented faces\n", rank); MPI_Abort(comm, 0); }

					if(initFACE[i].idE[0] < 0) initFACE[i].idE[0] = initFACE[i + 1].idE[0];
					else                       initFACE[i].idE[1] = initFACE[i + 1].idE[1];

					initFACE[i + 1].idV[0] = -1;
					initFACE[nIF] = initFACE[i];
					nIF += 1; i += 2;
				}
				else
				// грань не имеет пары
				// такая грань либо принадлежит границе расчетной области, либо ее пара в списке другого параллельного процесса
				{
					if(initFACE[i].idE[0] >= 0) nPF += 1;
					else                        nNF += 1;
					initFACE[nIF] = initFACE[i];
					nIF += 1; i += 1;
				}
			}
			else
			// анализ последней грани общего списка
			// проверка наличия пары в данном случае не производится
			{
				if(initFACE[i].idE[0] >= 0) nPF += 1;
				else                        nNF += 1;
				initFACE[nIF] = initFACE[i];
				nIF += 1; i += 1;
			}
		} // for i

		if(nIF*2 - nPF - nNF != nIFDUPL)
		{ fprintf(stderr, "csrHmeshDualMPI --> %4d: *** ERROR *** nIF*2 - nPF - nNF (%d) != nIFDUPL (%d)\n", rank, nIF*2 - nPF - nNF, nIFDUPL); MPI_Abort(comm, 0); }

		// выделение памяти под отдельный массив оставшихся без пары позитивно ориентированных граней общего списка
		if(nPF > 0)
		{
			posFACE = (t_Face *)malloc(nPF * sizeof(t_Face));
			if(posFACE == NULL) { fprintf(stderr, "csrHmeshDualMPI --> %4d: *** ERROR *** Can't allocate memory <%d * t_Face>\n", rank, nPF); MPI_Abort(comm, 0); }
		}
		// выделение памяти под отдельный массив оставшихся без пары негативно ориентированных граней общего списка
		if(nNF > 0)
		{
			negFACE = (t_Face *)malloc(nNF * sizeof(t_Face));
			if(negFACE == NULL) { fprintf(stderr, "csrHmeshDualMPI --> %4d: *** ERROR *** Can't allocate memory <%d * t_Face>\n", rank, nNF); MPI_Abort(comm, 0); }
		}

		// инициализация массивов парных, позитивных и негативных гарней
		for(i=0, j=0, k=0, m=0; i<nIF; i++)
		{
			// парная (инициализированная грань)
			if( (initFACE[i].idE[0] >= 0) && (initFACE[i].idE[1] >= 0) )
			{
				initFACE[j] = initFACE[i];
				j += 1;
			}
			else
			{
				// позитивно ориентированная грань
				if(initFACE[i].idE[0] >= 0)
				{
					posFACE[k] = initFACE[i];
					k += 1;
				}
				// отрицательно (негативно) ориентированная грань
				else
				{
					negFACE[m] = initFACE[i];
					m += 1;
				}
			}
		} // for i

		if( (k != nPF) || (m != nNF) || (j != nIF - k - m) )
		{ fprintf(stderr, "csrHmeshDualMPI --> %4d: *** ERROR *** nIF / nPF / nNF\n", rank); MPI_Abort(comm, 0); }

		// корректировка длины массива initFACE
		// число элементов массива сокращается за счет удаления повторений и перемещения непарных элементов в отдельные списки
		// в зависимости от ориентации
		{
			t_Face *PTR = NULL;
			
			nIF = j;

			// ??? следующие условие по логике алгоритма выполняется всегда ???
			if(nIF != nIFDUPL)
			{
				if(nIF != 0)
				{
					PTR = (t_Face *)realloc(initFACE, nIF * sizeof(t_Face));
					if(PTR == NULL) { fprintf(stderr, "csrHmeshDualMPI --> %4d: *** ERROR *** Can't reallocate memory <%d * t_Face>\n", rank, nIF); MPI_Abort(comm, 0); }
					if(PTR != initFACE) initFACE = PTR; 
				}
				else
				{
					free(initFACE); initFACE = NULL;
				}
			}			
		}
		// корректировка длины массива initFACE

#if(_DEBUG_OUTPUT_ == 1)
		if(outSTREAM != NULL)
		{ fprintf(outSTREAM, "csrHmeshDualMPI --> %4d: Faces lists before nIF = %8d nPF = %8d nNF = %8d\n", rank, nIF, nPF, nNF); }
#endif
	}
	// ДЕЛЕНИЕ ОБЩЕГО СПИСКА НА ИНИЦИАЛИЗИРОВАННЫЕ, ПОЗИТИВНЫЕ И НЕГАТИВНЫЕ ГРАНИ

	MPI_Barrier(comm);

	// ПАРАЛЛЕЛЬНАЯ ЧАСТЬ
	// позитивные грани остаются в памяти параллельного процесса
	// негативные грани совершают круговой обход по всем параллельным процессам и возвращаются назад
	// принятые на текущем шаге негативные грани другого процесса сравниваются с собственными позитивными гранями
	// если для негативной грани находится позитивная пара, то обе грани инициализируются
	// таким образом, по окончании процедуры (возвращении назад негативных граней) должны быть полностью инициализированны
	// все внутренние сеточные грани
	// не инициализированная грань принадлеждит границе расчетной области
	{
		t_Face *faceRECV, *faceSEND;
		int *negNUMBERS = (int *)malloc(size * sizeof(int));
		int lRECV, lSEND, idRECV;
		int rankRECV, rankSEND;

		// инициализация массива с числами негативных граней параллельных процессов
		// данный массив необходим для приема/передачи данных без запроса длин массивов и выделения
		// памяти под соответствующие буфера обмена
		if(negNUMBERS == NULL) { fprintf(stderr, "csrHmeshDualMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, size); MPI_Abort(comm, 0); }
		if(rank == 0)
		{
			negNUMBERS[0] = nNF;			
			for(i=1; i<size; i++)
			{
				MPI_Status status;
				MPI_Recv(&negNUMBERS[i], 1, MPI_INT, i, 1000 + i, comm, &status);
			} // for i
			MPI_Bcast(negNUMBERS, size, MPI_INT, 0, comm);
		}
		else
		{
			MPI_Send(&nNF, 1, MPI_INT, 0, 1000 + rank, comm);
			MPI_Bcast(negNUMBERS, size, MPI_INT, 0, comm);
		}

		// выделение памяти под буфера приема/передачи данных и инициализация соответствующих указателей
		{
			// определение максимальной длины массива отрицательно ориентированных граней среди всех параллельных процессов
			int nMAX = negNUMBERS[0];
			for(i=1; i<size; i++) if(nMAX < negNUMBERS[i]) nMAX = negNUMBERS[i];

			// выделение памяти под буфера приема/передачи данных
			// указатель на массив, выделенный под буфера совпадает с указателем на массив с параметрами негративных граней
			
			// выделение памяти, когда у параллельного процесса нет отрицательно ориентированных граней 
			if(negFACE == NULL)
			{
				negFACE = (t_Face *)malloc(nMAX * 2 * sizeof(t_Face));
				if(negFACE == NULL) { fprintf(stderr, "csrHmeshDualMPI --> %4d: *** ERROR *** Can't allocate memory <%d * t_Face>\n", rank, nMAX*2); MPI_Abort(comm, 0); }
			}
			else
			// выделение памяти, когда у параллельного процесса есть отрицательно ориентированные грани 
			{
				// перевод индексов элементов в списке негативных граней в глобальную нумерацию idE[1] -> [eDIST[rank]; eDIST[rank + 1])
				for(i=0; i<nNF; i++) negFACE[i].idE[1] += eDIST[rank];

				// корректировка длины массива negFACE
				{
					t_Face *PTR = (t_Face *)realloc(negFACE, nMAX * 2 * sizeof(t_Face));
					if(PTR == NULL) { fprintf(stderr, "csrHmeshDualMPI --> %4d: *** ERROR *** Can't reallocate memory <%d * t_Face>\n", rank, nMAX*2); MPI_Abort(comm, 0); }
					if(PTR != negFACE) negFACE = PTR; 
				}
			}

			// инициализация указателей на буфера приема и передачи данных
			faceRECV = negFACE; faceSEND = faceRECV + nMAX;

			MPI_Barrier(comm);
#if(_DEBUG_OUTPUT_ == 1)
			if( (rank == 0) && (outSTREAM != NULL) )
			{ fprintf(outSTREAM, "csrHmeshDualMPI --> MAIN: Negative faces buffer size = %d elements\n", nMAX); }
#endif
			MPI_Barrier(comm);
		}

		// инициализация номеров параллельных процессов "слева" (прием данных) и "справа" (передача данных)
		rankRECV = rank - 1; if(rankRECV ==   -1) rankRECV = size - 1;
		rankSEND = rank + 1; if(rankSEND == size) rankSEND =        0;

		// начальная инициализация номера параллельного процесса и длины соответствующего массива как собственных параметров 
		lRECV = nNF; idRECV = rank;

		// цикл по приему/обработке/передаче данных
		for(i=0; i<size; i++)
		{
			MPI_Request reqRECV, reqSEND;
			MPI_Status statRECV, statSEND;

			// "рокировка" массивов приема/передачи данных
			{ t_Face *TMP = faceSEND; faceSEND = faceRECV; faceRECV = TMP; lSEND = lRECV; }

			// определение идентификатора того процесса, чей массив потребуется обработать на данном шаге и длины этого массива
			idRECV -= 1; if(idRECV < 0) idRECV = size - 1;
			lRECV = negNUMBERS[idRECV];

			// запуск асинхронного приема и передачи данных
			if(lSEND > 0) MPI_Isend(faceSEND, lSEND * sizeof(t_Face), MPI_CHAR, rankSEND, 3000, comm, &reqSEND);
			if(lRECV > 0) MPI_Irecv(faceRECV, lRECV * sizeof(t_Face), MPI_CHAR, rankRECV, 3000, comm, &reqRECV);

			// ожидание завершения асинхронного обмена данными
			if(lSEND > 0) { MPI_Wait(&reqSEND, &statSEND); if(0) MPI_Request_free(&reqSEND); }
			if(lRECV > 0) { MPI_Wait(&reqRECV, &statRECV); if(0) MPI_Request_free(&reqRECV); }

			// обработка массива - поиск совпадений с собственным массивом позитивных граней
			// на последнем витке цикла не выполняется, поскольку здесь назад возвращается собственный массив негативных граней
			if( (i != size - 1) && (lRECV > 0) && (nPF > 0) )
			{
				// попытка поиска совпадений для еще не инициализированных граней из списка негативных граней другого процесса
				for(j=0; j<lRECV; j++) if(faceRECV[j].idE[0] == -1)
				{
					// поиск совпадения граней
					int get = FindFaceInOrderList(faceRECV[j], posFACE, nPF);

					// совпадение найдено
					if(get != -1)
					{
						// в негативную грань дописывается глобальный индекс элемента собственной позитивной грани
						faceRECV[j].idE[0] = eDIST[rank] + posFACE[get].idE[0];
						// в собственную позитивную грань дописывается глобальный индекс элемента чужой негативной грани
						posFACE[get].idE[1] = faceRECV[j].idE[1];
					}
				} // for j				
			} // if		
		} // for i

		if(nNF != lRECV)
		{ fprintf(stderr, "csrHmeshDualMPI --> %4d: *** ERROR *** nNF != lRECV\n", rank); MPI_Abort(comm, 0); }

		// копирование данных по нужному указателю (negFACE)
		if(faceRECV != negFACE) for(i=0; i<nNF; i++) negFACE[i] = faceRECV[i];

		// подсчет числа инициализированных собственных негативных граней со смещением данных (удаляются грани, принадлежащие границе расчетной области)
		// индекс "собственного" элемента процесса в ненативной грани возвращается на локальный числовой интервал idE[1] -> [0; Ne) 
		j=0; for(i=0; i<nNF; i++) if(negFACE[i].idE[0] >= 0) { negFACE[j] = negFACE[i]; negFACE[j].idE[1] -= eDIST[rank]; j += 1; } nNF = j;
		// подсчет числа инициализированных собственных позитивных граней
		j=0; for(i=0; i<nPF; i++) if(posFACE[i].idE[1] >= 0) { posFACE[j] = posFACE[i];                                   j += 1; } nPF = j;		

		// освобождение памяти, выделенной под массив с длинами негативных списков параллельных процессов
		free(negNUMBERS);

#if(_DEBUG_OUTPUT_ == 1)
		if(outSTREAM != NULL)
		{ fprintf(outSTREAM, "csrHmeshDualMPI --> %4d: Faces lists  after nIF = %8d nPF = %8d nNF = %8d\n", rank, nIF, nPF, nNF); }
#endif
	}
	// ПАРАЛЛЕЛЬНАЯ ЧАСТЬ

	MPI_Barrier(comm);

	// ИНИЦИАЛИЗАЦИЯ ИСКОМЫХ МАССИВОВ
	{
		int nE = eDIST[rank + 1] - eDIST[rank];
		int nF = nIF*2 + nNF + nPF;

		x[0] = (int *)malloc((nE + 1) * sizeof(int));
		if(x[0] == NULL) { fprintf(stderr, "csrHmeshDualMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, nE + 1); MPI_Abort(comm, 0); }
		a[0] = (int *)malloc( nF      * sizeof(int));
		if(a[0] == NULL) { fprintf(stderr, "csrHmeshDualMPI --> %4d: *** ERROR *** Can't allocate memory <%d * INT>\n", rank, nF    ); MPI_Abort(comm, 0); }

		for(i=0; i<=nE; i++) x[0][i] =  0;
		for(i=0; i< nF; i++) a[0][i] = -1;

		for(i=0; i<nIF; i++) x[0][initFACE[i].idE[0] + 1] += 1, x[0][initFACE[i].idE[1] + 1] += 1;
		for(i=0; i<nPF; i++) x[0][ posFACE[i].idE[0] + 1] += 1;
		for(i=0; i<nNF; i++) x[0][ negFACE[i].idE[1] + 1] += 1;
		for(i=2; i<=nE; i++) x[0][i] += x[0][i-1];

		if(x[0][nE] != nF)
		{ fprintf(stderr, "csrHmeshDualMPI --> %4d: *** ERROR *** nF != x[nE]\n", rank); MPI_Abort(comm, 0); }

#if(_DEBUG_OUTPUT_ == 1)
		if(outSTREAM != NULL)
		{ fprintf(outSTREAM, "csrHmeshDualMPI --> %4d: Local graph edges number = %8d\n", rank, nF); }
#endif

		for(i=0; i<nIF; i++)
		{
			a[0][x[0][initFACE[i].idE[0]]] = eDIST[rank] + initFACE[i].idE[1]; x[0][initFACE[i].idE[0]] += 1;
			a[0][x[0][initFACE[i].idE[1]]] = eDIST[rank] + initFACE[i].idE[0]; x[0][initFACE[i].idE[1]] += 1;
		} // for i
		free(initFACE);

		for(i=0; i<nPF; i++)
		{
			a[0][x[0][posFACE[i].idE[0]]] = posFACE[i].idE[1]; x[0][posFACE[i].idE[0]] += 1;
		} // for i
		free(posFACE);

		for(i=0; i<nNF; i++)
		{
			a[0][x[0][negFACE[i].idE[1]]] = negFACE[i].idE[0]; x[0][negFACE[i].idE[1]] += 1;
		} // for i
		free(negFACE);

		for(i=nE; i>0; i--) x[0][i] = x[0][i-1]; x[0][0] = 0;

		for(i=0; i<nE; i++) PuzSort(&a[0][x[0][i]], x[0][i+1] - x[0][i]);
	}
	// ИНИЦИАЛИЗАЦИЯ ИСКОМЫХ МАССИВОВ

	MPI_Barrier(comm);

	return 0;
}


// ===========================================================================================================
// инициализация структуры с параметрами грани элемента
static int eFaceStructureInit(int eTYPE, int *ePTR, int idELT, int idFACE, t_Face *face)
{
	int n = 4;
	int init = 0;

	face[0].idV[3] = -1;
	face[0].idE[0] = face[0].idE[1] = -1;

	if(eTYPE == 4)
	{
		if(idFACE == 0) { face[0].idV[0] = ePTR[1], face[0].idV[1] = ePTR[0], face[0].idV[2] = ePTR[2]; init = 1; }
		if(idFACE == 1) { face[0].idV[0] = ePTR[0], face[0].idV[1] = ePTR[1], face[0].idV[2] = ePTR[3]; init = 1; }
		if(idFACE == 2) { face[0].idV[0] = ePTR[1], face[0].idV[1] = ePTR[2], face[0].idV[2] = ePTR[3]; init = 1; }
		if(idFACE == 3) { face[0].idV[0] = ePTR[2], face[0].idV[1] = ePTR[0], face[0].idV[2] = ePTR[3]; init = 1; }
	} // тетраэдр

	if(eTYPE == 5)
	{
		if(idFACE == 0) { face[0].idV[0] = ePTR[0], face[0].idV[1] = ePTR[2], face[0].idV[2] = ePTR[3], face[0].idV[3] = ePTR[1]; init = 1; }
		if(idFACE == 1) { face[0].idV[0] = ePTR[0], face[0].idV[1] = ePTR[1], face[0].idV[2] = ePTR[4];                           init = 1; }
		if(idFACE == 2) { face[0].idV[0] = ePTR[1], face[0].idV[1] = ePTR[3], face[0].idV[2] = ePTR[4];                           init = 1; }
		if(idFACE == 3) { face[0].idV[0] = ePTR[3], face[0].idV[1] = ePTR[2], face[0].idV[2] = ePTR[4];                           init = 1; }
		if(idFACE == 4) { face[0].idV[0] = ePTR[2], face[0].idV[1] = ePTR[0], face[0].idV[2] = ePTR[4];                           init = 1; }						   
	} // пирамида

	if(eTYPE == 6)
	{
		if(idFACE == 0) { face[0].idV[0] = ePTR[0], face[0].idV[1] = ePTR[1], face[0].idV[2] = ePTR[4], face[0].idV[3] = ePTR[3]; init = 1; }
		if(idFACE == 1) { face[0].idV[0] = ePTR[1], face[0].idV[1] = ePTR[2], face[0].idV[2] = ePTR[5], face[0].idV[3] = ePTR[4]; init = 1; }
		if(idFACE == 2) { face[0].idV[0] = ePTR[2], face[0].idV[1] = ePTR[0], face[0].idV[2] = ePTR[3], face[0].idV[3] = ePTR[5]; init = 1; }
		if(idFACE == 3) { face[0].idV[0] = ePTR[0], face[0].idV[1] = ePTR[2], face[0].idV[2] = ePTR[1];                           init = 1; }
		if(idFACE == 4) { face[0].idV[0] = ePTR[3], face[0].idV[1] = ePTR[4], face[0].idV[2] = ePTR[5];                           init = 1; }						   
	} // призма

	if(eTYPE == 8)
	{
		if(idFACE == 0) { face[0].idV[0] = ePTR[0], face[0].idV[1] = ePTR[1], face[0].idV[2] = ePTR[5], face[0].idV[3] = ePTR[4]; init = 1; }
		if(idFACE == 1) { face[0].idV[0] = ePTR[1], face[0].idV[1] = ePTR[3], face[0].idV[2] = ePTR[7], face[0].idV[3] = ePTR[5]; init = 1; }
		if(idFACE == 2) { face[0].idV[0] = ePTR[3], face[0].idV[1] = ePTR[2], face[0].idV[2] = ePTR[6], face[0].idV[3] = ePTR[7]; init = 1; }
		if(idFACE == 3) { face[0].idV[0] = ePTR[2], face[0].idV[1] = ePTR[0], face[0].idV[2] = ePTR[4], face[0].idV[3] = ePTR[6]; init = 1; }
		if(idFACE == 4) { face[0].idV[0] = ePTR[1], face[0].idV[1] = ePTR[0], face[0].idV[2] = ePTR[2], face[0].idV[3] = ePTR[3]; init = 1; }
		if(idFACE == 5) { face[0].idV[0] = ePTR[4], face[0].idV[1] = ePTR[5], face[0].idV[2] = ePTR[7], face[0].idV[3] = ePTR[6]; init = 1; }						   
	} // гексаэдр

	if(init == 0) { fprintf(stderr, "eFaceStructureInit --> LOCAL: *** ERROR *** Can't create vertex list for face %d of %d-node element \n", idFACE, eTYPE); exit(0); }

	if(face[0].idV[3] == -1) n = 3;

	// инициализация ориентации: позитивная / негативная
	{
		int idMIN = 0;
		int i;

		for(i=1; i<n; i++) if(face[0].idV[idMIN] > face[0].idV[i]) idMIN = i;

		if(face[0].idV[(idMIN + 1) % n] < face[0].idV[(idMIN + n - 1) % n]) face[0].idE[0] = idELT;
		else                                                                face[0].idE[1] = idELT;
	}
	// инициализация ориентации: позитивная / негативная

	// сортировка целочисленных индексов вершин грани по возрастанию
	PuzSort(face[0].idV, n);

	return 0;
}
// ===========================================================================================================

// ===========================================================================================================
// пузырьковая сортировка целочисленного массива по возрастанию
static int PuzSort(int *P, int N)
{
	int i, j;
	
	for(i=0; i<N-1; i++)
		for(j=0; j<N-1-i; j++)
			if(P[j] > P[j+1])
			{
				int tmp;
				tmp = P[j+1];
				P[j+1] = P[j];
				P[j] = tmp;
			}

	return 0;
}
// ===========================================================================================================

// ===========================================================================================================
// сравнение двух граней для реализация процедур поиска и сортировки списка граней по возрастанию 
static int compareFACES(t_Face face1, t_Face face2)
{
	int i;

	for(i=0; i<4; i++)
	{
		if(face1.idV[i] > face2.idV[i]) return  1;
		if(face1.idV[i] < face2.idV[i]) return -1;
	} // for i

	return 0;
}
// ===========================================================================================================

// ===========================================================================================================
// сортировка структур с параметрами граней элементов
static void wsortFACES(t_Face *A, int l, int r)
{
	t_Face tmp;
	t_Face B = A[(l + r) / 2];
	int i = l, j = r;

	while( i <= j )
	{
		while( compareFACES(A[i], B) == (-1) ) i++;
		while( compareFACES(A[j], B) ==   1  ) j--;
		if( i <= j )
		{
			tmp = A[i];
			A[i] = A[j];
			A[j] = tmp;
			i++;
			j--;      
		} // if   
	} // while
  
	if( l < j ) wsortFACES(A, l, j);
	if( i < r ) wsortFACES(A, i, r);
}
// ===========================================================================================================

// ===========================================================================================================
static int FindFaceInOrderList(t_Face face, t_Face *FACES, int n)
{
	int L = 0;
	int R = n - 1;

	// проверка выхода значения за границы отсортированного массива
	if(-1 == compareFACES(face, FACES[L])) return -1;
	if( 1 == compareFACES(face, FACES[R])) return -1;
		
	// совпадение добавляемого элемента с первым или последним элементом списка
	if( 0 == compareFACES(face, FACES[L])) return L;
	if( 0 == compareFACES(face, FACES[R])) return R;
	
	while(R - L > 1)
	{
		int C = (R + L) / 2;
		int get = compareFACES(face, FACES[C]);
		if(0 == get) return C;
		if(get == 1) L = C; else R = C;
	} // while	

	return -1;
}
// ===========================================================================================================
