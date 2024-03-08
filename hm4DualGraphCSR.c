#include <stdio.h>
#include <stdlib.h>

#include "hm4DualGraphCSR.h"

typedef struct
{
	int iV;
	int E[4];
} t_Side;

#define _EEOFFS_ 24
#define _EeOFFS_ 4
#define _EtOFFS_ 6

static int _ELid_[9] = {-1, -1, -1, -1, 0, 1, 2, -1, 3};

static int _EdgeN_[4] = {4, 5, 5, 6};

static int _EdgeL_[24] = {3, 3, 3, 3, 0, 0,
                          4, 3, 3, 3, 3, 0,
						  4, 4, 4, 3, 3, 0,
						  4, 4, 4, 4, 4, 4};

static int _EdgeV_[96] = { 1,  0,  2, -1,
                           0,  1,  3, -1,
						   1,  2,  3, -1,
						   2,  0,  3, -1,
						  -1, -1, -1, -1,
						  -1, -1, -1, -1,

						   0,  2,  3,  1,
						   0,  1,  4, -1,
						   1,  3,  4, -1,
						   3,  2,  4, -1,
						   2,  0,  4, -1,
						  -1, -1, -1, -1,

						   0,  1,  4,  3,
						   1,  2,  5,  4,
						   2,  0,  3,  5,
						   0,  2,  1, -1,
						   3,  4,  5, -1,
						  -1, -1, -1, -1,

						   0,  1,  5,  4,
						   1,  3,  7,  5,
						   3,  2,  6,  7,
						   2,  0,  4,  6,
						   1,  0,  2,  3, 
						   4,  5,  7,  6};

static int PuzSort(int *P, int N);
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

// 0 --> a == b
// 1 --> a > b
// (-1) --> a < b
static int compare_tS(t_Side a, t_Side b);
static int compare_tS(t_Side a, t_Side b)
{
	if(a.E[0] > b.E[0]) return 1;
	if(a.E[0] < b.E[0]) return (-1);

	if(a.E[1] > b.E[1]) return 1;
	if(a.E[1] < b.E[1]) return (-1);

	if(a.E[2] > b.E[2]) return 1;
	if(a.E[2] < b.E[2]) return (-1);

	if(a.E[3] > b.E[3]) return 1;
	if(a.E[3] < b.E[3]) return (-1);

	return 0;
}

static void wsort_tS(t_Side *A, int l, int r);
static void wsort_tS(t_Side *A, int l, int r)
{
  t_Side tmp;
  t_Side B = A[(l + r) / 2];
  int i = l, j = r;

  while( i <= j )
  {
 	while( compare_tS(A[i], B) == (-1) ) i++;
    while( compare_tS(A[j], B) == 1 ) j--;
    if( i <= j )
    {
      tmp = A[i];
      A[i] = A[j];
      A[j] = tmp;
      i++;
      j--;      
    }   
  }
  
  if( l < j ) wsort_tS(A, l, j);
  if( i < r ) wsort_tS(A, i, r);
}

int hm4MeshDualGraphCSRSerial(int Ne, int *EX, int *EA, int **X, int **A)
{
	int Ne4 = 0; // число тетраэдров
	int Ne5 = 0; // число пирамид
	int Ne6 = 0; // число призм
	int Ne8 = 0; // число шестигранников
	int NGT;
	t_Side *tS;
	int *x, *a;

	int i;

	{
		int Nall;	

		for(i=0; i<Ne; i++)
		{
			int Nl = EX[i+1] - EX[i];
			if(Nl == 4) Ne4 += 1;
			if(Nl == 5) Ne5 += 1;
			if(Nl == 6) Ne6 += 1;
			if(Nl == 8) Ne8 += 1;
		} // for i

		Nall = Ne4 + Ne5 + Ne6 + Ne8;
		if(Nall != Ne) { fprintf(stderr, "ERROR 1\n"); exit(0); }

#if 0
		printf("DISTRIBUTION OF MESH ELEMENTS BY TYPE\n");
		printf("         TETRAHEDRON: %15d\n", Ne4);
		printf("             PYRAMID: %15d\n", Ne5);
		printf("               WEDGE: %15d\n", Ne6);
		printf("               BRICK: %15d\n", Ne8);		
		printf("\n");
#endif
	}

	{
		int Nloc = 0;

		NGT = Ne4*4 + Ne5*5 + Ne6*5 + Ne8*6;
		// printf("NGT = %d\n", NGT);

		tS = NULL; tS = (t_Side *)malloc(NGT * sizeof(t_Side));
		if(tS == NULL) { fprintf(stderr, "Can't allocate memory\n"); exit(0); }
		for(i=0; i<NGT; i++) { tS[i].iV = -1; tS[i].E[0] = tS[i].E[1] = tS[i].E[2] = tS[i].E[3] = -1; }
		// printf("Allocate EDGES\n");

		for(i=0; i<Ne; i++)
		{
			int ETYPE, *EPTR, ENUM;
			int j;

			EPTR = EA + EX[i];
			ETYPE = _ELid_[EX[i+1] - EX[i]];
			ENUM = _EdgeN_[ETYPE];

			for(j=0; j<ENUM; j++)
			{
				int k;
				int ELENGTH = _EdgeL_[ETYPE*_EtOFFS_ + j];
			
				tS[Nloc].iV = i;
				for(k=0; k<ELENGTH; k++) tS[Nloc].E[k] = EPTR[_EdgeV_[ETYPE*_EEOFFS_ + j*_EeOFFS_ + k]];
				Nloc += 1;			
			} // for j
		} // for i
		if(Nloc != NGT) { fprintf(stderr, "ERROR 2\n"); exit(0); }

		for(i=0; i<NGT; i++) 
		{
			int Nl;
			if(tS[i].E[3] >= 0) Nl = 4; else Nl = 3;
			PuzSort(tS[i].E, Nl);
		}
		// printf("Create all sides array\n");
	}

	{
		int Nloc = 0;

		wsort_tS(tS, 0, NGT-1);
		// printf("Sort\n");

		for(i=0; i<NGT; i++)
		{
			int get;

			if(i == NGT-1)
			{
				tS[i].iV = -1;
				break;
			}

			get = compare_tS(tS[i], tS[i+1]);
			if(get != 0)
			{
				tS[i].iV = -1;
			}
			else
			{
				i += 1;
			}
		} // for i
		// printf("Compare\n");

		for(i=0; i<NGT; i++)
		{
			if(tS[i].iV != -1)
			{
				tS[Nloc] = tS[i];
				Nloc += 1;
			}
		} // for i
		NGT = Nloc;
		if(Nloc%2 != 0) { fprintf(stderr, "ERROR 3\n"); exit(0); }
		// printf("Side number = %d\n", NGT/2);
	}

	{
		x = NULL; x = (int *)malloc((Ne + 1) * sizeof(int));
		if(x == NULL) { fprintf(stderr, "Can't allocate memory\n"); exit(0); }
		for(i=0; i<=Ne; i++) x[i] = 0;

		a = NULL; a = (int *)malloc(NGT * sizeof(int));
		if(a == NULL) { fprintf(stderr, "Can't allocate memory\n"); exit(0); }
		// printf("Allocate result memory\n");

		for(i=0; i<NGT; i++) x[tS[i].iV + 1] += 1;
		for(i=2; i<=Ne; i++) x[i] += x[i-1];
		if(x[Ne] != NGT) { fprintf(stderr, "ERROR 4\n"); exit(0); }

		for(i=0; i<NGT; i+=2)
		{
			int iL, iR;

			iL = tS[i  ].iV;
			iR = tS[i+1].iV;

			a[x[iL]] = iR; x[iL] += 1;
			a[x[iR]] = iL; x[iR] += 1;
		} // for i
		if(x[Ne] != x[Ne-1]) { fprintf(stderr, "ERROR 5\n"); exit(0); }
		for(i=Ne-1; i>0; i--) x[i] = x[i-1]; x[0] = 0;
		// printf("Create link list\n");

		for(i=0; i<Ne; i++)
		{
			int Nl = x[i+1] - x[i];
			int *ptr = a + x[i];
			
			PuzSort(ptr, Nl);
		} // for i
		// printf("Final sort\n");
	}

	free(tS);

	{
		int NN[7] = {0, 0, 0, 0, 0, 0, 0};

		for(i=0; i<Ne; i++)
		{
			int Nl = x[i+1] - x[i];
			NN[Nl] += 1;			
		} // for i

		// for(i=0; i<7; i++) printf("[%d]: %d\n", i, NN[i]);		
	}

	(*X) = x;
	(*A) = a;

	return 0;
}
