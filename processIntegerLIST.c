#define _INTEGER_ELEMENT_MAX_SIZE_ 5

#include <stdio.h>
#include <stdlib.h>

#include "processIntegerLIST.h"

// ÏÎÈÑÊ ÌÅÒÎÄÎÌ ÄÈÕÎÒÎÌÈÈ ÝËÅÌÅÍÒÀ TARGET Â ÎÒÑÎÐÒÈÐÎÂÀÍÍÎÌ ÏÎ
// ÂÎÇÐÀÑÒÀÍÈÞ ÌÀÑÑÈÂÅ ÖÅËÛÕ ×ÈÑÅË PTR[LENGTH]
// ÏÎÄÏÐÎÃÐÀÌÌÀ ÂÎÇÂÐÀÙÀÅÒ ÈÍÄÅÊÑ ÑÎÎÒÂÅÒÑÒÂÓÞÙÅÃÎ ÝËÅÌÅÍÒÀ, ÅÑËÈ ÎÍ ÍÀÉÄÅÍ,
// È ÇÍÀ×ÅÍÈÅÌ -1, ÅÑËÈ ÝËÅÌÅÍÒ TARGET Â ÌÀÑÑÈÂÅ ÎÒÑÓÒÑÒÂÓÅÒ
int DichotomySearchElementInOrderIntegerListUnit(int *PTR, int LENGTH, int TARGET)
{
	int L = 0;
	int R = LENGTH - 1;
	
	if(PTR[L] == TARGET) return L;
	if(PTR[R] == TARGET) return R;

	while(R - L > 1)
	{
		int C = (R + L) / 2;

		if(PTR[C] == TARGET) return C;

		if(PTR[C] < TARGET) L = C;
		else R = C;
	}

	return -1;
}

int DichotomySearchElementInOrderIntegerListSizeUnit(int *PTR, int LENGTH, int TARGET, int size)
{
	int L = 0;
	int R = LENGTH - 1;
	
	if(PTR[L * size] == TARGET) return L;
	if(PTR[R * size] == TARGET) return R;

	while(R - L > 1)
	{
		int C = (R + L) / 2;

		if(PTR[C * size] == TARGET) return C;

		if(PTR[C * size]  < TARGET) L = C;
		else R = C;
	}

	return -1;
}

int DichotomySearchElementInOrderIntegerListSizeFunc(int *PTR, int LENGTH, int *TARGET, int size, int (*cmp)(int *, int *))
{
	int L = 0;
	int R = LENGTH - 1;
	
	if((*cmp)(PTR + L*size, TARGET) == 0) return L;
	if((*cmp)(PTR + R*size, TARGET) == 0) return R;

	while(R - L > 1)
	{
		int C = (R + L) / 2;

		if((*cmp)(PTR + C*size, TARGET) ==  0) return C;

		if((*cmp)(PTR + C*size, TARGET) == -1) L = C;
		                                  else R = C;
	}

	return -1;
}

void wsortIntegerListUnit(int *A, int l, int r)
{
	int tmp;
	int B = A[(l + r) / 2];
	int i = l, j = r;

	while( i <= j )
	{
		while( A[i] < B ) i++;
		while( A[j] > B ) j--;
		if( i <= j )
		{
			tmp = A[i];
			A[i] = A[j];
			A[j] = tmp;
			i++;
			j--;      
		} // if   
	} // while
  
	if( l < j ) wsortIntegerListUnit(A, l, j);
	if( i < r ) wsortIntegerListUnit(A, i, r);
}

void wsortIntegerListSizeUnit(int *A, int l, int r, int size)
{
	int B[_INTEGER_ELEMENT_MAX_SIZE_];
	int i = l, j = r, k = 0, id = (l + r) / 2;

	for(k=0; k<size; k++) B[k] = A[id*size + k];
	
	while( i <= j )
	{
		while( A[i*size] < B[0] ) i++;
		while( A[j*size] > B[0] ) j--;
		if( i <= j )
		{
			int k;
			for(k=0; k<size; k++)
			{
				int tmp = A[i*size + k];
				A[i*size + k] = A[j*size + k];
				A[j*size + k] = tmp;
			} // for k
			i++;
			j--;      
		} // if   
	} // while
  
	if( l < j ) wsortIntegerListSizeUnit(A, l, j, size);
	if( i < r ) wsortIntegerListSizeUnit(A, i, r, size);
}

void wsortIntegerListSizeFunc(int *A, int l, int r, int size, int (*cmp)(int *, int *))
{
	int B[_INTEGER_ELEMENT_MAX_SIZE_];
	int i = l, j = r, k = 0, id = (l + r) / 2;

	for(k=0; k<size; k++) B[k] = A[id*size + k];
	
	while( i <= j )
	{
		while( (*cmp)(A + i*size, B) == -1) i++;
		while( (*cmp)(A + j*size, B) ==  1) j--;
		if( i <= j )
		{
			for(k=0; k<size; k++)
			{
				int tmp = A[i*size + k];
				A[i*size + k] = A[j*size + k];
				A[j*size + k] = tmp;
			} // for k
			i++;
			j--;      
		} // if   
	} // while
  
	if( l < j ) wsortIntegerListSizeFunc(A, l, j, size, cmp);
	if( i < r ) wsortIntegerListSizeFunc(A, i, r, size, cmp);
}

int SortAndClearIntegerListUnit(int oldN, int *list, int *newN)
{
	int N = 1;
	int i;

	if(oldN <= 1) { newN[0] = oldN; return 0; }

	wsortIntegerListUnit(list, 0, oldN - 1);

	for(i=1; i<oldN; i++)
		if(list[N - 1] != list[i])
		{
			list[N] = list[i];
			N += 1;
		} // if

	newN[0] = N;

	return 0;
}

int SortAndClearIntegerListSizeUnit(int oldN, int *list, int size, int *newN)
{
	int N = 1;
	int i, k;

	if(oldN <= 1) { newN[0] = oldN; return 0; }

	wsortIntegerListSizeUnit(list, 0, oldN - 1, size);

	for(i=1; i<oldN; i++)
		if(list[(N - 1) * size] != list[i * size])
		{
			for(k=0; k<size; k++) list[N*size + k] = list[i*size + k];
			N += 1;
		} // if

	newN[0] = N;

	return 0;
}

int SortAndClearIntegerListSizeFunc(int oldN, int *list, int size, int (*cmp)(int *, int *), int *newN)
{
	int N = 1;
	int i, k;

	if(oldN <= 1) { newN[0] = oldN; return 0; }

	wsortIntegerListSizeFunc(list, 0, oldN - 1, size, cmp);

	for(i=1; i<oldN; i++)
		if((*cmp)(list + (N-1)*size, list + i*size) != 0)
		{
			for(k=0; k<size; k++) list[N*size + k] = list[i*size + k];
			N += 1;
		} // if

	newN[0] = N;

	return 0;
}

int DeleteListFromOrderedList(int nE, int *list, int oldEO, int *listO, int *newEO)
{
	char *ids;
	int i, j;

	if(oldEO <= 0) { newEO[0] = oldEO; return 0; }

	ids = (char *)calloc(oldEO, sizeof(char));

	for(i=0; i<nE; i++)
	{
		int get = DichotomySearchElementInOrderIntegerListUnit(listO, oldEO, list[i]);
		if(get != -1) ids[get] = 1;		
	} // for i

	for(j=0, i=0; i<oldEO; i++)
	{
		if(ids[i] != 1)
		{
			listO[j] = listO[i];
			j += 1;
		} // if ids[i]		
	} // for i

	free(ids);
	newEO[0] = j;	

	return 0;
}
