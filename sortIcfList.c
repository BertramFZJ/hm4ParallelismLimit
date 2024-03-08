#include <stdio.h>
#include <stdlib.h>

#include "cfdHeader.h"
#include "sortIcfList.h"

// ========================================================================================================
// 0 --> a == b
// 1 --> a > b
// (-1) --> a < b
static int compareIcf(gtypeInteriorCellFace a, gtypeInteriorCellFace b, int *level);
static int compareIcf(gtypeInteriorCellFace a, gtypeInteriorCellFace b, int *level)
{
	int levelLA = level[a.eL];
	int levelRA = level[a.eR];

	int levelLB = level[b.eL];
	int levelRB = level[b.eR];
	
	if(levelLA > levelRA) { int B = levelLA; levelLA = levelRA; levelRA = B; }
	if(levelLB > levelRB) { int B = levelLB; levelLB = levelRB; levelRB = B; }

	if(levelLA > levelLB) return  1;
	if(levelLA < levelLB) return -1;

	if(levelRA > levelRB) return  1;
	if(levelRA < levelRB) return -1;	

	if(a.eL > b.eL) return  1;
	if(a.eL < b.eL) return -1;

	if(a.eR > b.eR) return  1;
	if(a.eR < b.eR) return -1;

	return 0;
}
// ========================================================================================================

void wsortIcfListLevel(gtypeInteriorCellFace *A, int l, int r, int *level)
{
	gtypeInteriorCellFace tmp;
	gtypeInteriorCellFace B = A[(l + r) / 2];
	int i = l, j = r;

	while( i <= j )
	{
		while( compareIcf(A[i], B, level) == -1 ) i++;
		while( compareIcf(A[j], B, level) ==  1 ) j--;
		if( i <= j )
		{
			tmp = A[i];
			A[i] = A[j];
			A[j] = tmp;
			i++;
			j--;      
		} // if   
	} // while
  
	if( l < j ) wsortIcfListLevel(A, l, j, level);
	if( i < r ) wsortIcfListLevel(A, i, r, level);
}

// ========================================================================================================
// 0 --> a == b
// 1 --> a > b
// (-1) --> a < b
static int compareBcf(gtypeBoundaryCellFace a, gtypeBoundaryCellFace b, int *level);
static int compareBcf(gtypeBoundaryCellFace a, gtypeBoundaryCellFace b, int *level)
{
	int plA = level[a.eL]; 
	int plB = level[b.eL];

	if(plA > plB) return  1;
	if(plA < plB) return -1;

	if(a.eL > b.eL) return  1;
	if(a.eL < b.eL) return -1;	

	return 0;
}
// ========================================================================================================

void wsortBcfListLevel(gtypeBoundaryCellFace *A, int l, int r, int *level)
{
	gtypeBoundaryCellFace tmp;
	gtypeBoundaryCellFace B = A[(l + r) / 2];
	int i = l, j = r;

	while( i <= j )
	{
		while( compareBcf(A[i], B, level) == -1 ) i++;
		while( compareBcf(A[j], B, level) ==  1 ) j--;
		if( i <= j )
		{
			tmp = A[i];
			A[i] = A[j];
			A[j] = tmp;
			i++;
			j--;      
		} // if   
	} // while
  
	if( l < j ) wsortBcfListLevel(A, l, j, level);
	if( i < r ) wsortBcfListLevel(A, i, r, level);
}
