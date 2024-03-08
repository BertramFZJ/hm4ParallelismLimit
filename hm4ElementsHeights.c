#define _CUBE_EDGE_HEIGHT_ 0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hm4ElementsFaces.h"

#include "hm4ElementsHeights.h"

#if _CUBE_EDGE_HEIGHT_
double hm4TetrahedronHeight(double *P1, double *P2, double *P3, double *P4)
{
	double xMinMax[2] = {P1[0], P1[0]};
	double yMinMax[2] = {P1[1], P1[1]};
	double zMinMax[2] = {P1[2], P1[2]};
	double h;

	if(xMinMax[0] > P2[0]) xMinMax[0] = P2[0]; if(xMinMax[1] < P2[0]) xMinMax[1] = P2[0];
	if(yMinMax[0] > P2[1]) yMinMax[0] = P2[1]; if(yMinMax[1] < P2[1]) yMinMax[1] = P2[1];
	if(zMinMax[0] > P2[2]) zMinMax[0] = P2[2]; if(zMinMax[1] < P2[2]) zMinMax[1] = P2[2];

	if(xMinMax[0] > P3[0]) xMinMax[0] = P3[0]; if(xMinMax[1] < P3[0]) xMinMax[1] = P3[0];
	if(yMinMax[0] > P3[1]) yMinMax[0] = P3[1]; if(yMinMax[1] < P3[1]) yMinMax[1] = P3[1];
	if(zMinMax[0] > P3[2]) zMinMax[0] = P3[2]; if(zMinMax[1] < P3[2]) zMinMax[1] = P3[2];

	if(xMinMax[0] > P4[0]) xMinMax[0] = P4[0]; if(xMinMax[1] < P4[0]) xMinMax[1] = P4[0];
	if(yMinMax[0] > P4[1]) yMinMax[0] = P4[1]; if(yMinMax[1] < P4[1]) yMinMax[1] = P4[1];
	if(zMinMax[0] > P4[2]) zMinMax[0] = P4[2]; if(zMinMax[1] < P4[2]) zMinMax[1] = P4[2];

	h = xMinMax[1] - xMinMax[0];
	yMinMax[1] -= yMinMax[0]; if(h > yMinMax[1]) h = yMinMax[1];
	zMinMax[1] -= zMinMax[0]; if(h > zMinMax[1]) h = zMinMax[1];	

	return h;	
}
#else
double hm4TetrahedronHeight(double *P1, double *P2, double *P3, double *P4)
{
	double *C[4] = {P1, P2, P3, P4};
	
	int index[4][3] = { {1,  0,  2} , {0,  1,  3}, {1,  2,  3}, {2,  0,  3} };
	int  vtxs[4]    = {3, 2, 0, 1};
	
	double h;

	int i;

	for(i=0; i<4; i++)
	{
		double TC[3], face[4], hLOC;

		hm4TriangleCenter(C[index[i][0]], C[index[i][1]], C[index[i][2]], TC  );
		  hm4TriangleFace(C[index[i][0]], C[index[i][1]], C[index[i][2]], face);

		face[3] = 1.0 / sqrt(face[0]*face[0] + face[1]*face[1] + face[2]*face[2]);
		
		face[0] *= face[3]; 
		face[1] *= face[3]; 
		face[2] *= face[3];

		face[3] = - (face[0]*TC[0] + face[1]*TC[1] + face[2]*TC[2]);

		hLOC = - face[0]*C[vtxs[i]][0] - face[1]*C[vtxs[i]][1] - face[2]*C[vtxs[i]][2] - face[3];
		if(i == 0) h = hLOC;
		else if(h > hLOC) h = hLOC;
	} // for i
		
	if(h <= 0.0) { fprintf(stderr, "ERROR hm4TetrahedronHeight\n"); exit(0); }
	
	return h;
}
#endif

#if _CUBE_EDGE_HEIGHT_
double hm4PyramidHeight(double *P0, double *P1, double *P2, double *P3, double *P4)
{
	double xMinMax[2] = {P0[0], P0[0]};
	double yMinMax[2] = {P0[1], P0[1]};
	double zMinMax[2] = {P0[2], P0[2]};
	double h;

	if(xMinMax[0] > P1[0]) xMinMax[0] = P1[0]; if(xMinMax[1] < P1[0]) xMinMax[1] = P1[0];
	if(yMinMax[0] > P1[1]) yMinMax[0] = P1[1]; if(yMinMax[1] < P1[1]) yMinMax[1] = P1[1];
	if(zMinMax[0] > P1[2]) zMinMax[0] = P1[2]; if(zMinMax[1] < P1[2]) zMinMax[1] = P1[2];

	if(xMinMax[0] > P2[0]) xMinMax[0] = P2[0]; if(xMinMax[1] < P2[0]) xMinMax[1] = P2[0];
	if(yMinMax[0] > P2[1]) yMinMax[0] = P2[1]; if(yMinMax[1] < P2[1]) yMinMax[1] = P2[1];
	if(zMinMax[0] > P2[2]) zMinMax[0] = P2[2]; if(zMinMax[1] < P2[2]) zMinMax[1] = P2[2];

	if(xMinMax[0] > P3[0]) xMinMax[0] = P3[0]; if(xMinMax[1] < P3[0]) xMinMax[1] = P3[0];
	if(yMinMax[0] > P3[1]) yMinMax[0] = P3[1]; if(yMinMax[1] < P3[1]) yMinMax[1] = P3[1];
	if(zMinMax[0] > P3[2]) zMinMax[0] = P3[2]; if(zMinMax[1] < P3[2]) zMinMax[1] = P3[2];

	if(xMinMax[0] > P4[0]) xMinMax[0] = P4[0]; if(xMinMax[1] < P4[0]) xMinMax[1] = P4[0];
	if(yMinMax[0] > P4[1]) yMinMax[0] = P4[1]; if(yMinMax[1] < P4[1]) yMinMax[1] = P4[1];
	if(zMinMax[0] > P4[2]) zMinMax[0] = P4[2]; if(zMinMax[1] < P4[2]) zMinMax[1] = P4[2];

	h = xMinMax[1] - xMinMax[0];
	yMinMax[1] -= yMinMax[0]; if(h > yMinMax[1]) h = yMinMax[1];
	zMinMax[1] -= zMinMax[0]; if(h > zMinMax[1]) h = zMinMax[1];	

	return h;	
}
#else
double hm4PyramidHeight(double *P0, double *P1, double *P2, double *P3, double *P4)
{
	double *C[5] = {P0, P1, P2, P3, P4};
	
	int index[5][4] = { {0,  2,  3,  1}, {0,  1,  4, -1}, {1,  3,  4, -1}, {3,  2,  4, -1}, {2,  0,  4, -1} };
	int  vtxs[4][2] = { {2, 3}, {0, 2}, {0, 1}, {1, 3} };
	
	double h;

	int i, j;

	{
		double QC[3], face[4];

		hm4QuadCenter(C[index[0][0]], C[index[0][1]], C[index[0][2]], C[index[0][3]], QC  );
		  hm4QuadFace(C[index[0][0]], C[index[0][1]], C[index[0][2]], C[index[0][3]], face);

		face[3] = 1.0 / sqrt(face[0]*face[0] + face[1]*face[1] + face[2]*face[2]);
		
		face[0] *= face[3]; 
		face[1] *= face[3]; 
		face[2] *= face[3];

		face[3] = - (face[0]*QC[0] + face[1]*QC[1] + face[2]*QC[2]);
				
		h = - face[0]*C[4][0] 
		    - face[1]*C[4][1] 
			- face[2]*C[4][2] 
			- face[3];		
	}

	for(i=1; i<5; i++)
	{
		double TC[3], face[4], hLOC;

		hm4TriangleCenter(C[index[i][0]], C[index[i][1]], C[index[i][2]], TC  );
		  hm4TriangleFace(C[index[i][0]], C[index[i][1]], C[index[i][2]], face);

		face[3] = 1.0 / sqrt(face[0]*face[0] + face[1]*face[1] + face[2]*face[2]);
		
		face[0] *= face[3]; 
		face[1] *= face[3]; 
		face[2] *= face[3];

		face[3] = - (face[0]*TC[0] + face[1]*TC[1] + face[2]*TC[2]);

		for(j=0; j<2; j++)
		{
			hLOC = - face[0]*C[vtxs[i-1][j]][0] 
			       - face[1]*C[vtxs[i-1][j]][1] 
				   - face[2]*C[vtxs[i-1][j]][2] 
				   - face[3];
			if(h > hLOC) h = hLOC;
		} // for j
	} // for i
		
	if(h <= 0.0) { fprintf(stderr, "ERROR hm4PyramidHeight\n"); exit(0); }
	
	return h;
}
#endif

#if _CUBE_EDGE_HEIGHT_
double hm4PrismHeight(double *P0, double *P1, double *P2, double *P3, double *P4, double *P5)
{
	double xMinMax[2] = {P0[0], P0[0]};
	double yMinMax[2] = {P0[1], P0[1]};
	double zMinMax[2] = {P0[2], P0[2]};
	double h;

	if(xMinMax[0] > P1[0]) xMinMax[0] = P1[0]; if(xMinMax[1] < P1[0]) xMinMax[1] = P1[0];
	if(yMinMax[0] > P1[1]) yMinMax[0] = P1[1]; if(yMinMax[1] < P1[1]) yMinMax[1] = P1[1];
	if(zMinMax[0] > P1[2]) zMinMax[0] = P1[2]; if(zMinMax[1] < P1[2]) zMinMax[1] = P1[2];

	if(xMinMax[0] > P2[0]) xMinMax[0] = P2[0]; if(xMinMax[1] < P2[0]) xMinMax[1] = P2[0];
	if(yMinMax[0] > P2[1]) yMinMax[0] = P2[1]; if(yMinMax[1] < P2[1]) yMinMax[1] = P2[1];
	if(zMinMax[0] > P2[2]) zMinMax[0] = P2[2]; if(zMinMax[1] < P2[2]) zMinMax[1] = P2[2];

	if(xMinMax[0] > P3[0]) xMinMax[0] = P3[0]; if(xMinMax[1] < P3[0]) xMinMax[1] = P3[0];
	if(yMinMax[0] > P3[1]) yMinMax[0] = P3[1]; if(yMinMax[1] < P3[1]) yMinMax[1] = P3[1];
	if(zMinMax[0] > P3[2]) zMinMax[0] = P3[2]; if(zMinMax[1] < P3[2]) zMinMax[1] = P3[2];

	if(xMinMax[0] > P4[0]) xMinMax[0] = P4[0]; if(xMinMax[1] < P4[0]) xMinMax[1] = P4[0];
	if(yMinMax[0] > P4[1]) yMinMax[0] = P4[1]; if(yMinMax[1] < P4[1]) yMinMax[1] = P4[1];
	if(zMinMax[0] > P4[2]) zMinMax[0] = P4[2]; if(zMinMax[1] < P4[2]) zMinMax[1] = P4[2];

	if(xMinMax[0] > P5[0]) xMinMax[0] = P5[0]; if(xMinMax[1] < P5[0]) xMinMax[1] = P5[0];
	if(yMinMax[0] > P5[1]) yMinMax[0] = P5[1]; if(yMinMax[1] < P5[1]) yMinMax[1] = P5[1];
	if(zMinMax[0] > P5[2]) zMinMax[0] = P5[2]; if(zMinMax[1] < P5[2]) zMinMax[1] = P5[2];

	h = xMinMax[1] - xMinMax[0];
	yMinMax[1] -= yMinMax[0]; if(h > yMinMax[1]) h = yMinMax[1];
	zMinMax[1] -= zMinMax[0]; if(h > zMinMax[1]) h = zMinMax[1];	

	return h;	
}
#else
double hm4PrismHeight(double *P0, double *P1, double *P2, double *P3, double *P4, double *P5)
{
	double *C[6] = {P0, P1, P2, P3, P4, P5};
	
	int index[5][4] = { {0,  1,  4,  3}, {1,  2,  5,  4}, {2,  0,  3,  5}, {0,  2,  1, -1}, {3,  4,  5, -1} };
	
	int vtxs2[3][2] = { {2, 5}, {0, 3}, {1, 4} };
	int vtxs3[2][3] = { {3, 4, 5}, {0, 1, 2} };
	
	double h;

	int i, j;

	for(i=0; i<3; i++)
	{
		double QC[3], face[4], hLOC;

		hm4QuadCenter(C[index[i][0]], C[index[i][1]], C[index[i][2]], C[index[i][3]], QC  );
		  hm4QuadFace(C[index[i][0]], C[index[i][1]], C[index[i][2]], C[index[i][3]], face);

		face[3] = 1.0 / sqrt(face[0]*face[0] + face[1]*face[1] + face[2]*face[2]);
		
		face[0] *= face[3]; 
		face[1] *= face[3]; 
		face[2] *= face[3];

		face[3] = - (face[0]*QC[0] + face[1]*QC[1] + face[2]*QC[2]);
		
		for(j=0; j<2; j++)
		{
			hLOC = - face[0]*C[vtxs2[i][j]][0] 
				   - face[1]*C[vtxs2[i][j]][1] 
				   - face[2]*C[vtxs2[i][j]][2] 
			       - face[3];
			if( (i == 0) && (j == 0) ) h = hLOC;
			else if(h > hLOC) h = hLOC;
		} // for j
	} // for i

	for(i=3; i<5; i++)
	{
		double TC[3], face[4], hLOC;

		hm4TriangleCenter(C[index[i][0]], C[index[i][1]], C[index[i][2]], TC  );
		  hm4TriangleFace(C[index[i][0]], C[index[i][1]], C[index[i][2]], face);

		face[3] = 1.0 / sqrt(face[0]*face[0] + face[1]*face[1] + face[2]*face[2]);
		
		face[0] *= face[3]; 
		face[1] *= face[3]; 
		face[2] *= face[3];

		face[3] = - (face[0]*TC[0] + face[1]*TC[1] + face[2]*TC[2]);

		for(j=0; j<3; j++)
		{
			hLOC = - face[0]*C[vtxs3[i-3][j]][0] 
			       - face[1]*C[vtxs3[i-3][j]][1] 
				   - face[2]*C[vtxs3[i-3][j]][2] 
				   - face[3];
			if(h > hLOC) h = hLOC;
		} // for j
	} // for i
		
	if(h <= 0.0) { fprintf(stderr, "ERROR hm4PrismHeight\n"); exit(0); }
	
	return h;
}
#endif

#if _CUBE_EDGE_HEIGHT_
double hm4HexahedronHeight(double *P0, double *P1, double *P2, double *P3, double *P4, double *P5, double *P6, double *P7)
{
	double xMinMax[2] = {P0[0], P0[0]};
	double yMinMax[2] = {P0[1], P0[1]};
	double zMinMax[2] = {P0[2], P0[2]};
	double h;

	if(xMinMax[0] > P1[0]) xMinMax[0] = P1[0]; if(xMinMax[1] < P1[0]) xMinMax[1] = P1[0];
	if(yMinMax[0] > P1[1]) yMinMax[0] = P1[1]; if(yMinMax[1] < P1[1]) yMinMax[1] = P1[1];
	if(zMinMax[0] > P1[2]) zMinMax[0] = P1[2]; if(zMinMax[1] < P1[2]) zMinMax[1] = P1[2];

	if(xMinMax[0] > P2[0]) xMinMax[0] = P2[0]; if(xMinMax[1] < P2[0]) xMinMax[1] = P2[0];
	if(yMinMax[0] > P2[1]) yMinMax[0] = P2[1]; if(yMinMax[1] < P2[1]) yMinMax[1] = P2[1];
	if(zMinMax[0] > P2[2]) zMinMax[0] = P2[2]; if(zMinMax[1] < P2[2]) zMinMax[1] = P2[2];

	if(xMinMax[0] > P3[0]) xMinMax[0] = P3[0]; if(xMinMax[1] < P3[0]) xMinMax[1] = P3[0];
	if(yMinMax[0] > P3[1]) yMinMax[0] = P3[1]; if(yMinMax[1] < P3[1]) yMinMax[1] = P3[1];
	if(zMinMax[0] > P3[2]) zMinMax[0] = P3[2]; if(zMinMax[1] < P3[2]) zMinMax[1] = P3[2];

	if(xMinMax[0] > P4[0]) xMinMax[0] = P4[0]; if(xMinMax[1] < P4[0]) xMinMax[1] = P4[0];
	if(yMinMax[0] > P4[1]) yMinMax[0] = P4[1]; if(yMinMax[1] < P4[1]) yMinMax[1] = P4[1];
	if(zMinMax[0] > P4[2]) zMinMax[0] = P4[2]; if(zMinMax[1] < P4[2]) zMinMax[1] = P4[2];

	if(xMinMax[0] > P5[0]) xMinMax[0] = P5[0]; if(xMinMax[1] < P5[0]) xMinMax[1] = P5[0];
	if(yMinMax[0] > P5[1]) yMinMax[0] = P5[1]; if(yMinMax[1] < P5[1]) yMinMax[1] = P5[1];
	if(zMinMax[0] > P5[2]) zMinMax[0] = P5[2]; if(zMinMax[1] < P5[2]) zMinMax[1] = P5[2];

	if(xMinMax[0] > P6[0]) xMinMax[0] = P6[0]; if(xMinMax[1] < P6[0]) xMinMax[1] = P6[0];
	if(yMinMax[0] > P6[1]) yMinMax[0] = P6[1]; if(yMinMax[1] < P6[1]) yMinMax[1] = P6[1];
	if(zMinMax[0] > P6[2]) zMinMax[0] = P6[2]; if(zMinMax[1] < P6[2]) zMinMax[1] = P6[2];

	if(xMinMax[0] > P7[0]) xMinMax[0] = P7[0]; if(xMinMax[1] < P7[0]) xMinMax[1] = P7[0];
	if(yMinMax[0] > P7[1]) yMinMax[0] = P7[1]; if(yMinMax[1] < P7[1]) yMinMax[1] = P7[1];
	if(zMinMax[0] > P7[2]) zMinMax[0] = P7[2]; if(zMinMax[1] < P7[2]) zMinMax[1] = P7[2];

	h = xMinMax[1] - xMinMax[0];
	yMinMax[1] -= yMinMax[0]; if(h > yMinMax[1]) h = yMinMax[1];
	zMinMax[1] -= zMinMax[0]; if(h > zMinMax[1]) h = zMinMax[1];	

	return h;	
}
#else
double hm4HexahedronHeight(double *P0, double *P1, double *P2, double *P3, double *P4, double *P5, double *P6, double *P7)
{
	double *C[8] = {P0, P1, P2, P3, P4, P5, P6, P7};
	
	int index[6][4] = { {0, 1, 5, 4}, {1, 3, 7, 5}, {3, 2, 6, 7}, {2, 0, 4, 6}, {1, 0, 2, 3}, {4, 5, 7, 6} };	
	int  vtxs[6][4] = { {2, 3, 6, 7}, {0, 2, 4, 6}, {0, 1, 4, 5}, {1, 3, 5, 7}, {4, 5, 6, 7}, {0, 1, 2, 3} };
		
	double h;

	int i, j;

	for(i=0; i<6; i++)
	{
		double QC[3], face[4], hLOC;

		hm4QuadCenter(C[index[i][0]], C[index[i][1]], C[index[i][2]], C[index[i][3]], QC  );
		  hm4QuadFace(C[index[i][0]], C[index[i][1]], C[index[i][2]], C[index[i][3]], face);

		face[3] = 1.0 / sqrt(face[0]*face[0] + face[1]*face[1] + face[2]*face[2]);
		
		face[0] *= face[3]; 
		face[1] *= face[3]; 
		face[2] *= face[3];

		face[3] = - (face[0]*QC[0] + face[1]*QC[1] + face[2]*QC[2]);
		
		for(j=0; j<4; j++)
		{
			hLOC = - face[0]*C[vtxs[i][j]][0] 
				   - face[1]*C[vtxs[i][j]][1] 
				   - face[2]*C[vtxs[i][j]][2] 
			       - face[3];
			if( (i == 0) && (j == 0) ) h = hLOC;
			else if(h > hLOC) h = hLOC;
		} // for j
	} // for i

	if(h <= 0.0) { fprintf(stderr, "ERROR hm4HexahedronHeight\n"); exit(0); }
	
	return h;
}
#endif

int hm4InitMeshElementsHeights(int nE, int *xE, int *aE, double *C, double *He)
{
	int i;

	for(i=0; i<nE; i++)
	{
		int n = xE[i + 1] - xE[i];
		int *iPTR = aE + xE[i];

		if(n == 4)
		{
			He[i] = hm4TetrahedronHeight(C + iPTR[0]*3, C + iPTR[1]*3, C + iPTR[2]*3, C + iPTR[3]*3);
		}

		if(n == 5)
		{
			He[i] = hm4PyramidHeight(C + iPTR[0]*3, C + iPTR[1]*3, C + iPTR[2]*3, C + iPTR[3]*3, C + iPTR[4]*3);
		}

		if(n == 6)
		{
			He[i] = hm4PrismHeight(C + iPTR[0]*3, C + iPTR[1]*3, C + iPTR[2]*3, C + iPTR[3]*3, C + iPTR[4]*3, C + iPTR[5]*3);
		}

		if(n == 8)
		{
			He[i] = hm4HexahedronHeight(C + iPTR[0]*3, C + iPTR[1]*3, C + iPTR[2]*3, C + iPTR[3]*3, 
				                        C + iPTR[4]*3, C + iPTR[5]*3, C + iPTR[6]*3, C + iPTR[7]*3);
		}
	} // for i

	return 0;
}
