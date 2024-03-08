#include <stdio.h>
#include <stdlib.h>

#include "hm4ElementsVolumes.h"

double hm4TetrahedronVolume(double *P1, double *P2, double *P3, double *P4)
{
	double V = 0.0;

	V =  P1[0] * ( P2[1]*(P3[2]-P4[2]) - P3[1]*(P2[2]-P4[2]) + P4[1]*(P2[2]-P3[2]) )
		-P1[1] * ( P2[0]*(P3[2]-P4[2]) - P3[0]*(P2[2]-P4[2]) + P4[0]*(P2[2]-P3[2]) )
		+P1[2] * ( P2[0]*(P3[1]-P4[1]) - P3[0]*(P2[1]-P4[1]) + P4[0]*(P2[1]-P3[1]) )
		-        ( P2[0]*(P3[1]*P4[2]-P4[1]*P3[2]) - P3[0]*(P2[1]*P4[2]-P4[1]*P2[2]) + P4[0]*(P2[1]*P3[2]-P3[1]*P2[2]) );

	V /= 6.0;

	if(V < 0.0)
		return -V;
	else
		return V;
}

double hm4PyramidVolume(double *P0, double *P1, double *P2, double *P3, double *P4)
{
	double  vFIRST = hm4TetrahedronVolume(P0, P1, P2, P4);
	double vSECOND = hm4TetrahedronVolume(P1, P2, P3, P4);
	
	return vFIRST + vSECOND;
}

double hm4PrismVolume(double *P0, double *P1, double *P2, double *P3, double *P4, double *P5)
{
	double  vFIRST = hm4TetrahedronVolume(P0, P1, P2, P4);
	double vSECOND = hm4TetrahedronVolume(P0, P2, P3, P4);
	double  vTHIRD = hm4TetrahedronVolume(P2, P5, P3, P4);

	return vFIRST + vSECOND + vTHIRD;
}

double hm4HexahedronVolume(double *P0, double *P1, double *P2, double *P3, double *P4, double *P5, double *P6, double *P7)
{
    double  vFIRST = hm4TetrahedronVolume(P4, P6, P5, P2);
	double vSECOND = hm4TetrahedronVolume(P5, P6, P7, P2);
	double  vTHIRD = hm4TetrahedronVolume(P4, P5, P0, P2);
	double vFOURTH = hm4TetrahedronVolume(P2, P3, P7, P5);
	double  vFIFTH = hm4TetrahedronVolume(P5, P3, P1, P2);
	double  vSIXTH = hm4TetrahedronVolume(P0, P1, P2, P5);

	return vFIRST + vSECOND + vTHIRD + vFOURTH + vFIFTH + vSIXTH;
}

int hm4InitMeshElementsVolumes(int nE, int *xE, int *aE, double *C, double *V)
{
	int i;

	for(i=0; i<nE; i++)
	{
		int n = xE[i + 1] - xE[i];
		int *iPTR = aE + xE[i];

		if(n == 4)
		{
			V[i] = hm4TetrahedronVolume(C + iPTR[0]*3, C + iPTR[1]*3, C + iPTR[2]*3, C + iPTR[3]*3);
		}

		if(n == 5)
		{
			V[i] = hm4PyramidVolume(C + iPTR[0]*3, C + iPTR[1]*3, C + iPTR[2]*3, C + iPTR[3]*3, C + iPTR[4]*3);
		}

		if(n == 6)
		{
			V[i] = hm4PrismVolume(C + iPTR[0]*3, C + iPTR[1]*3, C + iPTR[2]*3, C + iPTR[3]*3, C + iPTR[4]*3, C + iPTR[5]*3);
		}

		if(n == 8)
		{
			V[i] = hm4HexahedronVolume(C + iPTR[0]*3, C + iPTR[1]*3, C + iPTR[2]*3, C + iPTR[3]*3, 
				                       C + iPTR[4]*3, C + iPTR[5]*3, C + iPTR[6]*3, C + iPTR[7]*3);
		}
	} // for i

	return 0;
}
