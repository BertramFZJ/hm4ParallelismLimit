#include <stdio.h>
#include <stdlib.h>

#include "hm4ElementsVolumes.h"

#include "hm4ElementsCenters.h"

int hm4TetrahedronCenter(double *P1, double *P2, double *P3, double *P4, double *CEL)
{
	int i;

	for(i=0; i<3; i++) CEL[i] = 0.25 * (P1[i] + P2[i] + P3[i] + P4[i]);	

	return 0;
}

int hm4PyramidCenter(double *P0, double *P1, double *P2, double *P3, double *P4, double *CEL)
{
	double VOL[3];
	double TCENTER[6];

	int i;
		
	VOL[0] = hm4TetrahedronVolume(P0, P1, P2, P4);
	VOL[1] = hm4TetrahedronVolume(P1, P2, P3, P4);
	VOL[2] = VOL[0] + VOL[1];

	hm4TetrahedronCenter(P0, P1, P2, P4, TCENTER    );
	hm4TetrahedronCenter(P1, P2, P3, P4, TCENTER + 3);

	for(i=0; i<3; i++) CEL[i] = (TCENTER[i] * VOL[0] + TCENTER[i + 3] * VOL[1]) / VOL[2];

	return 0;
}

int hm4PrismCenter(double *P0, double *P1, double *P2, double *P3, double *P4, double *P5, double *CEL)
{
	double VOL[4];
	double TCENTER[9];

	int i;

	VOL[0] = hm4TetrahedronVolume(P0, P1, P2, P4);
	VOL[1] = hm4TetrahedronVolume(P0, P2, P3, P4); 
	VOL[2] = hm4TetrahedronVolume(P2, P5, P3, P4);
	VOL[3] = VOL[0] + VOL[1] + VOL[2];

	hm4TetrahedronCenter(P0, P1, P2, P4, TCENTER    );
	hm4TetrahedronCenter(P0, P2, P3, P4, TCENTER + 3); 
	hm4TetrahedronCenter(P2, P5, P3, P4, TCENTER + 6);

	for(i=0; i<3; i++) CEL[i] = (TCENTER[i] * VOL[0] + TCENTER[i + 3] * VOL[1] + TCENTER[i + 6] * VOL[2]) / VOL[3];

	return 0;
}

int hm4HexahedronCenter(double *P0, double *P1, double *P2, double *P3, double *P4, double *P5, double *P6, double *P7, double *CEL)
{
	double VOL[7];
	double TCENTER[18];

	int i;

	VOL[0] = hm4TetrahedronVolume(P4, P0, P5, P7);
	VOL[1] = hm4TetrahedronVolume(P0, P5, P1, P7);
	VOL[2] = hm4TetrahedronVolume(P0, P1, P7, P3);
	VOL[3] = hm4TetrahedronVolume(P0, P2, P3, P7);
	VOL[4] = hm4TetrahedronVolume(P0, P2, P6, P7);
	VOL[5] = hm4TetrahedronVolume(P0, P4, P6, P7);
	VOL[6] = VOL[0] + VOL[1] + VOL[2] + VOL[3] + VOL[4] + VOL[5];

	hm4TetrahedronCenter(P4, P0, P5, P7, TCENTER     );
	hm4TetrahedronCenter(P0, P5, P1, P7, TCENTER + 3 );
	hm4TetrahedronCenter(P0, P1, P7, P3, TCENTER + 6 );
	hm4TetrahedronCenter(P0, P2, P3, P7, TCENTER + 9 );
	hm4TetrahedronCenter(P0, P2, P6, P7, TCENTER + 12);
	hm4TetrahedronCenter(P0, P4, P6, P7, TCENTER + 15);

	for(i=0; i<3; i++) CEL[i] = ( TCENTER[i    ] * VOL[0] + TCENTER[i + 3 ] * VOL[1] + TCENTER[i + 6 ] * VOL[2]
	                            + TCENTER[i + 9] * VOL[3] + TCENTER[i + 12] * VOL[4] + TCENTER[i + 15] * VOL[5]) 
								/ VOL[6];

	return 0;
}

int hm4InitMeshElementsCenters(int nE, int *xE, int *aE, double *C, double *CEL)
{
	int i;

	for(i=0; i<nE; i++)
	{
		int n = xE[i + 1] - xE[i];
		int *iPTR = aE + xE[i];

		if(n == 4)
		{
			hm4TetrahedronCenter(C + iPTR[0]*3, C + iPTR[1]*3, C + iPTR[2]*3, C + iPTR[3]*3, CEL + i*3);
		}

		if(n == 5)
		{
			hm4PyramidCenter(C + iPTR[0]*3, C + iPTR[1]*3, C + iPTR[2]*3, C + iPTR[3]*3, C + iPTR[4]*3, CEL + i*3);
		}

		if(n == 6)
		{
			hm4PrismCenter(C + iPTR[0]*3, C + iPTR[1]*3, C + iPTR[2]*3, C + iPTR[3]*3, C + iPTR[4]*3, C + iPTR[5]*3, CEL + i*3);
		}

		if(n == 8)
		{
			hm4HexahedronCenter(C + iPTR[0]*3, C + iPTR[1]*3, C + iPTR[2]*3, C + iPTR[3]*3, 
				                C + iPTR[4]*3, C + iPTR[5]*3, C + iPTR[6]*3, C + iPTR[7]*3,
								CEL + i*3);
		}
	} // for i

	return 0;
}
