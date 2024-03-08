#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define _GAMMA_ 1.4
#define _GAMMA_1_ 0.4
#define _1_DIV_GAMMA_1_ 2.5

#include "solverEuler3D.h"

int solverRusanov3D(double *QphR, double *QphL,
					 double *flux,
					 double  nX, double  nY, double  nZ,
					 double  sX, double  sY, double  sZ,
					 double  tX, double  tY, double  tZ)
{
	double rhoL, uL, vL, wL, pL, eL, cL;
	double rhoR, uR, vR, wR, pR, eR, cR;
	double CC;
	double fluxEU[5];

	rhoL = QphL[0];                                      rhoR = QphR[0];
	uL = QphL[1] * nX + QphL[2] * nY + QphL[3] * nZ;     uR = QphR[1] * nX + QphR[2] * nY + QphR[3] * nZ;
	vL = QphL[1] * sX + QphL[2] * sY + QphL[3] * sZ;     vR = QphR[1] * sX + QphR[2] * sY + QphR[3] * sZ;
	wL = QphL[1] * tX + QphL[2] * tY + QphL[3] * tZ;     wR = QphR[1] * tX + QphR[2] * tY + QphR[3] * tZ;
	pL = QphL[4];                                        pR = QphR[4];

	eL = _1_DIV_GAMMA_1_ * pL + 0.5 * rhoL * (uL*uL + vL*vL + wL*wL);
	eR = _1_DIV_GAMMA_1_ * pR + 0.5 * rhoR * (uR*uR + vR*vR + wR*wR);
	cL = sqrt(_GAMMA_ * pL / rhoL);
	cR = sqrt(_GAMMA_ * pR / rhoR);

	fluxEU[0] = 0.5 * (rhoL * uL + rhoR * uR);
    fluxEU[1] = 0.5 * (rhoL * uL * uL + pL + rhoR * uR * uR + pR);
    fluxEU[2] = 0.5 * (rhoL * vL * uL +      rhoR * vR * uR);
    fluxEU[3] = 0.5 * (rhoL * wL * uL +      rhoR * wR * uR);
    fluxEU[4] = 0.5 * ((eL + pL) * uL + (eR + pR) * uR);

	                       CC = fabs(uR     );
	if(CC < fabs(uR + cR)) CC = fabs(uR + cR);
	if(CC < fabs(uR - cR)) CC = fabs(uR - cR);

	if(CC < fabs(uL     )) CC = fabs(uL     );
	if(CC < fabs(uL + cL)) CC = fabs(uL + cL);
	if(CC < fabs(uL - cL)) CC = fabs(uL - cL);

	fluxEU[0] += 0.5 * CC * (rhoL      - rhoR   );
	fluxEU[1] += 0.5 * CC * (rhoL * uL - rhoR*uR);
	fluxEU[2] += 0.5 * CC * (rhoL * vL - rhoR*vR);
	fluxEU[3] += 0.5 * CC * (rhoL * wL - rhoR*wR);
	fluxEU[4] += 0.5 * CC * (       eL -      eR);

	flux[0] =  fluxEU[0];
	flux[1] =  fluxEU[1] * nX + fluxEU[2] * sX + fluxEU[3] * tX;
    flux[2] =  fluxEU[1] * nY + fluxEU[2] * sY + fluxEU[3] * tY;
	flux[3] =  fluxEU[1] * nZ + fluxEU[2] * sZ + fluxEU[3] * tZ;
	flux[4] =  fluxEU[4];

	return 0;
}
