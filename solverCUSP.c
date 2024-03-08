#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "solverEuler3D.h"

#define _GAMMA_1_ 0.4
#define _1_DIV_GAMMA_1_ 2.5
#define _GAMMA_DIV_GAMMA_1_ 3.5
#define _GAMMA_ 1.4

int solverCUSP3D(double *QphR, double *QphL,
				 double *flux,
				 double    nX, double    nY, double    nZ,
				 double    sX, double    sY, double    sZ,
				 double    tX, double    tY, double    tZ)
{
	double rhoL, uL, vL, wL, pL, eL, HL, aL;
	double rhoR, uR, vR, wR, pR, eR, HR, aR;
	double rhoAV, uAV, vAV, wAV, HAV, cAV;
	double Mn, eigenPLUS, eigenMINUS, beta, alfaC;

	double fluxEU[5];

	// 1. œŒ¬Œ–Œ“ — Œ–Œ—“≈…
	rhoL = QphL[0];                                      rhoR = QphR[0];
	uL = QphL[1] * nX + QphL[2] * nY + QphL[3] * nZ;     uR = QphR[1] * nX + QphR[2] * nY + QphR[3] * nZ;
	vL = QphL[1] * sX + QphL[2] * sY + QphL[3] * sZ;     vR = QphR[1] * sX + QphR[2] * sY + QphR[3] * sZ;
	wL = QphL[1] * tX + QphL[2] * tY + QphL[3] * tZ;     wR = QphR[1] * tX + QphR[2] * tY + QphR[3] * tZ;
	pL = QphL[4];                                        pR = QphR[4];
	
	// 2. 
	eL = _1_DIV_GAMMA_1_ * pL + 0.5 * rhoL * (uL*uL + vL*vL + wL*wL);
	eR = _1_DIV_GAMMA_1_ * pR + 0.5 * rhoR * (uR*uR + vR*vR + wR*wR);

	HL = (eL + pL) / rhoL; aL = sqrt(_GAMMA_1_ * pL / rhoL);
	HR = (eR + pR) / rhoR; aR = sqrt(_GAMMA_1_ * pR / rhoR);

	// ÷≈Õ“–¿À‹Õ¿ﬂ –¿«ÕŒ—“‹
	fluxEU[0] = 0.5 * (rhoL * uL + rhoR * uR);
    fluxEU[1] = 0.5 * (rhoL * uL * uL + pL + rhoR * uR * uR + pR);
    fluxEU[2] = 0.5 * (rhoL * vL * uL +      rhoR * vR * uR);
    fluxEU[3] = 0.5 * (rhoL * wL * uL +      rhoR * wR * uR);
    fluxEU[4] = 0.5 * ((eL + pL) * uL + (eR + pR) * uR);

	// Œ—–≈ƒÕ≈Õ»≈ œŒ –Œ”
	{
		double rhoLsq = sqrt(rhoL);
		double rhoRsq = sqrt(rhoR);
		double kL = rhoLsq / (rhoLsq + rhoRsq);
		double kR = rhoRsq / (rhoLsq + rhoRsq);

		rhoAV = rhoLsq * rhoRsq;
		
		uAV = kL * uL + kR * uR;
		vAV = kL * vL + kR * vR;
		wAV = kL * wL + kR * wR;
		
		HAV = kL * HL + kR * HR;
		cAV = sqrt(_GAMMA_1_ * (HAV - 0.5 * (uAV*uAV + vAV*vAV + wAV*wAV)));
	}

	Mn = uAV / cAV;

	eigenPLUS  = 0.5 * (_GAMMA_ + 1.0) * uAV / _GAMMA_
		       + sqrt((uAV*0.5/_GAMMA_DIV_GAMMA_1_)*(uAV*0.5/_GAMMA_DIV_GAMMA_1_) + cAV*cAV / _GAMMA_);
	eigenMINUS = 0.5 * (_GAMMA_ + 1.0) * uAV / _GAMMA_
		       - sqrt((uAV*0.5/_GAMMA_DIV_GAMMA_1_)*(uAV*0.5/_GAMMA_DIV_GAMMA_1_) + cAV*cAV / _GAMMA_);

	if( (0.0 <= Mn) && (Mn < 1.0) )
	{
		beta = (uAV + eigenMINUS) / (uAV - eigenMINUS);
		if(beta < 0.0) beta = 0.0;
	}
	else
		if( (-1.0 <= Mn) && (Mn < 0.0) )
		{
			beta = (uAV + eigenPLUS) / (uAV - eigenPLUS);
			if(beta < 0.0) beta = 0.0;
			if(beta > 0.0) beta = - beta;
		}
		else
		{
			if(Mn < 0.0) beta = -1.0;
			else beta = 1.0;
		}

	if(beta == 0.0)
		alfaC = fabs(uAV);
	else
		if( (beta > 0.0) && (0.0 < Mn) && (Mn < 1.0) )
			alfaC = - (1.0 + beta) * eigenMINUS;
		else
			if( (beta < 0.0) && (-1.0 < Mn) && (Mn < 0.0) )
				alfaC = + (1.0 - beta) * eigenPLUS;
			else
				alfaC = 0.0;

	fluxEU[0] -= 0.5 * (alfaC * (rhoR      - rhoL     ) + beta * (rhoR * uR      - rhoL * uL               ));
	fluxEU[1] -= 0.5 * (alfaC * (rhoR * uR - rhoL * uL) + beta * (rhoR * uR * uR - rhoL * uL * uL + pR - pL));
	fluxEU[2] -= 0.5 * (alfaC * (rhoR * vR - rhoL * vL) + beta * (rhoR * vR * uR - rhoL * vL * uL          ));
	fluxEU[3] -= 0.5 * (alfaC * (rhoR * wR - rhoL * wL) + beta * (rhoR * wR * uR - rhoL * wL * uL          ));
	fluxEU[4] -= 0.5 * (alfaC * (rhoR * HR - rhoL * HL) + beta * (rhoR * HR * uR - rhoL * HL * uL          ));
	
	flux[0] =  fluxEU[0];
	flux[1] =  fluxEU[1] * nX + fluxEU[2] * sX + fluxEU[3] * tX;
    flux[2] =  fluxEU[1] * nY + fluxEU[2] * sY + fluxEU[3] * tY;
	flux[3] =  fluxEU[1] * nZ + fluxEU[2] * sZ + fluxEU[3] * tZ;
	flux[4] =  fluxEU[4];

	return 0;
}

#undef _GAMMA_1_
#undef _1_DIV_GAMMA_1_
#undef _GAMMA_DIV_GAMMA_1_
#undef _GAMMA_
