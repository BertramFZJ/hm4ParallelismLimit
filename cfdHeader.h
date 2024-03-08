#ifndef _gtypeInteriorCellFace_
typedef struct
{
	int eL, eR;
	double fMC[3];
	double n[3];
	double s[3];
	double t[3];
	double fS;	
} gtypeInteriorCellFace;
#define _gtypeInteriorCellFace_
#endif

#ifndef _gtypeBoundaryCellFace_
typedef struct
{
	int eL;
	int faceIndex;
	int surfaceIndex;
	int boundaryConditionType;
	double fMC[3];
	double n[3];
	double s[3];
	double t[3];
	double fS;	
} gtypeBoundaryCellFace;
#define _gtypeBoundaryCellFace_
#endif

#ifndef _gtypeHm4CfdHostCore_
typedef struct
{
	int nE;
	int nEinterior;
	double *cellCEL;
	double *cellVOL;
	double *cellH;
	double *cellQ;
	double *cellQph;

	int nICF;
	gtypeInteriorCellFace *ICF;
	int nBCF;
	gtypeBoundaryCellFace *BCF;
	double *faceCE;

	int *faceToElementsX;
	int *faceToElementsA;
	
	int *gaussX;
	int *gaussA;
	double *kXYZC;
	double *kXYZS;	
	double *celldQph;

	int levelMax;
	int *cellLevels;
	int *icfOffset;
	int *bcfOffset;

} gtypeHm4CfdHostCore;
#define _gtypeHm4CfdHostCore_
#endif
