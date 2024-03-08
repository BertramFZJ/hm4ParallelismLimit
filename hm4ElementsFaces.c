#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hm4ElementsFaces.h"

// люйпня бшвхякемхъ бейрнпмнцн опнхгбедемхъ дбсу бейрнпнб
// vc = va x vb
#define VectorProductMacro(va,vb,vc) { vc[0] = va[1]*vb[2] - va[2]*vb[1]; vc[1] = va[2]*vb[0] - va[0]*vb[2]; vc[2] = va[0]*vb[1] - va[1]*vb[0]; }

// нопедекемхе йннпдхмюр жемрпю люяя рпесцнкэмхйю
int hm4TriangleCenter(double *P0, double *P1, double *P2, double *TC)
{
	TC[0] = (P0[0] + P1[0] + P2[0]) / 3.0;
	TC[1] = (P0[1] + P1[1] + P2[1]) / 3.0;
	TC[2] = (P0[2] + P1[2] + P2[2]) / 3.0;

	return 0;
}

// бшвхякемхе окныюдх рпесцнкэмхйю
double hm4TriangleSquare(double *P1, double *P2, double *P3)
{
	double S = 0.0;
	double A, B, C;

	A = P1[1]*(P2[2]-P3[2]) - P2[1]*(P1[2]-P3[2]) + P3[1]*(P1[2]-P2[2]);
	A *= A;

	B = P1[2]*(P2[0]-P3[0]) - P2[2]*(P1[0]-P3[0]) + P3[2]*(P1[0]-P2[0]);
	B *= B;

	C = P1[0]*(P2[1]-P3[1]) - P2[0]*(P1[1]-P3[1]) + P3[0]*(P1[1]-P2[1]);
	C *= C;

	S = (A + B + C) / 4.0;
	S = sqrt(S);

	return S;
}

// нопедекемхе йннпдхмюр жемрпю люяя вершпеусцнкэмхйю
int hm4QuadCenter(double *P0, double *P1, double *P2, double *P3, double *QC)
{
	double dx, dy, dz;
	double v1, v2, v3, v4;

	v1 = hm4TriangleSquare(P0, P1, P2);
	dx = v1 * (P0[0] + P1[0] + P2[0]) / 3.0;
	dy = v1 * (P0[1] + P1[1] + P2[1]) / 3.0;
	dz = v1 * (P0[2] + P1[2] + P2[2]) / 3.0;

	v2 = hm4TriangleSquare(P0, P3, P2);
	dx += v2 * (P0[0] + P3[0] + P2[0]) / 3.0;
	dy += v2 * (P0[1] + P3[1] + P2[1]) / 3.0;
	dz += v2 * (P0[2] + P3[2] + P2[2]) / 3.0;

	v3 = hm4TriangleSquare(P1, P2, P3);
	dx += v3 * (P1[0] + P2[0] + P3[0]) / 3.0;
	dy += v3 * (P1[1] + P2[1] + P3[1]) / 3.0;
	dz += v3 * (P1[2] + P2[2] + P3[2]) / 3.0;

	v4 = hm4TriangleSquare(P1, P0, P3);
	dx += v4 * (P1[0] + P0[0] + P3[0]) / 3.0;
	dy += v4 * (P1[1] + P0[1] + P3[1]) / 3.0;
	dz += v4 * (P1[2] + P0[2] + P3[2]) / 3.0;

	QC[0] = dx / (v1 + v2 + v3 + v4);
	QC[1] = dy / (v1 + v2 + v3 + v4);
	QC[2] = dz / (v1 + v2 + v3 + v4);

	return 0;
}

// еярэ рпесцнкэмюъ цпюмэ мейнрнпнцн щкелемрю я хмдейянл iE
// мслепюжхъ бепьхм рпесцнкэмхйю - опнрхб вюянбни ярпекйх опх бгцкъде мю щкелемр "ямюпсфх"
// ондопнцпюллю бшвхякъер бейрнп окныюдх цпюмх face
// бейрнп окныюдх цпюмх - бейрнп бмеьмеи мнплюкх, рн еярэ, сякнбмн, хяундхр хг жемрпю люяя щкелемрю iE б
// мюопюбкемхх жемрпю люяя яняедмецн он пюяялюрпхбюелни цпюмх щкелемрю
int hm4TriangleFace(double *P0, double *P1, double *P2, double *face)
{
	int i;	
	double N[3], a[3], b[3];

	for(i=0; i<3; i++) { a[i] = P1[i] - P0[i]; b[i] = P2[i] - P0[i]; }
	VectorProductMacro(a,b,N)
	for(i=0; i<3; i++) face[i] = N[i] / 2.0;

	return 0;
}

// еярэ вершпеусцнкэмюъ цпюмэ мейнрнпнцн щкелемрю я хмдейянл iE
// мслепюжхъ бепьхм вершпеусцнкэмхйю - опнрхб вюянбни ярпекйх опх бгцкъде мю щкелемр "ямюпсфх"
// ондопнцпюллю бшвхякъер бейрнп окныюдх цпюмх face
// бейрнп окныюдх цпюмх - бейрнп бмеьмеи мнплюкх, рн еярэ, сякнбмн, хяундхр хг жемрпю люяя щкелемрю iE б
// мюопюбкемхх жемрпю люяя яняедмецн он пюяялюрпхбюелни цпюмх щкелемрю

// !!!!!
// б дюммнл яксвюе бейрнп окныюдх бшвхякъеряъ он бейрнпмнлс опнхгбедемхч дхюцнмюкэмшу бейрнпнб
// врн асдер, еякх вершпе рнвйх ме кефюр б ндмни окняйнярх?
int hm4QuadFace(double *P0, double *P1, double *P2, double *P3, double *face)
{
	int i;	
	double a[3], b[3], c[3];

	for(i=0; i<3; i++) { a[i] = P1[i] - P3[i]; b[i] = P2[i] - P0[i]; }
	VectorProductMacro(a,b,c)
	for(i=0; i<3; i++) face[i] = c[i] / 2.0;

	return 0;
}

// опюбюъ яхярелю йннпдхмюр, онярпнеммюъ он едхмхвмнлс бейрнпс бмеьмеи мнплюкх й цпюмх
int hm4FaceCoordinateSystem(double *n, double *s, double *t)
{
    if( (n[0] == 1.0) && (n[1] == 0.0) && (n[2] == 0.0) )
	{
		s[0] = 0.0; s[1] = 1.0; s[2] = 0.0;
		t[0] = 0.0; t[1] = 0.0; t[2] = 1.0;
	}
	else
        if( (n[0] == -1.0) && (n[1] == 0.0) && (n[2] == 0.0) )
		{
			s[0] = 0.0; s[1] = -1.0; s[2] = 0.0;
			t[0] = 0.0; t[1] = 0.0; t[2] = -1.0;
		}
		else
		{
			double zn = sqrt(n[1]*n[1] + n[2]*n[2]);

			s[0] = -zn;
			s[1] = n[0]*n[1] / zn;
			s[2] = n[0]*n[2] / zn;
			t[0] = 0.0;
			t[1] = (-n[2]) / zn;
			t[2] = n[1] / zn;
		}

	return 0;
}

// хмхжхюкхгюжхъ вхякю (*fN) х яохяйю бепьхм vF[(*fN)] наыеи цпюмх щкелемрнб
// рхонб (он вхякс бепьхм) eL х eR ян яохяйюлх бепьхм evL[eL] х evR[eR], яннрберярбеммн
// хмдейяш бепьхм наыеи цпюмх оепевхякъчряъ б онпъдйе опнрхб вюянбни ярпекйх опх
// бмеьмел бгцкъде мю щкелемр eL
// опедонкюцюеряъ, врн дбю щкелемрю наъгюрекэмн хлечр наысч цпюмэ
// нрясрярбхе наыеи цпюмх рпюйрсеряъ ондопнцпюллни йюй ньхайю
int hm4ElementsSameFace(int eL, int *evL, int eR, int *evR, int *fN, int *vF)
{
	int eLSameNodes[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	
	int i, j;

	for(i=0; i<eL; i++)
	{
		for(j=0; j<eR; j++)
		{
			if(evL[i] == evR[j])
			{
				eLSameNodes[i] = 1;
				break;
			} // if
		} // for j
	} // for i

	if(eL == 4)
	{
		int index[4][3] = { {1,  0,  2} , {0,  1,  3}, {1,  2,  3}, {2,  0,  3} };

		for(i=0; i<4; i++)
		{
			if(eLSameNodes[index[i][0]] + eLSameNodes[index[i][1]] + eLSameNodes[index[i][2]] == 3)
			{
				fN[0] = 3;
				for(j=0; j<3; j++) vF[j] = evL[index[i][j]];
				return 0;
			} // if
		} // for i
	}

	if(eL == 5)
	{
		int faceL[5]    =   {4, 3, 3, 3, 3};
		int index[5][4] = { {0,  2,  3,  1}, {0,  1,  4, -1}, {1,  3,  4, -1}, {3,  2,  4, -1}, {2,  0,  4, -1} };

		for(i=0; i<5; i++)
		{
			int nSAME = 0;

			for(j=0; j<faceL[i]; j++) nSAME += eLSameNodes[index[i][j]];
			
			if(nSAME == faceL[i])
			{
				fN[0] = faceL[i];
				for(j=0; j<faceL[i]; j++) vF[j] = evL[index[i][j]];
				return 0;
			} // if
		} // for i
	}

	if(eL == 6)
	{
		int faceL[5]    =   {4, 4, 4, 3, 3};
		int index[5][4] = { {0,  1,  4,  3}, {1,  2,  5,  4}, {2,  0,  3,  5}, {0,  2,  1, -1}, {3,  4,  5, -1} };

		for(i=0; i<5; i++)
		{
			int nSAME = 0;

			for(j=0; j<faceL[i]; j++) nSAME += eLSameNodes[index[i][j]];
			
			if(nSAME == faceL[i])
			{
				fN[0] = faceL[i];
				for(j=0; j<faceL[i]; j++) vF[j] = evL[index[i][j]];
				return 0;
			} // if
		} // for i
	}

	if(eL == 8)
	{
		int index[6][4] = { {0,  1,  5,  4}, {1,  3,  7,  5}, {3,  2,  6,  7}, {2,  0,  4,  6}, {1,  0,  2,  3}, {4,  5,  7,  6} };

		for(i=0; i<6; i++)
		{			
			if(eLSameNodes[index[i][0]] + eLSameNodes[index[i][1]] + eLSameNodes[index[i][2]] + eLSameNodes[index[i][3]] == 4)
			{
				fN[0] = 4;
				for(j=0; j<4; j++) vF[j] = evL[index[i][j]];
				return 0;
			} // if
		} // for i
	}

	fprintf(stderr, "hm4ElementsSameFace --> *** ERROR *** SAME FACE NOT FOUND\n");
	fprintf(stderr, " LEFT %d: ", eL); for(i=0; i<eL; i++) fprintf(stderr, "%d ", evL[i]); fprintf(stderr, "\n");
	fprintf(stderr, "RIGHT %d: ", eR); for(i=0; i<eR; i++) fprintf(stderr, "%d ", evR[i]); fprintf(stderr, "\n");
	exit(0);

	return -1;
}

// хмхжхюкхгюжхъ вхякю (*fN) х яохяйю бепьхм vF[(*fN)]
// дкъ цпюмх я мнлепнл indexFace щкелемрю рхою (он вхякс бепьхм) eL ян яохяйнл бепьхм evL[eL]
int hm4ElementFaceNodes(int eL, int *evL, int indexFace, int *fN, int *vF)
{
	int i;

	if(eL == 4)
	{
		int index[4][3] = { {1,  0,  2} , {0,  1,  3}, {1,  2,  3}, {2,  0,  3} };

		fN[0] = 3;
		for(i=0; i<3; i++) vF[i] = evL[index[indexFace][i]];
	}

	if(eL == 5)
	{
		int faceL[5]    =   {4, 3, 3, 3, 3};
		int index[5][4] = { {0,  2,  3,  1}, {0,  1,  4, -1}, {1,  3,  4, -1}, {3,  2,  4, -1}, {2,  0,  4, -1} };

		fN[0] = faceL[indexFace];
		for(i=0; i<faceL[indexFace]; i++) vF[i] = evL[index[indexFace][i]];		
	}

	if(eL == 6)
	{
		int faceL[5]    =   {4, 4, 4, 3, 3};
		int index[5][4] = { {0,  1,  4,  3}, {1,  2,  5,  4}, {2,  0,  3,  5}, {0,  2,  1, -1}, {3,  4,  5, -1} };

		fN[0] = faceL[indexFace];
		for(i=0; i<faceL[indexFace]; i++) vF[i] = evL[index[indexFace][i]];
	}

	if(eL == 8)
	{
		int index[6][4] = { {0,  1,  5,  4}, {1,  3,  7,  5}, {3,  2,  6,  7}, {2,  0,  4,  6}, {1,  0,  2,  3}, {4,  5,  7,  6} };
		
		fN[0] = 4;
		for(i=0; i<4; i++) vF[i] = evL[index[indexFace][i]];		
	}	

	return 0;
}

#undef VectorProductMacro
