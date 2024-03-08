#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ioHM3D.h"

// fname --> указатель на строку с именем создаваемого файла
// mname --> указатель на строку с названием сетки
// SP[0] --> Np
// SP[1] --> Ne
// SP[2] --> Ne4
// SP[3] --> Ne5
// SP[4] --> Ne6
// SP[5] --> Ne8
// SP[6] --> Nbc
// SP[7] --> Neb
// SP[8] --> Neb3
// SP[9] --> Neb4
// C --> координаты узлов
// ET --> типы элементов (массив смещений из Ne + 1 элемента)
// EV --> массив со списками вершин элементов. Ёлемент i содержит вершины от EV[ET[i]] до EV[ET[i+1]-1] включительно
// BNAME --> указатель на массив с описани€ми (названи€ми) граничных поверхностей
// BE --> описание граничных элементов, где 
// BE[i*3 + 0] - номер объемного элемента
// BE[i*3 + 1] - номер грани объемного элемента
// BE[i*3 + 2] - метка соответствующей поверхности
int WriteMesh_HMSH(char *fname, char *mname, int *SP, double *C, int *ET, int *EV, char *BNAME, int *BE)
{
	FILE *file = NULL;
	int i, j;

	file = fopen(fname, "w"); 
	if(file == NULL) { fprintf(stderr, "WriteMesh_HMSH: ERROR --> Can't create data file \"%s\"\n", fname); exit(0); }

	fprintf(file, "%s\n", mname);
	
	for(i=0; i<10; i++) fprintf(file, "%9d ", SP[i]); fprintf(file, "\n");
	
	for(i=0; i<SP[6]; i++) fprintf(file, "%s\n", BNAME + i*HMSH_STR_LENGTH);
	
	for(i=0; i<SP[0]; i++) fprintf(file, "%21.12e %21.12e %21.12e\n", C[i*3 + 0], C[i*3 + 1], C[i*3 + 2]);
	
	for(i=0; i<SP[1]; i++)
	{
		int TYPE = ET[i+1] - ET[i];
		fprintf(file, "%d ", TYPE);
		for(j=ET[i]; j<ET[i+1]; j++) fprintf(file, "%10d ", EV[j]);
		fprintf(file, "\n");
	} // for i

	for(i=0; i<SP[7]; i++) fprintf(file, "%10d %d %3d\n", BE[i*3 + 0], BE[i*3 + 1], BE[i*3 + 2]);

	fclose(file);

	return 0;
}

int ReadMesh_HMSH(char *fname, char **mname, int **SP, double **C, int **ET, int **EV, char **BNAME, int **BE)
{
	FILE *file = NULL;
	char strbuf[1024], *oldstr, *newstr, *getbuf;
	int i, j, k;

	file = fopen(fname, "r"); 
	if(file == NULL) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> Can't open data file \"%s\"\n", fname); exit(0); }

	(*mname) = NULL; (*mname) = (char *)malloc(HMSH_STR_LENGTH * sizeof(char));
	if((*mname) == NULL) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> Can't allocate memory\n"); exit(0); }
	for(i=0; i<HMSH_STR_LENGTH; i++) (*mname)[i] = '\0';

	if( feof(file) ) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> End of file\n"); exit(0); }
	getbuf = fgets(strbuf, 1024, file); if(getbuf == NULL) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> fgets return NULL pointer\n"); exit(0); }
	k = (int)strlen(strbuf); strbuf[k-1] = '\0'; k -= 1; if(k > HMSH_STR_LENGTH - 1) k = HMSH_STR_LENGTH - 1;
	for(i=0; i<k; i++) (*mname)[i] = strbuf[i];
	
	if(0) printf("MESH: %s\n", (*mname));

	if( feof(file) ) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> End of file\n"); exit(0); }
	getbuf = fgets(strbuf, 1024, file); if(getbuf == NULL) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> fgets return NULL pointer\n"); exit(0); }
	oldstr = strbuf; newstr = NULL;
	(*SP) = NULL; (*SP) = (int *)malloc(10 * sizeof(int));
	if((*SP) == NULL) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> Can't allocate memory\n"); exit(0); }
	for(i=0; i<10; i++) { (*SP)[i] = (int)strtol(oldstr, &newstr, 10); oldstr = newstr; }

	(*BNAME) = NULL; (*BNAME) = (char *)malloc(HMSH_STR_LENGTH * (*SP)[6] * sizeof(char));
	if((*BNAME) == NULL) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> Can't allocate memory\n"); exit(0); }
	for(i=0; i<HMSH_STR_LENGTH*(*SP)[6]; i++) (*BNAME)[i] = '\0';
	for(i=0; i<(*SP)[6]; i++)
	{
		if( feof(file) ) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> End of file\n"); exit(0); }
		getbuf = fgets(strbuf, 1024, file); if(getbuf == NULL) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> fgets return NULL pointer\n"); exit(0); }

		k = (int)strlen(strbuf); strbuf[k-1] = '\0'; k -= 1; if(k > HMSH_STR_LENGTH - 1) k = HMSH_STR_LENGTH - 1;
		for(j=0; j<k; j++) (*BNAME)[i*HMSH_STR_LENGTH + j] = strbuf[j];
	} // for i

	if(0) printf("\n");
	if(0) printf("BOUNDARY NAMES:\n");
	if(0) for(i=0; i<(*SP)[6]; i++) printf("%2d --> %s\n", i, (*BNAME) + i*HMSH_STR_LENGTH);

	(*C) = NULL; (*C) = (double *)malloc((*SP)[0] * 3 * sizeof(double));
	if((*C) == NULL) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> Can't allocate memory\n"); exit(0); }
	for(i=0; i<(*SP)[0]; i++)
	{
		if( feof(file) ) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> End of file\n"); exit(0); }
		getbuf = fgets(strbuf, 1024, file); if(getbuf == NULL) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> fgets return NULL pointer\n"); exit(0); }
		oldstr = strbuf; newstr = NULL;

		(*C)[i*3 + 0] = strtod(oldstr, &newstr); oldstr = newstr;
		(*C)[i*3 + 1] = strtod(oldstr, &newstr); oldstr = newstr;
		(*C)[i*3 + 2] = strtod(oldstr, &newstr);
	} // for i

	if(0) printf("\nREAD COORDINATES OF %d NODES\n", (*SP)[0]);

	k = (*SP)[2]*4 + (*SP)[3]*5 + (*SP)[4]*6 + (*SP)[5]*8;
	(*ET) = NULL; (*ET) = (int *)malloc(((*SP)[1] + 1) * sizeof(int));
	if((*ET) == NULL) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> Can't allocate memory\n"); exit(0); }
	(*ET)[0] = 0;
	(*EV) = NULL; (*EV) = (int *)malloc(k * sizeof(int));
	if((*EV) == NULL) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> Can't allocate memory\n"); exit(0); }
	k = 0;
	for(i=0; i<(*SP)[1]; i++)
	{
		if( feof(file) ) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> End of file\n"); exit(0); }
		getbuf = fgets(strbuf, 1024, file); if(getbuf == NULL) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> fgets return NULL pointer\n"); exit(0); }
		oldstr = strbuf; newstr = NULL;

		(*ET)[i+1] = (int)strtol(oldstr, &newstr, 10); oldstr = newstr;

		for(j=0; j<(*ET)[i+1]; j++)
		{
			(*EV)[k] = (int)strtol(oldstr, &newstr, 10); oldstr = newstr;
			k += 1;
		} // for j
	} // for i
	for(i=2; i<=(*SP)[1]; i++) (*ET)[i] += (*ET)[i-1];
	if(k != (*SP)[2]*4 + (*SP)[3]*5 + (*SP)[4]*6 + (*SP)[5]*8) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> Wrong final index EV\n"); exit(0); }
	if(k != (*ET)[(*SP)[1]]) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> Wrong final index ET\n"); exit(0); }

	if(0) printf("\n");
	if(0) printf("READ ET + EV ARRAYS FOR %d ELEMENTS\n", (*SP)[1]);

	(*BE) = NULL; (*BE) = (int *)malloc((*SP)[7] * 3 * sizeof(int));
	if((*BE) == NULL) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> Can't allocate memory\n"); exit(0); }
	for(i=0; i<(*SP)[7]; i++)
	{
		if( feof(file) ) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> End of file\n"); exit(0); }
		getbuf = fgets(strbuf, 1024, file); if(getbuf == NULL) { fprintf(stderr, "ReadMesh_HMSH: ERROR --> fgets return NULL pointer\n"); exit(0); }
		oldstr = strbuf; newstr = NULL;

		(*BE)[i*3 + 0] = (int)strtol(oldstr, &newstr, 10); oldstr = newstr;
		(*BE)[i*3 + 1] = (int)strtol(oldstr, &newstr, 10); oldstr = newstr;
		(*BE)[i*3 + 2] = (int)strtol(oldstr, &newstr, 10);
	} // for i

	if(0) printf("\n");
	if(0) printf("READ PARAMETERS OF %d BOUNDARY EDGES\n", (*SP)[7]);

	fclose(file);

	return 0;
}

int WriteMesh_HMSB(char *fname, char *mname, int *SP, double *C, int *ET, int *EV, char *BNAME, int *BE)
{
	FILE *file = NULL;
	int get;
	
	file = fopen(fname, "wb"); 
	if(file == NULL) { fprintf(stderr, "WriteMesh_HMSB: ERROR --> Can't create data file \"%s\"\n", fname); exit(0); }

	get = (int)fwrite(mname, sizeof(char), HMSH_STR_LENGTH, file);
	if(get != HMSH_STR_LENGTH) { fprintf(stderr, "WriteMesh_HMSB: ERROR --> Can't write data to file %s\n", fname); exit(0); }
	
	get = (int)fwrite(SP, sizeof(int), 10, file);
	if(get != 10) { fprintf(stderr, "WriteMesh_HMSB: ERROR --> Can't write data to file %s\n", fname); exit(0); }
	
	get = (int)fwrite(BNAME, sizeof(char), SP[6]*HMSH_STR_LENGTH, file);
	if(get != SP[6]*HMSH_STR_LENGTH) { fprintf(stderr, "WriteMesh_HMSB: ERROR --> Can't write data to file %s\n", fname); exit(0); }

	get = (int)fwrite(C, sizeof(double), SP[0]*3, file);
	if(get != SP[0]*3) { fprintf(stderr, "WriteMesh_HMSB: ERROR --> Can't write data to file %s\n", fname); exit(0); }

	get = (int)fwrite(ET, sizeof(int), SP[1]+1, file);
	if(get != SP[1]+1) { fprintf(stderr, "WriteMesh_HMSB: ERROR --> Can't write data to file %s\n", fname); exit(0); }

	get = (int)fwrite(EV, sizeof(int), ET[SP[1]], file);
	if(get != ET[SP[1]]) { fprintf(stderr, "WriteMesh_HMSB: ERROR --> Can't write data to file %s\n", fname); exit(0); }

	get = (int)fwrite(BE, sizeof(int), SP[7]*3, file);
	if(get != SP[7]*3) { fprintf(stderr, "WriteMesh_HMSB: ERROR --> Can't write data to file %s\n", fname); exit(0); }	
	
	fclose(file);

	return 0;
}

int ReadMesh_HMSB(char *fname, char **mname, int **SP, double **C, int **ET, int **EV, char **BNAME, int **BE)
{
	FILE *file = NULL;
	int i;
	int get;
	
	file = fopen(fname, "rb"); 
	if(file == NULL) { fprintf(stderr, "ReadMesh_HMSB: ERROR --> Can't open data file \"%s\"\n", fname); exit(0); }

	(*mname) = NULL; (*mname) = (char *)malloc(HMSH_STR_LENGTH * sizeof(char));
	if((*mname) == NULL) { fprintf(stderr, "ReadMesh_HMSB: ERROR --> Can't allocate memory\n"); exit(0); }
	for(i=0; i<HMSH_STR_LENGTH; i++) (*mname)[i] = '\0';
	get = (int)fread((*mname), sizeof(char), HMSH_STR_LENGTH, file);
	if(get != HMSH_STR_LENGTH) { fprintf(stderr, "ReadMesh_HMSB: ERROR --> Can't read data from file %s\n", fname); exit(0); }
	
	(*SP) = NULL; (*SP) = (int *)malloc(10 * sizeof(int));
	if((*SP) == NULL) { fprintf(stderr, "ReadMesh_HMSB: ERROR --> Can't allocate memory\n"); exit(0); }
	get = (int)fread((*SP), sizeof(int), 10, file);
	if(get != 10) { fprintf(stderr, "ReadMesh_HMSB: ERROR --> Can't read data from file %s\n", fname); exit(0); }

	(*BNAME) = NULL; (*BNAME) = (char *)malloc(HMSH_STR_LENGTH * (*SP)[6] * sizeof(char));
	if((*BNAME) == NULL) { fprintf(stderr, "ReadMesh_HMSB: ERROR --> Can't allocate memory\n"); exit(0); }
	for(i=0; i<HMSH_STR_LENGTH*(*SP)[6]; i++) (*BNAME)[i] = '\0';	
	get = (int)fread((*BNAME), sizeof(char), (*SP)[6]*HMSH_STR_LENGTH, file);
	if(get != (*SP)[6]*HMSH_STR_LENGTH) { fprintf(stderr, "ReadMesh_HMSB: ERROR --> Can't read data from file %s\n", fname); exit(0); }

	(*C) = NULL; (*C) = (double *)malloc((*SP)[0] * 3 * sizeof(double));
	if((*C) == NULL) { fprintf(stderr, "ReadMesh_HMSB: ERROR --> Can't allocate memory\n"); exit(0); }
	get = (int)fread((*C), sizeof(double), (*SP)[0]*3, file);
	if(get != (*SP)[0]*3) { fprintf(stderr, "ReadMesh_HMSB: ERROR --> Can't read data from file %s\n", fname); exit(0); }

	(*ET) = NULL; (*ET) = (int *)malloc(((*SP)[1] + 1) * sizeof(int));
	if((*ET) == NULL) { fprintf(stderr, "ReadMesh_HMSB: ERROR --> Can't allocate memory\n"); exit(0); }
	get = (int)fread((*ET), sizeof(int), (*SP)[1]+1, file);
	if(get != (*SP)[1]+1) { fprintf(stderr, "ReadMesh_HMSB: ERROR --> Can't read data from file %s\n", fname); exit(0); }

	(*EV) = NULL; (*EV) = (int *)malloc((*ET)[(*SP)[1]] * sizeof(int));
	if((*EV) == NULL) { fprintf(stderr, "ReadMesh_HMSB: ERROR --> Can't allocate memory\n"); exit(0); }
	get = (int)fread((*EV), sizeof(int), (*ET)[(*SP)[1]], file);
	if(get != (*ET)[(*SP)[1]]) { fprintf(stderr, "ReadMesh_HMSB: ERROR --> Can't read data from file %s\n", fname); exit(0); }

	(*BE) = NULL; (*BE) = (int *)malloc((*SP)[7] * 3 * sizeof(int));
	if((*BE) == NULL) { fprintf(stderr, "ReadMesh_HMSB: ERROR --> Can't allocate memory\n"); exit(0); }
	get = (int)fread((*BE), sizeof(int), (*SP)[7]*3, file);
	if(get != (*SP)[7]*3) { fprintf(stderr, "ReadMesh_HMSB: ERROR --> Can't read data from file %s\n", fname); exit(0); }

	fclose(file);

	return 0;
}

int IOHM3D_PlotHm4MeshPart(int nP, double *c, int nE, int *x, int *a, int *part, char *fname)
{
	FILE *file;
	int *id;

	int tetrahedraMesh = 1;
	int partMIN, partMAX;
	int i, j, k;

	partMIN = partMAX = part[0];
	if(x[1] - x[0] != 4) tetrahedraMesh = 0;
	for(i=1; i<nE; i++)
	{
		if(partMIN > part[i]) partMIN = part[i];
		if(partMAX < part[i]) partMAX = part[i];
		if(x[i + 1] - x[i] != 4) tetrahedraMesh = 0;
	} // for i

	id = (int *)malloc(nP * sizeof(int)); if(id == NULL) exit(0);
	file = fopen(fname, "w"); if(file == NULL) exit(0);

	fprintf(file, "VARIABLES= \"X\",\"Y\",\"Z\"\n");

	for(i=partMIN; i<=partMAX; i++)
	{
		int nPL = 0;
		int nEL = 0;

		for(j=0; j<nP; j++) id[j] = -1;
		for(j=0; j<nE; j++) if(part[j] == i)
		{
			nEL += 1;
			for(k=x[j]; k<x[j+1]; k++) id[a[k]] = 0;
		} // if
		for(j=0; j<nP; j++) if(id[j] == 0) id[j] = nPL, nPL += 1;

		if( (nPL != 0) && (nEL != 0) )
		{
			fprintf(file, "ZONE T=\"SUBDOMEN %d\"\n", i);
			if(tetrahedraMesh == 1) fprintf(file, "F=FEPOINT, ET=TETRAHEDRON, N=%d E=%d\n", nPL, nEL);
			if(tetrahedraMesh == 0) fprintf(file, "F=FEPOINT, ET=BRICK, N=%d E=%d\n", nPL, nEL);
			
			for(j=0; j<nP; j++) if(id[j] >= 0) fprintf(file, "%g %g %g\n", c[j*3 + 0], c[j*3 + 1], c[j*3 + 2]);			

			if(tetrahedraMesh == 0)
			for(j=0; j<nE; j++) if(part[j] == i)
			{
				int TYPE, *PTR, PE[8];

				TYPE = x[j+1] - x[j];
				PTR = a + x[j];
				for(k=0; k<TYPE; k++) PE[k] = id[PTR[k]];

				if(TYPE == 4)
				{
					int pos[8] = {0, 2, 3, 3, 1, 1, 1, 1};
					fprintf(file, "%d %d %d %d %d %d %d %d\n", PE[pos[0]]+1, PE[pos[1]]+1, PE[pos[2]]+1, PE[pos[3]]+1, PE[pos[4]]+1, PE[pos[5]]+1, PE[pos[6]]+1, PE[pos[7]]+1);
				}

				if(TYPE == 5)
				{
					int pos[8] = {0, 1, 3, 2, 4, 4, 4, 4};
					fprintf(file, "%d %d %d %d %d %d %d %d\n", PE[pos[0]]+1, PE[pos[1]]+1, PE[pos[2]]+1, PE[pos[3]]+1, PE[pos[4]]+1, PE[pos[5]]+1, PE[pos[6]]+1, PE[pos[7]]+1);
				}

				if(TYPE == 6)
				{
					int pos[8] = {0, 2, 5, 3, 1, 1, 4, 4};
					fprintf(file, "%d %d %d %d %d %d %d %d\n", PE[pos[0]]+1, PE[pos[1]]+1, PE[pos[2]]+1, PE[pos[3]]+1, PE[pos[4]]+1, PE[pos[5]]+1, PE[pos[6]]+1, PE[pos[7]]+1);
				}

				if(TYPE == 8)
				{
					int pos[8] = {4, 5, 1, 0, 6, 7, 3, 2};
					fprintf(file, "%d %d %d %d %d %d %d %d\n", PE[pos[0]]+1, PE[pos[1]]+1, PE[pos[2]]+1, PE[pos[3]]+1, PE[pos[4]]+1, PE[pos[5]]+1, PE[pos[6]]+1, PE[pos[7]]+1);
				}
			} // for j

			if(tetrahedraMesh == 1)
			for(j=0; j<nE; j++) if(part[j] == i)
			{
				int TYPE, *PTR, PE[8];

				TYPE = x[j+1] - x[j];
				PTR = a + x[j];
				for(k=0; k<TYPE; k++) PE[k] = id[PTR[k]];
				fprintf(file, "%d %d %d %d\n", PE[0]+1, PE[1]+1, PE[2]+1, PE[3]+1);
			} // for j
		} // if
	} // for i

	fclose(file);
	free(id);

	return 0;
}

int IOHM3D_PlotHm4Mesh(int nP, double *c, int nE, int *x, int *a, char *fname)
{
	FILE *file = fopen(fname, "w");
	
	int tetrahedraMesh = 1;
	int j;

	if(file == NULL) { fprintf(stderr, "Can't create file %s\n", fname); exit(0); }

	if(x[1] - x[0] != 4) tetrahedraMesh = 0;
	for(j=1; j<nE; j++) if(x[j + 1] - x[j] != 4) { tetrahedraMesh = 0; break; }
    tetrahedraMesh = 0;

	fprintf(file, "VARIABLES= \"X\",\"Y\",\"Z\"\n");
	fprintf(file, "ZONE T=\"HYBRID MESH\"\n");
	if(tetrahedraMesh == 1) fprintf(file, "F=FEPOINT, ET=TETRAHEDRON, N=%d E=%d\n", nP, nE);
	if(tetrahedraMesh == 0) fprintf(file, "F=FEPOINT, ET=BRICK, N=%d E=%d\n", nP, nE);

	for(j=0; j<nP; j++) fprintf(file, "%g %g %g\n", c[j*3 + 0], c[j*3 + 1], c[j*3 + 2]);

	if(tetrahedraMesh == 0)
	for(j=0; j<nE; j++)
	{
		int TYPE, *PTR;

		TYPE = x[j+1] - x[j];
		PTR = a + x[j];

		if(TYPE == 4)
		{
			int pos[8] = {0, 2, 3, 3, 1, 1, 1, 1};
			fprintf(file, "%d %d %d %d %d %d %d %d\n", PTR[pos[0]]+1, PTR[pos[1]]+1, PTR[pos[2]]+1, PTR[pos[3]]+1, PTR[pos[4]]+1, PTR[pos[5]]+1, PTR[pos[6]]+1, PTR[pos[7]]+1);
		}

		if(TYPE == 5)
		{
			int pos[8] = {0, 1, 3, 2, 4, 4, 4, 4};
			fprintf(file, "%d %d %d %d %d %d %d %d\n", PTR[pos[0]]+1, PTR[pos[1]]+1, PTR[pos[2]]+1, PTR[pos[3]]+1, PTR[pos[4]]+1, PTR[pos[5]]+1, PTR[pos[6]]+1, PTR[pos[7]]+1);
		}

		if(TYPE == 6)
		{
			int pos[8] = {0, 2, 5, 3, 1, 1, 4, 4};
			fprintf(file, "%d %d %d %d %d %d %d %d\n", PTR[pos[0]]+1, PTR[pos[1]]+1, PTR[pos[2]]+1, PTR[pos[3]]+1, PTR[pos[4]]+1, PTR[pos[5]]+1, PTR[pos[6]]+1, PTR[pos[7]]+1);
		}

		if(TYPE == 8)
		{
			int pos[8] = {4, 5, 1, 0, 6, 7, 3, 2};
			fprintf(file, "%d %d %d %d %d %d %d %d\n", PTR[pos[0]]+1, PTR[pos[1]]+1, PTR[pos[2]]+1, PTR[pos[3]]+1, PTR[pos[4]]+1, PTR[pos[5]]+1, PTR[pos[6]]+1, PTR[pos[7]]+1);
		}
	} // for j

	if(tetrahedraMesh == 1)
	{
        for(j=0; j<nE; j++)
        {
            fprintf(file, "%d %d %d %d\n", a[j*4 + 0] + 1, a[j*4 + 1] + 1, a[j*4 + 2] + 1, a[j*4 + 3] + 1);
        }		
	}

	fclose(file);
	
	return 0;
}

int IOHM3D_PlotHm4MeshFunc(int nP, double *c, 
						   int nE, int *x, int *a, 
						   char *fname,
						   int nF, double *func,
						   char *varNames)
{
	FILE *file = fopen(fname, "w");
	
	int tetrahedraMesh = 1;
	int i, j;

	if(file == NULL) { fprintf(stderr, "Can't create file %s\n", fname); exit(0); }

	if(x[1] - x[0] != 4) tetrahedraMesh = 0;
	for(j=1; j<nE; j++) if(x[j + 1] - x[j] != 4) { tetrahedraMesh = 0; break; }	

	fprintf(file, "VARIABLES= \"X\",\"Y\",\"Z\"");
	if(varNames == NULL) for(j=0; j<nF; j++) fprintf(file, ",\"VAR%d\"", j);
	else fprintf(file, "%s", varNames);
	fprintf(file, "\n");

	fprintf(file, "ZONE\n");
	
	if(tetrahedraMesh == 0) fprintf(file, "F=FEBLOCK, ET=BRICK, N=%d E=%d VARLOCATION=([4", nP, nE);
	if(tetrahedraMesh == 1) fprintf(file, "F=FEBLOCK, ET=TETRAHEDRON, N=%d E=%d VARLOCATION=([4", nP, nE);
	for(j=1; j<nF; j++) fprintf(file, ",%d", 4+j);
	fprintf(file, "]=CELLCENTERED)\n");

	for(i=0; i<nP; i++) { fprintf(file, "%g ", c[i*3 + 0]); if((i != 0) && (i%10 == 0) ) fprintf(file, "\n"); } if((i-1)%10 != 0) fprintf(file, "\n");
	for(i=0; i<nP; i++) { fprintf(file, "%g ", c[i*3 + 1]); if((i != 0) && (i%10 == 0) ) fprintf(file, "\n"); } if((i-1)%10 != 0) fprintf(file, "\n");
	for(i=0; i<nP; i++) { fprintf(file, "%g ", c[i*3 + 2]); if((i != 0) && (i%10 == 0) ) fprintf(file, "\n"); } if((i-1)%10 != 0) fprintf(file, "\n");
	for(j=0; j<nF; j++)
	for(i=0; i<nE; i++) { fprintf(file, "%g ", func[i*nF + j]); if((i != 0) && (i%10 == 0) ) fprintf(file, "\n"); } if((i-1)%10 != 0) fprintf(file, "\n");
	
	if(tetrahedraMesh == 0)
	for(j=0; j<nE; j++)
	{
		int TYPE, *PTR;

		TYPE = x[j+1] - x[j];
		PTR = a + x[j];

		if(TYPE == 4)
		{
			int pos[8] = {0, 2, 3, 3, 1, 1, 1, 1};
			fprintf(file, "%d %d %d %d %d %d %d %d\n", PTR[pos[0]]+1, PTR[pos[1]]+1, PTR[pos[2]]+1, PTR[pos[3]]+1, PTR[pos[4]]+1, PTR[pos[5]]+1, PTR[pos[6]]+1, PTR[pos[7]]+1);
		}

		if(TYPE == 5)
		{
			int pos[8] = {0, 1, 3, 2, 4, 4, 4, 4};
			fprintf(file, "%d %d %d %d %d %d %d %d\n", PTR[pos[0]]+1, PTR[pos[1]]+1, PTR[pos[2]]+1, PTR[pos[3]]+1, PTR[pos[4]]+1, PTR[pos[5]]+1, PTR[pos[6]]+1, PTR[pos[7]]+1);
		}

		if(TYPE == 6)
		{
			int pos[8] = {0, 2, 5, 3, 1, 1, 4, 4};
			fprintf(file, "%d %d %d %d %d %d %d %d\n", PTR[pos[0]]+1, PTR[pos[1]]+1, PTR[pos[2]]+1, PTR[pos[3]]+1, PTR[pos[4]]+1, PTR[pos[5]]+1, PTR[pos[6]]+1, PTR[pos[7]]+1);
		}

		if(TYPE == 8)
		{
			int pos[8] = {4, 5, 1, 0, 6, 7, 3, 2};
			fprintf(file, "%d %d %d %d %d %d %d %d\n", PTR[pos[0]]+1, PTR[pos[1]]+1, PTR[pos[2]]+1, PTR[pos[3]]+1, PTR[pos[4]]+1, PTR[pos[5]]+1, PTR[pos[6]]+1, PTR[pos[7]]+1);
		}
	} // for j

	if(tetrahedraMesh == 1)
	for(j=0; j<nE; j++)
	{
		int *PTR = a + j*4;
		fprintf(file, "%d %d %d %d\n", PTR[0] + 1, PTR[1] + 1, PTR[2] + 1, PTR[3] + 1);		
	} // for j

	fclose(file);
	
	return 0;
}

int IOHM3D_PlotHm4MeshFuncVertex(int nP, double *c, 
						         int nE, int *x, int *a, 
						         char *fname,
						         int nF, double *func,
						         char *varNames)
{
	FILE *file = fopen(fname, "w");
	
	int tetrahedraMesh = 1;
	int i, j;

	if(file == NULL) { fprintf(stderr, "Can't create file %s\n", fname); exit(0); }

	if(x[1] - x[0] != 4) tetrahedraMesh = 0;
	for(j=1; j<nE; j++) if(x[j + 1] - x[j] != 4) { tetrahedraMesh = 0; break; }	

	fprintf(file, "VARIABLES= \"X\",\"Y\",\"Z\"");
	if(varNames == NULL) for(j=0; j<nF; j++) fprintf(file, ",\"VAR%d\"", j);
	else fprintf(file, "%s", varNames);
	fprintf(file, "\n");

	fprintf(file, "ZONE\n");
	
	if(tetrahedraMesh == 0) fprintf(file, "F=FEPOINT, ET=BRICK, N=%d E=%d\n", nP, nE);
	if(tetrahedraMesh == 1) fprintf(file, "F=FEPOINT, ET=TETRAHEDRON, N=%d E=%d\n", nP, nE);

	for(i=0; i<nP; i++)
	{
		fprintf(file, "%g %g %g", c[i*3 + 0], c[i*3 + 1], c[i*3 + 2]);
		for(j=0; j<nF; j++) fprintf(file, " %g", func[i*nF + j]);
		fprintf(file, "\n");
	} // for i	
	
	if(tetrahedraMesh == 0)
	for(j=0; j<nE; j++)
	{
		int TYPE, *PTR;

		TYPE = x[j+1] - x[j];
		PTR = a + x[j];

		if(TYPE == 4)
		{
			int pos[8] = {0, 2, 3, 3, 1, 1, 1, 1};
			fprintf(file, "%d %d %d %d %d %d %d %d\n", PTR[pos[0]]+1, PTR[pos[1]]+1, PTR[pos[2]]+1, PTR[pos[3]]+1, PTR[pos[4]]+1, PTR[pos[5]]+1, PTR[pos[6]]+1, PTR[pos[7]]+1);
		}

		if(TYPE == 5)
		{
			int pos[8] = {0, 1, 3, 2, 4, 4, 4, 4};
			fprintf(file, "%d %d %d %d %d %d %d %d\n", PTR[pos[0]]+1, PTR[pos[1]]+1, PTR[pos[2]]+1, PTR[pos[3]]+1, PTR[pos[4]]+1, PTR[pos[5]]+1, PTR[pos[6]]+1, PTR[pos[7]]+1);
		}

		if(TYPE == 6)
		{
			int pos[8] = {0, 2, 5, 3, 1, 1, 4, 4};
			fprintf(file, "%d %d %d %d %d %d %d %d\n", PTR[pos[0]]+1, PTR[pos[1]]+1, PTR[pos[2]]+1, PTR[pos[3]]+1, PTR[pos[4]]+1, PTR[pos[5]]+1, PTR[pos[6]]+1, PTR[pos[7]]+1);
		}

		if(TYPE == 8)
		{
			int pos[8] = {4, 5, 1, 0, 6, 7, 3, 2};
			fprintf(file, "%d %d %d %d %d %d %d %d\n", PTR[pos[0]]+1, PTR[pos[1]]+1, PTR[pos[2]]+1, PTR[pos[3]]+1, PTR[pos[4]]+1, PTR[pos[5]]+1, PTR[pos[6]]+1, PTR[pos[7]]+1);
		}
	} // for j

	if(tetrahedraMesh == 1)
	for(j=0; j<nE; j++)
	{
		int *PTR = a + j*4;
		fprintf(file, "%d %d %d %d\n", PTR[0] + 1, PTR[1] + 1, PTR[2] + 1, PTR[3] + 1);		
	} // for j

	fclose(file);
	
	return 0;
}

