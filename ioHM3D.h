#define HMSH_STR_LENGTH 81

// fname --> ��������� �� ������ � ������ ������������ �����
// mname --> ��������� �� ������ � ��������� �����
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
// C --> ���������� �����
// ET --> ���� ��������� (������ �������� �� Ne + 1 ��������)
// EV --> ������ �� �������� ������ ���������. ������� i �������� ������� �� EV[ET[i]] �� EV[ET[i+1]-1] ������������
// BNAME --> ��������� �� ������ � ���������� (����������) ��������� ������������
// BE --> �������� ��������� ���������, ��� 
// BE[i*3 + 0] - ����� ��������� ��������
// BE[i*3 + 1] - ����� ����� ��������� ��������
// BE[i*3 + 2] - ����� ��������������� �����������
int WriteMesh_HMSH(char *fname, char *mname, int *SP, double *C, int *ET, int *EV, char *BNAME, int *BE);

int ReadMesh_HMSH(char *fname, char **mname, int **SP, double **C, int **ET, int **EV, char **BNAME, int **BE);

int WriteMesh_HMSB(char *fname, char *mname, int *SP, double *C, int *ET, int *EV, char *BNAME, int *BE);

int ReadMesh_HMSB(char *fname, char **mname, int **SP, double **C, int **ET, int **EV, char **BNAME, int **BE);

int IOHM3D_PlotHm4MeshPart(int nP, double *c, int nE, int *x, int *a, int *part, char *fname);

int IOHM3D_PlotHm4Mesh(int nP, double *c, int nE, int *x, int *a, char *fname);

int IOHM3D_PlotHm4MeshFunc(int nP, double *c, 
						   int nE, int *x, int *a, 
						   char *fname,
						   int nF, double *func,
						   char *varNames);

int IOHM3D_PlotHm4MeshFuncVertex(int nP, double *c, 
						         int nE, int *x, int *a, 
						         char *fname,
						         int nF, double *func,
						         char *varNames);

