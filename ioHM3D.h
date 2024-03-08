#define HMSH_STR_LENGTH 81

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

