// ������ ������������ ������������ ����� ��� ����� ������ � �����
int csrPartParMetisV4NoWeightsMPI(int *eDIST, int *csrX, int *csrA, int nPART, int *part, FILE *outSTREAM, MPI_Comm comm);
// ������ ������������ ������������ ����� � �������������� ������ ������, ������� ����� ����������� �� ������� �����
// ���� ����� ����� �� ���������������
int csrPartParMetisV4NodesWeightsByCsrEdgesMPI(int *eDIST, int *csrX, int *csrA, int nPART, int *part, FILE *outSTREAM, MPI_Comm comm);

// ������ ���������������� ������������ ����� � ������������ ��������� ������������� ����� ������ � �����
// wghtN[Ne] - ������ � ������ ������ ����� (wghtN == NULL - ����� ���)
// wghtE[csrX[Ne]] - ������ � ������ ����� ����� (wghtE == NULL - ����� ���)
int csrPartMetisNodesEdgesWeightsSerial(int nE, int *csrX, int *csrA, int *wghtN, int *wghtE, int nPART, int *part, FILE *outSTREAM);
