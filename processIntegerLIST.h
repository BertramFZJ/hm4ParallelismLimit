// ����� ������� ��������� �������� TARGET � ��������������� ��
// ����������� ������� ����� ����� PTR[LENGTH]
// ������������ ���������� ������ ���������������� ��������, ���� �� ������,
// � ��������� -1, ���� ������� TARGET � ������� �����������
int DichotomySearchElementInOrderIntegerListUnit(int *PTR, int LENGTH, int TARGET);
int DichotomySearchElementInOrderIntegerListSizeUnit(int *PTR, int LENGTH, int TARGET, int size);
/*NEW*/ int DichotomySearchElementInOrderIntegerListSizeFunc(int *PTR, int LENGTH, int *TARGET, int size, int (*cmp)(int *, int *));

void wsortIntegerListUnit(int *A, int l, int r);
void wsortIntegerListSizeUnit(int *A, int l, int r, int size);
void wsortIntegerListSizeFunc(int *A, int l, int r, int size, int (*cmp)(int *, int *));

int SortAndClearIntegerListUnit(int oldN, int *list, int *newN);
/*NEW*/ int SortAndClearIntegerListSizeUnit(int oldN, int *list, int size, int *newN);
/*NEW*/ int SortAndClearIntegerListSizeFunc(int oldN, int *list, int size, int (*cmp)(int *, int *), int *newN);

int DeleteListFromOrderedList(int nE, int *list, int oldEO, int *listO, int *newEO);
