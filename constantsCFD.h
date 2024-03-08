// 0 - ����� ��������
// 1 - ����� HLL / HLLE
// 2 - ����� ROEPIKE
// 3 - ����� HLLC
// 4 - ����� AUSM
// 6 - ����� CUSP
#define _EULER_SOLVER_TYPE_ 0

// ������������� ��� �������������� �� �������
// #define _FIX_TIME_STEP_ 0.00006

// ����� �������
   #define _CFL_ 0.1
// ����� �������������� �� �������
#define _RUNGE_TIME_SCHEME_ 1

// ����� ���� � ���������� ������
   #define _MACH_NUMBER_ 0.1
// �������� �� �������������
   #define _1_DIV_GAMMA_DIV_MACH2_ 71.428571428571428571428571428571
// ����� ����������
   #define _RE_ 500.0
// �������, �������� �� ����� ����������
   #define _MU_ 0.002
   #define _2_DIV_RE_DIV_PR_ 0.00555555555555555555555555555556

// ������ � ���������� ��������� � ������������ � ������� ������������
// 0 - ����/����� � ������������ ���������� �������������� ������
// 1 - ���������
// 2 - ����������
// 3 - ����� �� ������ �� ��������������
// 4 - ����� �� ������ ��������

// ���������� ��������� ����� - ��������� �����
static const int nBSF = 4;
// ��� ���� �����
// static const int BoundaryConditions[4] = {0, 3, 3, 1};
   static const int BoundaryConditions[4] = {0, 3, 3, 1};
// static const int BoundaryConditions[4] = {0, 0, 0, 0};

// ==================================== ��������� ��������� ====================================
#if(_EULER_SOLVER_TYPE_ == 0)
#define EULERSOLVER3D solverRusanov3D
#define EULERSOLVER3DCUDA solverRusanov3DCUDA
#define _USE_SHOCKFIX_ 0
#endif

#if(_EULER_SOLVER_TYPE_ == 6)
#define EULERSOLVER3D solverCUSP3D
#define EULERSOLVER3DCUDA solverCUSP3DCUDA
#define _USE_SHOCKFIX_ 0
#endif

#if(_EULER_SOLVER_TYPE_ == 1)
#define EULERSOLVER2D solverHLL2D
#define _USE_SHOCKFIX_ 0
#endif

#if(_EULER_SOLVER_TYPE_ == 2)
#define EULERSOLVER3D solverRoePike3D
#define EULERSOLVER3DCUDA solverRoePike3DCUDA
#define _USE_SHOCKFIX_ 0
#endif

#if(_EULER_SOLVER_TYPE_ == 3)
#define EULERSOLVER2D solverHLL2D
#define _USE_SHOCKFIX_ 0
#endif

#if(_EULER_SOLVER_TYPE_ == 4)
#define EULERSOLVER2D solverAUSMD2D
#define _USE_SHOCKFIX_ 1
#endif
// ==================================== ��������� ��������� ====================================

// ==================================== ��������������� ��������� ====================================
// �����
#define _GAMMA_ 1.4
// ����� ��� �������
#define _GAMMA_1_ 0.4
// �����, �������� �� ����� ��� �������
#define _GAMMA_DIV_GAMMA_1_ 3.5
// �������, �������� �� ����� ��� �������
#define _1_DIV_GAMMA_1_ 2.5
// ����� ��������
#define _PR_ 0.72
#define _GAMMA_DIV_PR_DIV_GAMMA_1_ 4.8611111111111111
#define _LAMBDA_ _GAMMA_DIV_PR_DIV_GAMMA_1_ 
// ==================================== ��������������� ��������� ====================================
