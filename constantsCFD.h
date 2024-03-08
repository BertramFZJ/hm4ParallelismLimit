// 0 - яуелю псяюмнбю
// 1 - яуелю HLL / HLLE
// 2 - яуелю ROEPIKE
// 3 - яуелю HLLC
// 4 - яуелю AUSM
// 6 - яуелю CUSP
#define _EULER_SOLVER_TYPE_ 0

// тхйяхпнбюммши ьюц хмрецпхпнбюмхъ он бпелемх
// #define _FIX_TIME_STEP_ 0.00006

// вхякн йспюмрю
   #define _CFL_ 0.1
// яуелю хмрецпхпнбюмхъ он бпелемх
#define _RUNGE_TIME_SCHEME_ 1

// вхякн люую б мюаецючыел онрнйе
   #define _MACH_NUMBER_ 0.1
// дюбкемхе мю аеяйнмевмнярх
   #define _1_DIV_GAMMA_DIV_MACH2_ 71.428571428571428571428571428571
// вхякн пеимнкэдяю
   #define _RE_ 500.0
// едхмхжю, декеммюъ мю вхякн пеимнкэдяю
   #define _MU_ 0.002
   #define _2_DIV_RE_DIV_PR_ 0.00555555555555555555555555555556

// люяяхб я цпюмхвмшлх сякнбхълх б яннрберярбхх я лерйюлх онбепумняреи
// 0 - бунд/бшунд я ондярюмнбйни оюпюлерпнб мебнглсыеммнцн онрнйю
// 1 - нрпюфемхе
// 2 - опхкхоюмхе
// 3 - бшунд ян ямнянл он уюпюйрепхярхйе
// 4 - бшунд ян ямнянл гмювемхъ

// днгбсйнбне нарейюмхе ятепш - цхапхдмюъ яерйю
static const int nBSF = 4;
// аег сцкю юрюйх
// static const int BoundaryConditions[4] = {0, 3, 3, 1};
   static const int BoundaryConditions[4] = {0, 3, 3, 1};
// static const int BoundaryConditions[4] = {0, 0, 0, 0};

// ==================================== мюярпнийх пеьюрекеи ====================================
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
// ==================================== мюярпнийх пеьюрекеи ====================================

// ==================================== тсмдюлемрюкэмше йнмярюмрш ====================================
// цюллю
#define _GAMMA_ 1.4
// цюллю аег едхмхжш
#define _GAMMA_1_ 0.4
// цюллю, декеммне мю цюллс аег едхмхжш
#define _GAMMA_DIV_GAMMA_1_ 3.5
// едхмхжю, декеммюъ мю цюллс аег едхмхжш
#define _1_DIV_GAMMA_1_ 2.5
// вхякн опюмдркъ
#define _PR_ 0.72
#define _GAMMA_DIV_PR_DIV_GAMMA_1_ 4.8611111111111111
#define _LAMBDA_ _GAMMA_DIV_PR_DIV_GAMMA_1_ 
// ==================================== тсмдюлемрюкэмше йнмярюмрш ====================================
