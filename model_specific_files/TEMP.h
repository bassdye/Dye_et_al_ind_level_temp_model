//*** TEMP.h   -- 

#ifndef BIFURCATION
#define BIFURCATION             0 
#endif

#if (BIFURCATION == 1)
#define LENGTHCURVES            0
#else
#define LENGTHCURVES            20
#endif

#define POPULATION_NR           1
#define I_STATE_DIM             4
#define I_CONST_DIM             14
#define ENVIRON_DIM             3
#define OUTPUT_VAR_NR           20	// 1st column is time so (1 + # outputs in .c)
#define PARAMETER_NR            59	
#define TIME_METHOD             RKCK
#define EVENT_NR                0
#define DYNAMIC_COHORTS         0
#define ALGEBRAIC_EQS           0
#define CHECK_EXTINCTION        1
	
