/*
FORCES_NLP_solver : A fast customized optimization solver.

Copyright (C) 2013-2020 EMBOTECH AG [info@embotech.com]. All rights reserved.


This software is intended for simulation and testing purposes only. 
Use of this software for any commercial purpose is prohibited.

This program is distributed in the hope that it will be useful.
EMBOTECH makes NO WARRANTIES with respect to the use of the software 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE. 

EMBOTECH shall not have any liability for any damage arising from the use
of the software.

This Agreement shall exclusively be governed by and interpreted in 
accordance with the laws of Switzerland, excluding its principles
of conflict of laws. The Courts of Zurich-City shall have exclusive 
jurisdiction in case of any dispute.

*/

/* Generated by FORCES PRO v3.1.0 on Tuesday, July 28, 2020 at 8:46:33 AM */

#ifndef SOLVER_STDIO_H
#define SOLVER_STDIO_H
#include <stdio.h>
#endif

#ifndef FORCES_NLP_solver_H
#define FORCES_NLP_solver_H

/* DATA TYPE ------------------------------------------------------------*/
typedef double FORCES_NLP_solver_float;

typedef double FORCES_NLP_solverinterface_float;

#ifndef SOLVER_STANDARD_TYPES
#define SOLVER_STANDARD_TYPES

typedef signed char solver_int8_signed;
typedef unsigned char solver_int8_unsigned;
typedef char solver_int8_default;
typedef signed short int solver_int16_signed;
typedef unsigned short int solver_int16_unsigned;
typedef short int solver_int16_default;
typedef signed int solver_int32_signed;
typedef unsigned int solver_int32_unsigned;
typedef int solver_int32_default;
typedef signed long long int solver_int64_signed;
typedef unsigned long long int solver_int64_unsigned;
typedef long long int solver_int64_default;

#endif


/* SOLVER SETTINGS ------------------------------------------------------*/

/* MISRA-C compliance */
#ifndef MISRA_C_FORCES_NLP_solver
#define MISRA_C_FORCES_NLP_solver (0)
#endif

/* restrict code */
#ifndef RESTRICT_CODE_FORCES_NLP_solver
#define RESTRICT_CODE_FORCES_NLP_solver (0)
#endif

/* print level */
#ifndef SET_PRINTLEVEL_FORCES_NLP_solver
#define SET_PRINTLEVEL_FORCES_NLP_solver    (2)
#endif

/* timing */
#ifndef SET_TIMING_FORCES_NLP_solver
#define SET_TIMING_FORCES_NLP_solver    (1)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define SET_MAXIT_FORCES_NLP_solver			(40)	

/* scaling factor of line search (FTB rule) */
#define SET_FLS_SCALE_FORCES_NLP_solver		(FORCES_NLP_solver_float)(0.99)      

/* maximum number of supported elements in the filter */
#define MAX_FILTER_SIZE_FORCES_NLP_solver	(40) 

/* maximum number of supported elements in the filter */
#define MAX_SOC_IT_FORCES_NLP_solver			(4) 

/* desired relative duality gap */
#define SET_ACC_RDGAP_FORCES_NLP_solver		(FORCES_NLP_solver_float)(0.0001)

/* desired maximum residual on equality constraints */
#define SET_ACC_RESEQ_FORCES_NLP_solver		(FORCES_NLP_solver_float)(1E-06)

/* desired maximum residual on inequality constraints */
#define SET_ACC_RESINEQ_FORCES_NLP_solver	(FORCES_NLP_solver_float)(1E-06)

/* desired maximum violation of complementarity */
#define SET_ACC_KKTCOMPL_FORCES_NLP_solver	(FORCES_NLP_solver_float)(1E-06)


/* RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define OPTIMAL_FORCES_NLP_solver      (1)

/* maximum number of iterations has been reached */
#define MAXITREACHED_FORCES_NLP_solver (0)

/* wrong number of inequalities error */
#define INVALID_NUM_INEQ_ERROR_FORCES_NLP_solver  (-4)

/* factorization error */
#define FACTORIZATION_ERROR_FORCES_NLP_solver   (-5)

/* NaN encountered in function evaluations */
#define BADFUNCEVAL_FORCES_NLP_solver  (-6)

/* no progress in method possible */
#define NOPROGRESS_FORCES_NLP_solver   (-7)

/* invalid values in parameters */
#define PARAM_VALUE_ERROR_FORCES_NLP_solver   (-11)

/* licensing error - solver not valid on this machine */
#define LICENSE_ERROR_FORCES_NLP_solver  (-100)

/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct
{
    /* vector of size 180 */
    FORCES_NLP_solver_float x0[180];

    /* vector of size 6 */
    FORCES_NLP_solver_float xinit[6];

    /* vector of size 240 */
    FORCES_NLP_solver_float all_parameters[240];


} FORCES_NLP_solver_params;


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct
{
    /* vector of size 9 */
    FORCES_NLP_solver_float x01[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x02[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x03[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x04[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x05[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x06[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x07[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x08[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x09[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x10[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x11[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x12[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x13[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x14[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x15[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x16[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x17[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x18[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x19[9];

    /* vector of size 9 */
    FORCES_NLP_solver_float x20[9];


} FORCES_NLP_solver_output;


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct
{
    /* iteration number */
    solver_int32_default it;

	/* number of iterations needed to optimality (branch-and-bound) */
	solver_int32_default it2opt;
	
    /* inf-norm of equality constraint residuals */
    FORCES_NLP_solver_float res_eq;
	
    /* inf-norm of inequality constraint residuals */
    FORCES_NLP_solver_float res_ineq;

	/* norm of stationarity condition */
    FORCES_NLP_solver_float rsnorm;

	/* max of all complementarity violations */
    FORCES_NLP_solver_float rcompnorm;

    /* primal objective */
    FORCES_NLP_solver_float pobj;	
	
    /* dual objective */
    FORCES_NLP_solver_float dobj;	

    /* duality gap := pobj - dobj */
    FORCES_NLP_solver_float dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    FORCES_NLP_solver_float rdgap;		

    /* duality measure */
    FORCES_NLP_solver_float mu;

	/* duality measure (after affine step) */
    FORCES_NLP_solver_float mu_aff;
	
    /* centering parameter */
    FORCES_NLP_solver_float sigma;
	
    /* number of backtracking line search steps (affine direction) */
    solver_int32_default lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    solver_int32_default lsit_cc;
    
    /* step size (affine direction) */
    FORCES_NLP_solver_float step_aff;
    
    /* step size (combined direction) */
    FORCES_NLP_solver_float step_cc;    

	/* solvertime */
	FORCES_NLP_solver_float solvetime;   

	/* time spent in function evaluations */
	FORCES_NLP_solver_float fevalstime;  

} FORCES_NLP_solver_info;







/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* Time of Solver Generation: (UTC) Tuesday, July 28, 2020 8:46:35 AM */
/* User License expires on: (UTC) Thursday, January 21, 2021 10:00:00 PM (approx.) (at the time of code generation) */
/* Solver Static License expires on: (UTC) Thursday, January 21, 2021 10:00:00 PM (approx.) */
/* Solver Generation Request Id: 9f99a5bf-57c6-43ae-a6dc-09ad4f36de30 */
/* examine exitflag before using the result! */
#ifdef __cplusplus
extern "C" {
#endif		

typedef void (*FORCES_NLP_solver_extfunc)(FORCES_NLP_solver_float* x, FORCES_NLP_solver_float* y, FORCES_NLP_solver_float* lambda, FORCES_NLP_solver_float* params, FORCES_NLP_solver_float* pobj, FORCES_NLP_solver_float* g, FORCES_NLP_solver_float* c, FORCES_NLP_solver_float* Jeq, FORCES_NLP_solver_float* h, FORCES_NLP_solver_float* Jineq, FORCES_NLP_solver_float* H, solver_int32_default stage, solver_int32_default iterations);

extern solver_int32_default FORCES_NLP_solver_solve(FORCES_NLP_solver_params *params, FORCES_NLP_solver_output *output, FORCES_NLP_solver_info *info, FILE *fs, FORCES_NLP_solver_extfunc evalextfunctions_FORCES_NLP_solver);	





#ifdef __cplusplus
}
#endif

#endif
