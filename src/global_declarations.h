#ifndef _Header_H
#define _Header_H

#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <string>	
#include <iosfwd>

/*****constants used throught the program******/
const double g=9.81, pi=4.0*atan(1.0), rad=180.0/pi;
const int ndvar=1, nconv=3;
extern int mxeqn;

/****global integer variables related to grid and domain*****/
extern int Nx, Ny, ib, jb, id1, id2, jd1, jd2, imax, jmax;

/****integers related to boundary conditions****/
extern int bind, sbind, ebind, pbind, spbind, epbind, idum1, idum2, jdum1, jdum2,
           iins1, iins2, jins1, jins2; 

/****other variables used throughout the code****/
extern double machinf, rhoinf, uinf, vinf, pinf, alpha, refrho, refvel, qinf, cp, cl, cd,
        minject, tinject, maxwchg, maxichg, ttinl, ptinl, pout, betainl, betaout;
           
/****global strings declared here, which will be used throughout the code*****/
extern std::string flow_type, flo_intr_extr, space_accuracy, flux_type, time_accuracy, test_case, grid_file, vort_corr, time_step;

/****global pointers declared here, which will be used throughout the code*****/
extern double cfl;
extern double *Icond, *Prim_var_L, *Prim_var_M, *Prim_var_R, *limref;
extern double ***cv, ***cvold, ***dv, ***rhs, ***dui, ***duj, ***si, ***sj, ***pbvar, ***gradfi, ***gradfj, ***diss, ***epsij;
extern double **x, **y, **xc, **yc, **area, **p, **sri, **srj, **tstep, **u, **v;
extern double **L_state, **R_state;


/*********************forward declaration of functions called in the main program **********************/
void Declaration_1d_array();
void Declaration_2d_array();
void Declaration_3d_array();
void Read_grid( int Nx, int id1, int id2, int Ny, int jd1, int jd2, double **&x, double **&y);
void Grid_Computations(int ib, int jb, int id1, int id2, int jd1, int jd2, double **&x, double **&y, double **&area, double ***&si, double ***&sj);
void Initialize_Circular_Dam_Break(int id1, int id2, int jd1, int jd2, int nconv, double **&x, double **&y, double ***&cv);
void Ic_Dam_Break(int id1, int id2, int jd1, int jd2, double **&x, double **&y, double ***&cv);
void Ic_Oblique_Hydraulic_Jump(int id2, int jd2, double ***&cv);
void Dependent_Variables(int id2, int jd2, double ***&cv, double ***&dv);
void Dependent_Variables_One(int i, int j, double ***&cv, double ***&dv);
void Bc_Transmitive(int ib, int id1, int id2, int jb, int jd1, int jd2, int bind, int sbind, int ebind, double ***&cv);
void Bc_Prescribed_Inflow(int ib, int id1, int id2, int jb, int jd1, int jd2, int bind, int sbind, int ebind, double ***&cv);
void BC_wall(int ib, int id1, int id2, int jb, int jd1, int jd2, int bind, int sbind, int ebind, double ***&cv, double ***&dv);
void BC_Interior_wall(int int_dum, int bind, int sbind, int ebind, double ***&cv, double ***&dv);
void Prim_Variables_Differences(int ib, int id2, int jb, int jd2, double ***cv, double ***dv, double ***&dui, double ***&duj);
void Conv_Variables_Differences(int ib, int id2, int jb, int jd2, double ***cv, double ***dv, double ***&dui, double ***&duj);
void Time_Step(int ib, int id1, int id2, int jb, int jd1, int jd2, int nconv, int ndepv, double ***&cv, double ***&dv, double **&area, double ***&si, double ***&sj, double **&tstep, double **&sri,  double **&srj, double ***&epsij);
void Avg_Conv_Flux1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double **&area, double ***&si, double ***&sj, double ***&diss, double ***&rhs);
void Avg_Conv_Flux2(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double **&area, double ***&si, double ***&sj, double ***&diss, double ***&rhs);
void Diss_MOVERS(int ib, int id1, int jb, int jd1, double ***&cv, double ***&dv, double ***&si, double ***&sj, double ***&diss);
void Diss_MOVERS_H1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si, double ***&sj, double ***&diss);
void Diss_MOVERS_LE1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si, double ***&sj, double ***&diss);
void Diss_MOVERS_H2_Prim(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si, double ***&sj, double ***&diss);
void Diss_MOVERS_H2_Conv(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si, double ***&sj, double ***&diss);
void Diss_MOVERS_LE2_Conv(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si, double ***&sj, double ***&diss);
void Diss_MOVERS_LE2_Prim(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si, double ***&sj, double ***&diss);
void MOVERS_H1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&dui, double ***&duj, double **&area, double ***&si, double ***&sj, double ***&diss, double ***&rhs);
void Flux_Boundary(int ib, int id1, int id2, int jb, int jd1, int jd2, int bind, int sbind, int ebind, double ***&cv, double ***&dv, double ***&si, double ***&sj, double ***&rhs);
void Solver_Explicit(int ib, int id1, int id2, int jb, int jd1, int jd2, int nconv, int ndepv, double **&x, double **&y, double ***&cv, double ***&dv, double ***&dui, double ***&duj, double **&area, double ***&si, double ***&sj, double **&tstep, double **&sri,  double **&srj, double ***&gradfi,  double ***&tgradfj, double ***&cvold, double ***&diss, double ***&rhs, double ***&epsij);
void Solver(int ib, int id1, int id2, int jb, int jd1, int jd2, int nconv, int ndepv, double **&x, double **&y, double ***&cv, double ***&dv, double ***&dui, double ***&duj, double **&area, double ***&si, double ***&sj, double **&tstep, double **&sri,  double **&srj, double ***&gradfi,  double ***&tgradfj, double ***&cvold, double ***&diss, double ***&rhs, double ***&epsij);
void Bc_Circular_Dam_Break(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double **&x, double **&y, double ***&si, double ***&sj);
void Bc_Dam_Break(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv);
void Bc_Oblique_Hydraulic_Jump(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double **&x, double **&y, double ***&si, double ***&sj);
void Write_Solution(int id1, int jd1, double **&x, double **&y, double ***&cv);

#endif  //#ifndef _Header_H
