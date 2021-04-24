#include "global_declarations.h"
#include "basic_functions.h"


/****global integer variables related to grid and domain*****/
//int Nx, Ny, ib, jb, id1, id2, jd1, jd2, imax, jmax;
int Nx=40, Ny=40, ib=Nx+1, jb=Ny+1, id1=Nx+2, id2=Nx+3, jd1=Ny+2, jd2=Ny+3, imax=Nx+4, jmax=Ny+4;
//int Nx=40, Ny=40, ib=Nx+1, jb=Ny+1, id1=Nx+2, id2=Nx+3, jd1=Ny+2, jd2=Ny+3, imax=Nx+4, jmax=Ny+4;
//int Nx=50, Ny=50, ib=51, jb=51, id1=52, id2=53, jd1=52, jd2=53, imax=1000, jmax=100;
//int Nx=40, Ny=30, ib=41, jb=31, id1=42, id2=43, jd1=32, jd2=53, imax=1000, jmax=100;
int maxeqn;

/****global strings declared here, which will be used throughout the code*****/
std::string flow_type, flo_intr_extr, space_accuracy, flux_type, time_accuracy, test_case, grid_file, vort_corr, time_step;
double cfl=0.8;

/****variables related to boundary conditons****/
int bind, sbind, ebind, pbind, spbind, epbind, idum1, idum2, jdum1, jdum2,
    iins1, iins2, jins1, jins2;

/****other variables used throughout the code****/
double machinf, rhoinf, uinf, vinf, pinf, alpha, refrho, refvel, qinf, cp, cl, cd,
        minject, tinject, maxwchg, maxichg, ttinl, ptinl, pout, betainl, betaout;

/****declaration of the global pointer variables used in the code****/
double *Icond, *Prim_var_L, *Prim_var_R, *limref;
double ***cv, ***cvold, ***dv, ***rhs, ***dui, ***duj, ***si, ***sj, ***pbvar, ***diss, ***epsij, ***gradfi, ***gradfj;
double **x, **y, **xc, **yc, **p, **area, **sri, **srj, **tstep, **u, **v;


/**************  to allocate 1d array***********/
void Declaration_1d_array(){
	Icond=new double[nconv]; Prim_var_L=new double[nconv]; Prim_var_R=new double[nconv], limref=new double[nconv];
}

/**************  to allocate 2d array***********/
void Declaration_2d_array(){
	Allocate_2D(x, imax, jmax); Allocate_2D(y, imax, jmax); Allocate_2D(xc, imax, jmax);
	Allocate_2D(yc, imax, jmax); Allocate_2D(area, imax, jmax); Allocate_2D(p, imax, jmax);  
	Allocate_2D(sri, imax, jmax); Allocate_2D(srj, imax, jmax); Allocate_2D(tstep, imax, jmax);
        Allocate_2D(u, imax, jmax); Allocate_2D(v, imax, jmax);
}

/**************  to allocate 3d array***********/
void Declaration_3d_array(){
	Allocate_3D(cv, nconv, imax, jmax); Allocate_3D(cvold, nconv, imax, jmax); Allocate_3D(dv, ndvar, imax, jmax);
	Allocate_3D(rhs, nconv, imax, jmax); Allocate_3D(dui, nconv, imax, jmax); Allocate_3D(duj, nconv, imax, jmax);
        Allocate_3D(si, 2, imax, jmax); Allocate_3D(sj, 2, imax, jmax); Allocate_3D(pbvar, 10, imax, jmax);
        Allocate_3D(gradfi, 6, imax, jmax); Allocate_3D(gradfj, 6, imax, jmax); Allocate_3D(diss, nconv, imax, jmax);
	Allocate_3D(epsij, 2, imax, jmax);
}
