#include "global_declarations.h"
#include "basic_functions.h"

int main(){

	std::cout << std::fixed;
        std::cout << std::setprecision(5);
    
	Declaration_1d_array();
	Declaration_2d_array();
	Declaration_3d_array();

	Read_grid( Nx, id1, id2, Ny, jd1, jd2, x, y);

        Grid_Computations(ib, jb, id1, id2, jd1, jd2, x, y, area, si, sj);
        
        //Initialize_Circular_Dam_Break(id1, id2, jd1, jd2, nconv, x, y, cv);
        //Ic_Oblique_Hydraulic_Jump(id2, jd2, cv);
        Ic_Dam_Break(id1, id2, jd1, jd2, x, y, cv);
     
        Dependent_Variables(id2, jd2, cv, dv);
        
        //Bc_Circular_Dam_Break(ib, id1, id2, jb, jd1, jd2, cv, dv, x, y, si, sj);
        //Bc_Oblique_Hydraulic_Jump(ib, id1, id2, jb, jd1, jd2, cv, dv, x, y, si, sj);
        Bc_Dam_Break(ib, id1, id2, jb, jd1, jd2, cv);
             
        
        int iter=0, maxiter=1000;
        double t=0, tmax=7.2;
        
        while(t<tmax){
            
          t=t+tstep[3][3];
          iter=iter+1;
          std::cout<<"iter number"<<iter<<"\t"<<"current time="<<tstep[2][2]<<"\t"<<"total time="<<t<<std::endl;
          
          //Solver_Explicit(ib, id1, id2, jb, jd1, jd2, nconv, ndvar, x, y, cv, dv, dui, duj, area, si, sj, tstep, sri,  srj, gradfi,  gradfj, cvold, diss, rhs, epsij);
          Solver(ib, id1, id2, jb, jd1, jd2, nconv, ndvar, x, y, cv, dv, dui, duj, area, si, sj, tstep, sri,  srj, gradfi,  gradfj, cvold, diss, rhs, epsij);
        }
        
        Write_Solution(id1, jd1, x, y, cv);
        
        return 0;
}
