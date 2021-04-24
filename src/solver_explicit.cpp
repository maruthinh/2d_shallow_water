#include "global_declarations.h"
#include "basic_functions.h"

void Solver_Explicit(int ib, int id1, int id2, int jb, int jd1, int jd2, int nconv, int ndepv, double **&x, double **&y, double ***&cv, double ***&dv, double ***&dui, double ***&duj, double **&area, double ***&si, double ***&sj, double **&tstep, double **&sri,  double **&srj, double ***&gradfi,  double ***&tgradfj, double ***&cvold, double ***&diss, double ***&rhs, double ***&epsij){
      
    double adtv; 
    
    for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            
            cvold[0][i][j] = cv[0][i][j];
            cvold[1][i][j] = cv[1][i][j];
            cvold[2][i][j] = cv[2][i][j];
            diss[0][i][j]  = 0.0;
            diss[1][i][j]  = 0.0;
            diss[2][i][j]  = 0.0;
            
        }
    }
    
    Time_Step(ib, id1, id2, jb, jd1, jd2, nconv, ndepv, cv, dv, area, si, sj, tstep, sri,  srj, epsij);
    Diss_MOVERS_H1(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
    //Prim_Variables_Differences(ib, id2, jb, jd2, cv, dv, dui, duj);
    //Conv_Variables_Differences(ib, id2, jb, jd2, cv, dv, dui, duj);
    //Diss_MOVERS_H2_Prim(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
    //Diss_MOVERS_H2_Conv(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
    //Diss_MOVERS_LE1(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
    //Diss_MOVERS_LE2_Conv(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
    //Diss_MOVERS_LE2_Prim(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
   // MOVERS_H1(ib, id1, id2, jb, jd1, jd2, cv, dv, dui, duj, area, si, sj, diss, rhs);
    Avg_Conv_Flux1(ib, id1, id2, jb, jd1, jd2, cv, dv, area, si, sj, diss, rhs);
    
    for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            adtv         = tstep[i][j]/area[i][j];
            rhs[0][i][j] = adtv*rhs[0][i][j];
            rhs[1][i][j] = adtv*rhs[1][i][j];
            rhs[2][i][j] = adtv*rhs[2][i][j];
            
        }
    }
    
    //update conserved variables
    for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            
            cv[0][i][j] = cvold[0][i][j] - rhs[0][i][j];
            cv[1][i][j] = cvold[1][i][j] - rhs[1][i][j];
            cv[2][i][j] = cvold[2][i][j] - rhs[2][i][j];
            
            if(cv[0][i][j]<0.0){
                std::cout<<"h became zero at i="<<i<<"\t"<<"and at j="<<j<<std::endl;
                exit(0);
            }
        }
    }
    
    Dependent_Variables(id2, jd2, cv, dv);
    //Bc_Circular_Dam_Break(ib, id1, id2, jb, jd1, jd2, cv, dv, x, y, si, sj);
    //Bc_Oblique_Hydraulic_Jump(ib, id1, id2, jb, jd1, jd2, cv, dv, x, y, si, sj);
    Bc_Dam_Break(ib, id1, id2, jb, jd1, jd2, cv);
}

