#include "global_declarations.h"
#include "basic_functions.h"

void MOVERS_H1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&dui, double ***&duj, double **&area, double ***&si, double ***&sj, double ***&diss, double ***&rhs){

    double *fc;
    double qsl, qsr, ghavg;
    
    fc = new double[nconv];
    
    for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            
            rhs[0][i][j]= -diss[0][i][j];
            rhs[1][i][j]= -diss[1][i][j];
            rhs[2][i][j]= -diss[2][i][j];
        }
    }
    
    for(int j=2;j<=jb;j++){
        //flux in i direction
        for(int i=3;i<=ib;i++){
            qsl = (cv[1][i-1][j]*si[0][i][j] + cv[2][i-1][j]*si[1][i][j])/cv[0][i-1][j];
            qsr = (cv[1][i][j]*si[0][i][j] + cv[2][i][j]*si[1][i][j])/cv[0][i][j];
            
            ghavg = 0.5*0.5*g*(cv[0][i-1][j]*cv[0][i-1][j]+cv[0][i][j]*cv[0][i][j]);
            
            fc[0] = 0.5*(qsl*cv[0][i-1][j]+qsr*cv[0][i][j]);
            fc[1] = 0.5*(qsl*cv[1][i-1][j]+qsr*cv[1][i][j]) + ghavg*si[0][i][j];
            fc[2] = 0.5*(qsl*cv[2][i-1][j]+qsr*cv[2][i][j]) + ghavg*si[1][i][j];
            
            rhs[0][i][j] = rhs[0][i][j] + fc[0];
            rhs[1][i][j] = rhs[1][i][j] + fc[1];
            rhs[2][i][j] = rhs[2][i][j] + fc[2];
            
            rhs[0][i-1][j] = rhs[0][i-1][j] - fc[0];
            rhs[1][i-1][j] = rhs[1][i-1][j] - fc[1];
            rhs[2][i-1][j] = rhs[2][i-1][j] - fc[2];
        }
        
        //flux in j-direction except at boundaries
        if(j>2){
            for(int i=2;i<=ib;i++){
                
                qsl = (cv[1][i][j-1]*sj[0][i][j] + cv[2][i][j-1]*sj[1][i][j])/cv[0][i][j-1];
                qsr = (cv[1][i][j]*sj[0][i][j] + cv[2][i][j]*sj[1][i][j])/cv[0][i][j];

                ghavg = 0.5*0.5*g*(cv[0][i][j-1]*cv[0][i][j-1]+cv[0][i][j]*cv[0][i][j]);

                fc[0] = 0.5*(qsl*cv[0][i][j-1]+qsr*cv[0][i][j]);
                fc[1] = 0.5*(qsl*cv[1][i][j-1]+qsr*cv[1][i][j]) + ghavg*sj[0][i][j];
                fc[2] = 0.5*(qsl*cv[2][i][j-1]+qsr*cv[2][i][j]) + ghavg*sj[1][i][j];

                rhs[0][i][j] = rhs[0][i][j] + fc[0];
                rhs[1][i][j] = rhs[1][i][j] + fc[1];
                rhs[2][i][j] = rhs[2][i][j] + fc[2];

                rhs[0][i][j-1] = rhs[0][i][j-1] - fc[0];
                rhs[1][i][j-1] = rhs[1][i][j-1] - fc[1];
                rhs[2][i][j-1] = rhs[2][i][j-1] - fc[2];
            }
        }
        
    }
    
    Flux_Boundary(ib, id1, id2, jb, jd1, jd2, 1, 2, ib, cv, dv, si, sj, rhs);
    Flux_Boundary(ib, id1, id2, jb, jd1, jd2, 2, 2, jb, cv, dv, si, sj, rhs);
    Flux_Boundary(ib, id1, id2, jb, jd1, jd2, 3, 2, ib, cv, dv, si, sj, rhs);
    Flux_Boundary(ib, id1, id2, jb, jd1, jd2, 4, 2, jb, cv, dv, si, sj, rhs);
    
    /*for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            std::cout<<"the computed rhs="<<rhs[1][i][j]<<std::endl;
        }
    }*/
    
    //at boundaries
    
    
    
    delete [] fc;
    
}
