#include "global_declarations.h"
#include "basic_functions.h"

void Avg_Conv_Flux1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double **&area, double ***&si, double ***&sj, double ***&diss, double ***&rhs){

    double *fc;
    double ha, hua, hva, vcont, ghavg;
    
    fc = new double[nconv];
    
    for(int j=2;j<=jb;j++){
        for(int i=2;i<=id1;i++){
            
            rhs[0][i][j]= -diss[0][i][j];
            rhs[1][i][j]= -diss[1][i][j];
            rhs[2][i][j]= -diss[2][i][j];
        }
    }
    
    for(int j=2;j<=jb;j++){
        //flux in i direction
        for(int i=2;i<=id1;i++){
            ha  = 0.5*(cv[0][i-1][j] + cv[0][i][j]);
            hua = 0.5*(cv[1][i-1][j] + cv[1][i][j]);
            hva = 0.5*(cv[2][i-1][j] + cv[2][i][j]);
            
            vcont = (hua*si[0][i][j] + hva*si[1][i][j])/ha;
            ghavg = 0.5*0.5*g*(cv[0][i-1][j]*cv[0][i-1][j]+cv[0][i][j]*cv[0][i][j]);
            
            fc[0] = vcont*ha;
            fc[1] = vcont*hua + ghavg*si[0][i][j];
            fc[2] = vcont*hva + ghavg*si[1][i][j];
            
            rhs[0][i][j] = rhs[0][i][j] + fc[0];
            rhs[1][i][j] = rhs[1][i][j] + fc[1];
            rhs[2][i][j] = rhs[2][i][j] + fc[2];
            
            rhs[0][i-1][j] = rhs[0][i-1][j] - fc[0];
            rhs[1][i-1][j] = rhs[1][i-1][j] - fc[1];
            rhs[2][i-1][j] = rhs[2][i-1][j] - fc[2];
        }
    }
        for(int j=2;j<=jd1;j++){
        //flux in j-direction except at boundaries
  
            for(int i=2;i<=ib;i++){
                
                ha  = 0.5*(cv[0][i][j-1] + cv[0][i][j]);
                hua = 0.5*(cv[1][i][j-1] + cv[1][i][j]);
                hva = 0.5*(cv[2][i][j-1] + cv[2][i][j]);

                vcont = (hua*sj[0][i][j] + hva*sj[1][i][j])/ha;
                ghavg = 0.5*0.5*g*(cv[0][i][j-1]*cv[0][i][j-1]+cv[0][i][j]*cv[0][i][j]);

                fc[0] = vcont*ha;
                fc[1] = vcont*hua + ghavg*sj[0][i][j];
                fc[2] = vcont*hva + ghavg*sj[1][i][j];

                rhs[0][i][j] = rhs[0][i][j] + fc[0];
                rhs[1][i][j] = rhs[1][i][j] + fc[1];
                rhs[2][i][j] = rhs[2][i][j] + fc[2];

                rhs[0][i][j-1] = rhs[0][i][j-1] - fc[0];
                rhs[1][i][j-1] = rhs[1][i][j-1] - fc[1];
                rhs[2][i][j-1] = rhs[2][i][j-1] - fc[2];
            }
        }
        
    
    
    //Flux_Boundary(ib, id1, id2, jb, jd1, jd2, 1, 2, ib, cv, dv, si, sj, rhs);
    //Flux_Boundary(ib, id1, id2, jb, jd1, jd2, 2, 2, jb, cv, dv, si, sj, rhs);
    //Flux_Boundary(ib, id1, id2, jb, jd1, jd2, 3, 2, ib, cv, dv, si, sj, rhs);
    //Flux_Boundary(ib, id1, id2, jb, jd1, jd2, 4, 2, jb, cv, dv, si, sj, rhs);
    
    /*for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            std::cout<<"the computed rhs="<<rhs[1][i][j]<<std::endl;
        }
    }*/
    
    //at boundaries
    
    
    
    delete [] fc;
    
}

