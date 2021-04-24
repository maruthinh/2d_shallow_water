#include "global_declarations.h"
#include "basic_functions.h"

void Time_Step(int ib, int id1, int id2, int jb, int jd1, int jd2, int nconv, int ndepv, double ***&cv, double ***&dv, double **&area, double ***&si, double ***&sj, double **&tstep, double **&sri,  double **&srj, double ***&epsij){
    
    double rh, u, v, sx, sy, vc, cs, tsmin;
    
    time_step=="global";
    
    for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            
            rh      = 1/cv[0][i][j];
            u         = cv[1][i][j]*rh;
            v         = cv[2][i][j]*rh;
            
            sx        = 0.5*(si[0][i][j]+si[0][i+1][j]);
            sy        = 0.5*(si[1][i][j]+si[1][i+1][j]);
            vc        = u*sx + v*sy;
            cs        = sqrt(g*cv[0][i][j])*sqrt(sx*sx+sy*sy);
            sri[i][j] = fabs(vc) + cs;
            
            sx        = 0.5*(sj[0][i][j]+sj[0][i][j+1]);
            sy        = 0.5*(sj[1][i][j]+sj[1][i][j+1]);
            vc        = u*sx + v*sy;
            cs        = sqrt(g*cv[0][i][j])*sqrt(sx*sx+sy*sy);
            srj[i][j] = fabs(vc) + cs;
            
            tstep[i][j] = cfl*area[i][j]/(sri[i][j]+srj[i][j]);
        }
    }
    
    /****values of sri and srj at dummy nodes****/
    
    for(int i=2;i<=ib;i++){
        
        sri[i][1]     = sri[i][2];
        sri[i][0]     = sri[i][2];
        sri[i][jdum1] = sri[i][jb];
        sri[i][jdum2] = sri[i][jb];
        srj[i][1]     = srj[i][2];
        srj[i][0]     = srj[i][2];
        srj[i][jdum1] = srj[i][jb];
        srj[i][jdum2] = srj[i][jb];
        
        }
    for(int j=2;j<=jb;j++){
        
        sri[1][j]     = sri[2][j];
        sri[0][j]     = sri[2][j];
        sri[idum1][j] = sri[ib][j];
        sri[idum2][j] = sri[ib][j];
        srj[1][j]     = srj[2][j];
        srj[0][j]     = srj[2][j];
        srj[idum1][j] = srj[ib][j];
        srj[idum2][j] = srj[ib][j];
    }
    
    /*BC_Periodic(ib, jb, id1, id2, jd1, jd2, bind, sbind, ebind, pbind, spbind, epbind, 1, pbvar);
    
    for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            sri[i][j]=pbvar[1][i][j];
        }
    }*/
    
    /*BC_Periodic(ib, jb, id1, id2, jd1, jd2, bind, sbind, ebind, pbind, spbind, epbind, 1, pbvar);
    
    for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            srj[i][j]=pbvar[1][i][j];
        }
    }*/
    
    /****for global time-step*/
    

        tsmin = 1e32;
        for(int j=2;j<=jb;j++){
            for(int i=2;i<=ib;i++){
                
                tsmin = std::min(tsmin, tstep[i][j]); 
            }
        }
        
        for(int j=2;j<=jb;j++){
            for(int i=2;i<=ib;i++){
                
                tstep[i][j]=tsmin;
                
            }
        }
        
  
    
        /*tsmin = tstep[2][2];
        for(int j=2;j<=jb;j++){
            for(int i=2;i<=ib;i++){
                
                if(tstep[i][j]<tsmin) tsmin=tstep[i][j]; 
            }
        }
        
        for(int j=2;j<=jb;j++){
            for(int i=2;i<=ib;i++){
                
                tstep[i][j]=tsmin;
                
            }
        }*/

    
    /*for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            std::cout<<"the computed time step values are="<<tstep[i][j]<<std::endl;
        }
    }*/
    
}
