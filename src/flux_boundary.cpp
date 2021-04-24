#include "global_declarations.h"
#include "basic_functions.h"

void Flux_Boundary(int ib, int id1, int id2, int jb, int jd1, int jd2, int bind, int sbind, int ebind, double ***&cv, double ***&dv, double ***&si, double ***&sj, double ***&rhs){
    
    //computes convective fluxes at for boundaries  
    double sx, sy, hinv, hl, ul, vl, hr, ur, vr, qsrl, qsrr, ghavg;
    double *fc;
    
    fc = new double[nconv];
    
    if(bind==1){

        jdum1 = 1;
        jdum2 = 0;
        jins1 = 2;
        jins2 = 3;

    }
    else if(bind==2){

        idum1 = id1;
        idum2 = id2;
        iins1 = ib;
        iins2 = ib-1;

    }
    else if(bind==3){

        jdum1 = jd1;
        jdum2 = jd2;
        jins1 = jb;
        jins2 = jb-1;

    }
    else if(bind==4){

        jdum1 = 1;
        jdum2 = 0;
        jins1 = 2;
        jins2 = 3;

    }
    
    if(bind==1 or bind==3){
        for(int i=sbind; i<=ebind; i++){
            if(bind==1){
                sx = sj[0][i][jins1];
                sy = sj[1][i][jins1];
            }
            else{
                sx = -sj[0][i][jdum1];
                sy = -sj[1][i][jdum1];
            }
            
            hinv = 1.0/cv[0][i][jdum1];
            hl   = cv[0][i][jdum1];
            ul   = cv[1][i][jdum1]*hinv;
            vl   = cv[2][i][jdum1]*hinv;
            
            hinv = 1.0/cv[0][i][jins1];
            hr   = cv[0][i][jins1];
            ur   = cv[1][i][jins1]*hinv;
            vr   = cv[2][i][jins1]*hinv;
            
            qsrl  = (ul*sx + vl*sy)*hl;
            qsrr  = (ur*sx + vr*sy)*hr;
            
            ghavg = 0.5*0.5*g*(hl*hl+hr*hr);
            
            fc[0] = 0.5*(qsrl + qsrr);
            fc[1] = 0.5*(qsrl*ul + qsrr*ur)+ghavg*sx;
            fc[2] = 0.5*(qsrl*vl + qsrr*vr)+ghavg*sy;
            
            rhs[0][i][jins1] = rhs[0][i][jins1] + fc[0];
            rhs[1][i][jins1] = rhs[1][i][jins1] + fc[1];
            rhs[2][i][jins1] = rhs[2][i][jins1] + fc[2];
            
        }
        
    }
    else if(bind==2 or bind==4){
        for(int j=sbind; j<=ebind; j++){
            if(bind==4){
                sx = si[0][iins1][j];
                sy = si[1][iins1][j];
            }
            else{
                sx = -si[0][idum1][j];
                sy = -si[1][idum1][j];
            }
            
            hinv = 1.0/cv[0][idum1][j];
            hl   = cv[0][idum1][j];
            ul   = cv[1][idum1][j]*hinv;
            vl   = cv[2][idum1][j]*hinv;
            
            hinv = 1.0/cv[0][iins1][j];
            hr   = cv[0][iins1][j];
            ur   = cv[1][iins1][j]*hinv;
            vr   = cv[2][iins1][j]*hinv;
            
            qsrl  = (ul*sx + vl*sy)*hl;
            qsrr  = (ur*sx + vr*sy)*hr;
            
            ghavg = 0.5*0.5*g*(hl*hl+hr*hr);
            
            fc[0] = 0.5*(qsrl + qsrr);
            fc[1] = 0.5*(qsrl*ul + qsrr*ur)+ghavg*sx;
            fc[2] = 0.5*(qsrl*vl + qsrr*vr)+ghavg*sy;
            
            rhs[0][iins1][j] = rhs[0][iins1][j] + fc[0];
            rhs[1][iins1][j] = rhs[1][iins1][j] + fc[1];
            rhs[2][iins1][j] = rhs[2][iins1][j] + fc[2];
            
        }
        
    }
    
   /* for(int j=sbind;j<=ebind;j++){
        std::cout<<"Flux boundary values are="<<rhs[0][iins1][j]<<std::endl;
    }*/
    
    delete[] fc;
}
