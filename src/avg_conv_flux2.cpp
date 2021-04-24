#include "global_declarations.h"
#include "basic_functions.h"

template <typename T>
T Minmod_Flux(T del_p, T del_m);

void Avg_Conv_Flux2(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double **&area, double ***&si, double ***&sj, double ***&diss, double ***&rhs){

    double *fc, *deltl, *deltr;
    double ghavg;
    int im1, jm1, ip1, jp1;
    double hinv, hl, ul, vl, hr, ur, vr, qslr, qsll;
    
    
    
    fc = new double[nconv]; deltl = new double[nconv]; deltr = new double[nconv];
    
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
            
            im1=i-1;
            ip1=i+1;
            
            deltl[0] = 0.5*dui[0][im1][j]*Minmod_Flux(dui[0][i][j], dui[0][im1][j]);
            deltr[0] = 0.5*dui[0][i][j]*Minmod_Flux(dui[0][ip1][j], dui[0][i][j]);
            deltl[1] = 0.5*dui[1][im1][j]*Minmod_Flux(dui[1][i][j], dui[1][im1][j]);
            deltr[1] = 0.5*dui[1][i][j]*Minmod_Flux(dui[1][ip1][j], dui[1][i][j]);
            deltl[2] = 0.5*dui[2][im1][j]*Minmod_Flux(dui[2][i][j], dui[2][im1][j]);
            deltr[2] = 0.5*dui[2][i][j]*Minmod_Flux(dui[2][ip1][j], dui[2][i][j]);
            
            hinv = 1.0/cv[0][im1][j];
            hl   = cv[0][im1][j] + deltl[0];
            ul   = cv[1][im1][j]*hinv + deltl[1];
            vl   = cv[2][im1][j]*hinv + deltl[2];
            qsll  = (ul*si[0][i][j]+vl*si[1][i][j])*hl;
            
            hinv = 1.0/cv[0][i][j];
            hr   = cv[0][i][j] - deltr[0];
            ur   = cv[1][i][j]*hinv - deltr[1];
            vr   = cv[2][i][j]*hinv - deltr[2];
            qslr  = (ur*si[0][i][j]+vr*si[1][i][j])*hr;
            
            ghavg = 0.5*0.5*g*(hl*hl+hr*hr);
            
            fc[0] = 0.5*(qsll + qslr);
            fc[1] = 0.5*(qsll*ul + qslr*ur) + ghavg*si[0][i][j];
            fc[2] = 0.5*(qsll*vl + qslr*vr) + ghavg*si[1][i][j];
            
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
                
                jm1=j-1;
                jp1=j+1;
                
                deltl[0] = 0.5*duj[0][i][jm1]*Minmod_Flux(duj[0][i][j], duj[0][i][jm1]);
                deltr[0] = 0.5*duj[0][i][j]*Minmod_Flux(duj[0][i][jp1], duj[0][i][j]);
                deltl[1] = 0.5*duj[1][i][jm1]*Minmod_Flux(duj[1][i][j], duj[1][i][jm1]);
                deltr[1] = 0.5*duj[1][i][j]*Minmod_Flux(duj[1][i][jp1], duj[1][i][j]);
                deltl[2] = 0.5*duj[2][i][jm1]*Minmod_Flux(duj[2][i][j], duj[2][i][jm1]);
                deltr[2] = 0.5*duj[2][i][j]*Minmod_Flux(duj[2][i][jp1], duj[2][i][j]);

                hinv = 1.0/cv[0][i][jm1];
                hl   = cv[0][i][jm1] + deltl[0];
                ul   = cv[1][i][jm1]*hinv + deltl[1];
                vl   = cv[2][i][jm1]*hinv + deltl[2];
                qsll = (ul*sj[0][i][j]+vl*sj[1][i][j])*hl;

                hinv = 1.0/cv[0][i][j];
                hr   = cv[0][i][j] - deltr[0];
                ur   = cv[1][i][j]*hinv - deltr[1];
                vr   = cv[2][i][j]*hinv - deltr[2];
                qslr = (ur*sj[0][i][j]+vr*sj[1][i][j])*hr;
                
                ghavg = 0.5*0.5*g*(hl*hl+hr*hr);

                fc[0] = 0.5*(qsll + qslr);
                fc[1] = 0.5*(qsll*ul + qslr*ur) + ghavg*sj[0][i][j];
                fc[2] = 0.5*(qsll*vl + qslr*vr) + ghavg*sj[1][i][j];

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
    delete [] deltl;
    delete [] deltr;
    
}


template <typename T>
T Minmod_Flux(T del_p, T del_m){
    
    double r;
    
    r=0.0;
    
    r=del_p/(del_m+1e-5);
    
    return std::max(0.0,std::min(1.0,r));
    
}