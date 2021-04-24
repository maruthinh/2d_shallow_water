#include "global_declarations.h"
#include "basic_functions.h"

void Gradient_FaceI(int ib, int id1, int id2, int jb, int jd1, int jd2, int nconv, int ndvar, double ***&cv, double ***&dv, double **&u, double **&v, double **&area, double ***&si, double ***&sj, double ***&gradfi){
    
    //computes gradients of u, v, Temperature
    
    double sx, sy, uav, vav, tav, rvol;
    double *fgx, *fgy;
    
    fgx=new double[3]; fgy=new double[3];
    
    /****right face of the auxiliary control volume****/
    
    for(int j=2; j<=jd1; j++){
        for(int i=2;i<=ib;i++){
            
            sx     = 0.5*(si[0][i][j]+si[0][i+1][j]);
            sy     = 0.5*(si[1][i][j]+si[1][i+1][j]);
            fgx[0] = u[i][j]*sx;
            fgx[1] = v[i][j]*sx;
            fgx[2] = dv[1][i][j]*sx;
            fgy[0] = u[i][j]*sy;
            fgy[1] = v[i][j]*sy;
            fgy[2] = dv[1][i][j]*sy;
            
            gradfi[0][i][j] = gradfi[0][i][j] - fgx[0];
            gradfi[1][i][j] = gradfi[1][i][j] - fgy[0];
            gradfi[2][i][j] = gradfi[2][i][j] - fgx[1];
            gradfi[3][i][j] = gradfi[3][i][j] - fgy[1];
            gradfi[4][i][j] = gradfi[4][i][j] - fgx[2];
            gradfi[5][i][j] = gradfi[5][i][j] - fgy[2];
            
            gradfi[0][i+1][j] = gradfi[0][i][j] + fgx[0];
            gradfi[1][i+1][j] = gradfi[1][i][j] + fgy[0];
            gradfi[2][i+1][j] = gradfi[2][i][j] + fgx[1];
            gradfi[3][i+1][j] = gradfi[3][i][j] + fgy[1];
            gradfi[4][i+1][j] = gradfi[4][i][j] + fgx[2];
            gradfi[5][i+1][j] = gradfi[5][i][j] + fgy[2];
 
            /****bottom face of auxiliary cv****/
            
            sx  = 0.5*(sj[0][i][j]+sj[0][i-1][j]);
            sy  = 0.5*(sj[1][i][j]+sj[1][i-1][j]);
            uav = 0.25*(u[i][j]+u[i-1][j]+u[i][j-1]+u[i-1][j-1]);
            vav = 0.25*(v[i][j]+v[i-1][j]+v[i][j-1]+v[i-1][j-1]);
            tav = 0.25*(dv[1][i][j]++dv[1][i-1][j]+dv[1][i][j-1]+dv[1][i-1][j-1]);
            
            fgx[0] = uav*sx;
            fgx[1] = vav*sx;
            fgx[2] = tav*sx;
            
            fgy[0] = uav*sy;
            fgy[1] = vav*sy;
            fgy[2] = tav*sy;
            
            gradfi[0][i][j] = gradfi[0][i][j] + fgx[0];
            gradfi[1][i][j] = gradfi[1][i][j] + fgy[0];
            gradfi[2][i][j] = gradfi[2][i][j] + fgx[1];
            gradfi[3][i][j] = gradfi[3][i][j] + fgy[1];
            gradfi[4][i][j] = gradfi[4][i][j] + fgx[2];
            gradfi[5][i][j] = gradfi[4][i][j] + fgy[2]; 
           
        }
    }
}

