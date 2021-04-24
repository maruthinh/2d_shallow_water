#include "global_declarations.h"

void Grid_Computations(int ib, int jb, int id1, int id2, int jd1, int jd2, double **&x, double **&y, double **&area, double ***&si, double ***&sj){
    
    double sumx, sumy, smax, sum, area_min, area_max;
    
    /*normal vectors in i-direction*/
    for(int j=0;j<=jd1;j++){
        for(int i=0;i<=id2;i++){
            
            si[0][i][j] = y[i][j] - y[i][j+1];
            si[1][i][j] = x[i][j+1] - x[i][j];
        }
    }
    
    /*normal vectors in j-direction*/
    for(int j=0;j<=jd2;j++){
        for(int i=0;i<=id1;i++){
            sj[0][i][j] = y[i+1][j] - y[i][j];
            sj[1][i][j] = x[i][j] - x[i+1][j];
        }
    }
    
    /*to compute cell areas*/
    for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            
            area[i][j] = 0.5*((x[i][j]-x[i+1][j+1])*(y[i+1][j]-y[i][j+1])+(x[i][j+1]-x[i+1][j])*(y[i][j]-y[i+1][j+1]));
            
        }
    }
    
    for(int i=2;i<=ib;i++){
        area[i][1] = area[i][2];
        area[i][0] = area[i][2];
        area[i][jd1] = area[i][jb];
        area[i][jd2] = area[i][jb];
    }
    
    for(int j=2;j<=jb;j++){
        area[1][j] = area[2][j];
        area[0][j] = area[2][j];
        area[id1][j] = area[ib][j];
        area[id2][j] = area[ib][j];
    }
    
   // BC_Periodic(ib, jb, id1, id2, jd1, jd2, 1, 1, 32, 2, 32, 1, 1, pbvar);
    
    /*for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            area[i][j]=pbvar[1][i][j];
        }
    }*/
    
        /*to get the sum of face vectors and min/max area of the cells*/
    smax = -1e32;
    area_min = 1e32;
    area_max = -1e32;
    
    for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            
            sumx     = si[0][i][j]-si[0][i+1][j]+sj[0][i][j]-sj[0][i][j+1];
            sumy     = si[1][i][j]-si[1][i+1][j]+sj[1][i][j]-sj[1][i][j+1];
            sum      = sqrt(sumx*sumx+sumy*sumy);
            smax     = std::max(sum, smax);
            area_min = std::min(area_min, area[i][j]);
            area_max = std::max(area_min, area[i][j]);
        }
    }
    
    /*for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            std::cout<<"the computed normals are="<<si[0][i][j]<<std::endl;
        }
    }*/
    
    std::cout<<"maximum of directional normals and min and max area"<<"\t"<<"smax="<<smax<<"\t"<<"area_min="<<area_min<<"\t"<<"area_max="<<area_max<<std::endl;
    
}

