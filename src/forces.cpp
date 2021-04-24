#include "global_declarations.h"
#include "basic_functions.h"

void Forces(int ib, int id1, int id2, int jb, int jd1, int jd2, int ndvar, double **&x, double **&y, double ***&dv, double ***&si, double ***&sj){
    
    double cx, cy, cm, sx, sy, pwall, xa, ya, dcy, dcx;
    double xref, yref, cref;
    
    cx = 0; cy = 0; cm = 0; 
    xref = 0.25; yref = 0.0; cref = 1.0;
    
    if(bind==1){
        idum1 = 1;
        iins1 = 2;
    }
    else if(bind==2){
        idum1 = id1;
        iins1 = ib;
    }
    else if(bind==3){
        jdum1 = 1;
        jins1 = 2;
    }
    else if(bind==4){
        jdum1 = jd1;
        jins1 = jb;
    }    
    
    if(bind==3 or bind==4){
        
        //then sbind, ebind will be i
        
        for(int i=sbind; i<=ebind; i++){
            
            if(bind==3){
                
                sx = sj[0][i][jins1];
                sy = sj[1][i][jins1];
            }
            else{
                sx = -sj[0][i][jdum1];
                sy = -sj[1][i][jdum1];
            }
            
            pwall = 0.5*(dv[0][i][jdum1]+dv[0][i][jins1]);
            cp    = 2.0*(pwall-pinf)/(rhoinf*qinf*qinf);
            xa    = (0.5*(x[i][jins1]+x[i+1][jins1])-xref)/cref;
            ya    = (0.5*(y[i][jins1]+y[i+1][jins1])-yref)/cref;
            dcy   = sy*cp;
            dcx   = sx*cp;
            cy    = cy + dcy;
            cx    = cx + dcx;
            cm    = cm + dcx*ya - dcy*xa;
            
        }
    }
    if(bind==1 or bind==2){
        
        //then sbind, ebind will be j
        
        for(int j=sbind; j<=ebind; j++){
            
            if(bind==1){
                
                sx = sj[0][iins1][j];
                sy = sj[1][iins1][j];
            }
            else{
                sx = -sj[0][idum1][j];
                sy = -sj[1][idum1][j];
            }
            
            pwall = 0.5*(dv[0][idum1][j]+dv[0][iins1][j]);
            cp    = 2.0*(pwall-pinf)/(rhoinf*qinf*qinf);
            xa    = (0.5*(x[iins1][j]+x[iins1][j+1])-xref)/cref;
            ya    = (0.5*(y[iins1][j]+y[iins1][j+1])-yref)/cref;
            dcy   = sy*cp;
            dcx   = sx*cp;
            cy    = cy + dcy;
            cx    = cx + dcx;
            cm    = cm + dcx*ya - dcy*xa;
            
        }
    }
    
    cl = cy*cos(alpha) - cx*sin(alpha);
    cd = cy*sin(alpha) - cx*cos(alpha);
    
    std::cout<<"the values of cl and cd in forces"<<"\t"<<cl<<"\t"<<cd<<"\t"<<std::endl;
}

