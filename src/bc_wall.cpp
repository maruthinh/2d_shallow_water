#include "global_declarations.h"
#include "basic_functions.h"

void BC_wall(int ib, int id1, int id2, int jb, int jd1, int jd2, int bind, int sbind, int ebind, double ***&cv, double ***&dv){
    
    double nx, ny, ds, vnor, vtan, u, v;
    
    maxwchg = 0.0;
    
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

        idum1 = 1;
        idum2 = 0;
        iins1 = 2;
        iins2 = 3;

    }
    
        if(bind==1 or bind==3){
        
           // if(sbind==2) sbind = 1;
           // if(ebind==ib) ebind = id1;
        
            for(int i=sbind; i<=ebind; i++){
                if(bind==1){
                    ds = (sj[0][i][jins1]*sj[0][i][jins1]+sj[1][i][jins1]*sj[1][i][jins1]);
                    nx = sj[0][i][jins1]/ds;
                    ny = sj[1][i][jins1]/ds;
                }
                else{
                    ds = (sj[0][i][jdum1]*sj[0][i][jdum1]+sj[1][i][jdum1]*sj[1][i][jdum1]);
                    nx = -sj[0][i][jdum1]/ds;
                    ny = -sj[1][i][jdum1]/ds;
                }
                
                vnor = (cv[1][i][jins1]*nx + cv[2][i][jins1]*ny)/cv[0][i][jins1];
                vtan = (cv[1][i][jins1]*ny - cv[2][i][jins1]*nx)/cv[0][i][jins1];
                vnor = -vnor;
                u    = vnor*nx + vtan*ny;
                v    = vnor*ny - vtan*nx;
                
                cv[0][i][jdum1] = cv[0][i][jins1];
                cv[1][i][jdum1] = cv[0][i][jins1]*u;
                cv[2][i][jdum1] = cv[0][i][jins1]*v;
                    
                cv[0][i][jdum2] = 2.0*cv[0][i][jdum1] - cv[0][i][jins1];
                cv[1][i][jdum2] = 2.0*cv[1][i][jdum1] - cv[1][i][jins1];
                cv[2][i][jdum2] = 2.0*cv[2][i][jdum1] - cv[2][i][jins1];
                
                Dependent_Variables_One(i, jdum1, cv, dv);
                Dependent_Variables_One(i, jdum2, cv, dv);
            
            }
    
        }
    
        if(bind==2 or bind==4){
        
           // if(sbind==2) sbind = 1;
           // if(ebind==jb) ebind = jd1;
        
            for(int j=sbind; j<=ebind; j++){
                    
                if(bind==4){
                    ds = (si[0][iins1][j]*si[0][iins1][j]+si[1][iins1][j]*si[1][iins1][j]);
                    nx = si[0][iins1][j]/ds;
                    ny = si[1][iins1][j]/ds;
                }
                else{
                    ds = (si[0][idum1][j]*si[0][idum1][j]+si[1][idum1][j]*si[1][idum1][j]);
                    nx = -si[0][idum1][j]/ds;
                    ny = -si[1][idum1][j]/ds;
                }
                
                vnor = (cv[1][iins1][j]*nx + cv[2][iins1][j]*ny)/cv[0][iins1][j];
                vtan = (cv[1][iins1][j]*ny - cv[2][iins1][j]*nx)/cv[0][iins1][j];
                vnor = -vnor;
                u    = vnor*nx + vtan*ny;
                v    = vnor*ny - vtan*nx;
                
                cv[0][idum1][j] = cv[0][iins1][j];
                cv[1][idum1][j] = cv[0][iins1][j]*u;
                cv[2][idum1][j] = cv[0][iins1][j]*v;
           
                cv[0][idum2][j] = 2.0*cv[0][idum1][j] - cv[0][iins1][j];
                cv[1][idum2][j] = 2.0*cv[1][idum1][j] - cv[1][iins1][j];
                cv[2][idum2][j] = 2.0*cv[2][idum1][j] - cv[2][iins1][j];
           
                Dependent_Variables_One(idum1, j, cv, dv);
                Dependent_Variables_One(idum2, j, cv, dv);
                
            }
        }
    
    /*for(int j=sbind;j<=ebind;j++){
        std::cout<<"the EULER_wall boundary values are="<<cv[3][idum1][j]<<std::endl;
    }*/
}

void BC_Interior_wall(int int_dum, int bind, int sbind, int ebind, double ***&cv, double ***&dv){
    
    double nx, ny, ds, vnor, vtan, u, v;
    int dum1, dum2, ins1;
    //int ins2;
    
    maxwchg = 0.0;
    
    if(bind==1 or bind==4 ){

        dum1 = int_dum;
        dum2 = int_dum-1;
        ins1 = int_dum+1;
        //ins2 = int_dum+2;
    }
    else if(bind==2 or bind==3){

        dum1 = int_dum;
        dum2 = int_dum+1;
        ins1 = int_dum-1;
        //ins2 = int_dum-2;

    }
    else{
        std::cout<<"boundary index error at interior of the grid"<<std::endl;
        exit(0);
    }
    
        if(bind==1 or bind==3){
        
           // if(sbind==2) sbind = 1;
           // if(ebind==ib) ebind = id1;
        
            for(int i=sbind; i<=ebind; i++){
                if(bind==1){
                    ds = (sj[0][i][ins1]*sj[0][i][ins1]+sj[1][i][ins1]*sj[1][i][ins1]);
                    nx = sj[0][i][ins1]/ds;
                    ny = sj[1][i][ins1]/ds;
                }
                else{
                    nx = -sj[0][i][jdum1]/ds;
                    ny = -sj[1][i][jdum1]/ds;
                }
                
                vnor = (cv[1][i][ins1]*nx + cv[2][i][ins1]*ny)/cv[0][i][ins1];
                vtan = (cv[1][i][ins1]*ny - cv[2][i][ins1]*nx)/cv[0][i][ins1];
                vnor = -vnor;
                u    = vnor*nx + vtan*ny;
                v    = vnor*ny - vtan*nx;
                
                cv[0][i][dum1] = cv[0][i][ins1];
                cv[1][i][dum1] = cv[0][i][ins1]*u;
                cv[2][i][dum1] = cv[0][i][ins1]*v;
                    
                //cv[0][i][dum2] = 2.0*cv[0][i][dum1] - cv[0][i][ins1];
                //cv[1][i][dum2] = 2.0*cv[1][i][dum1] - cv[1][i][ins1];
                //cv[2][i][dum2] = 2.0*cv[2][i][dum1] - cv[2][i][ins1];
                
                Dependent_Variables_One(i, dum1, cv, dv);
                //Dependent_Variables_One(i, dum2, cv, dv);
            
            }
    
        }
    
        if(bind==2 or bind==4){
        
           // if(sbind==2) sbind = 1;
           // if(ebind==jb) ebind = jd1;
        
            for(int j=sbind; j<=ebind; j++){
                    
                if(bind==4){
                    ds = (si[0][iins1][j]*si[0][ins1][j]+si[1][ins1][j]*si[1][ins1][j]);
                    nx = si[0][ins1][j]/ds;
                    ny = si[1][ins1][j]/ds;
                }
                else{
                    
                    nx = -si[0][ins1][j]/ds;
                    ny = -si[1][ins1][j]/ds;
                }
                
                vnor = (cv[1][ins1][j]*nx + cv[2][ins1][j]*ny)/cv[0][ins1][j];
                vtan = (cv[1][ins1][j]*ny - cv[2][ins1][j]*nx)/cv[0][ins1][j];
                vnor = -vnor;
                u    = vnor*nx + vtan*ny;
                v    = vnor*ny - vtan*nx;
                
                cv[0][dum1][j] = cv[0][ins1][j];
                cv[1][dum1][j] = cv[0][ins1][j]*u;
                cv[2][dum1][j] = cv[0][ins1][j]*v;
           
                //cv[0][dum2][j] = 2.0*cv[0][dum1][j] - cv[0][ins1][j];
                //cv[1][dum2][j] = 2.0*cv[1][dum1][j] - cv[1][ins1][j];
                //cv[2][dum2][j] = 2.0*cv[2][dum1][j] - cv[2][ins1][j];
           
                Dependent_Variables_One(dum1, j, cv, dv);
                //Dependent_Variables_One(dum2, j, cv, dv);
                
            }
        }
    
    /*for(int j=sbind;j<=ebind;j++){
        std::cout<<"the EULER_wall boundary values are="<<cv[3][idum1][j]<<std::endl;
    }*/
}

