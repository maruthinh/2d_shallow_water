#include "global_declarations.h"

void Dependent_Variables(int id2, int jd2, double ***&cv, double ***&dv){
    
    double h;
    
    for (int j=0; j<=jd2; j++){
        for(int i=0; i<=id2; i++){
            
            h=cv[0][i][j];
            dv[0][i][j]=sqrt(g*h);            //sound speed
            
        }
    }
    
    /*for (int j=0; j<=jd2; j++){
        for(int i=0; i<=id2; i++){
            
            std::cout<<"the initialized dependent variables are="<<"\t"<<dv[0][i][j]<<std::endl;
            
        }
    }*/
}


void Dependent_Variables_One(int i, int j, double ***&cv, double ***&dv){
    
    double h;
    
    h=cv[0][i][j];
    dv[0][i][j]=sqrt(g*h);            //sound speed
    
}