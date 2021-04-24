#include "global_declarations.h"

void Initialize_Circular_Dam_Break(int id1, int id2, int jd1, int jd2, int nconv, double **&x, double **&y, double ***&cv){
    
    
    for (int j = 2; j <= jd1; j++) {
        for (int i = 2; i <= jd2; i++) {
            
            if(((x[i][j]-25.0)*(x[i][j]-25.0)+(y[i][j]-25.0)*(y[i][j]-25.0))<=(11.0*11.0)){
                
                cv[0][i][j] = 10.0;
                cv[1][i][j] = 0.0;
                cv[2][i][j] = 0.0;
                
            }
            else{
                cv[0][i][j] = 1.0;
                cv[1][i][j] = 0.0;
                cv[2][i][j] = 0.0;
            }
        }
    }
    
    for (int i = 2; i <= id1; i++) {
            
        cv[0][i][0]=cv[0][i][2];
        cv[1][i][0]=cv[1][i][2];
        cv[2][i][0]=cv[2][i][2];
        
        cv[0][i][1]=cv[0][i][2];
        cv[1][i][1]=cv[1][i][2];
        cv[2][i][1]=cv[2][i][2];
        
        cv[0][i][jd2]=cv[0][i][jd1];
        cv[1][i][jd2]=cv[1][i][jd1];
        cv[2][i][jd2]=cv[2][i][jd1];
        
    }
    
    for(int j = 0; j <= jd2; j++){
        
        cv[0][0][j]=cv[0][2][j];
        cv[1][0][j]=cv[1][2][j];
        cv[2][0][j]=cv[2][2][j];
        
        cv[0][1][j]=cv[0][2][j];
        cv[1][1][j]=cv[1][2][j];
        cv[2][1][j]=cv[2][2][j];
        
        cv[0][id2][j]=cv[0][id1][j];
        cv[1][id2][j]=cv[1][id1][j];
        cv[2][id2][j]=cv[2][id1][j];
    }

    /*for (int j=0; j<=jd2; j++){
        for(int i=0; i<=id2; i++){
            
            std::cout<<"the initialized conserved variables are="<<"\t"<<cv[0][i][j]<<std::endl;
            
        }
    }*/
    
}



