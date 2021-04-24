#include "global_declarations.h"
#include "basic_functions.h"

void Bc_Transmitive(int ib, int id1, int id2, int jb, int jd1, int jd2, int bind, int sbind, int ebind, double ***&cv){
    
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
        for(int i=sbind; i<=ebind; i++){
            
                /*cv[0][i][jdum1] = 2.0*cv[0][i][jins1] - cv[0][i][jins2];
                cv[1][i][jdum1] = 2.0*cv[1][i][jins1] - cv[1][i][jins2];
                cv[2][i][jdum1] = 2.0*cv[2][i][jins1] - cv[2][i][jins2];
                
                cv[0][i][jdum2] = 2.0*cv[0][i][jdum1] - cv[0][i][jins1];
                cv[1][i][jdum2] = 2.0*cv[1][i][jdum1] - cv[1][i][jins1];
                cv[2][i][jdum2] = 2.0*cv[2][i][jdum1] - cv[2][i][jins1];*/
            
                cv[0][i][jdum1] = cv[0][i][jins1];
                cv[1][i][jdum1] = cv[1][i][jins1];
                cv[2][i][jdum1] = cv[2][i][jins1];
                
                Dependent_Variables_One(i, jdum1, cv, dv);
                Dependent_Variables_One(i, jdum2, cv, dv);
                
        }
    }
    
    else if(bind==2 or bind==4){
        for(int j=sbind; j<=ebind; j++){
            
                cv[0][idum1][j] = cv[0][iins1][j];
                cv[1][idum1][j] = cv[1][iins1][j];
                cv[2][idum1][j] = cv[2][iins1][j];
                
                /*cv[0][idum1][j] = 2.0*cv[0][iins1][j] - cv[0][iins2][j];
                cv[1][idum1][j] = 2.0*cv[1][iins1][j] - cv[1][iins2][j];
                cv[2][idum1][j] = 2.0*cv[2][iins1][j] - cv[2][iins2][j];
                
                cv[0][idum2][j] = 2.0*cv[0][idum1][j] - cv[0][iins1][j];
                cv[1][idum2][j] = 2.0*cv[1][idum1][j] - cv[1][iins1][j];
                cv[2][idum2][j] = 2.0*cv[2][idum1][j] - cv[2][iins1][j];*/
                
                Dependent_Variables_One(idum1, j, cv, dv);
                Dependent_Variables_One(idum2, j, cv, dv);
        }
    }
}
