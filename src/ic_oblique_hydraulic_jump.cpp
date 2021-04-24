#include "global_declarations.h"
#include "basic_functions.h"

void Ic_Oblique_Hydraulic_Jump(int id2, int jd2, double ***&cv){
    
    double const h=1.0, u=8.57, v=0.0;
                
    for(int j=0;j<=jd2;j++){
        for(int i=0;i<=id2;i++){
            
            cv[0][i][j] = h;
            cv[1][i][j] = h*u;
            cv[2][i][j] = h*v;
            
        }
    }
    
}
