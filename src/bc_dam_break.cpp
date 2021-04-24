#include "global_declarations.h"

void Bc_Dam_Break(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv){
    
    /****fill dummy nodes at the corners****/
    for(int neqn=0;neqn<nconv;neqn++){
        
        cv[neqn][1][1]     = 0.5*(cv[neqn][1][2]+cv[neqn][2][1]);
        cv[neqn][id1][1]   = 0.5*(cv[neqn][ib][1]+cv[neqn][id1][2]);
        cv[neqn][1][jd1]   = 0.5*(cv[neqn][2][jd1]+cv[neqn][1][jb]);
        cv[neqn][id1][jd1] = 0.5*(cv[neqn][id1][jb]+cv[neqn][ib][jd1]);
        
    }
    
    /*for(int neqn=0;neqn<nconv;neqn++){
        std::cout<<cv[neqn][1][1]<<std::endl;
    }*/
    
    Dependent_Variables_One(1, 1, cv, dv);
    Dependent_Variables_One(id1, 1, cv, dv);
    Dependent_Variables_One(1, jd1, cv, dv);
    Dependent_Variables_One(id1, jd1, cv, dv);
    
    /*for(int neqn=0;neqn<ndvar;neqn++){
        std::cout<<dv[neqn][1][1]<<std::endl;
    }*/
    /*Bc_Transmitive(ib, id1, id2, jb, jd1, jd2, 4, 2, jb, cv);
    Bc_Transmitive(ib, id1, id2, jb, jd1, jd2, 2, 2, jb, cv);
    BC_wall(ib, id1, id2, jb, jd1, jd2, 1, 2, ib, cv, dv);
    BC_wall(ib, id1, id2, jb, jd1, jd2, 3, 2, ib, cv, dv);
    BC_Interior_wall(20, 2, 2, 20, cv, dv);
    BC_Interior_wall(20, 2, 36, jb, cv, dv);
    BC_Interior_wall(21, 4, 2, 20, cv, dv);
    BC_Interior_wall(21, 4, 36, jb, cv, dv);*/
    
    Bc_Transmitive(ib, id1, id2, jb, jd1, jd2, 4, 2, jb, cv);
    Bc_Transmitive(ib, id1, id2, jb, jd1, jd2, 2, 2, jb, cv);
    BC_wall(ib, id1, id2, jb, jd1, jd2, 1, 2, ib, cv, dv);
    BC_wall(ib, id1, id2, jb, jd1, jd2, 3, 2, ib, cv, dv);
    //BC_Interior_wall(21, 2, 2, 21, cv, dv);
    //BC_Interior_wall(21, 2, 35, jb, cv, dv);
    BC_Interior_wall(22, 4, 2, 21, cv, dv);
    BC_Interior_wall(22, 4, 35, jb, cv, dv);
}


