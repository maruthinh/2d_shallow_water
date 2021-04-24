#include "global_declarations.h"
#include "basic_functions.h"

void Write_Solution(int id1, int jd1, double **&x, double **&y, double ***&cv){

    double h,u,v;
    
    std::ofstream other("./output/Result.dat");
    other.flags(std::ios::dec | std::ios::scientific);
    other.precision(5);
       
    if(!other){
    	std::cerr<<"File couldn't be opened to write the solution"<<std::endl;
	exit(1);
    }

    other << "TITLE = Flow" << std::endl << "VARIABLES = xc, yc, h, u, v" << std::endl;
    other << "Zone T = Omega I = " << Nx+1 << " J = " << Ny+1 << std::endl ;


    for (int j=2;j<=jd1;j++){
        for (int i=2;i<=id1;i++){
                            
            h = 0.25*(cv[0][i][j]+cv[0][i-1][j]+cv[0][i-1][j-1]+cv[0][i][j-1]);
            u = 0.25*(cv[1][i][j]/cv[0][i][j]+cv[1][i-1][j]/cv[0][i-1][j]+cv[1][i-1][j-1]/cv[0][i-1][j-1]+cv[1][i][j-1]/cv[0][i][j-1]);
            v = 0.25*(cv[2][i][j]/cv[0][i][j]+cv[2][i-1][j]/cv[0][i-1][j]+cv[2][i-1][j-1]/cv[0][i-1][j-1]+cv[2][i][j-1]/cv[0][i][j-1]);
            
            other<<x[i][j]<<"\t"<<y[i][j]<<"\t"<<h<<"\t"<<u<<"\t"<<v<<std::endl;
                
        }
    }

    other.close();
}
