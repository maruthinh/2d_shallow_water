#include "global_declarations.h"

void Read_grid(int Nx, int id1, int id2, int Ny, int jd1, int jd2, double **&x, double **&y) {

    int k;
    std::ifstream infile;
    //infile.open("./CircularDamBreakGrid/CircularDamBreakGrid5050.dat");
    //infile.open("./ObliqueHydraulicJumpGrid/ObliqueHydraulicJump4030.dat");
    infile.open("./DamBreakGrid/DamBreakGrid4040.dat");



    if (infile.fail()) {
        std::cerr << "File couldn't be opened to read the grid points" << std::endl;
        std::exit(1);
    }

    std::cout << "---------------------------------------------------------" << std::endl
            << "Reading grid points from the file" << std::endl
            << "---------------------------------------------------------" << std::endl;


    while (!infile.eof()) {
        for (int j = 2; j <= Ny + 2; j++) {
            for (int i = 2; i <= Nx + 2; i++) {

                infile >> k >> x[i][j] >> y[i][j];
                //std::cout << k << "\t" << x[i][j] << "\t" << y[i][j] << std::endl;
            }
        }
    }

    infile.close();

    //grid points for dummy cells

    
    for (int i = 2; i <= id1; i++) {
            
        x[i][0]=x[i][2];
        y[i][0]=y[i][2];
        x[i][1]=x[i][2];
        y[i][1]=y[i][2];
        x[i][jd2]=x[i][jd1];
        y[i][jd2]=y[i][jd1];
        
    }
    
    for(int j = 0; j <= jd2; j++){
        x[0][j] = x[2][j];
        y[0][j] = y[2][j];
        x[1][j] = x[2][j];
        y[1][j] = y[2][j];
        x[id2][j] = x[id1][j];
        y[id2][j] = y[id1][j];
    }

    /*for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {

            std::cout << k << "\t" << x[i][j] << "\t" << y[i][j] << std::endl;
        }
    }*/
}
    
