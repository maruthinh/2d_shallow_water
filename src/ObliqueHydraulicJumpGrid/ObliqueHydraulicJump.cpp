#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>


using namespace std;

//Global variable declaration
int const Nx1=10, Nx2=30, Ny=30;
double **x, **y, *yvar;
double ax=0.0, bx=10.0, cx=40, ay=0.0, by=30;
double dx1, dx2, dy1, dy2;

//Dynamic memory allocation
template <typename T>
T **Allocate_2D(T ** &m, int t1, int t2) {
	m=new T* [t1];
        for (int i=0; i<t1; ++i) {
                m[i]=new T [t2];
                for (int j=0; j<t2; ++j)
                        m[i][j]=0.0;
        }
        return m;
}

//Funtion to write the grid points to file 
void WriteToFile(){

	ofstream GridObliqueHydraulicJump("ObliqueHydraulicJump4030.dat");
	ofstream GridObliqueHydraulicJumpTecPlot("ObliqueHydraulicJumpTecPlot4030.dat");
	//GridObliqueHydraulicJump12040.flags(ios::dec | ios::scientific);
	//GridObliqueHydraulicJump12040.precision(5);

	if(!GridObliqueHydraulicJump){
		cerr<<"File couldn't be opened to write the solution"<<endl;
		exit(1);
	}
	
	GridObliqueHydraulicJumpTecPlot << "TITLE = Flow" << endl << "VARIABLES = xc, yc, rho" << endl;
	GridObliqueHydraulicJumpTecPlot << "Zone T = Omega I = " << Nx1+Nx2+1 << " J = " << Ny+1 << endl ;
	double count;
	for(int j=1;j<=Ny+1;j++){
		for(int i=1;i<=Nx1+Nx2+1;i++){
			count=count+1;
			GridObliqueHydraulicJump<<count<<"\t"<<x[i][j]<<"\t"<<y[i][j]<<"\t"<<endl;
			GridObliqueHydraulicJumpTecPlot<<x[i][j]<<"\t"<<y[i][j]<<"\t"<<0.0<<endl;

		}
	}

	GridObliqueHydraulicJump.close();
	GridObliqueHydraulicJumpTecPlot.close();


}

//Main function
int main(){
	
	Allocate_2D(x, Nx1+Nx2+2, Ny+2); Allocate_2D(y, Nx1+Nx2+2, Ny+2); yvar=new double[Nx1+Nx2+1];
	dx1=(bx-ax)/(Nx1); dx2=(cx-bx)/(Nx2); dy1=(by-ay)/(Ny);
	double h2, opp;
	
	opp=(30*tan((M_PI*8.95)/180.0));
	h2=(30-opp);
	dy2=h2/Ny;
	
	for (int j=1;j<=Ny+1;j++){
	    for(int i=1;i<=Nx1+1;i++){
	    
		x[i][j]=ax+(i-1)*dx1;
		y[i][j]=ay+(j-1)*dy1;  
		
	    }	  
	  
	}
/*	
	int k1=Nx1+Nx2+1;
	for (int j=1;j<=Ny+1;j++){
	    for(int i=Nx1+Nx2+1;i<=Nx1+Nx2+Nx3+1;i++){
	    
		x[i][j]=cx+(i-k1)*dx3;
		y[i][j]=opp+(j-1)*dy3;  
		
	    }	  
	  
	}
*/	
	int k2=Nx1+1;
	for (int j=1;j<=Ny+1;j++){
	    for(int i=Nx1+1;i<=Nx1+Nx2+1;i++){
	    
		x[i][j]=bx+(i-k2)*dx2;
		//yvar[i]=(1.0-(x[i][1]-x[k2][1])*tan(M_PI/7.0))/Ny;
		dy2=(30-(x[i][1]-x[k2][1])*tan((M_PI*8.95)/180.0))/Ny;
		y[i][j]=(x[i][1]-x[k2][1])*tan((M_PI*8.95)/180.0)+(j-1)*dy2;
				
	    }	  
	  
	}

	WriteToFile();
	return 0;
}
