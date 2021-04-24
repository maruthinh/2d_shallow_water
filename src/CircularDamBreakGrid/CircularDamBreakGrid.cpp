#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>


using namespace std;

int const Nx=100, Ny=100;
double **x, **y;
double dx,dy,xll=0.0,xul=50.0,yll=0.0,yul=50.0;

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


void WriteToFile(){

	ofstream GridCircularDamBreakGrid("CircularDamBreakGrid100100.dat");
	ofstream GridCircularDamBreakGridTecPlot("CircularDamBreakGridTecPlot100100.dat");
	//GridCircularDamBreakGrid12040.flags(ios::dec | ios::scientific);
	//GridCircularDamBreakGrid12040.precision(5);

	if(!GridCircularDamBreakGrid){
		cerr<<"File couldn't be opened to write the solution"<<endl;
		exit(1);
	}
	
	GridCircularDamBreakGridTecPlot << "TITLE = Flow" << endl << "VARIABLES = xc, yc, rho" << endl;
	GridCircularDamBreakGridTecPlot << "Zone T = Omega I = " << Nx+1 << " J = " << Ny+1 << endl ;
	int count=0;

	for(int j=1;j<=Ny+1;j++){
		for(int i=1;i<=Nx+1;i++){
			count=count+1;
			GridCircularDamBreakGrid<<count<<"\t"<<x[i][j]<<"\t"<<y[i][j]<<"\t"<<endl;
			GridCircularDamBreakGridTecPlot<<x[i][j]<<"\t"<<y[i][j]<<"\t"<<0.0<<endl;

		}
	}

	GridCircularDamBreakGrid.close();
	GridCircularDamBreakGridTecPlot.close();

}

int main(){
	
	Allocate_2D(x,Nx+2, Ny+2); Allocate_2D(y,Nx+2, Ny+2);
	dx=(xul-xll)/(Nx); dy=(yul-yll)/(Ny);
	
	
	
	for (int i=1;i<=Nx+1;i++){
	    for(int j=1;j<=Ny+1;j++){
	    
		x[i][j]=xll+(i-1)*dx;
		y[i][j]=yll+(j-1)*dy;  
		
	    }	  
	  
	}
	WriteToFile();
	return 0;
}
