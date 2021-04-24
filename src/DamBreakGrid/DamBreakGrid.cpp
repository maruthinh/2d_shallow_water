#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>


using namespace std;

int const Nx=40, Ny=40;
double **x, **y;
double dx,dy,xll=0.0,xul=200.0,yll=0.0,yul=200.0;

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

	ofstream GridDamBreak("DamBreakGrid4040.dat");
	ofstream GridDamBreakTecPlot("DamBreakTecPlot4040.dat");
	//GridDamBreak12040.flags(ios::dec | ios::scientific);
	//GridDamBreak12040.precision(5);

	if(!GridDamBreak){
		cerr<<"File couldn't be opened to write the solution"<<endl;
		exit(1);
	}
	
	GridDamBreakTecPlot << "TITLE = Flow" << endl << "VARIABLES = xc, yc, rho" << endl;
	GridDamBreakTecPlot << "Zone T = Omega I = " << Nx+1 << " J = " << Ny+1 << endl ;
	int count=0;

	for(int j=1;j<=Ny+1;j++){
		for(int i=1;i<=Nx+1;i++){
			count=count+1;
			GridDamBreak<<count<<"\t"<<x[i][j]<<"\t"<<y[i][j]<<"\t"<<endl;
			GridDamBreakTecPlot<<x[i][j]<<"\t"<<y[i][j]<<"\t"<<0.0<<endl;

		}
	}

	GridDamBreak.close();
	GridDamBreakTecPlot.close();

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
