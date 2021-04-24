#include "global_declarations.h"
#include "basic_functions.h"

template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min);

template <typename T>
T MaxEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny);

template <typename T>
T MinEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny);


void Diss_MOVERS(int ib, int id1, int jb, int jd1, double ***&cv, double ***&dv, double ***&si, double ***&sj, double ***&diss){
    
    int im1, jm1;
    double ds, nx, ny, hinv, hl, ul, vl, hr, ur, vr;
    double ar, al, max_eig, min_eig;
    double *fd, *fr, *fl, *Ur, *Ul;
    
   fd = new double [nconv]; fr = new double [nconv]; fl = new double [nconv]; Ur = new double [nconv]; Ul = new double [nconv]; 
    
    //compute in i-direction
    for(int j=2;j<=jb;j++){
        for(int i=2;i<=id1;i++){
            
            im1=i-1;
            
            ds = sqrt(si[0][i][j]*si[0][i][j]+si[1][i][j]*si[1][i][j]);
            nx = si[0][i][j]/ds;
            ny = si[1][i][j]/ds;
            
            //left and right states
            
            hinv = 1.0/cv[0][im1][j];
            hl   = cv[0][im1][j];
            ul   = cv[1][im1][j]*hinv;
            vl   = cv[2][im1][j]*hinv;
            al   = sqrt(g*hl);
            
            hinv = 1.0/cv[0][i][j];
            hr   = cv[0][i][j];
            ur   = cv[1][i][j]*hinv;
            vr   = cv[2][i][j]*hinv;
            ar   = sqrt(g*hr);
            
            max_eig = MaxEigVal(ur, ul, vr, vl, ar, al, nx, ny);
            min_eig = MinEigVal(ur, ul, vr, vl, ar, al, nx, ny);
            
            fr[0] = (hr*ur)*nx+(hr*vr)*ny;
            fr[1] = (hr*ur*ur+0.5*g*hr*hr)*nx + (hr*ur*vr)*ny;
            fr[2] = (hr*ur*vr)*nx + (hr*vr*vr+0.5*g*hr*hr)*ny;
            
            fl[0] = (hl*ul)*nx+(hl*vl)*ny;
            fl[1] = (hl*ul*ul+0.5*g*hl*hl)*nx + (hl*ul*vl)*ny;
            fl[2] = (hl*ul*vl)*nx + (hl*vl*vl+0.5*g*hl*hl)*ny;
            
            Ur[0] = hr;
            Ur[1] = hr*ur;
            Ur[2] = hr*vr;
            
            Ul[0] = hl;
            Ul[1] = hl*ul;
            Ul[2] = hl*vl;
            
            fd[0] = 0.5*Movers(fr[0], fl[0], Ur[0], Ul[0], max_eig, min_eig)*(Ur[0]-Ul[0]);
            fd[1] = 0.5*Movers(fr[1], fl[1], Ur[1], Ul[1], max_eig, min_eig)*(Ur[1]-Ul[1]);
            fd[2] = 0.5*Movers(fr[2], fl[2], Ur[2], Ul[2], max_eig, min_eig)*(Ur[2]-Ul[2]);
            
            //final dissipation terms
            
            diss[0][i][j] = diss[0][i][j] - fd[0]*ds;
            diss[1][i][j] = diss[1][i][j] - fd[1]*ds;
            diss[2][i][j] = diss[2][i][j] - fd[2]*ds;
            
            diss[0][im1][j] = diss[0][im1][j] + fd[0]*ds;
            diss[1][im1][j] = diss[1][im1][j] + fd[1]*ds;
            diss[2][im1][j] = diss[2][im1][j] + fd[2]*ds;                                 
        }
    }
    
       //compute in j-direction
    for(int i=2;i<=ib;i++){
        for(int j=2;j<=jd1;j++){
            
            jm1=j-1;
            
            ds = sqrt(sj[0][i][j]*sj[0][i][j]+sj[1][i][j]*sj[1][i][j]);
            nx = sj[0][i][j]/ds;
            ny = sj[1][i][j]/ds;
            
            //left and right states
            
            hinv = 1.0/cv[0][i][jm1];
            hl   = cv[0][i][jm1];
            ul   = cv[1][i][jm1]*hinv;
            vl   = cv[2][i][jm1]*hinv;
            al   = sqrt(g*hl);
            
            hinv = 1.0/cv[0][i][j];
            hr   = cv[0][i][j];
            ur   = cv[1][i][j]*hinv;
            vr   = cv[2][i][j]*hinv;
            ar   = sqrt(g*hr);
            
            max_eig = MaxEigVal(ur, ul, vr, vl, ar, al, nx, ny);
            min_eig = MinEigVal(ur, ul, vr, vl, ar, al, nx, ny);
            
            fr[0] = (hr*ur)*nx+(hr*vr)*ny;
            fr[1] = (hr*ur*ur+0.5*g*hr*hr)*nx + (hr*ur*vr)*ny;
            fr[2] = (hr*ur*vr)*nx + (hr*vr*vr+0.5*g*hr*hr)*ny;
            
            fl[0] = (hl*ul)*nx+(hl*vl)*ny;
            fl[1] = (hl*ul*ul+0.5*g*hl*hl)*nx + (hl*ul*vl)*ny;
            fl[2] = (hl*ul*vl)*nx + (hl*vl*vl+0.5*g*hl*hl)*ny;
            
            Ur[0] = hr;
            Ur[1] = hr*ur;
            Ur[2] = hr*vr;
            
            Ul[0] = hl;
            Ul[1] = hl*ul;
            Ul[2] = hl*vl;
            
            fd[0] = 0.5*Movers(fr[0], fl[0], Ur[0], Ul[0], max_eig, min_eig)*(Ur[0]-Ul[0]);
            fd[1] = 0.5*Movers(fr[1], fl[1], Ur[1], Ul[1], max_eig, min_eig)*(Ur[1]-Ul[1]);
            fd[2] = 0.5*Movers(fr[2], fl[2], Ur[2], Ul[2], max_eig, min_eig)*(Ur[2]-Ul[2]);
            
            //final dissipation terms
            
            diss[0][i][j] = diss[0][i][j] - fd[0]*ds;
            diss[1][i][j] = diss[1][i][j] - fd[1]*ds;
            diss[2][i][j] = diss[2][i][j] - fd[2]*ds;
            
            diss[0][i][jm1] = diss[0][i][jm1] + fd[0]*ds;
            diss[1][i][jm1] = diss[1][i][jm1] + fd[1]*ds;
            diss[2][i][jm1] = diss[2][i][jm1] + fd[2]*ds;       
        }
    }
    
    /*for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            std::cout<<"the computed dissipation values="<<diss[2][i][j]<<std::endl;
        }
    }*/
    
    delete [] fd;
    delete [] fr;
    delete [] fl;
    delete [] Ur;
    delete [] Ul;
}


template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){

	double S;
	const double epsilon=1e-4;

	if (fabs(Fr-Fl)<epsilon) return 0.0;
	else if (fabs(Ur-Ul)<epsilon) return L_min;
	else if (fabs(Ur-Ul)>epsilon and fabs(Fr-Fl)>epsilon) S=fabs(((Fr-Fl)/(Ur-Ul)));
	else S=L_min;
	
	if (S<epsilon)	return 0;
	else if ((S)>=L_max) return (L_max);
	else if((S)<=L_min) return (L_min);
	else return (S);
	return S;
}

template <typename T>
T MaxEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny){
	double L1r, L2r, L3r, L1l, L2l, L3l;
		L1r=fabs(ur*nx+vr*ny+ar);  L1l=fabs(ul*nx+vl*ny+al);
		L2r=fabs(ur*nx+vr*ny); 	   L2l=fabs(ul*nx+vl*ny);
		L3r=fabs(ur*nx+vr*ny-ar);  L3l=fabs(ul*nx+vl*ny-al);

	//return Max3(max(L1l,L1r), max(L2l,L2r), max(L3l,L3r));
	return Max2(Max3(L1r, L2r, L3r), Max3(L1l, L2l, L3l));
}

template <typename T>
T MinEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny){
	double L1r, L2r, L3r, L1l, L2l, L3l;
		L1r=fabs(ur*nx+vr*ny+ar);  L1l=fabs(ul*nx+vl*ny+al);
		L2r=fabs(ur*nx+vr*ny); 	   L2l=fabs(ul*nx+vl*ny);
		L3r=fabs(ur*nx+vr*ny-ar);  L3l=fabs(ul*nx+vl*ny-al);

	//return Max3(min(L1l,L1r), min(L2l,L2r), min(L3l,L3r));
	return Max2(Min3(L1r, L2r, L3r), Min3(L1l, L2l, L3l));
}

