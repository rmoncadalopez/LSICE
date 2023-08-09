
#ifndef LSEVOLVE_H_
#define LSEVOLVE_H_

#include "definitions.h"
#include "levelsetElementaryFunctions.h"

ArrayXd LSEvolve( const ArrayXd & g,
						const array<ArrayXd, 2> & lsIn,
						const int & xd,
						const int & yd,
						const int & zd,
						const int & maxiter,
						const double & dt,
						const double & mu,
						const double & lambda,
						const double & alpha,
						const double & beta) {
	
	array<ArrayXd,2> ls = lsIn;
	ArrayXd px(g.size());
	ArrayXd py(g.size());
	ArrayXd pz(g.size());
	ArrayXd pgmag(g.size());
	ArrayXd npx(g.size());
	ArrayXd npy(g.size());
	ArrayXd npz(g.size());
	ArrayXd dps;
	array<ArrayXd, 2> d;
	array<ArrayXd, 2> h;
	ArrayXd f1, f2, f3, f41, f42;
	f1.resize(g.size());
	f2.resize(g.size());
	
	for (size_t i = 0; i < (size_t)maxiter; i++) {
//		cout << i << endl;
		for (int j = 0; j < 2; j++) {
			gradient(ls[j],xd,yd,zd,px,py,pz);
			pgmag = sqrt(px*px + py*py + pz*pz) + 1e-12;
			npx = px/pgmag;
			npy = py/pgmag;
			npz = pz/pgmag;
			d[j] = sDirac(ls[j]);
			h[j] = sHeavi(-ls[j]);
			
			dps = dp(pgmag);
			f1 = mu*divergence(dps*px, dps*py, dps*pz, xd, yd, zd);
			f2 = lambda*d[j]*divergence(g*npx, g*npy, g*npz, xd, yd, zd);
			f3 = alpha*g*d[j];
			
			ls[j] += dt*(f1+f2+f3);
		}
		f41 = beta*d[0]*h[1];
		f42 = beta*d[1]*h[0];
		ls[0] += dt*f41;
		ls[1] += dt*f42;
		
	}
	
	return ls[0];
}


#endif
