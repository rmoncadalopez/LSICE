/*
 * levelsetGrainFunctions.h
 *
 *  Created on: Sep 24, 2014
 *      Author: Reid
 */

#ifndef LEVELSETGRAINFUNCTIONS_H_
#define LEVELSETGRAINFUNCTIONS_H_

// contains grain methods that operate on level sets such as center of mass, discretization, and moment of inertia


// discretizes the grain into points in the current configuration of the grain with respect to the level set
vector<Vector3d> discretizeGrain(const Levelset & lset, const size_t & npoints, const size_t & niter, const double & damping, const size_t & id) {
	
	vector<Vector3d> points(npoints);

	// place points
	vector<Vector3d> zeroCrossings;
	for (size_t k = 0; k < lset.getZdim()-1; k++) {
		for (size_t j = 0; j < lset.getYdim()-1; j++) {
			for (size_t i = 0; i < lset.getXdim()-1; i++) {
				double val = lset(i,j,k);
				if (val*lset(i,j,k+1) < 0) {
					zeroCrossings.push_back(Vector3d(double(i), double(j), double(k)+fabs(val)/(fabs(val)+fabs(lset(i,j,k+1)))));
				}
				if (val*lset(i,j+1,k) < 0) {
					zeroCrossings.push_back(Vector3d(double(i), double(j)+fabs(val)/(fabs(val) + fabs(lset(i,j+1,k))) , double(k)));
				}
				if (val*lset(i+1,j,k) < 0) {
					zeroCrossings.push_back(Vector3d(double(i)+fabs(val)/(fabs(val) + fabs(lset(i+1,j,k))), double(j), double(k)));
				}
			}
		}
	}
	
	if (npoints > zeroCrossings.size()) {
		cout << "grain " << id << " is a degenerate grain, number of points less than possible points" << endl;
		return points;
	}
	
	double step = double(zeroCrossings.size())/double(npoints);
	
	for (size_t i = 0; i < npoints; i++) {
		points[i] = zeroCrossings[size_t(double(i)*step)];
	}
	
	// set up temp vars for iteration
	vector<Vector3d> glist(npoints); vector<double> sdlist(npoints);
	Vector3d grad; double sd;
	Vector3d fij; double dij; Vector3d nij; Vector3d nji; double magp;
	vector<Vector3d> forces(npoints);
	vector<Vector3d> v(npoints);
	for (size_t i = 0; i < npoints; i++) {
		v[i] << 0., 0., 0.;
	}
	
	double sa = sDirac(lset.getGrid()).sum();
	double dt = sqrt(sa/npoints)*0.1;
	// begin iteration
	for (size_t t = 0; t < niter; t++ ) { // niter
		// store all of the sds and gradients first, also zero out the force
		for (size_t i = 0; i < npoints; i++) {
			if (lset.findValueGradient( points[i], sd, grad)) {
				glist[i] = grad;
				sdlist[i] = sd;
			}
			else {
				glist[i] << 0,0,0;
				sdlist[i] = -1;
				points[i] << 5+(double)rand()/RAND_MAX,5+(double)rand()/RAND_MAX,5+(double)rand()/RAND_MAX;
			}
			forces[i] << 0., 0., 0.;
		}
		// iterate through points
		for (size_t i = 0; i < npoints; i++) {
			// drag points closer to surface using gradient/signed distance
			points[i] -= .8*sdlist[i]*glist[i];
			// push points that are close away from each other
			for (size_t j = i+1; j < npoints; j++) {
				if (i == j) { break; }
				fij = points[i] - points[j];
				dij = fij.norm();
				if (dij < 1.1*sqrt(sa/npoints)) {
					nij =  fij - glist[i]*fij.dot(glist[i]);
					nji = -fij + glist[j]*fij.dot(glist[j]);
					magp = exp(-dij/sqrt(sa/npoints));
					
					forces[i] += nij/(nij.norm()+DBL_MIN)*magp;
					forces[j] += nji/(nji.norm()+DBL_MIN)*magp;
				}
			}
		}
		for (size_t i = 0; i < npoints; i++) {
			v[i] = 1/(1+damping*dt/2)*((1-damping*dt/2)*v[i] + dt*forces[i]);
			points[i] += dt*v[i];
		}
	}
	
	// additional iterations to move points to zero contour
	for (int tt = 0; tt < 200; tt++) {
		for (size_t i = 0; i < npoints; i++) {
			if (lset.findValueGradient( points[i], sd, grad)) {
				points[i] -= .8*sd*grad;
			}
		}
	}
	
//	double energy = 0;
//	for (size_t i = 0; i < npoints; i++) {
//		if (lset.findValueGradient( points[i], sd, grad)) {
//			if (!isnan(sd)) {
//				energy += fabs(sd);
//			}
//		}
//	}
//	cout << energy << endl;
	
	return points;
}



struct GrainInfo {
	double _mass;
	double _radius;
	Vector3d _moin;
	Vector3d _cm;
	Vector3d _cmLset;
	Vector3d _cmGlobal;
	Vector4d _quat;
	Matrix3d _R;
	vector<Vector3d> _points; // points in the principal frame
	Levelset _lset;
//	Vector3d _bboxShift; // how much trimming the level set moves its bounding box
};

// computes the fields in the above struct definition.  also moves the points into the principal frame from the current frame
GrainInfo computeGrainInfoNoTrimNoPoints( const Levelset & lset, const Vector3d & bboxCorner) {
	
	
	GrainInfo info;
// compute mass and center of mass (assuming a density of 1)
	double vol = 0;
	Vector3d cm(0., 0., 0.);
	ArrayXd heavi = sHeavi(-lset.getGrid());
	size_t iter = 0;
	for (size_t k = 0; k < lset.getZdim(); k++) {
		for (size_t j = 0; j < lset.getYdim(); j++) {
			for (size_t i = 0; i < lset.getXdim(); i++) {
				vol += heavi(iter);
				cm(0) += heavi(iter)*(double)i;
				cm(1) += heavi(iter)*(double)j;
				cm(2) += heavi(iter)*(double)k;
				iter++;
			}
		}
	}
	cm /= vol;
	info._mass = vol;
	info._cm = cm;
	
// compute moment of inertia and rotation from principal frame -> current frame
	Matrix3d I = Matrix3d::Zero();

	double rx, ry, rz;
//	ArrayXd heavi2 = pow(heavi, 1.082);
	ArrayXd heavi2 = heavi;
	
	//	ArrayXd rx(heavi.size()); ArrayXd ry(heavi.size()); ArrayXd rz(heavi.size());
	iter = 0;
	for (size_t k = 0; k < lset.getZdim(); k++) {
		for (size_t j = 0; j < lset.getYdim(); j++) {
			for (size_t i = 0; i < lset.getXdim(); i++) {
				rx = (double)i - cm(0);
				ry = (double)j - cm(1);
				rz = (double)k - cm(2);
				I(0,0) += heavi2(iter)*(ry*ry + rz*rz);
				I(1,1) += heavi2(iter)*(rx*rx + rz*rz);
				I(2,2) += heavi2(iter)*(rx*rx + ry*ry);
				I(0,1) += -heavi2(iter)*rx*ry;
				I(0,2) += -heavi2(iter)*rx*rz;
				I(1,2) += -heavi2(iter)*ry*rz;
				iter++;
			}
		}
	}
	I(1,0) = I(0,1);
	I(2,0) = I(0,2);
	I(2,1) = I(1,2);
	// find eigenvalues and eigenvectors of moment of inertia matrix
	// which are the principle moment of inertia and the rotation from principal -> current respectively
	EigenSolver<Matrix3d> eigs(I);
	Matrix3d evec; Vector3d eval;
	evec = eigs.eigenvectors().real();
	info._moin = eigs.eigenvalues().real();
	if (evec.determinant() < 0) {
		evec = -evec;
	}
	
	
	info._R = evec;
	// convert rotation matrix into quaternion
	double tr = evec.trace();
	info._quat(3) = sqrt(tr+1.)/2.;
	info._quat(0) = (evec(0,2)-evec(2,0))/(4*info._quat(3));
	info._quat(1) = (evec(1,2)-evec(2,1))/(4*info._quat(3));
	info._quat(2) = (evec(0,1)-evec(1,0))/(4*info._quat(3));
	
	// TODO: rotate the level set into the principal frame
	size_t maxDiam = size_t(ceil(sqrt(double((lset.getXdim()+1)*(lset.getXdim()+1) + 
														  (lset.getYdim()+1)*(lset.getYdim()+1) +
														  (lset.getZdim()+1)*(lset.getZdim()+1)))));
	ArrayXd grid(maxDiam*maxDiam*maxDiam);
	Levelset lsetPrin(grid, maxDiam, maxDiam, maxDiam);
	Vector3d cmPrin( double(maxDiam-1)/2., double(maxDiam-1)/2.,double(maxDiam-1)/2.);
	Vector3d pt;
	for (size_t k = 0; k < maxDiam; k++) {
		for (size_t j = 0; j < maxDiam; j++) {
			for (size_t i = 0; i < maxDiam; i++) {
				pt << (double)i, (double)j, (double)k;
				lsetPrin(i,j,k) = lset.findValueOnly( evec*(pt-cmPrin)+ cm);
			}
		}
	}
//	lset = lsetPrin;
	info._cmLset = cmPrin;
	info._cmGlobal = bboxCorner + cm;
	
//	computes the maximum radius of the grain
	double maxradius = 0;
	info._radius = maxradius;
	info._lset = lsetPrin;
	return info;
}

// finds the maximum radius given a set of points
double computeRadius(vector<Vector3d> points) {
	
	double maxradius = 0;
	for (size_t i = 0; i < points.size(); i++) {
		if (points[i].norm() > maxradius) {
			maxradius = points[i].norm();
//			cout << maxradius << endl;
		}
	}
	return maxradius;
	
}

GrainInfo computeGrainInfo( const Levelset & lset, const vector<Vector3d> & points, const Vector3d & bboxCorner) {
	
	
	GrainInfo info;
// compute mass and center of mass (assuming a density of 1)
	double vol = 0;
	Vector3d cm(0., 0., 0.);
	ArrayXd heavi = sHeavi(-lset.getGrid());
	size_t iter = 0;
	for (size_t k = 0; k < lset.getZdim(); k++) {
		for (size_t j = 0; j < lset.getYdim(); j++) {
			for (size_t i = 0; i < lset.getXdim(); i++) {
				vol += heavi(iter);
				cm(0) += heavi(iter)*(double)i;
				cm(1) += heavi(iter)*(double)j;
				cm(2) += heavi(iter)*(double)k;
				iter++;
			}
		}
	}
	cm /= vol;
	info._mass = vol;
//	info._cm = cm;
	
// compute moment of inertia and rotation from principal frame -> current frame
	Matrix3d I = Matrix3d::Zero();

	double rx, ry, rz;
//	ArrayXd heavi2 = pow(heavi, 1.082);
	ArrayXd heavi2 = heavi;
	
	//	ArrayXd rx(heavi.size()); ArrayXd ry(heavi.size()); ArrayXd rz(heavi.size());
	iter = 0;
	for (size_t k = 0; k < lset.getZdim(); k++) {
		for (size_t j = 0; j < lset.getYdim(); j++) {
			for (size_t i = 0; i < lset.getXdim(); i++) {
				rx = (double)i - cm(0);
				ry = (double)j - cm(1);
				rz = (double)k - cm(2);
				I(0,0) += heavi2(iter)*(ry*ry + rz*rz);
				I(1,1) += heavi2(iter)*(rx*rx + rz*rz);
				I(2,2) += heavi2(iter)*(rx*rx + ry*ry);
				I(0,1) += -heavi2(iter)*rx*ry;
				I(0,2) += -heavi2(iter)*rx*rz;
				I(1,2) += -heavi2(iter)*ry*rz;
				iter++;
			}
		}
	}
	I(1,0) = I(0,1);
	I(2,0) = I(0,2);
	I(2,1) = I(1,2);
	// find eigenvalues and eigenvectors of moment of inertia matrix
	// which are the principle moment of inertia and the rotation from principal -> current respectively
	EigenSolver<Matrix3d> eigs(I);
	Matrix3d evec; Vector3d eval;
	evec = eigs.eigenvectors().real();
	info._moin = eigs.eigenvalues().real();
	if (evec.determinant() < 0) {
		evec = -evec;
	}
	
	
	info._R = evec;
	// convert rotation matrix into quaternion
	double tr = evec.trace();
	info._quat(3) = sqrt(tr+1.)/2.;
	info._quat(0) = (evec(0,2)-evec(2,0))/(4*info._quat(3));
	info._quat(1) = (evec(1,2)-evec(2,1))/(4*info._quat(3));
	info._quat(2) = (evec(0,1)-evec(1,0))/(4*info._quat(3));
	
	// TODO: rotate the level set into the principal frame
	size_t maxDiam = size_t(ceil(sqrt(double((lset.getXdim()+1)*(lset.getXdim()+1) + 
														  (lset.getYdim()+1)*(lset.getYdim()+1) +
														  (lset.getZdim()+1)*(lset.getZdim()+1)))));
	ArrayXd grid(maxDiam*maxDiam*maxDiam);
	Levelset lsetPrin(grid, maxDiam, maxDiam, maxDiam);
	Vector3d cmPrin( double(maxDiam-1)/2., double(maxDiam-1)/2.,double(maxDiam-1)/2.);
	Vector3d pt;
	for (size_t k = 0; k < maxDiam; k++) {
		for (size_t j = 0; j < maxDiam; j++) {
			for (size_t i = 0; i < maxDiam; i++) {
				pt << (double)i, (double)j, (double)k;
				lsetPrin(i,j,k) = lset.findValueOnly( evec*(pt-cmPrin)+ cm);
			}
		}
	}
//	lset = lsetPrin;
	info._cmLset = cmPrin;
	info._cmGlobal = bboxCorner + cm;
	
	
//	rotate the points into the principal frame 
	vector<Vector3d> pointsPrin(points.size());
	for (size_t i = 0; i < points.size(); i++) {
		pointsPrin[i] = evec.transpose()*(points[i] - cm);
	}
	info._points = pointsPrin;
	
//	computes the maximum radius of the grain
	double maxradius = 0;
	for (size_t i = 0; i < points.size(); i++) {
		if (pointsPrin[i].norm() > maxradius) {
			maxradius = pointsPrin[i].norm();
//			cout << maxradius << endl;
		}
	}
	info._radius = maxradius;
	

	// finally, trim the level set (this changes the center of mass and the level set)
	int zmin = 0; int ymin = 0; int xmin = 0;
	// z min
	for (size_t k = 0; k < lsetPrin.getZdim(); k++) {
		for (size_t j = 0; j < lsetPrin.getYdim(); j++) {
			for (size_t i = 0; i < lsetPrin.getXdim(); i++) {
				if (lsetPrin(i,j,k) < 0) {
					zmin = k-1;
					goto finishzmin;
				}
			}
		}
	}
	finishzmin:
	// y min
	for (size_t j = 0; j < lsetPrin.getYdim(); j++) {
		for (size_t k = 0; k < lsetPrin.getZdim(); k++) {
			for (size_t i = 0; i < lsetPrin.getXdim(); i++) {
				if (lsetPrin(i,j,k) < 0) {
					ymin = j-1;
					goto finishymin;
				}
			}
		}
	}
	finishymin:
	// x min
	for (size_t i = 0; i < lsetPrin.getXdim(); i++) {
		for (size_t k = 0; k < lsetPrin.getZdim(); k++) {
			for (size_t j = 0; j < lsetPrin.getYdim(); j++) {
				if (lsetPrin(i,j,k) < 0) {
					xmin = i-1;
					goto finishxmin;
				}
			}
		}
	}
	finishxmin:
//	cout << "hello" << endl;
	// find maxes
	size_t zmax = lsetPrin.getZdim()-1; size_t ymax = lsetPrin.getYdim()-1; size_t xmax = lsetPrin.getXdim()-1;
	// zmax
	for (size_t k = lsetPrin.getZdim()-1; k > 0; k--) {
		for (size_t j = 0; j < lsetPrin.getYdim(); j++) {
			for (size_t i = 0; i < lsetPrin.getXdim(); i++) {
//				cout << i << " " << j << " " << k << endl;
				if (lsetPrin(i,j,k) < 0) {
					zmax = k+1;
					goto finishzmax;
				}
			}
		}
	}
	finishzmax:
	// ymax
	for (size_t j = lsetPrin.getYdim()-1; j > 0; j--) {
		for (size_t k = 0; k < lsetPrin.getZdim(); k++) {
			for (size_t i = 0; i < lsetPrin.getXdim(); i++) {
				if (lsetPrin(i,j,k) < 0) {
					ymax = j+1;
					goto finishymax;
				}
			}
		}
	}
	finishymax:
	// xmax
	for (size_t i = lsetPrin.getXdim()-1; i > 0; i--) {
		for (size_t k = 0; k < lsetPrin.getZdim(); k++) {
			for (size_t j = 0; j < lsetPrin.getYdim(); j++) {
				if (lsetPrin(i,j,k) < 0) {
					xmax = i+1;
					goto finishxmax;
				}
			}
		}
	}
	finishxmax:
	
	
	
	zmin = max(zmin,0);
	ymin = max(ymin,0);
	xmin = max(xmin,0);
	
	zmax = min(zmax, lsetPrin.getZdim()-1);
	ymax = min(ymax, lsetPrin.getYdim()-1);
	xmax = min(xmax, lsetPrin.getXdim()-1);
	
	size_t zdim = zmax-zmin+1;
	size_t ydim = ymax-ymin+1;
	size_t xdim = xmax-xmin+1;
	
	ArrayXd trim(1); trim(0) = 5;
	if (xdim > 0 && ydim > 0 && zdim > 0) {
		trim.resize(zdim*ydim*xdim);
		iter = 0;
		for (size_t k = zmin; k < zmax+1; k++) {
			for (size_t j = ymin; j < ymax+1; j++) {
				for (size_t i = xmin; i < xmax+1; i++) {
					trim(iter) = lsetPrin(i,j,k);
					iter++;
				}
			}
		}
	}
	
	Levelset lstrim(trim, xdim, ydim, zdim);
	info._cmLset -= Vector3d(xmin,ymin,zmin);
	info._lset = lstrim;
//	info._lset = lset;
//	info._bboxShift = Vector3d(xmin,ymin,zmin);
	
	return info;
}

GrainInfo trim( const GrainInfo & info0) { // const Levelset & lset, const vector<Vector3d> & points, const Vector3d & bboxCorner) {
	
	size_t iter = 0;
	
	Levelset lset = info0._lset;
	vector<Vector3d> points = info0._points;
//	Vector3d bboxCorner = info0._cmGlobal - info0._cmLset;
	
	GrainInfo info = info0;

	// finally, trim the level set (this changes the center of mass and the level set)
	int zmin = 0; int ymin = 0; int xmin = 0;
	// z min
	for (size_t k = 0; k < lset.getZdim(); k++) {
		for (size_t j = 0; j < lset.getYdim(); j++) {
			for (size_t i = 0; i < lset.getXdim(); i++) {
				if (lset(i,j,k) < 0) {
					zmin = k-1;
					goto finishzmin;
				}
			}
		}
	}
	finishzmin:
	// y min
	for (size_t j = 0; j < lset.getYdim(); j++) {
		for (size_t k = 0; k < lset.getZdim(); k++) {
			for (size_t i = 0; i < lset.getXdim(); i++) {
				if (lset(i,j,k) < 0) {
					ymin = j-1;
					goto finishymin;
				}
			}
		}
	}
	finishymin:
	// x min
	for (size_t i = 0; i < lset.getXdim(); i++) {
		for (size_t k = 0; k < lset.getZdim(); k++) {
			for (size_t j = 0; j < lset.getYdim(); j++) {
				if (lset(i,j,k) < 0) {
					xmin = i-1;
					goto finishxmin;
				}
			}
		}
	}
	finishxmin:
//	cout << "hello" << endl;
	// find maxes
	size_t zmax = lset.getZdim()-1; size_t ymax = lset.getYdim()-1; size_t xmax = lset.getXdim()-1;
	// zmax
	for (size_t k = lset.getZdim()-1; k > 0; k--) {
		for (size_t j = 0; j < lset.getYdim(); j++) {
			for (size_t i = 0; i < lset.getXdim(); i++) {
//				cout << i << " " << j << " " << k << endl;
				if (lset(i,j,k) < 0) {
					zmax = k+1;
					goto finishzmax;
				}
			}
		}
	}
	finishzmax:
	// ymax
	for (size_t j = lset.getYdim()-1; j > 0; j--) {
		for (size_t k = 0; k < lset.getZdim(); k++) {
			for (size_t i = 0; i < lset.getXdim(); i++) {
				if (lset(i,j,k) < 0) {
					ymax = j+1;
					goto finishymax;
				}
			}
		}
	}
	finishymax:
	// xmax
	for (size_t i = lset.getXdim()-1; i > 0; i--) {
		for (size_t k = 0; k < lset.getZdim(); k++) {
			for (size_t j = 0; j < lset.getYdim(); j++) {
				if (lset(i,j,k) < 0) {
					xmax = i+1;
					goto finishxmax;
				}
			}
		}
	}
	finishxmax:
	
	
	
	zmin = max(zmin,0);
	ymin = max(ymin,0);
	xmin = max(xmin,0);
	
	zmax = min(zmax, lset.getZdim()-1);
	ymax = min(ymax, lset.getYdim()-1);
	xmax = min(xmax, lset.getXdim()-1);
	
	size_t zdim = zmax-zmin+1;
	size_t ydim = ymax-ymin+1;
	size_t xdim = xmax-xmin+1;
	
	ArrayXd trim(1); trim(0) = 5;
	if (xdim > 0 && ydim > 0 && zdim > 0) {
		trim.resize(zdim*ydim*xdim);
		iter = 0;
		for (size_t k = zmin; k < zmax+1; k++) {
			for (size_t j = ymin; j < ymax+1; j++) {
				for (size_t i = xmin; i < xmax+1; i++) {
					trim(iter) = lset(i,j,k);
					iter++;
				}
			}
		}
	}
	
	Levelset lstrim(trim, xdim, ydim, zdim);
	info._cmLset -= Vector3d(xmin,ymin,zmin);
	info._lset = lstrim;
//	info._lset = lset;
//	info._bboxShift = Vector3d(xmin,ymin,zmin);
	
	return info;
}

bool removeSingles(Levelset & ls) {
	bool removed = false;
	for (size_t k = 1; k < ls.getZdim()-1; k++) {
		for (size_t j = 1; j < ls.getYdim()-1; j++) {
			for (size_t i = 1; i < ls.getXdim()-1; i++) {
				if (ls(i,j,k) < 0 ){
					if ( ls(i+1,j,k)>0 && ls(i-1,j,k)>0 && ls(i,j+1,k)>0 && ls(i,j-1,k)>0 && ls(i,j,k-1)>0 && ls(i,j,k+1)>0){
						ls(i,j,k) = 1;
						removed = true;
					}
				}
			}
		}
	}
	return removed;
}




#endif /* LEVELSETGRAINFUNCTIONS_H_ */








