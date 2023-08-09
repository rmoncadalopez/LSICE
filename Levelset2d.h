/*
 * Levelset2d.h
 *
 *  Created on: May 31, 2014
 * Author: Reid Kawamoto - Konstantinos Karapiperis - Liuchi Li
 */

#ifndef LEVELSET2D_H_
#define LEVELSET2D_H_

#include "definitions.h"

class Levelset2d {

public:

	// default constructor
	Levelset2d() {
		_xdim = 0; _ydim = 0;
	}
	// actual constructor
	Levelset2d(const vector<double> & levelset, const size_t & xdim, const size_t & ydim):
		_levelset(levelset), _xdim(xdim), _ydim(ydim) {
		if (_xdim*_ydim != _levelset.size()) {
			cout << "ERROR: levelset size not consistent with dimensions" << endl;
		}
	}
	
	// checks if there is penetration, if there is, finds the penetration amount and the normalized gradient
	// and stores them in input fields value and gradient.  If not, returns false.
	bool isPenetration(const Vector2d & point, double & value, Vector2d & gradient) const {
		double x = point(0);
		double y = point(1);

		// check if the point exists in the level set, if not, return false
		if (x+1 > double(_xdim) || y+1 > double(_ydim) || x < 0 || y < 0){
			return false;
		}
		size_t xr = (size_t)round(x);
		size_t yr = (size_t)round(y);

		// check if the point is close to the surface, if not, return false
		if (getGridValue(xr,yr) > 1) {
			return false;
		}

		// if inside the level set, do bilinear interpolation
		size_t xf 	= floor(x);
		size_t yf 	= floor(y);
		size_t xc 	= ceil(x);
		size_t yc 	= ceil(y);
		double dx 	= x - xf;
		double dy 	= y - yf;
		double b1 	= getGridValue(xf, yf);
		double b2 	= getGridValue(xc, yf) - b1;
		double b3 	= getGridValue(xf, yc) - b1;
		double b4 	= -b2 - getGridValue(xf, yc) + getGridValue(xc, yc);

		value = b1 + b2*dx + b3*dy + b4*dx*dy;
		if (value > 0) {
			return false;
		}
		gradient << b2 + b4*dy, b3 + b4*dx;
		gradient /= gradient.norm();
		return true;
	}

	// Finds if point is outside the grain but within cohesive distance to the surface
	// Returns actual distance and gradient 
	bool isCohesion(const Vector2d & point, const double & cohesiveDistance, 
					double & value, Vector2d & gradient) const {

		double x = point(0);
		double y = point(1);

		// check if the point exists in the level set, if not, return false 
		// since the level sets are already constructed to have a border
		// larger than the maximum cohesive distance 
		if (x+1 > double(_xdim) || y+1 > double(_ydim) || x < 0 || y < 0){
			return false;
		}

		// if inside the level set, do bilinear interpolation
		size_t xf 	= floor(x);
		size_t yf 	= floor(y);
		size_t xc 	= ceil(x);
		size_t yc 	= ceil(y);
		double dx 	= x - xf;
		double dy 	= y - yf;
		double b1 	= getGridValue(xf, yf);
		double b2 	= getGridValue(xc, yf) - b1;
		double b3 	= getGridValue(xf, yc) - b1;
		double b4 	= -b2 - getGridValue(xf, yc) + getGridValue(xc, yc);

		value = b1 + b2*dx + b3*dy + b4*dx*dy;
		if (value < 0 || value > cohesiveDistance) { 
			return false;
		}
		gradient << b2 + b4*dy, b3 + b4*dx;
		gradient /= gradient.norm();
		return true;
	}
	
	//Avoid penetration bond breakage
	bool getNormal(const Vector2d & point, const double & cohesiveDistance, 
					double & value, Vector2d & gradient) const {

		double x = point(0);
		double y = point(1);

		// check if the point exists in the level set, if not, return false 
		// since the level sets are already constructed to have a border
		// larger than the maximum cohesive distance 
		if (x+1 > double(_xdim) || y+1 > double(_ydim) || x < 0 || y < 0){
			return false;
		}

		// if inside the level set, do bilinear interpolation
		size_t xf 	= floor(x);
		size_t yf 	= floor(y);
		size_t xc 	= ceil(x);
		size_t yc 	= ceil(y);
		double dx 	= x - xf;
		double dy 	= y - yf;
		double b1 	= getGridValue(xf, yf);
		double b2 	= getGridValue(xc, yf) - b1;
		double b3 	= getGridValue(xf, yc) - b1;
		double b4 	= -b2 - getGridValue(xf, yc) + getGridValue(xc, yc);

		value = b1 + b2*dx + b3*dy + b4*dx*dy;
		// //if (value < 0 || value > cohesiveDistance) { 
		//if (value < 0 ) {  //Turn ON to break bonds if points enter other grain Level Set
		//	return false;
		//}
		gradient << b2 + b4*dy, b3 + b4*dx;
		gradient /= gradient.norm();
		return true;
	}

	// Computes actual distance if point is in level set, else INT_MAX 
	double findCohesiveDistance(const Vector2d & point, Vector2d & gradient, bool & outOfLset) const {

		outOfLset = false;
		double x = point(0);
		double y = point(1);

		// check if the point exists in the level set, if not, return false 
		// since the level sets are already constructed to have a border
		// larger than the maximum cohesive distance 
		if (x+1 > double(_xdim) || y+1 > double(_ydim) || x < 0 || y < 0){
			outOfLset = true;
		}

		// if inside the level set, do bilinear interpolation
		size_t xf 	= floor(x);
		size_t yf 	= floor(y);
		size_t xc 	= ceil(x);
		size_t yc 	= ceil(y);
		double dx 	= x - xf;
		double dy 	= y - yf;
		double b1 	= getGridValue(xf, yf);
		double b2 	= getGridValue(xc, yf) - b1;
		double b3 	= getGridValue(xf, yc) - b1;
		double b4 	= -b2 - getGridValue(xf, yc) + getGridValue(xc, yc);

		double value = b1 + b2*dx + b3*dy + b4*dx*dy;
		gradient << b2 + b4*dy, b3 + b4*dx;
		gradient /= gradient.norm();
		return value;
	}

	// Methods to deform the level sets
	void shearVertical(const double & gamma) {
		for (size_t j = 0; j < _ydim; j++) {
			for (size_t i = 0; i < _xdim; i++) {
//				_levelset[j*_xdim + i] += ((double)j - (double)_ydim/2.)*gamma;
				_levelset[j*_xdim + i] += ((double)j - 45.)*gamma;
			}
		}
	}
	
	void shearHorizontal(const double & gamma) {
		for (size_t j = 0; j < _ydim; j++) {
			for (size_t i = 0; i < _xdim; i++) {
//				_levelset[j*_xdim + i] += ((double)i - (double)_xdim/2.)*gamma;
				_levelset[j*_xdim + i] += ((double)i - 45.)*gamma;
			}
		}
	}

   bool findValueGradient(const Vector2d & point, double & value, Vector2d & gradient) const {
		double x = point(0);
		double y = point(1);
		if (x+1 > (double)_xdim || y+1 > (double)_ydim || x < 0 || y < 0 ) {
			return false;
        }
        // find the signed distance and the gradient by performing trilinear interpolation
		size_t x0 	= (size_t)floor(x);
		size_t y0 	= (size_t)floor(y);
		size_t x1 	= (size_t)ceil(x);
		size_t y1 	= (size_t)ceil(y);
		double p000 = getGridValue2(x0, y0);
		double xm 	= getGridValue2(x1, y0) - p000;
		double ym 	= getGridValue2(x0, y1) - p000;
		double xym	= -xm - getGridValue2(x0,y1) + getGridValue2(x1,y1);
		//double xyzm = -xym + getGridValue2(x0,y0) - getGridValue2(x1,y0) - getGridValue2(x0,y1) + getGridValue2(x1,y1);
		double dx 	= x - double(x0);
		double dy 	= y - double(y0);
		value = p000 + xm*dx + ym*dy + xym*dx*dy;
		gradient << xm + xym*dy, ym + xym*dx; 
		gradient /= gradient.norm(); // do we need to normalize the gradient?
		return true;
	}

	vector<double> makeFracSurfaceLine(const Vector2d & pt1, const Vector2d & pt2){
		vector<double> splitset(_levelset.size());
		for(size_t y=0;y<_ydim;y++){
			for(size_t x=0;x<_xdim;x++){//same size as grain to be split
				splitset[y*_xdim+x] = ((pt2(1)-pt1(1)) * double(x) - (pt2(0)-pt1(0)) * double(y)  + pt2(0) * pt1(1) - pt2(1) * pt1(0)) / (pt2-pt1).norm();
			}
		}
		return splitset;
	}



    //New arbitrary shape line function
	vector<double> makeFracSurfaceLine_arb(vector<Vector2d> & pointL, size_t & xdom){
		
		//Define all level set vector objects
		vector<double> splitset(_levelset.size()); //Final vector to output
		size_t nptsL = pointL.size(); //We need to loop every single point pair
		
		vector<double> splitset_temp(_levelset.size()); //Final vector to output
		vector<vector<double>> LS_Split; //Save value for each point pair, then filter based on linear location
		
		//For vector of pairwise points
		for(size_t z=0;z<nptsL-1;z++){
		    //Choose pair wise points up to nptsL - 1 
		    Vector2d pt1 = pointL[z]; 
		    Vector2d pt2 = pointL[z+1]; 
		    
    		for(size_t y=0;y<_ydim;y++){
    			for(size_t x=0;x<_xdim;x++){//same size as grain to be split
    			    if (xdom == 1){
    				    splitset_temp[y*_xdim+x] = ((pt2(1)-pt1(1)) * double(x) - (pt2(0)-pt1(0)) * double(y)  + pt2(0) * pt1(1) - pt2(1) * pt1(0)) / (pt2-pt1).norm();
    			    }
    			    else{
    			        splitset_temp[y*_xdim+x] = ((pt2(0)-pt1(0)) * double(y) - (pt2(1)-pt1(1)) * double(x)  + pt2(1) * pt1(0) - pt2(0) * pt1(1)) / (pt2-pt1).norm();
    			    }
    			}
    		}
    		LS_Split.push_back(splitset_temp);
		}

        //Filter according to point location, we need to inspect each cell individually and see what point pair they fullfil.
        for(size_t y=0;y<_ydim;y++){
    		for(size_t x=0;x<_xdim;x++){//same size as grain to be split
                for(size_t z=0;z<nptsL-1;z++){
                    Vector2d pt1 = pointL[z]; 
		            Vector2d pt2 = pointL[z+1]; 
                    vector<double> LS_Split_temp = LS_Split[z];
                    if (xdom == 1)
                    {
                        if (z == 0){
                            if ( double(x) < pt2(0) ){
                                splitset[y*_xdim+x] = LS_Split_temp[y*_xdim+x];
                            }
                        }
                        else if ( double(z) > 0 && double(z) < (nptsL -2)  )
                        {
                            if ( double(x) > pt1(0) && double(x) < pt2(0) ){
                                splitset[y*_xdim+x] = LS_Split_temp[y*_xdim+x];
                            }
                        }
                        else{
                            if ( double(x) > pt1(0) ){
                               splitset[y*_xdim+x] = LS_Split_temp[y*_xdim+x];
                            }
                        }
                    }
                    else
                    {
                        if (z == 0){
                            if ( double(y) < pt2(1) ){
                                splitset[y*_xdim+x] = LS_Split_temp[y*_xdim+x];
                            }
                        }
                        else if ( double(z) > 0 && double(z) < (nptsL -2)  )
                        {
                            if ( double(y) > pt1(1) && double(y) < pt2(1) ){
                                splitset[y*_xdim+x] = LS_Split_temp[y*_xdim+x];
                            }
                        }
                        else{
                            if ( double(y) > pt1(1) ){
                               splitset[y*_xdim+x] = LS_Split_temp[y*_xdim+x];
                            }
                        }
                    }
                }
    		}
        }

		return splitset;
	}




    double getGridValue2(const size_t & x, const size_t & y) const {
		return _levelset[y*_xdim + x];
	}
	
//  public get methods for debugging
	int getXdim() const {
		return _xdim;
	}
	int getYdim() const {
		return _ydim;
	}
	vector<double> getLevelset() const {
		return _levelset;
	}

    void changeLevelset(const vector<double> & NewLSet) {
		_levelset = NewLSet;
	}



	
private:

	double getGridValue(size_t & x, size_t & y) const {
		return _levelset[y*_xdim + x];
	}
	
	vector<double> _levelset;
	size_t _xdim;
	size_t _ydim;
};

#endif /* LEVELSET2D_H_ */
