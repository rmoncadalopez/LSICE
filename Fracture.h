/*
 * Fracture.h
 *
 *  Created on: Jan 10, 2018
 *      Author: john
 *      Editted: rigo
 */

#ifndef FRACTURE_H_
#define FRACTURE_H_

#include "definitions.h"

//#include "Levelset2d.h"
//#include "Grain2d.h"
//#include "World2d.h"

struct PointInfo {
	Vector2d point;
	double shear;
	size_t contact;

	//bool operator< ( const PointsInfo& rhs ) const
	//{ return pointsort(points,rhs.points); }
};

struct CForce {
	double force;
	Vector2d cpoint;
	size_t cid;

	bool operator<(const CForce & rhs) const {
		return force > rhs.force;
	}

};

class FracProps2D {
public:
	// Quick Constructor
	FracProps2D() {
		_RoMratio 	= 0.05;
		_maxRoM		= 0.5;
		_minpts  	= 2;
		_crumbles 	= 0;
	}
	// Full Constructor
	FracProps2D(const string & dir, const double & RoMratio, const double & maxRoM, const size_t & minpts):
		_dir(dir), _RoMratio(RoMratio), _maxRoM(maxRoM), _minpts(minpts){
		_crumbles = 0;
	}


    //FOR ICE DENSITY
	// save grain to file
	void SaveGrain(const Grain2d & g, const vector<Vector2d> & gp0){
        double rhoice = 0.910e-6; // or 0.910
		
		cout << "Start Creation" << endl;
		FILE * newgrain1;
		stringstream morphfile1;
		morphfile1 << _dir << "Grain_" << g.getId() << ".dat";
		//cout<< _maxId+1 << endl;
		remove(morphfile1.str().c_str());
		newgrain1 = fopen (morphfile1.str().c_str(),"w");
		//ngrains
		fprintf(newgrain1, "%d\n", 1);
		//mass
		fprintf(newgrain1, "%f\n", g.getMass()/rhoice);
		//position/velocity unnecessary
		fprintf(newgrain1, "0\n0\n");
		//MoI
		fprintf(newgrain1, "%f\n", g.getMoI());
		//theta/omega unnecessary
		fprintf(newgrain1, "0\n0\n");
		//cmLset
		fprintf(newgrain1,"%f %f\n",g.getCmLset()(0),g.getCmLset()(1));
		//npoints
		vector<Vector2d> pts = gp0;
		fprintf(newgrain1,"%lu\n", pts.size());
		//points
		for (size_t i=0; i<pts.size() ; i++){
			fprintf(newgrain1,"%f %f ",pts[i](0),pts[i](1));
		}
		//radius
		fprintf(newgrain1,"\n%f\n", g.getRadius());
		cout << "Save Level Set!" << endl;
		//Lset dims
		fprintf(newgrain1,"%lu %lu\n", g.getLset().getXdim(),g.getLset().getYdim());
		//Lset
		cout << "Level Set Size: " << g.getLset().getLevelset().size() << endl;
		for (size_t i=0;i<g.getLset().getLevelset().size();i++){
			fprintf(newgrain1,"%f ", g.getLset().getLevelset()[i]);
		}
		fclose(newgrain1);
	}

	//save splitter to file
	void SaveSplitter(const vector<Vector2d> & splpts0, const vector<double> & splitset, const size_t & Xdim, const size_t & Ydim) {
		//output splitset
		FILE * splitvis;
		stringstream splitfile;
		splitfile << _dir << "split.dat";
		splitvis = fopen (splitfile.str().c_str(),"w");
		// npts
		fprintf(splitvis,"%lu\n",splpts0.size());
		//points
		for(size_t i=0;i<splpts0.size();i++) {
			fprintf(splitvis,"%f %f ",splpts0[i](0),splpts0[i](1));
		}
		//lset dim
		fprintf(splitvis,"\n%lu %lu\n",Xdim,Ydim);
		//lset
		for(size_t i=0;i<splitset.size();i++) {
			fprintf(splitvis,"%f ",splitset[i]);
		}
	}

	// find grain properties function

	//build splitter points
	vector<Vector2d> MakeSplitter(const Grain2d & g, const Vector2d & splitpt1, const Vector2d & splitpt2) {
		Vector2d g1cm = g.getCmLset();
		Levelset2d g1lset = g.getLset();
		vector<double> g1Lset = g1lset.getLevelset();
		size_t Xdim = g1lset.getXdim();
		size_t Ydim = g1lset.getYdim();

		//create splitter with 2 points found (2 steps: 1) make surface points 2) make level set)
		// 1) make points
		Vector2d pDiff = splitpt2-splitpt1;
		//create rotation matrix
		double cosR = (pDiff(0))/pDiff.norm();
		double sinR = (pDiff(1))/pDiff.norm();

		size_t splitptnum = 2*Xdim;

		//make points initially horizontal line
		vector<Vector2d> splpts(splitptnum);
		vector<Vector2d> splpts0(splitptnum);

		for (size_t i=0;i<splitptnum;i++){
			splpts[i] = Vector2d(-5+2.*double(i)*max(Xdim,Ydim)/splitptnum,splitpt1(1));
		}
		//rotate points
		for (size_t i=0;i<splitptnum;i++){
			splpts[i] = Vector2d((splpts[i](0)-splitpt1(0))*cosR - (splpts[i](1)-splitpt1(1))*sinR,
								 (splpts[i](0)-splitpt1(0))*sinR + (splpts[i](1)-splitpt1(1))*cosR) + splitpt1;
			splpts0[i] = splpts[i]-g1cm;
		}

		return splpts0;
	}

	// determine new surface points from splitter
	void KeepSplitterPoints(const vector<Vector2d> & splpts0,const Levelset2d & g1lset, const Vector2d & g1cm,
												vector<PointInfo> & g2i, vector<PointInfo> & g3i){

		//decides which splitter points stay
		for ( size_t i=0; i<splpts0.size();i++){
			double value; Vector2d gradient; //make to pass function, but actually unneeded
			if (g1lset.isPenetration(splpts0[i]+g1cm,value,gradient)){

				PointInfo addpoint;
				addpoint.point = splpts0[i];
				addpoint.shear = 0.;
				addpoint.contact = 0;

				g2i.push_back(addpoint);
				g3i.push_back(addpoint);
			}
		}
	}
	
	//New function for arbitrary line
	// determine new surface points from splitter
	void KeepSplitterPoints_arb(const vector<Vector2d> & splpts0,const Levelset2d & g1lset, const Vector2d & g1cm,
												vector<PointInfo> & g2i, vector<PointInfo> & g3i){

		//decides which splitter points stay
		for ( size_t i=0; i<splpts0.size();i++){
			double value; Vector2d gradient; //make to pass function, but actually unneeded
			if (g1lset.isPenetration(splpts0[i]+g1cm,value,gradient)){

				PointInfo addpoint;
				addpoint.point = splpts0[i];
				addpoint.shear = 0.;
				addpoint.contact = 0;

				g2i.push_back(addpoint); //Add line splitter points inside level set
				g3i.push_back(addpoint);
			}
		}
	}
	
	
	// determine new surface points from grain
	void KeepGrainPoints(const Vector2d & splitpt1,const Vector2d & splitpt2, const vector<PointInfo> & g1i, const Vector2d & g1cm,
												vector<PointInfo> & g2i, vector<PointInfo> & g3i){
		for (size_t i=0; i<g1i.size(); i++){
			//easier variables for split line points
			double x1 = splitpt1(0);
			double y1 = splitpt1(1);
			double x2 = splitpt2(0);
			double y2 = splitpt2(1);
			//easier variables for g1 surface points
			double x0 = g1i[i].point(0)+g1cm(0);
			double y0 = g1i[i].point(1)+g1cm(1);
			if (((y2-y1)*x0-(x2-x1)*y0+x2*y1-y2*x1)<0){
				//we save gxp in g1 origin reference frame
				g2i.push_back(g1i[i]);
			}
			else{
				g3i.push_back(g1i[i]);
			}
		}
	}
	
	
	//New function for arbitrary shaped line
	// determine new surface points from grain
	void KeepGrainPoints_arb(vector<Vector2d> & pointL, const vector<PointInfo> & g1i, const Vector2d & g1cm,
												vector<PointInfo> & g2i, vector<PointInfo> & g3i, size_t & xdom){
		//Loop point list to decide
		for (size_t i=0; i<g1i.size(); i++){
		    //easier variables for g1 surface points
			double x0 = g1i[i].point(0)+g1cm(0);
			double y0 = g1i[i].point(1)+g1cm(1);
		    
		    //Loop point pairs
		    for (size_t k=0; k<pointL.size()-1; k++){
    			//easier variables for split line points
    			Vector2d splitpt1 = pointL[k];
    			Vector2d splitpt2 = pointL[k+1];
    			double x1 = splitpt1(0);
    			double y1 = splitpt1(1);
    			double x2 = splitpt2(0);
    			double y2 = splitpt2(1);
                //XDOM
                if (xdom == 1)
                {
                    if (k == 0){
                        if (x0 < x2){
                			if (((y2-y1)*x0-(x2-x1)*y0+x2*y1-y2*x1)<0){
                				//we save gxp in g1 origin reference frame
                				g2i.push_back(g1i[i]);
                			}
                			else{
                				g3i.push_back(g1i[i]);
                			}
                        }
                    }
                    else if (k > 0 && k < pointL.size()-2){
                        if (x0 > x1 && x0 < x2){
                			if (((y2-y1)*x0-(x2-x1)*y0+x2*y1-y2*x1)<0){
                				//we save gxp in g1 origin reference frame
                				g2i.push_back(g1i[i]);
                			}
                			else{
                				g3i.push_back(g1i[i]);
                			}
                        }
                    }
                   else{
                        if (x0 > x1){
                			if (((y2-y1)*x0-(x2-x1)*y0+x2*y1-y2*x1)<0){
                				//we save gxp in g1 origin reference frame
                				g2i.push_back(g1i[i]);
                			}
                			else{
                				g3i.push_back(g1i[i]);
                			}
                        }
                    }
                }
                //YDOM
                else
                {
                    if (k == 0){
                        if (y0 < y2){
                			if (((x2-x1)*y0-(y2-y1)*x0+y2*x1-x2*y1)<0){
                				//we save gxp in g1 origin reference frame
                				g2i.push_back(g1i[i]);
                			}
                			else{
                				g3i.push_back(g1i[i]);
                			}
                        }
                    }
                    else if (k > 0 && k < pointL.size()-2){
                        if (y0 > y1 && y0 < y2){
                			if (((x2-x1)*y0-(y2-y1)*x0+y2*x1-x2*y1)<0){
                				//we save gxp in g1 origin reference frame
                				g2i.push_back(g1i[i]);
                			}
                			else{
                				g3i.push_back(g1i[i]);
                			}
                        }
                    }
                   else{
                        if (y0 > y1){
                			if (((x2-x1)*y0-(y2-y1)*x0+y2*x1-x2*y1)<0){
                				//we save gxp in g1 origin reference frame
                				g2i.push_back(g1i[i]);
                			}
                			else{
                				g3i.push_back(g1i[i]);
                			}
                        }
                    }
                }
		    } //end splitter loop
		} //end pointList loop

	}
	
	
	
	
	
	

	void findGrainProps0(const Levelset2d & glset, double & mass, Vector2d & gcm, double & I) {
		//mass and cm
		mass = 0;
		gcm = Vector2d(0., 0.);

		vector<double> gLset = glset.getLevelset();
		size_t Xdim = glset.getXdim();
		size_t Ydim = glset.getYdim();

		vector<double> heavi = sHeaviOpp(gLset);
		size_t iter = 0;

		for (size_t j = 0; j < Ydim; j++) {
			for (size_t i = 0; i < Xdim; i++) {
				mass += heavi[iter];
				gcm(0) += heavi[iter] * (double) i;
				gcm(1) += heavi[iter] * (double) j;
				iter++;
			}
		}
		gcm /= mass;

		//MoI
		double rx, ry;
		I = 0;
		iter = 0;
		for (size_t j = 0; j < Ydim; j++) {
			for (size_t i = 0; i < Xdim; i++) {
				rx = (double) i - gcm(0);
				ry = (double) j - gcm(1);
				I += heavi[iter] * (rx * rx + ry * ry);
				iter++;
			}
		}
	}
	
    void findGrainProps0D(const Levelset2d & glset, double & mass, Vector2d & gcm, double & I, const double & density) {
		//mass and cm
		mass = 0;
		gcm = Vector2d(0., 0.);

		vector<double> gLset = glset.getLevelset();
		size_t Xdim = glset.getXdim();
		size_t Ydim = glset.getYdim();

		vector<double> heavi = sHeaviOpp(gLset);
		size_t iter = 0;

		for (size_t j = 0; j < Ydim; j++) {
			for (size_t i = 0; i < Xdim; i++) {
				mass += heavi[iter];
				gcm(0) += heavi[iter] * (double) i;
				gcm(1) += heavi[iter] * (double) j;
				iter++;
			}
		}
		gcm /= mass;
		
		//Same as in melt function
		mass *= density; 
		//mass *= thick*0.001; //CHECK EFFECT!!!! Later in FracRoutine_Nova
		

		//MoI
		double rx, ry;
		I = 0;
		iter = 0;
		for (size_t j = 0; j < Ydim; j++) {
			for (size_t i = 0; i < Xdim; i++) {
				rx = (double) i - gcm(0);
				ry = (double) j - gcm(1);
				I += heavi[iter] * (rx * rx + ry * ry);
				iter++;
			}
		}
		
		//Same as melt function
		I *= density; 
		//I *= thick;  //*0.001; Later in FracRoutine_Nova
	}
	

	void findGrainProps(const Levelset2d & glset,
									double & mass, Vector2d & gcm, double & I, const double & density) {
		//mass and cm
		mass = 0;
		gcm = Vector2d(0., 0.);

		vector<double> gLset = glset.getLevelset();
		size_t Xdim = glset.getXdim();
		size_t Ydim = glset.getYdim();

		vector<double> heavi = gLset;

        double eps = 1.5; //Check if it is okay  //MODIFY FOR MELT TEMP 
		for (size_t j = 0; j < glset.getYdim(); j++) {
			for (size_t i = 0; i < glset.getXdim(); i++) {   //Substitute of the Heaviside Function to Identify Inside and Outside of Grain
				if (-heavi[(j*glset.getXdim())+i] > eps) {  
                   heavi[(j*glset.getXdim())+i] = 1; 
                }
				else if (-heavi[(j*glset.getXdim())+i] < -eps) {
                   heavi[(j*glset.getXdim())+i] = 0;
				}	
				else{
                   heavi[(j*glset.getXdim())+i] = 0.5*(1. + heavi[(j*glset.getXdim())+i]/eps + sin(M_PI*heavi[(j*glset.getXdim())+i]/eps)/M_PI); //0.5*(1. + vec(i)/eps + sin(M_PI*vec(i)/eps)/M_PI);
				}	
			}
        }


		size_t iter = 0;

		for (size_t j = 0; j < Ydim; j++) {
			for (size_t i = 0; i < Xdim; i++) {
				mass += heavi[iter];
				gcm(0) += heavi[iter] * (double) i;
				gcm(1) += heavi[iter] * (double) j;
				iter++;
			}
		}
		gcm /= mass;

		mass = mass*density*0.89/(1);  //////

		//MoI
		double rx, ry;
		I = 0;
		iter = 0;
		for (size_t j = 0; j < Ydim; j++) {
			for (size_t i = 0; i < Xdim; i++) {
				rx = (double) i - gcm(0);
				ry = (double) j - gcm(1);
				I += heavi[iter] * (rx * rx + ry * ry);
				iter++;
			}
		}
	}


	// Filter out problematic grains. True = Passes Filter, False = Fails Filer

	//ONLY FOR ICE!!!!!!
	bool ApplyFilter (const Grain2d & g) {
//		size_t npts = g.getPointList().size();
//		double R    = g.getRadius();
//		double bigV = M_PI*R*R;

 
		double rhoice = 0.910e-6; // or would it be 0.910

		double gV   = g.getMass()/rhoice;
//		if (npts > _minpts && bigV/gV < _maxRoM){
		if (gV > 15.){
			return true;
		}
		else {
			_crumbleid.push_back(g.getId());
			_crumbles += g.getMass();
			return false;
		}
	}

	bool ApplyFilter_2 (const Grain2d & g) {
//		size_t npts = g.getPointList().size();
//		double R    = g.getRadius();
//		double bigV = M_PI*R*R;

 
		double rhoice = 0.910e-6; // or would it be 0.910

		double gV   = g.getMass()/rhoice;
//		if (npts > _minpts && bigV/gV < _maxRoM){
		if (gV > 1.){
			return true;
		}
		else {
			_crumbleid.push_back(g.getId());
			_crumbles += g.getMass();
			return false;
		}
	}


    //Geometric Level Set Splitting Function
    void VertSplitLSETSG(const Grain2d & g, const int & pbk, vector<double> & g2LV, vector<double> & g3LV){
		
		for (size_t i = 0; i<g.getLset().getXdim(); i++) 
   		{            
            for (size_t j = 0; j<g.getLset().getYdim(); j++)
            {
            	if (i>pbk-1){  //Perfectly vertical fracture (for now)  //*****FOR VERTICAL EDGE******** Left grain
                        g2LV[(j*g.getLset().getXdim())+i] = 1;
                }        
            	if (i<pbk+1){  //Perfectly vertical fracture (for now)  //*****FOR VERTICAL EDGE******** Right grain
                        g3LV[(j*g.getLset().getXdim())+i] = 1;
                }                   
            }
    	} 
    }

    //Thermal Level Set Splitting Function
    void VertSplitLSETS(const Grain2d & g, const int & pbk, vector<double> & g2LV, vector<double> & g3LV){
	double Tw = 20*1.5; //Must be equal to that of compute world!!!!!!!!

		for (size_t i = 0; i<g.getgrainTemp2D().getXdim(); i++) 
   		{            
            for (size_t j = 0; j<g.getgrainTemp2D().getYdim(); j++)
            {
            	if (i>pbk-1){  //Perfectly vertical fracture (for now)  //*****FOR VERTICAL EDGE******** Left grain
                        g2LV[(j*g.getgrainTemp2D().getXdim())+i] = Tw;
                }        
            	if (i<pbk+1){  //Perfectly vertical fracture (for now)  //*****FOR VERTICAL EDGE******** Right grain
                        g3LV[(j*g.getgrainTemp2D().getXdim())+i] = Tw;
                }                   
            }
    	} 
    }


	//Point Generating Function using a Level Set
    void PointsGen(const Levelset2d & glset, vector<Vector2d> & gpts, const int & npts, const Vector2d & newcm, bool & fail_ind){
    	bool fail_indi = false; //Assume it's okay
    	vector<double> gvec = glset.getLevelset();		
        double damping= 0.3;
		vector<Vector2d> pointsnew(npts); 

		//cout<<"PrintGVec"<<endl;	
        //for (size_t j = 0; j < glset.getYdim()-1; j++) 
        //{
		//	for (size_t i = 0; i < glset.getXdim()-1; i++) 
		//	{
        //         //cout<<"Value GVec at i: "<< i <<" j: " << j <<" is: " << gvec[(j*glset.getXdim())+i] << endl;
		//	}
		//}

        //cout<<"Zero Cross"<<endl;
    	//  place points interpolating around 0 contour or meltTemp
		vector<Vector2d> zeroCrossings;
		for (size_t j = 0; j < glset.getYdim()-1; j++) {
			for (size_t i = 0; i < glset.getXdim()-1; i++) {
				double val = gvec[(j*glset.getXdim())+i];
				if (val*gvec[((j+1)*glset.getXdim())+i] < (0)   ) { 
					zeroCrossings.push_back(Vector2d( double(i), double(j)+fabs(val)/(fabs(val) + fabs(gvec[((j+1)*glset.getXdim())+i]))  ));
				}
				if (val*gvec[(j*glset.getXdim())+(i+1)] < (0)   ) {	
					zeroCrossings.push_back(Vector2d(  double(i)+fabs(val)/(fabs(val) + fabs(gvec[(j*glset.getXdim())+(i+1)])), double(j)  ));
				}
			}
		}	

	    double step = double(zeroCrossings.size())/double(npts);

         //cout<<"Temporal pointsnew"<<endl;
         //cout<<"Numberpoints: "<< npts <<endl;
         //cout<<"ZeroC size: "<< zeroCrossings.size() <<endl;
		for (size_t i = 0; i < npts; i++) {
            pointsnew[i] = zeroCrossings[size_t(double(i)*step)];  
		}

	    // set up temp vars for iteration
		vector<Vector2d> glist(npts); vector<double> sdlist(npts);
		Vector2d grad; double sd;
		Vector2d fij; double dij; Vector2d nij; Vector2d nji; double magp;
		vector<Vector2d> forces(npts);
		vector<Vector2d> v(npts);
		for (size_t i = 0; i < npts; i++) {
			v[i] << 0., 0.;
		}
		
		//cout<<"Dirac"<<endl;
	    //double sa = sDirac(lset.getGrid()).sum();      for loop perhaps double sa = sDirac(Tvec[(j*_lset.getXdim())+i])  .sum();
	    vector<double> gvec0=glset.getLevelset(); //for Geometric Level Set -->
	    double eps = 0.015; //Check if it is okay   
	    double sa = 0;
        for (size_t j = 0; j < glset.getYdim(); j++) {
				for (size_t i = 0; i < glset.getXdim(); i++) {   //Substitute of the sDirac Function to Identify Inside and Outside of Grain
					if (fabs(gvec0[(j*glset.getXdim())+i]) > eps) {  
                       gvec0[(j*glset.getXdim())+i] = 0; 
                    }
					else{

                       gvec0[(j*glset.getXdim())+i] = (1.+cos(M_PI*gvec0[(j*glset.getXdim())+i]/eps))/(2.*eps); //(1.+cos(M_PI*vec(i)/eps))/(2.*eps);
					   
					}
					sa = sa + gvec0[(j*glset.getXdim())+i];
				}
		}

		double ddt = sqrt(sa/npts)*0.1;
        size_t niter=300;  //ADJUST

        //cout<<"Iterations"<<endl;
	// begin iteration
	 	for (size_t t = 0; t < niter; t++ ) { // niter
	        // store all of the sds and gradients first, also zero out the force
			for (size_t i = 0; i < npts; i++) {
				if (glset.findValueGradient(pointsnew[i], sd, grad)) {   //change _lset for Geometric Level Set, use _grainTemp2D for Heat Level Set -->
					glist[i] = grad;
					sdlist[i] = sd;
				}
				else {
					glist[i] << 0,0;
					sdlist[i] = -1;
					//pointsnew[i] << 5+(double)rand()/RAND_MAX,5+(double)rand()/RAND_MAX; //or << //PSOE
				}
				forces[i] << 0., 0.; // or <<
			}
	// 		// iterate through points
			for (size_t i = 0; i < npts; i++) {
 			//drag points closer to surface using gradient/signed distance
	        //pointsnew[i] -= .8*sdlist[i]*glist[i];                                          //PSOE
				// push points that are close away from each other
				for (size_t j = i+1; j < npts; j++) {
					if (i == j) { break; }
					fij = pointsnew[i] - pointsnew[j];
					dij = fij.norm();
					if (dij < 1.1*sqrt(sa/npts)) {
						nij =  fij - glist[i]*fij.dot(glist[i]);
						nji = -fij + glist[j]*fij.dot(glist[j]);
						magp = exp(-dij/sqrt(sa/npts));

						forces[i] += nij/(nij.norm()+DBL_MIN)*magp;
						forces[j] += nji/(nji.norm()+DBL_MIN)*magp;
					}
				}
			}
			for (size_t i = 0; i < npts; i++) {
				v[i] = 1/(1+damping*ddt/2)*((1-damping*ddt/2)*v[i] + ddt*forces[i]);
		        //pointsnew[i] += ddt*v[i];                                                 //PSOE
			}
		}
		
		// additional iterations to move points to zero contour
		for (int tt = 0; tt < 200; tt++) {   //200
			for (size_t i = 0; i < npts; i++) {
				if (glset.findValueGradient(pointsnew[i], sd, grad)) {                   //change _lset for Geometric Level Set, use _grainTemp2D for Heat Level Set -->   //Does it work for Binary Lset mmmmm unlikely
					  //pointsnew[i] -= .8*sd*grad;                                        //PSOE
				}
			}
		}

	        //Check if pointsnew were done ok
		bool nanp_before = false;
		cout << "Check Nan Points in Point Gen Frac!!!!" << endl;
		for (size_t ip = 0; ip < pointsnew.size(); ip++)
		{
			if (isnan(pointsnew[ip](0)) || isnan(pointsnew[ip](1)))
			{
			  nanp_before = true;
			  fail_indi = true;
			  break;
			}
		}  

		if (nanp_before)
		{
			cout << "CONFIRMED: Nan Points in Point Gen Frac!!!!" << endl;
			return;  //Skip to next healthy grain
		}	


        //REORDER POINTS CLOCKWISELY

        //FIND CLOCKWISE ANGLE FROM CENTROID TO EACH POINT, THEN ORDER THE POINTS FROM LOWEST TO HIGHEST ANGLE (ALTERNATIVE 2)
        
        //cout<<"Sorting CW for npoints ph1: "<< npts <<endl;

        vector<Vector2d> pneworder=pointsnew;

        // Angle , Point Number pairs
        vector<Vector2d> pangle=pointsnew;  //To pick from and delete from pool
        vector<Vector2d> pangle2=pointsnew;  //To reorder on and overwrite

        vector<int> indicesorder(npts);
        size_t minindex=0;
        const double pii=atan(1)*4;
        Vector2d savior = Vector2d (0.,0.);

        for (size_t i = 0; i < npts; i++) {       //Loop for Setting Up Angles
        	
             pangle[i](0) = (atan (  (  pneworder[i](1)-newcm(1) ) /  (  pneworder[i](0)-newcm(0) ) ) )*(180/pii); //Get CW angle in degrees
             pangle[i](1) = i; //Point index storage for convenience
             
             //Get CW angle in degrees and in the right quadrant measured as CW from theta=0

             //MODIFICATION 1: ADD 90 TO START AT LOWER BORDER WHICH IS USUALLY LESS ALTERED THAN RIGHT OR ZERO; OR EVEN 180 WHICH IS RIGHT BORDER

             if (   (pneworder[i](0)-newcm(0)) >= 0 && (pneworder[i](1)-newcm(1)) <= 0  ) {
             	   pangle[i](0)=-pangle[i](0) + 90; 
             } else if (   (pneworder[i](0)-newcm(0)) < 0 && (pneworder[i](1)-newcm(1)) < 0  ) {
             	   pangle[i](0)=180-pangle[i](0) + 90;

             }  else if (   (pneworder[i](0)-newcm(0)) < 0 && (pneworder[i](1)-newcm(1)) > 0  ) {
             	   pangle[i](0)=180-pangle[i](0) + 90;

             }  else { //if (   (pneworder[i](0)-_cmLset(0)) > 0 && (pneworder[i](1)-_cmLset(1)) > 0  ) {
                   pangle[i](0)=360-pangle[i](0) + 90;
             } //else {
                 //cout << "grain " << _id << " has a problem at its newpointagle" << endl; 
             //}
             pangle2[i](0)=pangle[i](0);
             pangle2[i](1)=pangle[i](1);
        }
        
        //cout<<"Sorting CW for npoints ph2: "<< npts <<endl;
        //cout<<"Sorting CW Part 2"<<endl;
        //Sorting Loop CW
        double closeness=0.0;
        for (size_t i = 0; i < npts; i++) {          //sort(v.begin(), v.end()); 
        	minindex=0;
        	if (i == npts-1){      //Check Closing Process for Points
                pangle2[i] = pangle[0];   //Should be only pangle2[i] = pangle[1]; //This is only to get the last remaining point
        	}
        	else if (i<3){   // We can map n points clockwise to set up right clockwise direction and then go with nearest distance instead of angle
                for (size_t j = 0; j < pangle.size() ; j++) {	 //Here we explore all points to find the three of minimum angle
                	if ( pangle[minindex](0)>pangle[j](0) ) {    //Improve criterion 
                   		minindex=j;
                	}
                }	
               savior=pangle[minindex];
               pangle.erase(pangle.begin() + minindex);
               pangle2[i]= savior;
        	}
            else{   //This is to sort all the other points using nearest distance to next instead of angle ordering, given regular spacing this should be fine, but still!!! 
            	
                //Should a competition for nearest distance be established for improving quality?????? Given preference to closer on Y???
 
            	closeness=(pneworder[pangle2[i-1](1)]-pneworder[pangle[minindex](1)]).norm();
            	for (size_t j = 1; j < pangle.size() ; j++) //pangle.size()
            	{     
            		                                              

            		if (  ((pneworder[pangle2[i-1](1)]-pneworder[pangle[j](1)]).norm())<closeness  )
            		{
                           
                           closeness=(pneworder[pangle2[i-1](1)]-pneworder[pangle[j](1)]).norm();
                           minindex=j;
                    }
                    //If the distance using NN is the same, use angle to solve the conflict
                    else if (  ((pneworder[pangle2[i-1](1)]-pneworder[pangle[j](1)]).norm()) == closeness  )
                    {
                    	if ( pangle[minindex](0)>pangle[j](0) )  //Improve criterion 
                    	{    
                   			minindex=j;
                		}	 

                    }	
                }
               
               savior=pangle[minindex];
               pangle.erase(pangle.begin() + minindex);
               pangle2[i]= savior;
               //cout << "pangle " << pangle2[i](0) << " order" << endl; 
            } 
            //cout << "pangle " << pangle2[i](0) << " order" << endl; 
        }
        
        //pangle2[i] = pangle[0];
        //cout<<"Export"<<endl;
        //cout<<"Sorting CW for npoints ph3: "<< npts <<endl;
        for (size_t i = 0; i < npts; i++) {
        	
        	minindex=pangle2[i](1);
        	
        	//minindex=indicesorder[i];
        	//cout << "pangle " << pangle2[i](0) << " grain no." << _id << endl;

        	gpts[i]=pointsnew[minindex];  //Using Order pangle2 array, pointsnew generated from ZeroCrossing are ideally fed CW to export gpts

        }	

    }



	// get methods
	const vector<size_t> getCrumbleId() const {
		return _crumbleid;
	}
	const double & getCrumbles() const {
		return _crumbles;
	}
	const double & getRoMratio() const {
		return _RoMratio;
	}
	const double & getMaxRoM() const {
		return _maxRoM;
	}


private:
	string 			_dir;      	//save directory for new grains
	double  			_RoMratio; 	//cutoff radius/mass value for grains
	double			_maxRoM;
	double  			_minpts;  	//minimum allowable mass before sent to crumbles
	double			_crumbles; 	//mass of crumbles
	vector<size_t> 	_crumbleid;
};






#endif /* FRACTURE_H_ */
