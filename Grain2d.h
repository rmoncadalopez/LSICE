/*
 * Grain2d.h
 * Bond contact version
 * Author: Reid Kawamoto - Konstantinos Karapiperis - Liuchi Li
 */

#ifndef GRAIN2D_H_
#define GRAIN2D_H_

#include "definitions.h"
#include "Levelset2d.h"
#include "WorldStates.h"

#include "SmoothReinit.h"

typedef vector<vector<double> > Matrixx;
typedef vector<double> Row;


class Grain2d {
public:
	// constructors
	Grain2d() {
		_radius = 0; _mass = 0; _momentInertia = 0;      
		_mass0 = 0; _temper = 0; _thick = 0; _temper0 = 0; _thick0 = 0; //_Utemper=0.; //_Utemper0={};                                                                //Thermal
		_omega = 0; _theta = 0; _id = 0; _density = 1;
		_morphologyID = 0; _kn = 0; _ks = 0; _mu = 0; _ncontacts = 0; _critStress = 0; _remove = false;_d =0;
	}


	Grain2d(const double & mass, const Vector2d & position, const Vector2d & velocity, 
		      const double & mass0, const double & temper, const double & thick,  const Levelset2d & Mthick, const double & temper0, const double & thick0, const vector<double> & Utemper, const vector<double> & Utemper0, 
		      const Levelset2d & grainTemp2D, const Levelset2d & grainTemp2D0, //Thermal
			  const double & momentInertia, const double & theta, const double & omega,
			  const Vector2d & cmLset, const vector<Vector2d> & pointList, const int & npoints,
			  const double & radius, const Levelset2d & lset, const Levelset2d & lset0, const int & id,
			  const int & morphologyID, const double & kn, const double & ks,
			  const double & mu, const bool & fracFlag, const int & fracLoc,
			  const Levelset2d & grainStress2D, const vector<Vector3d> & grainDamage, const MatrixXd & grainThickness, const vector<Vector2d> & refpointList, const double & time_fail, const double & origin_time_fail):
			  _mass(mass), _position(position), _velocity(velocity),   
              _mass0(mass0), _temper(temper), _thick(thick), _Mthick(Mthick) ,_temper0(temper0), _thick0(thick0), _Utemper(Utemper), _Utemper0(Utemper0), _grainTemp2D(grainTemp2D), _grainTemp2D0(grainTemp2D0),    //Thermal
			  _momentInertia(momentInertia), _npoints(npoints),
			  _theta(theta), _omega(omega), _cmLset(cmLset), _radius(radius), _lset(lset), _lset0(lset0), _id(id),
			  _morphologyID(morphologyID), _kn(kn), _ks(ks), _mu(mu), _fracFlag(fracFlag), _fracLoc(fracLoc),
			  _grainStress2D(grainStress2D), _grainDamage(grainDamage), _grainThickness(grainThickness), _time_fail(time_fail), _origin_time_fail(origin_time_fail)  {

		_pointList = pointList;
		_refpointList = refpointList;
		_ncontacts = 0;
		_density = 1;
		_vr = 4. * _radius;
		_nodeShears.resize(_pointList.size());
		_nodeContact.resize(_pointList.size());
		// _nodeNormals.resize(_pointList.size());

		// Increase bbox radius by cohesive distance to capture cohesive interaction check

		// Compute rotation matrix required for updating the pointlist from ref. to actual config. below
		Matrix2d rotMatrix;
		rotMatrix << cos(_theta), -sin(_theta), sin(_theta), cos(_theta);

		for (size_t i = 0; i < _pointList.size(); i++) {
			_nodeShears[i] = 0;
			_nodeContact[i] = 0;
			_pointList[i] = rotMatrix*_pointList[i] + _position;
			// _nodeNormals[i] << 0,0;
		}
		_scalePixToMM = 15/1000.;
  //       double cosd = cos(_theta);
		// double sind = sin(_theta);
		// for (size_t ptid = 0; ptid < _refpointList.size(); ptid++) {
		// 	_pointList[ptid] << _refpointList[ptid](0)*cosd - _refpointList[ptid](1)*sind,
		// 									  _refpointList[ptid](0)*sind + _refpointList[ptid](1)*cosd;
		// 	_pointList[ptid] += _position;
		// }
        
       	_rsq = _radius*_radius;
		_fracFlag = false;
		_yield = 1e10 * 1e-9; //TODO Reduce for more break
		_critStress = 0;
		_remove = false;

		size_t closestpt = 0; double closestdist = 1e10;
		for (size_t i=0;i<_pointList.size();i++){
			double squareddist = _refpointList[i](0)*_refpointList[i](0)+_refpointList[i](1)*_refpointList[i](1);
			if (squareddist < closestdist) {
				closestdist = squareddist;
				closestpt = i;
			}
		}
		size_t otherpt = 0; double otherdist = 1e10;
		for (size_t i=0;i<_pointList.size();i++){
			double squareddist = fabs(M_PI-fabs(atan2(_refpointList[i](0),_refpointList[i](1))-atan2(_refpointList[closestpt](0),_refpointList[closestpt](1))));
			if (squareddist < otherdist) {
				otherdist = squareddist;
				otherpt = i;
			}
		}
		_d = sqrt(pow(_refpointList[otherpt](0)-_refpointList[closestpt](0),2)+pow(_refpointList[otherpt](1)-_refpointList[closestpt](1),2));

	}

	// 1st level check
	bool bcircleGrainIntersection(const Grain2d & other) const {
		if ((other.getPosition()(0)-_position(0))*(other.getPosition()(0)-_position(0)) +
			 (other.getPosition()(1)-_position(1))*(other.getPosition()(1)-_position(1)) <
			 (other.getRadius()+_radius)*(other.getRadius()+ _radius) ) {
			return true;
		}
		return false;
	}
	
	//For regular DEM
	bool vRadiusCheck(const Grain2d & other) const {
		Vector2d d = other.getPosition() - _position;
		return d.norm() < _vr*_vr;
	}
	
	bool vRadiusCheck2(const Grain2d & other) const {
		Vector2d d = other.getPosition() - _position;
		return d.norm() < _vr*0.5 + _cohesiveDistance*1.2;
	}

	// 1st level check for Wall
	// bool bcircleGrainIntersection(const Grain2d & other) const {
	// 	if ((other.getPosition()(0)-_position(0))*(other.getPosition()(0)-_position(0)) +
	// 		 (other.getPosition()(1)-_position(1))*(other.getPosition()(1)-_position(1)) <
	// 		 (other.getRadius()+_radius)*(other.getRadius()+ _radius) ) {
	// 		return true;
	// 	}
	// 	return false;
	// }
	
	//Ice contact function (give grain and give point)
	//Find if Ocean is exposed to ice or atmosphere (improve for PBC)
    bool under_ice(Vector2d & ptDetect)
    {
        bool            result_cover = false; //We assume no ice cover by default and prove the opposite.
        double          penetration;    // penetration amount (is negative by convention of the level set)
        Vector2d        normal;         // surface normal
        Vector2d        ptOtherCM;      // point wrt the center of mass of other in real space
        Vector2d        ptOtherLset;    // point in the reference config of other's level set
    
        //Start with a simple BBox radius criterion for simplicity (improve for PBC)
        if (  (ptDetect-_position).norm() <= _radius  )  //No need for PBC for initialization, but needs for shifting grains
        {
            
            //You need to shift ptDetect to the right coordinate system before checking penetration
            const double cosdet = cos(_theta);
            const double sindet = sin(_theta);
            ptOtherCM = ptDetect - _position;
            ptOtherLset(0) =  ptOtherCM(0)*cosdet + ptOtherCM(1)*sindet;
            ptOtherLset(1) = -ptOtherCM(0)*sindet + ptOtherCM(1)*cosdet;
            ptOtherLset += _cmLset;
    
    
            //If within BBox radius now check if grid point penetrates ice level set
            if ( _lset.isPenetration(ptOtherLset, penetration, normal) )  
            {
                result_cover = true;
            }
        }
    
        return result_cover;
    }
	
	bool under_iceDEM(Vector2d & ptDetect)
    {
        bool            result_cover = false; //We assume no ice cover by default and prove the opposite.
        
        //Start with a simple BBox radius criterion for simplicity (improve for PBC)
        if (  (ptDetect-_position).norm() <= _radius  )  //No need for PBC for initialization, but needs for shifting grains
        {
                result_cover = true;
        }
    
        return result_cover;
    }
	
	

	// 1st level check for periodic bcs with a particular offset along x
	bool bcircleGrainIntersectionXOffset(const Grain2d & other, const double & offset) const {
		if ((other.getPosition()(0)-_position(0)+offset)*(other.getPosition()(0)-_position(0)+offset) +
			 (other.getPosition()(1)-_position(1))*(other.getPosition()(1)-_position(1)) <
			 (other.getRadius()+_radius)*(other.getRadius()+ _radius) ) {
			return true;
		}
		else if ((other.getPosition()(0)-_position(0)-offset)*(other.getPosition()(0)-_position(0)-offset) +
			 (other.getPosition()(1)-_position(1))*(other.getPosition()(1)-_position(1)) <
			 (other.getRadius()+_radius)*(other.getRadius()+ _radius) ) {
			return true;
		}
		return false;
	}
	// 1st level check for periodic bcs with a particular offset along x
	bool bcircleGrainIntersectionYOffset(const Grain2d & other, const double & offset) const {
		if ((other.getPosition()(0)-_position(0))*(other.getPosition()(0)-_position(0)) +
			 (other.getPosition()(1)-_position(1)+offset)*(other.getPosition()(1)-_position(1)+offset) <
			 (other.getRadius()+_radius)*(other.getRadius()+ _radius) ) {
			return true;
		}
		else if ((other.getPosition()(0)-_position(0))*(other.getPosition()(0)-_position(0)) +
			 (other.getPosition()(1)-_position(1)-offset)*(other.getPosition()(1)-_position(1)-offset) <
			 (other.getRadius()+_radius)*(other.getRadius()+ _radius) ) {
			return true;
		}
		return false;
	}

	// Actual (2nd level) contact check between *this and other. If there is no contact, returns false.
	// Compares points of *this to the level set of other.
	// If there is contact, returns true and updates force, which is the force on *this,
	// thisMoment, which is the moment on *this, and otherMoment, the moment on other.
	bool findInterparticleForceMoment(Grain2d & other, const double & dt, Vector2d & force, double & thisMoment, double & otherMoment, size_t & nContacts, const Vector2d & offset, const size_t & numGrainWalls) {

		// Declare temporary variables
	    force << 0., 0.;	    								 // force at a contact point (shared by the two grains in contact)
	    thisMoment = 0.; 										 // moment on this particle due to a contact
	    otherMoment = 0.;										 // moment on the other particle due to a contact
	    nContacts = 0;											 // number of contacts for this interaction
		Vector2d		   ptOtherCM; 		                     // point wrt the center of mass of other in real space
		Vector2d		   ptOtherLset; 	                     // point in the reference config of other's level set
		double		       penetration;	                         // penetration amount (is negative by convention of the level set)
		Vector2d		   normal; 			                     // surface normal pointing from other to *this
		const double       cos2 = cos(other.getTheta());
		const double       sin2 = sin(other.getTheta());
		Vector2d		   ptThisCM; 		                     // point wrt the center of mass of this in real space
		Vector2d		   df; 				                     // force increment from a single point
		Vector2d		   tangent; 		                     // surface tangent
		double		       sdot; 			                     // relative velocity of a point of contact
		double 			   ndot;
		double		       ds;									 // relative displacement of a point of contact
		double 			   dn;
		Vector2d		   Fs; 				                     // vector of frictional/shear force
		Vector2d           relativeVel;							 // vector of relative velocity between the two grains

		Vector2d		   forceNormal;
		double 			   sigMoment;

		double Fsmag;											 // magnitude of shear force
		bool isContact = false;									 // return flag (assigned to true if contact exists)
		
		//Validate for large difference in mass (NEW MODIF)
		double diff_factor = 100;
		double adjust_redux = 10;  //(NEW MODIF)
		double kn_redux = 1;
		if ( _mass > other.getMass() * diff_factor || _mass * diff_factor < other.getMass() ) //Identify large mass difference
		{
		    kn_redux = adjust_redux*(1/diff_factor); //Change reduction factor from identity (unity) to an actual reduction
		    //cout << "GRAIN MASS DIFFERENCE HUGE!*!*!" << endl;
		}
		
		double cres = .30;										 // auxilliary variable for computing contact damping coeff below   .3
		double GamaN = -2*sqrt(_kn*_mass*other.getMass()/(_mass+other.getMass()))*log(cres)/sqrt(M_PI*M_PI + log(cres)*log(cres));

		Vector2d offsetVec = Vector2d(0,0);
		if (other.getPosition()(0) > _position(0)) {
			offsetVec += Vector2d( -offset(0), 0.);
		}
		else {
			offsetVec += Vector2d( offset(0), 0.);
		}
		if (other.getPosition()(1) > _position(1)) {
			offsetVec += Vector2d( 0, -offset(1));
		}
		else {
			offsetVec += Vector2d( 0, offset(1));
		}


		// Frictional interaction
        for (size_t ptidx = 0; ptidx < _pointList.size(); ptidx++) {
			ptOtherCM = _pointList[ptidx] - other.getPosition() - offsetVec;
			ptOtherLset(0) =  ptOtherCM(0)*cos2 + ptOtherCM(1)*sin2;
			ptOtherLset(1) = -ptOtherCM(0)*sin2 + ptOtherCM(1)*cos2;
			ptOtherLset += other.getCmLset();

			//	if there is penetration, finds forces and moments due to each contact
			if ( other.getLset().isPenetration(ptOtherLset, penetration, normal) ) {

				isContact = true;
				nContacts++;
                _ncontacts++;
                


                //Needed for breakage, specially with Wall contact, optional for rest. //----------------------------------->> !!!!
                //if (numGrainWalls > 0)
                //{
                	ptThisCM = _pointList[ptidx] - _position;
                //}
				
				

				// rotate the normal from the reference config of 2's level set to real space
				normal << normal(0)*cos2 - normal(1)*sin2, normal(0)*sin2 + normal(1)*cos2;
				// find the tangent
				tangent << -normal(1), normal(0);
				// update force: normal force contribution
				relativeVel << other.getVelocity()(0) - other.getOmega()*ptOtherCM(1) - (_velocity(0) - _omega*ptThisCM(1)),
				               other.getVelocity()(1) + other.getOmega()*ptOtherCM(0) - (_velocity(1) + _omega*ptThisCM(0));
				// viscoelastic contact law
				df = penetration*normal*_kn - GamaN*normal.dot(relativeVel)*normal;
				
				//Validate for large difference in mass (NEW MODIF)
				df *= kn_redux; //Reduce force for large mass difference
				
				if (isnan(df(0)) == 1 || isnan(df(1)) == 1 )
				{
					//cout << "penetration: " << penetration << endl;
					//cout << "normal: " << normal << endl;
					//cout << "GamaN: " << GamaN << endl;
					//cout << "RelativeVel X: " << other.getVelocity()(0) - other.getOmega()*ptOtherCM(1) - (_velocity(0) - _omega*ptThisCM(1)) << endl;
					//cout << "Term 1: " << _velocity(0) << endl;
					//cout << "Term 2: " << _omega << endl;
					//cout << "Term 3: " << ptThisCM(1) << endl;
					//cout << "Term 3a: " << ptThisCM(0) << endl;
					//cout << "Term 4: " << _omega*ptThisCM(1) << endl;
				}
				force -= df;
				// update moments: eccentric loading contribution
				otherMoment -= df(0)*ptOtherCM(1) - df(1)*ptOtherCM(0);
				thisMoment += df(0)*ptThisCM(1) - df(1)*ptThisCM(0);
				// force/moment calculations based on friction
				//sdot = tangent.dot(_velocity - other.getVelocity()) - (_omega*ptThisCM.norm() + other.getOmega()*ptOtherCM.norm());
				sdot = tangent.dot(-relativeVel);
				ds = sdot*dt;
				_nodeContact[ptidx] = other.getId();
				// Static friction law
				_nodeShears[ptidx] += _ks*ds;	// elastic predictor
				if (_nodeShears[ptidx] > 0) {
					Fsmag = std::min(_nodeShears[ptidx], df.norm()*_mu );
				}
				else {
					Fsmag = std::max(_nodeShears[ptidx], -df.norm()*_mu );
				}
				Fs = -tangent*Fsmag;
				_nodeShears[ptidx] = Fsmag;
				force += Fs;
				thisMoment -= Fs(0)*ptThisCM(1) - Fs(1)*ptThisCM(0);
				otherMoment += Fs(0)*ptOtherCM(1) - Fs(1)*ptOtherCM(0);
			}
			// if there is no contact between the point and other, reset the shear force if other was the last contact
			else if (_nodeContact[ptidx] == other.getId() ){
				_nodeContact[ptidx] = 0;
				_nodeShears[ptidx] = 0;
				// _nodeNormals[ptidx] << 0,0;
			}
		}
		
		//Bonded Particle Method Properties (May 7, 2023)
        // bond interactions
        //cout << "Do bonding if needed" << endl;
        if (_bondInfo.size() > 0){
            for (size_t bondidx = 0; bondidx < _bondInfo.size(); bondidx++) {
    			if (size_t(_bondInfo[bondidx][0]) == 1 && size_t(_bondInfo[bondidx][2]) == other.getId()) { // intact bonds
    				isContact = true;
    			//	ncontacts++;
    				size_t ptidx = _bondpts[bondidx];
    				ptThisCM = _pointList[ptidx] - _position;
    				ptOtherCM = _pointList[ptidx] - other.getPosition();
    				ptOtherLset(0) =  ptOtherCM(0)*cos2 + ptOtherCM(1)*sin2;
    				ptOtherLset(1) = -ptOtherCM(0)*sin2 + ptOtherCM(1)*cos2;
    				ptOtherLset += other.getCmLset();
    				double p2pdist = 0; // not used
    				if (other.getLset().getNormal(ptOtherLset, _cohesiveDistance, p2pdist, normal)){ //Change function to getNormal, use isCohesion only for bond formation
    					// relative v of this wrt the other
    					relativeVel << -other.getVelocity()(0) + other.getOmega()*ptOtherCM(1) + (_velocity(0) - _omega*ptThisCM(1)),
    				               	   -other.getVelocity()(1) - other.getOmega()*ptOtherCM(0) + (_velocity(1) + _omega*ptThisCM(0));
    					
    					//normal point out of *this
    					normal << -normal(0)*cos2 + normal(1)*sin2, -normal(0)*sin2 - normal(1)*cos2; 
    					normal = normal/(normal.norm()+DBL_MIN);
    
    	               	// normal and shear displacement increments
    					Vector2d sdot_bond = (relativeVel - relativeVel.dot(normal)*normal);
    					Vector2d ndot_bond = normal.dot(-relativeVel)*normal;
    					Vector2d ds_bond = sdot_bond*dt;
    					Vector2d dn_bond = ndot_bond*dt;
    
    					// normal and shear rotation increments
    					double w = _omega - other.getOmega();
    					double sOmega = w;
    					double sTheta = sOmega*dt;
    
    					//rotate shear components
    					double k = -normal(0)*_bondNormals[bondidx](1) + normal(1)*_bondNormals[bondidx](0); // axis of rotation from n_old to n_cur
    					double sint = abs(k); // compute the sin and cos of the rotation //double sint = abs(k); or Magnitude
    					double cost = normal.dot(_bondNormals[bondidx]);
    					k = k/(sint+DBL_MIN); // normalize k to unit magnitude, also don't divide by 0  //SIgned
    					_bondForceShear[bondidx] = _bondForceShear[bondidx]*cost + Vector2d(-_bondForceShear[bondidx](1), _bondForceShear[bondidx](0))*k*sint;
    					_bondForceNormal[bondidx] = _bondForceNormal[bondidx]*cost + Vector2d(-_bondForceNormal[bondidx](1), _bondForceNormal[bondidx](0))*k*sint;
    					_bondThisMomentNormal[bondidx] = _bondThisMomentNormal[bondidx]*cost + k*k*_bondThisMomentNormal[bondidx]*(1.-cost);
    					_bondNormals[bondidx] = normal;
    
    					// bond forces and moments
    					_bondForceNormal[bondidx] += _kn_bond*dn_bond*_bondAreas;
    					_bondForceShear[bondidx] += -_ks_bond*ds_bond*_bondAreas;
    
    					// Check for yield
    					_fSig[bondidx] = fabs(normal.dot(_bondForceNormal[bondidx])/_bondAreas); // magnitute of compressive or tensile stress
    					_fTau[bondidx] = _bondForceShear[bondidx].norm()/_bondAreas;
    					
    					// apply your selected failure criteria. For example:
    					//if (_fSig[bondidx] > _sigC || _fTau[bondidx] > _tauC) {
    					if (_fSig[bondidx] > _sigC * _sigrF[bondidx] || _fTau[bondidx] > _tauC * _taurF[bondidx]) {
    					    cout << "Bond broken due to exceeding bond strength!!!" << endl;
    					    cout << "fSig[bondidx]: " << _fSig[bondidx] << " sigC: " << _sigC << " fTau[bondidx]: " << _fTau[bondidx] << " tauC: " << _tauC << endl;
    						// bonds broken
    						_bondInfo[bondidx][0] = 0;
    						_bondForceNormal[bondidx] = Vector2d(0,0);
    						_bondForceShear[bondidx] = Vector2d(0,0);
    						_bondThisMomentNormal[bondidx] = 0;
    						_fSig[bondidx] = 0;
    						_fTau[bondidx] = 0;
    						_bondNormals[bondidx] = Vector2d(0,0);
    					}else{
    						df = _bondForceNormal[bondidx] + _bondForceShear[bondidx];
    						df += ndot_bond*2.*sqrt(_kn_bond*_bondAreas*(_mass*other.getMass()/(_mass+other.getMass()))) - sdot_bond*2.*sqrt(_ks_bond*_bondAreas*(_mass*other.getMass()/(_mass+other.getMass()))); // Damping
    						force += df;
    						thisMoment += _bondThisMomentNormal[bondidx] + (ptThisCM(0)*df(1) - ptThisCM(1)*df(0));//moment of beam + from beam
    						otherMoment += -_bondThisMomentNormal[bondidx] + (-ptOtherCM(0)*df(1) + ptOtherCM(1)*df(0));
    
    						// double nMOIThis = (_momentInertia * normal).dot(normal) + _mass*(ptThisCM(0)*normal(1) - ptThisCM(1)*normal(0)).squaredNorm();  //For off-plane moment ?
    						// double nMOIOther = (other.getMomentInertia() * -normal).dot(-normal) + _mass*(-ptOtherCM(0)*normal(1) + ptOtherCM(1)*normal(0)).squaredNorm(); //For twist moment ?
    
    						//Additional damping that might be redundant.
    						// if (sOmega.norm() > 0) {
    						// 	Vector2d sheardir = sOmega/sOmega.norm(); //SIgn
    						// 	double shearMOIthis = (_momentInertia * sheardir).dot(sheardir) + _mass*(ptThisCM(0)*sheardir(1) - ptThisCM(1)*sheardir(0)).squaredNorm();
    						// 	double shearMOIother = (other.getMomentInertia() * -sheardir).dot(-sheardir) + _mass*(-ptOtherCM(0)*sheardir(1) + ptOtherCM(1)*sheardir(0)).squaredNorm();
    
    						// 	thisMoment  -= sOmega*2.*sqrt(_kn_bond*_bondMomentsInertia*shearMOIthis);
    						// 	otherMoment += sOmega*2.*sqrt(_kn_bond*_bondMomentsInertia*shearMOIother);
    						// }
    					}
    				}else{
    					// When in this loop, it means that bonds have not failed yet, but this surface node
    					// has fallen out of the other's level set, or the cohesive distance, or has penetrated
    					// *this. Depend on different applications, you can consider different approaches.
    					// Maybe re-adjust failure parameter to make sure bonds break before this happens.
    					// Or, you can simply set the bonds to be broken like the following
    					cout << "Bond broken due to falling out of level set" << endl;
    					_bondInfo[bondidx][0] = 0;
    					_bondForceNormal[bondidx] = Vector2d(0,0);
    					_bondForceShear[bondidx] = Vector2d(0,0);
    					_bondThisMomentNormal[bondidx] = 0;
    					_fSig[bondidx] = 0;
    					_fTau[bondidx] = 0;
    					_bondNormals[bondidx] = Vector2d(0,0);
    				}
    			}
    		} // compute for all bond reaction force	
        }
		return isContact;
	} // end findInterparticleForceMoment

    //Same as above, but for regular DEM
	bool findInterparticleForceMomentDEM(Grain2d & other, const double & dt, Vector2d & force, double & thisMoment, double & otherMoment, size_t & nContacts, const Vector2d & offset, const size_t & numGrainWalls, size_t & ntension, size_t & nshear) {

		// Declare temporary variables
	    force << 0., 0.;	    								 // force at a contact point (shared by the two grains in contact)
	    thisMoment = 0.; 										 // moment on this particle due to a contact
	    otherMoment = 0.;										 // moment on the other particle due to a contact

		Vector2d ptThisCM; 		// point wrt the center of mass of *this in real space
		Vector2d ptOtherCM; 		// point wrt the center of mass of other in real space
		Vector2d normal; 			// surface normal pointing out of other in the reference config (initially) and then pointing out of *this in real space (after)
		Vector2d Fn; 				// normal force from a single point
		Vector2d v; 				// velocity of *this relative to other at a point of contact
		Vector2d relativeVel; 
		Vector2d Fs; 				// vector of frictional/shear force
		Vector2d ds, dn;				// tangential increment (projection of relative velocity v in tangengial direction)
		double Fsmag;
		double sint, cost;
		Vector2d Qs, Pn;
		Vector2d sdot, ndot; 
		double w, k, sOmega, nOmega, sTheta, nTheta, Kn;
		double bondKn, bondKs;
		double sigC, tauC;
		Vector2d dFn,dFs,df;
		double dMn,dMs;
		
		double distmag;
		Vector2d dist;
		
		Vector2d offsetVec = Vector2d(0,0);
		if (other.getPosition()(0) > _position(0)) {
			offsetVec += Vector2d( -offset(0), 0.);
		}
		else {
			offsetVec += Vector2d( offset(0), 0.);
		}
		if (other.getPosition()(1) > _position(1)) {
			offsetVec += Vector2d( 0, -offset(1));
		}
		else {
			offsetVec += Vector2d( 0, offset(1));
		}

		Vector2d dr = _position - other.getPosition() - offsetVec; //Test PBCs, we need to check for this both for penetration and ptOtherCM
		double drMag = dr.norm();
		double sigma = _radius + other.getRadius();
		double penetration = drMag - sigma ; //negative if contact 
		size_t ptidx = other.getId();

		bool isContact = false;									 // return flag (assigned to true if contact exists)
		
		//Validate for large difference in mass (NEW MODIF)
		double diff_factor = 100;
		double adjust_redux = 10;  //(NEW MODIF)
		double kn_redux = 1;
		if ( _mass > other.getMass() * diff_factor || _mass * diff_factor < other.getMass() ) //Identify large mass difference
		{
		    kn_redux = adjust_redux*(1/diff_factor); //Change reduction factor from identity (unity) to an actual reduction
		    cout << "GRAIN MASS DIFFERENCE HUGE!*!*!" << endl;
		}
		
		double cres = .30;										 // auxiliary variable for computing contact damping coeff below   .3
		double GamaN = -2*sqrt(_kn*_mass*other.getMass()/(_mass+other.getMass()))*log(cres)/sqrt(M_PI*M_PI + log(cres)*log(cres));

		// Frictional interaction
		
		// iterate through all of the points of *this and check for contact for each one
		if (penetration <= 0.){
			isContact = true;
			// point out of *this to match KWL's 3d GEM paper
			normal = -dr/drMag;
			ptThisCM = _radius * normal;
			ptOtherCM = ptThisCM + _position - other.getPosition() - offsetVec;
			// update force: normal force contribution
			// note: penetration is negative which negates the negative sign in eq (2) of KWL
			//v = _velocity - other.getVelocity() + _omega*(ptThisCM) - other.getOmegaGlobal().cross(ptOtherCM); // eq (6) //3D
			v << -other.getVelocity()(0) + other.getOmega()*ptOtherCM(1) + (_velocity(0) - _omega*ptThisCM(1)),
		         -other.getVelocity()(1) - other.getOmega()*ptOtherCM(0) + (_velocity(1) + _omega*ptThisCM(0));
			
			Fn = penetration*normal*_kn - GamaN*normal.dot(v)*normal;
			Fn *= kn_redux;
			force += Fn;
			
			//Control
			if ( isnan(force(0)) || isnan(force(1))  || isnan(thisMoment)  || isnan(otherMoment)  ){
			    cout << "Force Contact: " << force(0) << " " << force(1) << endl;
			    cout << "Moments Contact: " << thisMoment << " " << otherMoment << endl;
			    cout << "penetration: " << penetration << " dr: " << dr << " v: " << v(0) << " " << v(1) << endl;
			}
			
			thisMoment += (ptThisCM(0)*Fn(1) - ptThisCM(1)*Fn(0)); 
    		otherMoment += (-ptOtherCM(0)*Fn(1) + ptOtherCM(1)*Fn(0));

			ds = (v - v.dot(normal)*normal)*dt; // eq (7) and (9)

            //Rotate the shear force into the new tangent direction (Danilo Style)
            if(_nodeNormals.find(ptidx) != _nodeNormals.end()){
                Qs = _nodeShearsv2[ptidx];
                Pn = _nodeNormals[ptidx];
                //Kn = Pn.cross(normal); //Rob's amendment 3D
                //2D
                Kn = (Pn(0)*normal(1) - Pn(1)*normal(0)); 

                sint = abs(Kn);
                cost = normal.dot(normal); //CHECK ASSUMPTION!!! sqrt(1.0 - sint*sint);

                Kn = Kn/(sint + 10.0*DBL_MIN);
                _nodeShearsv2[ptidx] = cost*Qs + sint*Kn*(Qs) + (1.0 - cost)*Kn*(Qs)*Kn; //CHECK!!!
            }
            else{
                _nodeShearsv2[ptidx] = Eigen::Vector2d::Zero(2);
            }
            
            _nodeNormals[ptidx] = normal;
            _nodeShearsv2[ptidx] -= _ks*ds; // eq (8) and (12)

			Fsmag = min(Fn.norm()*_mu, _nodeShearsv2[ptidx].norm() ); // eq (14)
			if (Fsmag > 0) {
				Fs = Fsmag*_nodeShearsv2[ptidx]/_nodeShearsv2[ptidx].norm(); // eq (13)
				_nodeShearsv2[ptidx] = Fs;
				force += Fs;
			    thisMoment += (ptThisCM(0)*Fs(1) - ptThisCM(1)*Fs(0)); 
    		    otherMoment += (-ptOtherCM(0)*Fs(1) + ptOtherCM(1)*Fs(0));

			}
			//Control
			if ( isnan(force(0)) || isnan(force(1))  || isnan(thisMoment)  || isnan(otherMoment)  ){
                cout << "Force Contact: " << force(0) << " " << force(1) << endl;
			    cout << "Moments Contact: " << thisMoment << " " << otherMoment<< endl;
			    cout << "Fs: " << Fs(0) << " " << Fs(1) << " _nodeShearsv2[ptidx]: " << _nodeShearsv2[ptidx](0)  << " " << _nodeShearsv2[ptidx](1) << endl;
			}
			
		}else if (_nodeNormals.find(ptidx) !=  _nodeNormals.end()){
			// if there is no contact between the point and other, reset the shear force if other was the last contact
			_nodeNormals.erase(ptidx);
            _nodeShearsv2.erase(ptidx);
		}
		
		// Bond Interaction v1
// 		if (_bondInfo.size() > 0){
//     		for (size_t bondidx = 0; bondidx < _bondInfo.size(); bondidx++) {
//     			if (size_t(_bondInfo[bondidx][0]) == 1 && size_t(_bondInfo[bondidx][2]) == other.getId()) { // intact bonds
//     				isContact = true;
//     				size_t bpt1 = _bondpts[bondidx];
//     				size_t bpt2 = _bondOtherpts[bondidx]; // point to point bond
    				
//     				normal = (other.getPointList()[bpt2] - _pointList[bpt1])/((other.getPointList()[bpt2] - _pointList[bpt1]).norm()+DBL_MIN); //normal pointing out of this and into other
    
//     				ptThisCM = _pointList[bpt1] - _position;
//     				ptOtherCM = _pointList[bpt1] - other.getPosition();
    				 
// 					// relative v of this wrt the other
// 					v <<  -other.getVelocity()(0) + other.getOmega()*ptOtherCM(1) + (_velocity(0) - _omega*ptThisCM(1)),
// 				          -other.getVelocity()(1) - other.getOmega()*ptOtherCM(0) + (_velocity(1) + _omega*ptThisCM(0));
    
//     				// normal and shear displacement increments
//     				sdot = (v - v.dot(normal)*normal);
//     				ndot = normal.dot(-v)*normal;
//     				ds = sdot*dt;
//     				dn = ndot*dt;
    
//     				// normal and shear rotation increments
//     				w = _omega - other.getOmega();
// 					sOmega = w; //Check!!!
// 					sTheta = sOmega*dt;
// 					nOmega = w; 
// 					nTheta = nOmega*dt;
    
//     				//rotate shear components
// 					k = -normal(0)*_bondNormals[bondidx](1) + normal(1)*_bondNormals[bondidx](0); // axis of rotation from n_old to n_cur
// 					sint = abs(k); // compute the sin and cos of the rotation //double sint = abs(k); or Magnitude
// 					cost = normal.dot(_bondNormals[bondidx]);
// 					k = k/(sint+DBL_MIN); // normalize k to unit magnitude, also don't divide by 0  //SIgned
// 					_bondForceShear[bondidx] = _bondForceShear[bondidx]*cost + Vector2d(-_bondForceShear[bondidx](1), _bondForceShear[bondidx](0))*k*sint;
// 					_bondForceNormal[bondidx] = _bondForceNormal[bondidx]*cost + Vector2d(-_bondForceNormal[bondidx](1), _bondForceNormal[bondidx](0))*k*sint;
// 					_bondThisMomentNormal[bondidx] = _bondThisMomentNormal[bondidx]*cost + k*k*_bondThisMomentNormal[bondidx]*(1.-cost);
// 					_bondThisMomentShear[bondidx] = _bondThisMomentShear[bondidx]*cost + k*k*_bondThisMomentShear[bondidx]*(1.-cost);
// 					_bondNormals[bondidx] = normal;
					
// 					//Control
// 					if ( isnan(_bondForceShear[bondidx](0)) || isnan(_bondForceShear[bondidx](1)) || isnan(_bondForceNormal[bondidx](0))  || isnan(_bondForceNormal[bondidx](1))  || isnan(_bondThisMomentNormal[bondidx]) ||  isnan(_bondThisMomentShear[bondidx]) ){
// 					    cout << "Bond Force Shear Contact: " << _bondForceShear[bondidx](0) << " " << _bondForceShear[bondidx](1) << endl;
// 					    cout << "Bond Force Normal Contact: " << _bondForceNormal[bondidx](0) << " " << _bondForceNormal[bondidx](1) << endl;
// 			            cout << "Bond Moments Contact: " << _bondThisMomentNormal[bondidx] << " " << _bondThisMomentShear[bondidx]<< endl;
// 			            cout << "Bond Normals: " <<  _bondNormals[bondidx](0) << " " << _bondNormals[bondidx](1) << endl;
// 			            cout << "sint/ cost: " <<  sint << " " << cost << endl;
// 					}
    				
//     				//Optional stiffen outer bonds (special applications)
//     				// if (_isSurface || other.isSurface()){
//     				// 	sigC = _sigC * _stiffVar;
//     				// 	tauC = _tauC * _stiffVar;
//     				// 	bondKn = _kn_bond * _stiffVar;
//     				// 	bondKs = _ks_bond * _stiffVar;
//     				//}else{
//     				sigC = _sigC;
//     				tauC = _tauC;
//     				bondKn = _kn_bond;
//     				bondKs = _ks_bond;
//     				//}
    				
// 				    // bond forces and moments
// 					_bondForceNormal[bondidx] += _kn_bond*dn*_bondAreas;
// 					_bondForceShear[bondidx] += -_ks_bond*ds*_bondAreas;
					
// 					dMn = -bondKs*_bondMomentsInertia*nTheta; //CHECK!!!
//     				dMs = -bondKn*_bondMomentsInertia*sTheta;
//     				_bondThisMomentNormal[bondidx] += dMn;
//     				_bondThisMomentShear[bondidx]  += dMs;
    				
//     				//Control
//     				if ( isnan(_bondForceShear[bondidx](0)) || isnan(_bondForceShear[bondidx](1)) || isnan(_bondForceNormal[bondidx](0))  || isnan(_bondForceNormal[bondidx](1))  || isnan(_bondThisMomentNormal[bondidx]) ||  isnan(_bondThisMomentShear[bondidx]) ){
//     				    cout << "Bond Force Shear Contact: " << _bondForceShear[bondidx](0) << " " << _bondForceShear[bondidx](1) << endl;
// 					    cout << "Bond Force Normal Contact: " << _bondForceNormal[bondidx](0) << " " << _bondForceNormal[bondidx](1) << endl;
// 			            cout << "Bond Moments Contact: " << _bondThisMomentNormal[bondidx] << " " << _bondThisMomentShear[bondidx]<< endl;
// 			            cout << "Bond Normals: " <<  _bondNormals[bondidx](0) << " " << _bondNormals[bondidx](1) << endl;
// 			            cout << "dn / ds: " << dn << " " << ds << endl;
// 			            cout << "Normal: " << normal(0) << " " << normal(1) << endl;
// 			            cout << "v: " << v(0) << " " << v(1) << endl;
//     				}

// 					// Check for yield
// 				    _fSig[bondidx] = fabs(normal.dot(_bondForceNormal[bondidx])/_bondAreas); // magnitude of compressive or tensile stress
// 				    _fTau[bondidx] = _bondForceShear[bondidx].norm()/_bondAreas;
// 				    // Check for yield (Potyondy and Cundall)
//     				//_fSig[bondidx] = fabs(normal.dot(_bondForceNormal[bondidx])/_bondAreas + fabs(_bondThisMomentShear[bondidx])*_bondRadius/_bondMomentsInertia; // magnitute of compressive or tensile stress
//     				//_fTau[bondidx] = _bondForceShear[bondidx].norm()/(_bondAreas) + fabs(_bondThisMomentNormal[bondidx])*_bondRadius/_bondMomentsInertia;
					
// 					// apply your selected failure criteria. For example:
// 					if (_fSig[bondidx] > _sigC || _fTau[bondidx] > _tauC) {
// 					    cout << "Bond broken due to exceeding bond strength!!!" << endl;
// 					    cout << "fSig[bondidx]: " << _fSig[bondidx] << " sigC: " << _sigC << " fTau[bondidx]: " << _fTau[bondidx] << " tauC: " << _tauC << endl;
// 						// bonds broken
// 						_bondInfo[bondidx][0] = 0;
// 						_bondForceNormal[bondidx] = Vector2d(0,0);
// 						_bondForceShear[bondidx] = Vector2d(0,0);
// 						_bondThisMomentNormal[bondidx] = 0;
// 						_bondThisMomentShear[bondidx] = 0;
// 						_fSig[bondidx] = 0;
// 						_fTau[bondidx] = 0;
// 						_bondNormals[bondidx] = Vector2d(0,0);
// 					}else{
// 						Vector2d Fbond = _bondForceNormal[bondidx] + _bondForceShear[bondidx];
// 						Fbond += ndot*2.*sqrt(_kn_bond*_bondAreas*(_mass*other.getMass()/(_mass+other.getMass()))) - sdot*2.*sqrt(_ks_bond*_bondAreas*(_mass*other.getMass()/(_mass+other.getMass()))); // Damping
// 						force += Fbond;
// 						thisMoment += _bondThisMomentNormal[bondidx] + (ptThisCM(0)*df(1) - ptThisCM(1)*df(0));//moment of beam + from beam
// 						otherMoment += -_bondThisMomentNormal[bondidx] + (-ptOtherCM(0)*df(1) + ptOtherCM(1)*df(0));


// 			            //Control
// 			            if ( isnan(force(0)) || isnan(force(1))  || isnan(thisMoment)  || isnan(otherMoment)  ){
//             			    cout << "BG Force Contact: " << force(0) << " " << force(1) << endl;
//             			    cout << "BG Moments Contact: " << thisMoment << " " << otherMoment<< endl;
//             			    cout << "BG Fbond: " << Fbond(0) << " " << Fbond(1) << endl;
//             			    cout << "ndot / sdot: " << ndot << " " << sdot << endl;
// 			            }
			            
// 						// double nMOIThis = (_momentInertia * normal).dot(normal) + _mass*(ptThisCM(0)*normal(1) - ptThisCM(1)*normal(0)).squaredNorm();  //For off-plane moment ?
// 						// double nMOIOther = (other.getMomentInertia() * -normal).dot(-normal) + _mass*(-ptOtherCM(0)*normal(1) + ptOtherCM(1)*normal(0)).squaredNorm(); //For twist moment ?

// 						//Additional damping that might be redundant.
// 						// if (sOmega.norm() > 0) {
// 						// 	Vector2d sheardir = sOmega/sOmega.norm(); //SIgn
// 						// 	double shearMOIthis = (_momentInertia * sheardir).dot(sheardir) + _mass*(ptThisCM(0)*sheardir(1) - ptThisCM(1)*sheardir(0)).squaredNorm();
// 						// 	double shearMOIother = (other.getMomentInertia() * -sheardir).dot(-sheardir) + _mass*(-ptOtherCM(0)*sheardir(1) + ptOtherCM(1)*sheardir(0)).squaredNorm();

// 						// 	thisMoment  -= sOmega*2.*sqrt(_kn_bond*_bondMomentsInertia*shearMOIthis);
// 						// 	otherMoment += sOmega*2.*sqrt(_kn_bond*_bondMomentsInertia*shearMOIother);
// 						// }
// 					}
//     			}
//     		} // compute for all bond reaction force
// 		}

        // Bond Interaction v2  
        //Bonded Particle Method Properties (May 7, 2023)
        // bond interactions
        //cout << "Do bonding if needed" << endl;
        size_t ptidxf; //For follower only point, if needed //BB
        if (_bondInfo.size() > 0){
            for (size_t bondidx = 0; bondidx < _bondInfo.size(); bondidx++) {
                //cout << "Bond: " << _bondInfo[bondidx][0] << " other: " << _bondInfo[bondidx][2] << " other2: " << other.getId() << endl;
    			if (size_t(_bondInfo[bondidx][0]) == 1 && size_t(_bondInfo[bondidx][2]) == other.getId()) { // intact bonds
    				isContact = true;
    			//	ncontacts++;
    				ptidx = _bondpts[bondidx];
    				ptidxf = _bondptsf[bondidx]; //BB
    				ptThisCM = _pointList[ptidx] - _position;
    				ptOtherCM = _pointList[ptidx] - other.getPosition() - offsetVec; //Test PBCs, to fix PBCs adequately for bonds.
    				
    				// ptOtherLset(0) =  ptOtherCM(0)*cos2 + ptOtherCM(1)*sin2;
    				// ptOtherLset(1) = -ptOtherCM(0)*sin2 + ptOtherCM(1)*cos2;
    				// ptOtherLset += other.getCmLset();
    				// double p2pdist = 0; // not used
    				
    				//normal = (other.getPointList()[ptidxf] - _pointList[ptidx])/((other.getPointList()[ptidxf] - _pointList[ptidx]).norm()+DBL_MIN);  //NEW: Redundant BB
    				
    				dist = _pointList[ptidx] - other.getPointList()[ptidxf] - offsetVec; //Single bond proposal for surface BB  //Test PBCs, to fix PBCs adequately for bonds.
    				//dist = _pointList[ptidx] - other.getPosition(); //Multiple bond default method AA
		        	
		        	distmag = dist.norm();
		        	//normal = dist / distmag;
		        	normal = -dist / distmag;  //Warning: check if Test PBCs affects normal from dist  //EXTRA WARNING: CHECK THE EFFECT OF THIS NEG. NORMAL GIVEN SIGN CONVENTION UP
		        	 
		        	sigma = distmag;   //Single bond proposal for surface BB
		        	//sigma = distmag - other.getRadius(); //Multiple bond default method AA
		        	
                    //cout << "Sigma: " << sigma << " cohesiveDist: " << _cohesiveDistance << endl;
                    //relax to check
                    double cohM = 1.0; //Default 1.0 //Larger, more lenient with separation
			        if (sigma <= _cohesiveDistance * cohM) {  //Change function to getNormal, use isCohesion only for bond formation, analogous using radial distance, normal before
    					// relative v of this wrt the other
    					relativeVel << -other.getVelocity()(0) + other.getOmega()*ptOtherCM(1) + (_velocity(0) - _omega*ptThisCM(1)),
    				               	   -other.getVelocity()(1) - other.getOmega()*ptOtherCM(0) + (_velocity(1) + _omega*ptThisCM(0));
    					//cout << "relativeVel: " << relativeVel(0) << " " << relativeVel(1) << endl;
    	               	// normal and shear displacement increments
    					Vector2d sdot_bond = (relativeVel - relativeVel.dot(normal)*normal);
    					Vector2d ndot_bond = normal.dot(-relativeVel)*normal;
    					Vector2d ds_bond = sdot_bond*dt;
    					Vector2d dn_bond = ndot_bond*dt;
    
    					// normal and shear rotation increments
    					double w = _omega - other.getOmega();
    					double sOmega = w;
    					double sTheta = sOmega*dt;
    
    					//rotate shear components
    					double k = -normal(0)*_bondNormals[bondidx](1) + normal(1)*_bondNormals[bondidx](0); // axis of rotation from n_old to n_cur
    					double sint = abs(k); // compute the sin and cos of the rotation //double sint = abs(k); or Magnitude
    					double cost = normal.dot(_bondNormals[bondidx]);
    					k = k/(sint+DBL_MIN); // normalize k to unit magnitude, also don't divide by 0  //SIgned
    					_bondForceShear[bondidx] = _bondForceShear[bondidx]*cost + Vector2d(-_bondForceShear[bondidx](1), _bondForceShear[bondidx](0))*k*sint;
    					_bondForceNormal[bondidx] = _bondForceNormal[bondidx]*cost + Vector2d(-_bondForceNormal[bondidx](1), _bondForceNormal[bondidx](0))*k*sint;
    					_bondThisMomentNormal[bondidx] = _bondThisMomentNormal[bondidx]*cost + k*k*_bondThisMomentNormal[bondidx]*(1.-cost);
    					_bondNormals[bondidx] = normal;
    
    					// bond forces and moments
    					_bondForceNormal[bondidx] += _kn_bond*dn_bond*_bondAreas;
    					_bondForceShear[bondidx] += -_ks_bond*ds_bond*_bondAreas;
    
    					// Check for yield
    					_fSig[bondidx] = fabs(normal.dot(_bondForceNormal[bondidx])/_bondAreas); // magnitute of compressive or tensile stress
    					_fTau[bondidx] = _bondForceShear[bondidx].norm()/_bondAreas;
    					
    					// apply your selected failure criteria. For example:
    					//if (_fSig[bondidx] > _sigC || _fTau[bondidx] > _tauC) {
    					if (_fSig[bondidx] > _sigC * _sigrF[bondidx] || _fTau[bondidx] > _tauC * _taurF[bondidx]) {
    					    cout << "Exceeding bond strength -> fSig[bondidx]: " << _fSig[bondidx] << " sigC: " << _sigC * _sigrF[bondidx] << " fTau[bondidx]: " << _fTau[bondidx] << " tauC: " << _tauC * _taurF[bondidx] << endl;
    					    if (   ( _fSig[bondidx] / (_sigC * _sigrF[bondidx]) ) >  ( _fTau[bondidx]/(_tauC * _taurF[bondidx]) )  ){
    					        cout << "Tension Failure" << endl;
    					        ntension++;
    					    }
    					    else{
    					        cout << "Shear Failure" << endl;
    					        nshear++;
    					    }
    						// bonds broken
    						_bondInfo[bondidx][0] = 0;
    						_bondForceNormal[bondidx] = Vector2d(0,0);
    						_bondForceShear[bondidx] = Vector2d(0,0);
    						_bondThisMomentNormal[bondidx] = 0;
    						_bondThisMomentShear[bondidx] = 0;
    						_fSig[bondidx] = 0;
    						_fTau[bondidx] = 0;
    						_bondNormals[bondidx] = Vector2d(0,0);
    					}else{
    						df = _bondForceNormal[bondidx] + _bondForceShear[bondidx];
    						df += ndot_bond*2.*sqrt(_kn_bond*_bondAreas*(_mass*other.getMass()/(_mass+other.getMass()))) - sdot_bond*2.*sqrt(_ks_bond*_bondAreas*(_mass*other.getMass()/(_mass+other.getMass()))); // Damping
    						force += df;
    						
    						thisMoment += _bondThisMomentNormal[bondidx] + (ptThisCM(0)*df(1) - ptThisCM(1)*df(0));//moment of beam + from beam  //Check effect of removing moment for TEST sheet
    						otherMoment += -_bondThisMomentNormal[bondidx] + (-ptOtherCM(0)*df(1) + ptOtherCM(1)*df(0));                         //Check effect of removing moment for TEST sheet
    
    						// double nMOIThis = (_momentInertia * normal).dot(normal) + _mass*(ptThisCM(0)*normal(1) - ptThisCM(1)*normal(0)).squaredNorm();  //For off-plane moment ?
    						// double nMOIOther = (other.getMomentInertia() * -normal).dot(-normal) + _mass*(-ptOtherCM(0)*normal(1) + ptOtherCM(1)*normal(0)).squaredNorm(); //For twist moment ?
    
    						//Additional damping that might be redundant.
    						// if (sOmega.norm() > 0) {
    						// 	Vector2d sheardir = sOmega/sOmega.norm(); //SIgn
    						// 	double shearMOIthis = (_momentInertia * sheardir).dot(sheardir) + _mass*(ptThisCM(0)*sheardir(1) - ptThisCM(1)*sheardir(0)).squaredNorm();
    						// 	double shearMOIother = (other.getMomentInertia() * -sheardir).dot(-sheardir) + _mass*(-ptOtherCM(0)*sheardir(1) + ptOtherCM(1)*sheardir(0)).squaredNorm();
    
    						// 	thisMoment  -= sOmega*2.*sqrt(_kn_bond*_bondMomentsInertia*shearMOIthis);
    						// 	otherMoment += sOmega*2.*sqrt(_kn_bond*_bondMomentsInertia*shearMOIother);
    						// }
    						
    						//Control
    			            //if ( 1 > 0  ){
    			            if ( isnan(force(0)) || isnan(force(1))  || isnan(thisMoment)  || isnan(otherMoment)  ){
                			    cout << "BG Force Contact: " << force(0) << " " << force(1) << endl;
                			    cout << "BG Moments Contact: " << thisMoment << " " << otherMoment<< endl;
                			    cout << "normal: " << normal(0) << " " << normal(1) << endl;
                			    cout << "ndot / sdot: " << ndot << " " << sdot << endl;
    			            }
    						
    						
    					}
    				}else{
    					// When in this loop, it means that bonds have not failed yet, but this surface node
    					// has fallen out of the other's level set, or the cohesive distance, or has penetrated
    					// *this. Depend on different applications, you can consider different approaches.
    					// Maybe re-adjust failure parameter to make sure bonds break before this happens.
    					// Or, you can simply set the bonds to be broken like the following
    					
    // 					//Just trying
    // 					// relative v of this wrt the other
    // 					relativeVel << -other.getVelocity()(0) + other.getOmega()*ptOtherCM(1) + (_velocity(0) - _omega*ptThisCM(1)),
    // 				               	   -other.getVelocity()(1) - other.getOmega()*ptOtherCM(0) + (_velocity(1) + _omega*ptThisCM(0));
    					
    // 	               	// normal and shear displacement increments
    // 					Vector2d sdot_bond = (relativeVel - relativeVel.dot(normal)*normal);
    // 					Vector2d ndot_bond = normal.dot(-relativeVel)*normal;
    // 					Vector2d ds_bond = sdot_bond*dt;
    // 					Vector2d dn_bond = ndot_bond*dt;
    
    // 					// normal and shear rotation increments
    // 					double w = _omega - other.getOmega();
    // 					double sOmega = w;
    // 					double sTheta = sOmega*dt;
    
    // 					//rotate shear components
    // 					double k = -normal(0)*_bondNormals[bondidx](1) + normal(1)*_bondNormals[bondidx](0); // axis of rotation from n_old to n_cur
    // 					double sint = abs(k); // compute the sin and cos of the rotation //double sint = abs(k); or Magnitude
    // 					double cost = normal.dot(_bondNormals[bondidx]);
    // 					k = k/(sint+DBL_MIN); // normalize k to unit magnitude, also don't divide by 0  //SIgned
    // 					_bondForceShear[bondidx] = _bondForceShear[bondidx]*cost + Vector2d(-_bondForceShear[bondidx](1), _bondForceShear[bondidx](0))*k*sint;
    // 					_bondForceNormal[bondidx] = _bondForceNormal[bondidx]*cost + Vector2d(-_bondForceNormal[bondidx](1), _bondForceNormal[bondidx](0))*k*sint;
    // 					_bondThisMomentNormal[bondidx] = _bondThisMomentNormal[bondidx]*cost + k*k*_bondThisMomentNormal[bondidx]*(1.-cost);
    // 					_bondNormals[bondidx] = normal;
    
    // 					// bond forces and moments
    // 					_bondForceNormal[bondidx] += _kn_bond*dn_bond*_bondAreas;
    // 					_bondForceShear[bondidx] += -_ks_bond*ds_bond*_bondAreas;
    					
    					
				// 	    df = _bondForceNormal[bondidx] + _bondForceShear[bondidx];
				// 		df += ndot_bond*2.*sqrt(_kn_bond*_bondAreas*(_mass*other.getMass()/(_mass+other.getMass()))) - sdot_bond*2.*sqrt(_ks_bond*_bondAreas*(_mass*other.getMass()/(_mass+other.getMass()))); // Damping
				// 		force += df;
				// 		thisMoment += _bondThisMomentNormal[bondidx] + (ptThisCM(0)*df(1) - ptThisCM(1)*df(0));//moment of beam + from beam
				// 		otherMoment += -_bondThisMomentNormal[bondidx] + (-ptOtherCM(0)*df(1) + ptOtherCM(1)*df(0));
    // 					//End trying
    					bool sep_fail = false; //Turn off to remove separation failure (this needs stronger higher kN though)
    					
    					
    				// 	//Control
    				    if (sep_fail){
        					cout << "Bond broken due to falling out of level set, separation distance: " << sigma << " cohesive distance: " << _cohesiveDistance << endl;
        					//cout << "Bond broken due to falling out of level set, separation distance: " << sigma << " cohesive distance: " << _cohesiveDistance << endl;
        				// 	cout << "Bond forces for reference: " << endl;
    				    //     cout << "BG Force Contact: " << force(0) << " " << force(1) << endl;
            // 			    cout << "BG Moments Contact: " << thisMoment << " " << otherMoment<< endl;
            // 			    cout << "normal: " << normal(0) << " " << normal(1) << endl;
            // 			    cout << "ndot / sdot: " << ndot << " " << sdot << endl;
            			    
            			    //Turn back on
        					_bondInfo[bondidx][0] = 0;
        					_bondForceNormal[bondidx] = Vector2d(0,0);
        					_bondForceShear[bondidx] = Vector2d(0,0);
        					_bondThisMomentNormal[bondidx] = 0;
        					_bondThisMomentShear[bondidx] = 0;
        					_fSig[bondidx] = 0;
        					_fTau[bondidx] = 0;
        					_bondNormals[bondidx] = Vector2d(0,0);
    				    }
    				}
    			}
    		} // compute for all bond reaction force	
        }
        
        
		return isContact;
	} // end findInterparticleForceMomentDEM

	size_t findInterparticleForceforFrac(const Grain2d & other, double & force){
		const double cos2 = cos(other.getTheta());
		const double sin2 = sin(other.getTheta());

		Vector2d		ptOtherCM; 		// point wrt the center of mass of other in real space
		Vector2d		ptOtherLset; 	// point in the reference config of other's level set
		double		penetration;	// penetration amount (is negative by convention of the level set)
		Vector2d		normal; 			// surface normal
		Vector2d		df; 				// force increment from a single point
		Vector2d		tangent; 		// surface tangent
		Vector2d		forcevector;    // force as a vector
		size_t		contactpoint=1;
		force = 0;			// interparticle force


		for (size_t p = 0; p < _pointList.size(); p++){
			double forcei;
			ptOtherCM = _pointList[p] - other.getPosition();
			ptOtherLset(0) =  ptOtherCM(0)*cos2 + ptOtherCM(1)*sin2;
			ptOtherLset(1) = -ptOtherCM(0)*sin2 + ptOtherCM(1)*cos2;
			ptOtherLset += other.getCmLset();

			if (other.getLset().isPenetration(ptOtherLset, penetration, normal)){;
				// rotate the normal from the reference config of 2's level set to real space
				normal << normal(0)*cos2 - normal(1)*sin2, normal(0)*sin2 + normal(1)*cos2;
				// find the tangent
				tangent << -normal(1), normal(0);
				df = penetration*normal*_kn;
				forcevector = -df - tangent*_nodeShears[p];
				forcei = forcevector.norm();
				if (forcei > force) {
					force = forcei;
					contactpoint = p;
				}
			}
		}
		return contactpoint;
	}



	//Start Fluid Modification (Interaction)
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------


    //Simpler function
		//Constant water and air velocities        //Uses only a simple drag force, no skin drag moment assuming center of mass application only
    void fluidInteractionSimple (const double & Cha, const double & Cva, const double & Chw, const double & Cvw,
                           const double & rhoice, const double & rhoair, const double & rhowater,
                           const double & hice, double & hair, double & hwater, Vector2d & Ua, Vector2d & Uw,
                           Vector2d & fluidForceha, Vector2d & fluidForcehw, Vector2d & fluidForceh,
                           Vector2d & fluidForceva, Vector2d & fluidForcevw, Vector2d & fluidForcev,
                           Vector2d & ppv, Vector2d & ppvn, Vector2d & midp, Vector2d & fluidForce, double & fluidMoment, size_t & tstep, double & slope_dir,  double & flowangle,  double & flowforce, const Vector2d & offset) 
    {

		//USEFUL IN ABSOLUTE COORDS
		Vector2d floeCentroid = findCentroid(_pointList); //Use points for area

		//VORTEX-like field  
		// double field_vel = 2.0e12;  //8.0e10 * 0.0 //8.0e1; //2.0e12
		// //double field_vel = 8.0e11; //2.0e12
		// //double field_vel = 1.0 * 0.00006; //8 //300 //2500.0;  //8.0 //400 bk  //TODO: Adjust to scale
		// double limit = 2000;
		// if (tstep > 0)
		// {
		// 	if (_position(1)<limit && _position(0)<limit)
		// 	{
		// 	  Uw << -field_vel , -field_vel  ; 
		// 	  //Uw << field_vel , field_vel  ;  //CONVERGE
		// 	    //Uw << field_vel , 0. ;   //7800 stable w/o 1D heat, 5000 stable for 1D heat (depending on grain morphology) , 2500 +/- under melting
				
		// 		//Uw << -900.,900. ;
		// 	}
		// 	else
		// 	{    if(_position(1)>=limit && _position(0)<limit)
		// 	    {
		// 	     Uw << -field_vel , -field_vel  ; 
		// 	     //Uw << field_vel , -field_vel ;  //CONVERGE
		// 	        //Uw << 0.,-field_vel ;

		// 	        //Uw << -900.,900. ;
		// 	    }
		// 	    else
		// 	    {   if(_position(1)>=limit && _position(0)>limit)
		// 	        {
		// 	       Uw << -field_vel , -field_vel  ; 
		// 	       //Uw << -field_vel , -field_vel ;  //CONVERGE
		// 	          //Uw << -field_vel,0. ;
			         
		// 	          //Uw << -900.,900. ;
		// 	        }
		// 	        else
		// 	        {
		// 	          Uw << -field_vel , -field_vel  ; 
		// 	          //Uw << -field_vel , field_vel;  //CONVERGE
		// 	            //Uw << 0.,field_vel ;
			         
		// 	            //Uw << -900.,900. ;
		// 	        }
		// 	    }
		// 	}
		// }	
        
  //       if (slope_dir > 0)
  //       {
  //           Uw << field_vel , field_vel  ; 
  //           //Uw << field_vel , -field_vel  ; 
  //       }
  //       else if (slope_dir < 0)
  //       {
  //           Uw << -field_vel , -field_vel  ; 
  //       }
  //       else
  //       {
  //          cout << "WARNING SLOPE ZERO, INSPECT!!!" << endl;
  //       }
  
//         //1. Convergent field for debugging (change for other fields) (NEW MODIF)
// 		if (_position(1)<offset(1)*0.5 && _position(0)<offset(0)*0.5)
// 		{
// 		  Uw << flowforce, flowforce ;
// 		}
// 		else
// 		{    if(_position(1)>offset(1)*0.5 && _position(0)<offset(0)*0.5)
// 		    {
// 		        Uw << flowforce, -flowforce ;
// 		    }
// 		    else
// 		    {   if(_position(1)>offset(1)*0.5 && _position(0)>offset(0)*0.5)
// 		        {
// 		          Uw << -flowforce, -flowforce ;
// 		        }
// 		        else
// 		        {
// 		          Uw << -flowforce, flowforce ;
// 		        }
// 		    }
// 		}
		
		//2. Vortex field with cycle change (cyclone, anticyclone) (NEW MODIF)
// 		if (_position(1)<offset(1)*0.5 && _position(0)<offset(0)*0.5) //DL
// 		{   
// 		    if (slope_dir >= 0){
// 		        Uw << flowforce, 0.0000 ;
// 		    }
// 		    else{
// 		        Uw << 0.0000, flowforce;
// 		    }
// 		}    
// 		else
// 		{   if(_position(1)>offset(1)*0.5 && _position(0)<offset(0)*0.5) //UL
// 		    {
//     		    if (slope_dir >= 0){
// 	    	        Uw << 0.000, -flowforce ;
// 		        }
// 		        else{
// 		            Uw << flowforce, 0.0000 ;
// 		        }
// 		    }
// 		    else
// 		    {   if(_position(1)>offset(1)*0.5 && _position(0)>offset(0)*0.5) //UR
// 		        {
// 		          if (slope_dir >= 0){
// 		            Uw << -flowforce, 0.0000 ;
// 		          }
// 		          else{
// 		            Uw << 0.000, -flowforce ;  
// 		          }
// 		        }
// 		        else  //DR
// 		        {
// 		          if (slope_dir >= 0){
// 		              Uw << 0.000, flowforce ;
// 		          }
// 		          else{
// 		              Uw << -flowforce, 0.0000 ; 
// 		          }
// 		        }
// 		    }
// 		}

		//3. Random direction field  (all move the same using data of magnitude and angle) (NEW MODIF)

        double field_force = flowforce; //4.0e11  //8.0e10 * 0.0 //8.0e1; //2.0e12    //-- 8.0e11; //2.0e12   (NEW MODIF)
        ////double field_vel = 1.0 * 0.00006; //8 //300 //2500.0;  //8.0 //400 bk  //TODO: Adjust to scale
        Uw << field_force * cos(flowangle) , field_force * sin(flowangle) ;   //(NEW MODIF)

        //No surf. area adjust
        fluidForce = Uw;
        
        //Aug 18, 2022 Change
        //BEGIN
        //Adjust for thickness (less force)
        fluidForce *= 1.0; //Aug 22, 2022 //Due to adjust on ice density. Check if it must be reduced. More massive ice needs a stronger force to accelarated sufficiently. Less massive needs a weaker force for stablility.
        //fluidForce *= 0.0;  //JUST FOR NOW!!!!
        //fluidForce *= 0.001;
        //END
        
		               
		//Generate a more arbitrary function
		//        Ua << 3000*_position(1), 3000*(-_position(0)-_position(1));  //What type of randomization be given to wind
		//        Ua << 3000*sqrt(pow((_position(0)-500),2) + pow((_position(1)-500),2)), 3000*sqrt(pow((_position(0)-500),2) + pow((_position(1)-500),2));  //What type of randomization be given to wind
		//        Uw << _position(0)*2 + _position(1)*2 , _position(0)*2 + _position(1)*2;
		//        Ua << _position(0)*2 + _position(1)*2 , _position(0)*2 + _position(1)*2;  //What type of randomization be given to wind
		//        Uw << _position(0)*2 + _position(1)*2 , _position(0)*2 + _position(1)*2;

		double adjust_factor = 1; //100 //1000000
		double moment_adjust_factor = 0.00000; //1000000
		//SKIN DRAG SIMPLE
		fluidForcehw = rhowater*Chw*((_mass*0.91/adjust_factor)/_density)*(Uw-_velocity).norm()*(Uw-_velocity);		
		//fluidForceha = rhoair*Cha*((floeArea*_thick)/_density)*Ua.norm()*Ua; //use mass0/density to get volumne or in this case to get 2D area
		//fluidForcehw = rhowater*Chw*((floeArea*_thick)/_density)*(Uw-_velocity).norm()*(Uw-_velocity);  //???????????
		fluidForceha = rhoair*Cha*((_mass*0.91/adjust_factor)/_density)*Ua.norm()*Ua; //use mass0/density to get volumne or in this case to get 2D area
		//fluidForcehw = rhowater*Chw*((_mass*0.91/1000000)/_density)*(Uw-_velocity).norm()*(Uw-_velocity);  //???????????
		fluidForceh = fluidForceha + fluidForcehw;
        //cout << "Fluid Moment after Skin: " << fluidMoment << endl;
		//A constant wind and ocean velocity will not induce moment due to skin drag if applied at center. But using cells that can change.



		//FORM DRAG
		//Input all points of grain to get each normal vector and compare each one to the Ua or Uw and get total Perimeter vertically pressured

		//Adjusted height based on area proportionality
		//hice=function of ( (_mass/_density) ) * hice; But applied on hair and hwater

		hair = _thick*(rhowater-rhoice)/rhowater;   //Let _thick replace hice as a grain property
		hwater = _thick-hair;

		//hair=hice*(rhowater-rhoice)/rhowater;   //Let _thick replace hice as a grain property
		//hwater=hice-hair;

		for (size_t ptidx = 0; ptidx < _pointList.size(); ptidx++) {
		  
			if (ptidx == _pointList.size()-1)
			{
			  ppv=_pointList[ptidx]-_pointList[0];
			  midp=ppv*0.5+_pointList[0];  //Define midpoint for better calculation of moments
			}
			else
			{
			  ppv=_pointList[ptidx]-_pointList[ptidx+1];  //Dist between 2 points
			  midp=ppv*0.5+_pointList[ptidx+1];
			}
			   		 
			ppvn << ppv(1), -1*ppv(0);  //Rotation of 90 degrees for a ccw ppvector s.t. normal faces outwardm, not normal, has a norm ~= 1
			Vector2d temporal_v = Vector2d(0.0 , 0.0);
			
			if(ppvn.dot(Ua) < 0)       //if(ppvn(0)*Ua(0)+ppvn(1)*Ua(1) < 0)
			{
			  temporal_v = rhoair*Cva*hair*abs(ppvn.dot(Ua))*Ua;
			  fluidForceva += temporal_v;   //Accumulation of forces if Ua crashes against face, otherwise, no body drag occurs, it occurs on the projection
			  //fluidForceva+=rhoair*Cva*hair*ppv.norm()*Ua.norm()*Ua;   //Accumulation of forces if Ua crashes against face, otherwise, no body drag occurs
			  //fluidMoment += fluidForceva(0)*(midp(0)-_cmLset(0))+fluidForceva(1)*(midp(1)-_cmLset(1));  //Is this CM of Mass Shifted as we go
			 
			  //fluidMoment += - temporal_v(0)*(midp(0)-floeCentroid(0)) + temporal_v(1)*(midp(1)-floeCentroid(1));  //?? +
			}
			   
			//cout << "Form Drag CHECK" << endl;
			//cout << "ppv: " << ppv(0) << " , " << ppv(1) << endl;
			//cout << "ppvn: " << ppvn(0) << " , " << ppvn(1) << endl;
			//cout << "Rel Vel: " << Uw(0)-_velocity(0) << " , " << Uw(1)-_velocity(1) << endl;

	        //Update to defined grid

            //cout << "MIDP: " << midp(0) << " , " << midp(1) << endl; 
	        //Find Uw at midp
	        //cout << "Form Drag Interpolation" << endl;
			//cout << "Midp: " << midp(0) << "  " << midp(1) << endl;
	        Vector2d midUwg = Uw; //Interpolate at point using grid info
	        //cout << "MIDUWG: " << midUwg(0) << " , " << midUwg(1) << endl;
	        //Vector2d midUwg = Uw;

	        if(ppvn.dot(midUwg - _velocity) < 0) 
	        { 
	        	//cout << "Water Form Drag" << endl;
	        	//cout << "Form Drag Water" << endl;
	        	temporal_v = rhowater * Cvw * hwater * abs(ppvn.dot(midUwg-_velocity))*(midUwg-_velocity);  //ppvn norm scales perimeter of floe going against current
	        	fluidForcevw += temporal_v;
	        	//fluidMoment += fluidForcevw(0)*(midp(0)-_cmLset(0)) + fluidForcevw(1)*(midp(1)-_cmLset(1));  //More works
	        	
	        	//fluidMoment += - temporal_v(0)*(midp(0)-floeCentroid(0)) + temporal_v(1)*(midp(1)-floeCentroid(1));  //?? +
	        }	

			// if(ppvn.dot(Uw-_velocity) < 0) // if(ppvn(0)*(Uw(0)-_velocity(0))+ppvn(1)*(Uw(1)-_velocity(1)) < 0)
			// {
			//   cout << "Water Form Drag" << endl;
			//   fluidForcevw+=rhowater*Cvw*hwater*abs(ppvn.dot(Uw-_velocity))*(Uw-_velocity);
			//   //fluidForcevw+=rhowater*Cvw*hwater*ppv.norm()*(Uw-_velocity).norm()*(Uw-_velocity);
			//   fluidMoment+=fluidForcevw(0)*(midp(0)-_cmLset(0)) + fluidForcevw(1)*(midp(1)-_cmLset(1)); //Should be + for ccw positive moment
			//   //cout << "X: " << fluidForcevw(0) << " Y: " <<  fluidForcevw(1) << " Moment: " << fluidForcevw(0)*(midp(0)-_cmLset(0))+fluidForcevw(1)*(midp(1)-_cmLset(1)) << endl;

			//   //fluidMoment*=0.0; //Inspect Moment
			// }
		}

		//fluidForceva=rhoair*Cva*hair*Perair*Ua.norm()*Ua;     //use Perimeter*exposed_height to get 2D area
		//fluidForcevw= rhowater*Cvw*hwater*Perwater*(Uw-_velocity).norm()*(Uw-_velocity);
		//fluidForcev = fluidForceva + fluidForcevw;
		//fluidForce = fluidForcev + fluidForceh;

        //Try out
		//fluidForce << 100.0, 100.0;

		fluidMoment = moment_adjust_factor*0.0;

		// cout << "Skin Drag Force"  << endl;
		// cout << "X: " << fluidForceh(0) << " Y: " <<  fluidForceh(1) << endl;
		// cout << "Form Drag Force"  << endl;
		// cout << "X: " << fluidForcev(0) << " Y: " <<  fluidForcev(1) << endl;
		// cout << "Fluid Moment after Form: " << fluidMoment << endl;
		//cout << "Total Force"  << endl;
		//cout << "X: " << fluidForce(0) << " Y: " <<  fluidForce(1) << endl;
       
    }


    //Simpler function
	//Constant water and air velocities        //Uses only a simple drag force, no skin drag moment assuming center of mass application only
    void fluidInteraction_Nova (const double & Cha, const double & Cva, const double & Chw, const double & Cvw,
                           const double & rhoice, const double & rhoair, const double & rhowater,
                           const double & hice, double & hair, double & hwater, Vector2d & Ua, Vector2d & Uw,
                           Vector2d & fluidForceha, Vector2d & fluidForcehw, Vector2d & fluidForceh,
                           Vector2d & fluidForceva, Vector2d & fluidForcevw, Vector2d & fluidForcev,
                           Vector2d & ppv, Vector2d & ppvn, Vector2d & midp, Vector2d & fluidForce, double & fluidMoment, 
                           size_t & tstep, double & slope_dir,  double & flowangle,  double & flowforce, const Vector2d & offset, size_t & cell_sizex, size_t & cell_sizey, vector<Vector2d> & Uwg, const size_t & x_cells, const size_t & y_cells, const vector<Vector2d> & fluid_coord, size_t & index_i) 
    {

		//Initialize forces and moments
		Vector2d force_out;
		force_out(0) = 0.00;
		force_out(1) = 0.00;
		double moment_out = 0.00;
		Vector2d fluid_drag_force;
		
		//Find grain centroid in terms of global coordinates (not local)
		Vector2d floeLocation = _position; //Use points for area
		
		double adjust_fff = 1e-4 * (1/0.000007) * 1; //1e-4 WORKS!!!! with 1/dt   ####1e3 few, 1e7 too much, even 5e3 too much //Initially 1e2 seems okay but still tooooo much , 5e-2 still too high converge, 5e-4 still a bitty too high but improving
		double adjust_mmm = 10.0 * 1; //1000 is excessive rot //100 is bit too much but less //Try 10  //Default 1
		
		//Find overlapping points and get forces and moments for each one (and accumulate)
		
		double floeArea = PointsArea(_pointList);
		bool small_floe = false;
		bool nocontact = true; //For fringe sized floes
        
        //A. If floe is smaller than cell size, add an if and use nearest grid point function round_Ocean to find force proportional to floe area 
		if (floeArea < (cell_sizex * cell_sizey)){
		    small_floe = true;
            //Assign current velocity
            Uw = round_OceanV(_position, Uwg, x_cells, y_cells, offset);
            //SKIN DRAG SIMPLE
            fluid_drag_force = rhowater*Chw*(Uw-_velocity).norm()*(Uw-_velocity) * floeArea;	
            
            force_out(0) += fluid_drag_force(0) * adjust_fff;
            force_out(1) += fluid_drag_force(1) * adjust_fff;
            
            //Based on centroid location
            //double dist_x = pt_compare(0) - floeLocation(0);
            //double dist_y = pt_compare(1) - floeLocation(1);
            //moment_out += dist_x * force_out(1) * adjust_mmm;  //Moment due to y force using x dist.
            //moment_out += -1.0 * dist_y * force_out(0) * adjust_mmm;  //Moment due to x force using y dist.
            moment_out += 0.00; //Too small to add moment
            cout << "CASE A: Veeeery small floe compared to grid!!!" << endl;
		}
		
		//B. Find all drag cell forces for larger floes (usual case)
		if (small_floe == false){
    		for (size_t i = 0; i < y_cells; i++) {
                for (size_t j = 0; j < x_cells; j++) {
                   
                    Vector2d pt_compare = fluid_coord[j+i*x_cells];
                    bool contact_p = under_ice(pt_compare);
    
                    
                    if (contact_p)
                    {
                        //Assign current velocity
                        Uw = Uwg[j+i*x_cells];
                        //SKIN DRAG SIMPLE
    		            fluid_drag_force = rhowater*Chw*(Uw-_velocity).norm()*(Uw-_velocity) * (cell_sizex * cell_sizey);	
                        
                        force_out(0) += fluid_drag_force(0) * adjust_fff;
                        force_out(1) += fluid_drag_force(1) * adjust_fff;
                        
                        //Based on centroid location
                        double dist_x = pt_compare(0) - floeLocation(0);
                        double dist_y = pt_compare(1) - floeLocation(1);
                        moment_out += dist_x * force_out(1) * adjust_mmm;  //Moment due to y force using x dist.
                        moment_out += -1.0 * dist_y * force_out(0) * adjust_mmm;  //Moment due to x force using y dist.
                        nocontact = false;
                    }
                }    
            } 
		}
		
		//C. Extreme case, for fringe floe larger than size resolution but not intersecting any grid points
		if (small_floe == false && nocontact == true){
            //Assign current velocity
            Uw = round_OceanV(_position, Uwg, x_cells, y_cells, offset);
            //SKIN DRAG SIMPLE
            fluid_drag_force = rhowater*Chw*(Uw-_velocity).norm()*(Uw-_velocity) * floeArea;	
            
            force_out(0) += fluid_drag_force(0) * adjust_fff;
            force_out(1) += fluid_drag_force(1) * adjust_fff;
            
            //Based on centroid location
            //double dist_x = pt_compare(0) - floeLocation(0);
            //double dist_y = pt_compare(1) - floeLocation(1);
            //moment_out += dist_x * force_out(1) * adjust_mmm;  //Moment due to y force using x dist.
            //moment_out += -1.0 * dist_y * force_out(0) * adjust_mmm;  //Moment due to x force using y dist.
            moment_out += 0.00; //Still too small to add moment
		}
		
		//Save values for output
		fluidForce = force_out;
		fluidMoment = moment_out;
		
        //Aug 18, 2022 Change
        //BEGIN
        //Adjust for thickness (less force) (this is done when mass is adjusted)
        //double adjust_fff = 1; 
        //fluidForce(0) = fluidForce(0) * adjust_fff; //0.1 too low //Aug 22, 2022 //Due to adjust on ice density. Check if it must be reduced. More massive ice needs a stronger force to accelarated sufficiently. Less massive needs a weaker force for stablility.
        //fluidForce(1) = fluidForce(1) * adjust_fff;
        //fluidMoment = fluidMoment; // 0.00001 too low

		//SKIN DRAG SIMPLE
		//fluidForcehw = rhowater*Chw*((_mass*0.91/adjust_factor)/_density)*(Uw-_velocity).norm()*(Uw-_velocity);		
		//fluidForceha = rhoair*Cha*((floeArea*_thick)/_density)*Ua.norm()*Ua; //use mass0/density to get volumne or in this case to get 2D area
		//fluidForcehw = rhowater*Chw*((floeArea*_thick)/_density)*(Uw-_velocity).norm()*(Uw-_velocity);  //???????????
		//fluidForceha = rhoair*Cha*((_mass*0.91/adjust_factor)/_density)*Ua.norm()*Ua; //use mass0/density to get volumne or in this case to get 2D area
		//fluidForcehw = rhowater*Chw*((_mass*0.91/1000000)/_density)*(Uw-_velocity).norm()*(Uw-_velocity);  //???????????
		//fluidForceh = fluidForceha + fluidForcehw;
		//FORM DRAG //Check other functions
        //Assume negligible for now
		//fluidMoment = moment_adjust_factor*0.0;

        //Verify
        if (index_i == 0){
    // 		cout << "Position" << endl;
    // 		cout << "X: " << _position(0) << " Y: " << _position(1) << endl;
    // 		cout << "Velocity" << endl;
    // 		cout << "X: " << _velocity(0) << " Y: " << _velocity(1) << endl;    		
    // 		cout << "Skin Drag Force"  << endl;
    // 		cout << "Fluid Moment: " << fluidMoment << endl;
    // 		cout << "Total Force"  << endl;
    // 		cout << "X: " << fluidForce(0) << " Y: " <<  fluidForce(1) << endl;
        }
        return;
    }

    //Simpler function for regular
	//Constant water and air velocities        //Uses only a simple drag force, no skin drag moment assuming center of mass application only
    void fluidInteraction_NovaDEM (const double & Cha, const double & Cva, const double & Chw, const double & Cvw,
                           const double & rhoice, const double & rhoair, const double & rhowater,
                           const double & hice, double & hair, double & hwater, Vector2d & Ua, Vector2d & Uw,
                           Vector2d & fluidForceha, Vector2d & fluidForcehw, Vector2d & fluidForceh,
                           Vector2d & fluidForceva, Vector2d & fluidForcevw, Vector2d & fluidForcev,
                           Vector2d & ppv, Vector2d & ppvn, Vector2d & midp, Vector2d & fluidForce, double & fluidMoment, 
                           size_t & tstep, double & slope_dir,  double & flowangle,  double & flowforce, const Vector2d & offset, size_t & cell_sizex, size_t & cell_sizey, vector<Vector2d> & Uwg, const size_t & x_cells, const size_t & y_cells, const vector<Vector2d> & fluid_coord, size_t & index_i, double & dt) 
    {

		//Initialize forces and moments
		Vector2d force_out;
		force_out(0) = 0.00;
		force_out(1) = 0.00;
		double moment_out = 0.00;
		Vector2d fluid_drag_force;
		
		//Find grain centroid in terms of global coordinates (not local)
		Vector2d floeLocation = _position; //Use points for area
		
		//double adjust_fff = 1e-4 * (1/0.000007) * 1; //0.01 too weak for mini alone //1e-4 WORKS!!!! with 1/dt   ####1e3 few, 1e7 too much, even 5e3 too much //Initially 1e2 seems okay but still tooooo much , 5e-2 still too high converge, 5e-4 still a bitty too high but improving
		double adjust_fff = 1e-4 * (1/dt) * 1 * 100; //300 aggressive //100 is a bit slow  //Gotta go fast * 20, add factor of 10 for more aggressive loading //Remove for other cases
		double adjust_mmm = 10.0 * 1; //1000 is excessive rot //100 is bit too much but less //Try 10  //Default 1
		
		//Find overlapping points and get forces and moments for each one (and accumulate)
		double floeArea = 0;
    	size_t n = _pointList.size();

    	for (size_t i = 0; i < n-1; i++)
    	{
    		floeArea += ( _pointList[i](0) * _pointList[i+1](1) -  _pointList[i](1) * _pointList[i+1](0) ); 
    	}
    	floeArea += (_pointList[n-1](0) * _pointList[0](1) -  _pointList[n-1](1) * _pointList[0](0) ); 

    	floeArea = 0.5*abs(floeArea);
        
        //Assume regular DEM is very small (no shape effects)
        //A. If floe is smaller than cell size, add an if and use nearest grid point function round_Ocean to find force proportional to floe area 
		
		if (isnan(floeArea) == false){
    		if (floeArea < (cell_sizex * cell_sizey)){
                //Assign current velocity
                Uw = round_OceanV(_position, Uwg, x_cells, y_cells, offset);
                //SKIN DRAG SIMPLE
                fluid_drag_force = rhowater*Chw*(Uw-_velocity).norm()*(Uw-_velocity) * floeArea;	
                
                force_out(0) += fluid_drag_force(0) * adjust_fff;
                force_out(1) += fluid_drag_force(1) * adjust_fff;
                
                //Based on centroid location
                //double dist_x = pt_compare(0) - floeLocation(0);
                //double dist_y = pt_compare(1) - floeLocation(1);
                //moment_out += dist_x * force_out(1) * adjust_mmm;  //Moment due to y force using x dist.
                //moment_out += -1.0 * dist_y * force_out(0) * adjust_mmm;  //Moment due to x force using y dist.
                moment_out += 0.00; //Too small to add moment //CHECK THIS!!!!
                //cout << "CASE A: Veeeery small floe compared to grid!!!" << endl;
    		}
    		else{
    		    cout << "Error: DO NOT USE NORMAL DEM!!!" << endl;
    		    cout << "Floe area: " << floeArea << " cell size: " << cell_sizex * cell_sizey << endl;
    		    exit(1);
    		}
		}
		else
		{
		    cout << "NAN points, check FORCES!!!!!" << endl;
		    cout << "Position X: " << _position(0) << " Position Y: " << _position(1) << " id: " << index_i << endl;
		}
		
		//Save values for output
		fluidForce = force_out;
		fluidMoment = moment_out;
		
		//Control
		//cout << "Fluid force: " << fluidForce(0) << " " << fluidForce(1) << endl;
		//cout << "Fluid moment: " << fluidMoment << endl;
		
        return;
    }
    
    void fluidInteraction_NovaDEM_Damp_small (const double & Cha, const double & Cva, const double & Chw, const double & Cvw,
                           const double & rhoice, const double & rhoair, const double & rhowater,
                           const double & hice, double & hair, double & hwater, Vector2d & Ua, Vector2d & Uw,
                           Vector2d & fluidForceha, Vector2d & fluidForcehw, Vector2d & fluidForceh,
                           Vector2d & fluidForceva, Vector2d & fluidForcevw, Vector2d & fluidForcev,
                           Vector2d & ppv, Vector2d & ppvn, Vector2d & midp, Vector2d & fluidForce, double & fluidMoment, 
                           size_t & tstep, double & slope_dir,  double & flowangle,  double & flowforce, const Vector2d & offset, size_t & cell_sizex, size_t & cell_sizey, vector<Vector2d> & Uwg, const size_t & x_cells, const size_t & y_cells, const vector<Vector2d> & fluid_coord, size_t & index_i, double & dt, vector<double> & dampMat) 
    {

		//Initialize forces and moments
		Vector2d force_out;
		force_out(0) = 0.00;
		force_out(1) = 0.00;
		double moment_out = 0.00;
		Vector2d fluid_drag_force;
		
		//Find grain centroid in terms of global coordinates (not local)
		Vector2d floeLocation = _position; //Use points for area
		
		//double adjust_fff = 1e-4 * (1/0.000007) * 1; //0.01 too weak for mini alone //1e-4 WORKS!!!! with 1/dt   ####1e3 few, 1e7 too much, even 5e3 too much //Initially 1e2 seems okay but still tooooo much , 5e-2 still too high converge, 5e-4 still a bitty too high but improving
		double adjust_fff = 1e-4 * (1/dt) * 1 * 100; //300 aggressive //100 is a bit slow  //Gotta go fast * 20, add factor of 10 for more aggressive loading //Remove for other cases
		double adjust_mmm = 10.0 * 1; //1000 is excessive rot //100 is bit too much but less //Try 10  //Default 1
		
		//Find overlapping points and get forces and moments for each one (and accumulate)
		double floeArea = 0;
    	size_t n = _pointList.size();

    	for (size_t i = 0; i < n-1; i++)
    	{
    		floeArea += ( _pointList[i](0) * _pointList[i+1](1) -  _pointList[i](1) * _pointList[i+1](0) ); 
    	}
    	floeArea += (_pointList[n-1](0) * _pointList[0](1) -  _pointList[n-1](1) * _pointList[0](0) ); 

    	floeArea = 0.5*abs(floeArea);
        
        //Assume regular DEM is very small (no shape effects)
        //A. If floe is smaller than cell size, add an if and use nearest grid point function round_Ocean to find force proportional to floe area 
		double dampG;
		if (isnan(floeArea) == false){
    		if (floeArea < (cell_sizex * cell_sizey)){
                //Assign current velocity
                Uw = round_OceanV(_position, Uwg, x_cells, y_cells, offset);
                dampG = round_Ocean(_position, dampMat, x_cells, y_cells, offset);
                //SKIN DRAG SIMPLE
                fluid_drag_force = rhowater*Chw*(Uw-_velocity).norm()*(Uw-_velocity) * floeArea;	
                
                force_out(0) += fluid_drag_force(0) * adjust_fff * dampG;
                force_out(1) += fluid_drag_force(1) * adjust_fff * dampG;
                
                //Based on centroid location
                //double dist_x = pt_compare(0) - floeLocation(0);
                //double dist_y = pt_compare(1) - floeLocation(1);
                //moment_out += dist_x * force_out(1) * adjust_mmm;  //Moment due to y force using x dist.
                //moment_out += -1.0 * dist_y * force_out(0) * adjust_mmm;  //Moment due to x force using y dist.
                moment_out += 0.00; //Too small to add moment //CHECK THIS!!!!
                //cout << "CASE A: Veeeery small floe compared to grid!!!" << endl;
    		}
    		else{
    		    cout << "Error: DO NOT USE NORMAL DEM!!!" << endl;
    		    cout << "Floe area: " << floeArea << " cell size: " << cell_sizex * cell_sizey << endl;
    		    exit(1);
    		}
		}
		else
		{
		    cout << "NAN points, check FORCES!!!!!" << endl;
		    cout << "Position X: " << _position(0) << " Position Y: " << _position(1) << " id: " << index_i << endl;
		}
		
		//Save values for output
		fluidForce = force_out;
		fluidMoment = moment_out;
		
		//Control
		//cout << "Fluid force: " << fluidForce(0) << " " << fluidForce(1) << endl;
		//cout << "Fluid moment: " << fluidMoment << endl;
		
        return;
    }
    
    //Constant water and air velocities        //Uses only a simple drag force, no skin drag moment assuming center of mass application only
    void fluidInteraction_NovaDEM_Damp_big (const double & Cha, const double & Cva, const double & Chw, const double & Cvw,
                           const double & rhoice, const double & rhoair, const double & rhowater,
                           const double & hice, double & hair, double & hwater, Vector2d & Ua, Vector2d & Uw,
                           Vector2d & fluidForceha, Vector2d & fluidForcehw, Vector2d & fluidForceh,
                           Vector2d & fluidForceva, Vector2d & fluidForcevw, Vector2d & fluidForcev,
                           Vector2d & ppv, Vector2d & ppvn, Vector2d & midp, Vector2d & fluidForce, double & fluidMoment, 
                           size_t & tstep, double & slope_dir,  double & flowangle,  double & flowforce, const Vector2d & offset, size_t & cell_sizex, size_t & cell_sizey, vector<Vector2d> & Uwg, const size_t & x_cells, const size_t & y_cells, const vector<Vector2d> & fluid_coord, size_t & index_i, double & dt, vector<double> & dampMat, double & curr_factor) 
    {

		//Initialize forces and moments
		Vector2d force_out;
		force_out(0) = 0.00;
		force_out(1) = 0.00;
		double moment_out = 0.00;
		Vector2d fluid_drag_force;
		double dampG;
		
		//Find grain centroid in terms of global coordinates (not local)
		Vector2d floeLocation = _position; //Use points for area
		
		//double adjust_fff = 1e-4 * (1/dt) * 1 * 100 * 1.0 * 0.1; //1007 0.01 //1006 0.1 //1008 0.001   //NEW: //0.001 small RT, larger RT try 0.0001 for 10x step is to adjust strength of velocity drag in real time   1/3 is to make arbitrary similar to idealized periodic currents //300 aggressive //100 is a bit slow  //Gotta go fast * 20, add factor of 10 for more aggressive loading //Remove for other cases
		//double current_factor = 1.0; //Default 100
		double current_factor = curr_factor; //Default 100
		double adjust_fff = 1e-4 * (1/dt) * 1.0 * current_factor * 1.0; //100
		double adjust_mmm = 10.0 * 1; //1000 is excessive rot //100 is bit too much but less //Try 10  //Default 1    
		
		//Find overlapping points and get forces and moments for each one (and accumulate)
		double floeArea = 0;
    	size_t n = _pointList.size();

    	for (size_t i = 0; i < n-1; i++)
    	{
    		floeArea += ( _pointList[i](0) * _pointList[i+1](1) -  _pointList[i](1) * _pointList[i+1](0) ); 
    	}
    	floeArea += (_pointList[n-1](0) * _pointList[0](1) -  _pointList[n-1](1) * _pointList[0](0) ); 

    	floeArea = 0.5*abs(floeArea);
		
		
		bool small_floe = false;
		bool nocontact = true; //For fringe sized floes
		
		if (isnan(floeArea)){
		    cout << "ERROR: check position update, NAN area, check forces!!!" << endl;
		}
        
        //A. If floe is smaller than cell size, add an if and use nearest grid point function round_Ocean to find force proportional to floe area 
		if (floeArea < (cell_sizex * cell_sizey)){
		    small_floe = true;
            //Assign current velocity
            Uw = round_OceanV(_position, Uwg, x_cells, y_cells, offset);
            dampG = round_Ocean(_position, dampMat, x_cells, y_cells, offset);
            //SKIN DRAG SIMPLE
            fluid_drag_force = rhowater*Chw*(Uw-_velocity).norm()*(Uw-_velocity) * floeArea;	
            
            force_out(0) += fluid_drag_force(0) * adjust_fff * dampG;
            force_out(1) += fluid_drag_force(1) * adjust_fff * dampG;
            
            //Based on centroid location
            //double dist_x = pt_compare(0) - floeLocation(0);
            //double dist_y = pt_compare(1) - floeLocation(1);
            //moment_out += dist_x * force_out(1) * adjust_mmm;  //Moment due to y force using x dist.
            //moment_out += -1.0 * dist_y * force_out(0) * adjust_mmm;  //Moment due to x force using y dist.
            moment_out += 0.00; //Too small to add moment
            //cout << "CASE A: Veeeery small floe compared to grid!!!" << endl;
		}
		
		//B. Find all drag cell forces for larger floes (usual case)
		if (small_floe == false){
    		for (size_t i = 0; i < y_cells; i++) {
                for (size_t j = 0; j < x_cells; j++) {
                   
                    Vector2d pt_compare = fluid_coord[j+i*x_cells];
                    bool contact_p = under_iceDEM(pt_compare);
    
                    
                    if (contact_p)
                    {
                        //Assign current velocity
                        Uw = Uwg[j+i*x_cells];
                        dampG = dampMat[j+i*x_cells];
                        //SKIN DRAG SIMPLE
    		            fluid_drag_force = rhowater*Chw*(Uw-_velocity).norm()*(Uw-_velocity) * (cell_sizex * cell_sizey);	
                        
                        force_out(0) += fluid_drag_force(0) * adjust_fff * dampG;
                        force_out(1) += fluid_drag_force(1) * adjust_fff * dampG;
                        
                        //Based on centroid location
                        double dist_x = pt_compare(0) - floeLocation(0);
                        double dist_y = pt_compare(1) - floeLocation(1);
                        moment_out += dist_x * force_out(1) * adjust_mmm;  //Moment due to y force using x dist.
                        moment_out += -1.0 * dist_y * force_out(0) * adjust_mmm;  //Moment due to x force using y dist.
                        nocontact = false;
                    }
                }    
            } 
		}
		
		//C. Extreme case, for fringe floe larger than size resolution but not intersecting any grid points
		if (small_floe == false && nocontact == true){
            //Assign current velocity
            Uw = round_OceanV(_position, Uwg, x_cells, y_cells, offset);
            dampG = round_Ocean(_position, dampMat, x_cells, y_cells, offset);
            //SKIN DRAG SIMPLE
            fluid_drag_force = rhowater*Chw*(Uw-_velocity).norm()*(Uw-_velocity) * floeArea;	
            
            force_out(0) += fluid_drag_force(0) * adjust_fff * dampG;
            force_out(1) += fluid_drag_force(1) * adjust_fff * dampG;
            
            moment_out += 0.00; //Still too small to add moment
		}
		
		//Save values for output
		fluidForce = force_out;
		fluidMoment = moment_out;
		
        return;
    }

    
    //For loading BCs
    //Simpler function for regular
	//Constant water and air velocities        //Uses only a simple drag force, no skin drag moment assuming center of mass application only
    void fluidInteractionLoad_NovaDEM (const double & Cha, const double & Cva, const double & Chw, const double & Cvw,
                           const double & rhoice, const double & rhoair, const double & rhowater,
                           const double & hice, double & hair, double & hwater, Vector2d & Ua, Vector2d & Uw,
                           Vector2d & fluidForceha, Vector2d & fluidForcehw, Vector2d & fluidForceh,
                           Vector2d & fluidForceva, Vector2d & fluidForcevw, Vector2d & fluidForcev,
                           Vector2d & ppv, Vector2d & ppvn, Vector2d & midp, Vector2d & fluidForce, double & fluidMoment, 
                           size_t & tstep, double & slope_dir,  double & flowangle,  double & flowforce, const Vector2d & offset, size_t & cell_sizex, size_t & cell_sizey, vector<Vector2d> & Uwg, const size_t & x_cells, const size_t & y_cells, const vector<Vector2d> & fluid_coord, size_t & index_i) 
    {
        
        bool force_based = false; //True use an escalating force over time, False: use a constant velocity if the grain has to be loaded
        Vector2d tryvel;
        
        double nTref = 4233600.0; //Scale based on nT
        double dTref = 0.000007;  //Time step for dt in kinematics and contacts
        double desired_disp = 1e-2; //Desired displacement to attain at the end //1e-1
        double vmag = desired_disp/(dTref * nTref);
        //Depends on how the displacements are applied
        tryvel(0) = 0.000; //Y tension
        tryvel(1) = vmag;  //Y tension
        // tryvel(0) = vmag;  //X shear
        // tryvel(1) = 0.000; //X shear 
        
		//Initialize forces and moments
		if (force_based){
    		
    		Vector2d force_out;
    		force_out(0) = 0.00;
    		force_out(1) = 0.00;
    		double moment_out = 0.00;
    		Vector2d fluid_drag_force;
    		
    		//Find grain centroid in terms of global coordinates (not local)
    		Vector2d floeLocation = _position; //Use points for area
    		
    		double adjust_fff = 1e-4 * (1/0.000007) * 1; //0.01 too weak for mini alone //1e-4 WORKS!!!! with 1/dt   ####1e3 few, 1e7 too much, even 5e3 too much //Initially 1e2 seems okay but still tooooo much , 5e-2 still too high converge, 5e-4 still a bitty too high but improving
    		double adjust_mmm = 10.0 * 1; //1000 is excessive rot //100 is bit too much but less //Try 10  //Default 1
    		
    		//Find overlapping points and get forces and moments for each one (and accumulate)
    		double floeArea = 0;
        	size_t n = _pointList.size();
    
        	for (size_t i = 0; i < n-1; i++)
        	{
        		floeArea += ( _pointList[i](0) * _pointList[i+1](1) -  _pointList[i](1) * _pointList[i+1](0) ); 
        	}
        	floeArea += (_pointList[n-1](0) * _pointList[0](1) -  _pointList[n-1](1) * _pointList[0](0) ); 
    
        	floeArea = 0.5*abs(floeArea);
            
            //Assume regular DEM is very small (no shape effects)
            //A. If floe is smaller than cell size, add an if and use nearest grid point function round_Ocean to find force proportional to floe area
            
            //Time factor
            double nTs = 4233600.0; //As reference (just for 2018!!!)
            double fBase = 1; //Power of fixed force will go from 0 to fBase
            
            double ffactor = (double(tstep)/nTs)*fBase; //Upper load (Positive)
            //double ffactor = -(double(tstep)/nTs)*fBase; //Down load (Negative)
            
            //Round force to get average, always positive
            double ave_forceUwX = 0;
            double ave_forceUwY = 0;
            for (size_t i = 0; i < Uwg.size(); i++)
            {
                ave_forceUwX += abs(Uwg[i](0));
                ave_forceUwY += abs(Uwg[i](1));
            }
            ave_forceUwX /= Uwg.size();
            ave_forceUwY /= Uwg.size();
    		
    		if (isnan(floeArea) == false){
        		if (floeArea < (cell_sizex * cell_sizey)){
                    //Assign current velocity
                    //Uw = round_OceanV(_position, Uwg, x_cells, y_cells, offset);
                    Uw(0) = ave_forceUwX;
                    Uw(1) = ave_forceUwY;
                    //SKIN DRAG SIMPLE
                    //fluid_drag_force = rhowater*Chw*(Uw-_velocity).norm()*(Uw-_velocity) * floeArea;
                    fluid_drag_force = rhowater*Chw*(Uw).norm()*(Uw) * floeArea;
                    
                    force_out(0) += fluid_drag_force(0) * adjust_fff * 0.0 * ffactor; //Only Vertical Loading for testing
                    force_out(1) += fluid_drag_force(1) * adjust_fff * ffactor;
                    
                    //Based on centroid location
                    //double dist_x = pt_compare(0) - floeLocation(0);
                    //double dist_y = pt_compare(1) - floeLocation(1);
                    //moment_out += dist_x * force_out(1) * adjust_mmm;  //Moment due to y force using x dist.
                    //moment_out += -1.0 * dist_y * force_out(0) * adjust_mmm;  //Moment due to x force using y dist.
                    moment_out += 0.00; //Too small to add moment //CHECK THIS!!!!
                    //cout << "CASE A: Veeeery small floe compared to grid!!!" << endl;
        		}
        		else{
        		    cout << "Error: DO NOT USE NORMAL DEM!!!" << endl;
        		    cout << "Floe area: " << floeArea << " cell size: " << cell_sizex * cell_sizey << endl;
        		    exit(1);
        		}
    		}
    		else
    		{
    		    cout << "NAN points, check FORCES!!!!!" << endl;
    		    cout << "Position X: " << _position(0) << " Position Y: " << _position(1) << " id: " << index_i << endl;
    		}
    		
    		//Save values for output
    		fluidForce = force_out;
    		fluidMoment = moment_out;
    		
    		//Control
    		//cout << "Fluid force: " << fluidForce(0) << " " << fluidForce(1) << endl;
    		//cout << "Fluid moment: " << fluidMoment << endl;
		}
		else{
		    _velocity = tryvel;
		    //_position += tryvel; //Try this out
		}
		
        return;
    }


	//Constant water and air velocities        //add later   vector<Vector2d> & ppos, if needed
    void fluidInteraction (const double & Cha, const double & Cva, double & Chw, double & Cvw,
                           const double & rhoice, const double & rhoair, const double & rhowater,
                           const double & hice, double & hair, double & hwater, Vector2d & Ua, Vector2d & Uw,
                           Vector2d & fluidForceha, Vector2d & fluidForcehw, Vector2d & fluidForceh,
                           Vector2d & fluidForceva, Vector2d & fluidForcevw, Vector2d & fluidForcev,
                           Vector2d & ppv, Vector2d & ppvn, Vector2d & midp, Vector2d & fluidForce, double & fluidMoment, size_t & tstep, const Vector2d & offset, 
                           const vector<Vector2d> fluid_coord, vector<Vector2d> Uwg, const size_t x_cells, const size_t y_cells) 
    {

		//VORTEX-like field  
		if (_position(1)<3000 && _position(0)<3000)
		{
		 //Uw << 500.,0. ;   //7800 stable w/o 1D heat, 5000 stable for 1D heat (depending on grain morphology) , 2500 +/- under melting
		 Uw << 0,1000. ;
		}
		else
		{    if(_position(1)>3000 && _position(0)<3000)
		    {
		        Uw << 1000.,0. ;
		    }
		    else
		    {   if(_position(1)>3000 && _position(0)>3000)
		        {
		          //Uw << -500,0 ;
		          Uw << 0.,-1000. ;
		        }
		        else
		        {
		          //Uw << 0,500 ;
		            Uw << -1000.,0. ;
		        }
		    }
		}
		               
        double floeArea = PointsArea(_pointList); //Use points for area
        Vector2d floeCentroid = findCentroid(_pointList); //Use points for area
        //cout << "Total area: " << floeArea  << endl;
        //cout << "Ref. Point: " << _pointList[0](0)  << " , " << _pointList[0](1) << endl;
		//Generate a more arbitrary function
		//        Ua << 3000*_position(1), 3000*(-_position(0)-_position(1));  //What type of randomization be given to wind
		//        Ua << 3000*sqrt(pow((_position(0)-500),2) + pow((_position(1)-500),2)), 3000*sqrt(pow((_position(0)-500),2) + pow((_position(1)-500),2));  //What type of randomization be given to wind
		//        Uw << _position(0)*2 + _position(1)*2 , _position(0)*2 + _position(1)*2;
		//        Ua << _position(0)*2 + _position(1)*2 , _position(0)*2 + _position(1)*2;  //What type of randomization be given to wind
		//        Uw << _position(0)*2 + _position(1)*2 , _position(0)*2 + _position(1)*2;

		
		//SKIN DRAG
	    //Discretize Skin Drag
	    Vector2d UwSkin = Vector2d (0.,0.); //For each cell contribution
	    double Vxx, Vyy;

	    //1. Find extents of grains to use less cells for inspection and get extent of Xmax, Xmin, Ymax and Ymin to scan cells
	    double xmin = offset(0)*1000000.000;
	    double xmax = 0.0;
	    double ymin = offset(1)*1000000.000;
	    double ymax = 0.0;
	    size_t idx1, idx2, idy1, idy2;  // j are x cols, i are y rows
	    Vector2d fQ11; Vector2d fQ12;  Vector2d fQ21; Vector2d fQ22;
	    bool per_ind = false; //Consider PBC
	    double x1, x2, y1, y2;

	    //PointsBounds(xmin, xmax, ymin, ymax, idx1, idx2, idy1, idy2, x_cells, y_cells, offset, per_ind, fluid_coord); //For all points to find x y bounds of grain
	    //cout << "Find Centroid Bounds" << endl;
	    CentroidBounds(xmin, xmax, ymin, ymax, idx1, idx2, idy1, idy2, x_cells, y_cells, offset, per_ind, fluid_coord, x1, x2, y1, y2);  //For  applying grid force based on position of grain (grid with relatively gradual fluid velocity changes)
        //cout << "Cell Data x1: " << x1 << " x2: " << x2 << " y1: " << y1 << " y2: " << y2 << endl;
        bool simple_vel = true; //Assuming the criterion above

        //This is for all bounds (TODO: MODIFY BOUNDS AND WAY OF READING IF PBC GRAIN) //-------------->>>
	 //    double x1, x2, y1, y2;
		// x1 = fluid_coord[idx1 + 0*x_cells](0);
		// x2 = fluid_coord[idx2 + 0*x_cells](0);
		// y1 = fluid_coord[0 + idy1*x_cells](1);
		// y2 = fluid_coord[0 + idy2*x_cells](1);


	    //2. Go through all cells in rectangle (in 4 point fashion)
	    //3. Find if No Points inside (No force). 4 Points Inside (Full Force, bilin interp center).  1-3 Points inside (AreaPropForce Use fx.) (Find if it contains centroid or average other points)

	    
	    Vector2d temporal_h = Vector2d(0.0 , 0.0);
	    //Cells involved in skin drag
	    size_t n_cells = (idx2 - idx1) * (idy2 - idy1);
	    size_t n_cellsx = (idx2 - idx1);
	    size_t n_cellsy = (idy2 - idy1);
	    //Floe completely inside the same cell (Interpolate Velocity at center)	(Single Cell Contribution)
		
		//if (n_cells <= 1)
	    

	    //Floe vel at centroid using grid 
	    if (simple_vel == true) //Simple Skin grad even for multiple cells
		{
			//cout <<"Floe Smaller than Cell" << endl;
			//Find function values based on location indices
			//fQ11 = Uwg[idx1 + idy1*x_cells];
			//fQ12 = Uwg[idx1 + idy2*x_cells];
			//fQ21 = Uwg[idx2 + idy1*x_cells];
			//fQ22 = Uwg[idx2 + idy2*x_cells];
			fQ11 = Uwg[idx1];
			fQ12 = Uwg[idy1];
			fQ21 = Uwg[idx2];
			fQ22 = Uwg[idy2];

			//Find Value using bilnear interpolation
			//cout << "Skin Drag Interpolation" << endl;

			//cout << "idx1: " << idx1  << endl;
			//cout << "idx2: " << idx2 << endl;
			//cout << "idy1: " << idy1 << endl;
			//cout << "idy2: " << idy2 << endl;
			//cout << "XDir: fQ11: " << fQ11(0) << " fQ12: " << fQ12(0) << " fQ21: " << fQ21(0) << " fQ22: " << fQ22(0) << endl;
			//cout << "YDir: fQ11: " << fQ11(1) << " fQ12: " << fQ12(1) << " fQ21: " << fQ21(1) << " fQ22: " << fQ22(1) << endl;
			//cout << "x1: " << x1 << " x2: " << x2 << " y1: " << y1 << " y2: " << y2 << endl;
			//cout << "*FCxp: " << floeCentroid(0) << " *Fcyp: " << floeCentroid(1) << endl;

	       
	        Vxx = bilinear(fQ11(0), fQ12(0), fQ21(0), fQ22(0), x1, x2, y1, y2, floeCentroid(0), floeCentroid(1)); //Centroid should update in world
	        Vyy = bilinear(fQ11(1), fQ12(1), fQ21(1), fQ22(1), x1, x2, y1, y2, floeCentroid(0), floeCentroid(1)); //Centroid should update in world

	        UwSkin << Vxx, Vyy;
	        //cout << "Vxx: " << Vxx << " Vyy: " << Vyy << endl;

	        //SKIN DRAG COMPLEX (USE FULL AREA, _mass*1)
	        fluidForcehw = rhowater*Chw*((floeArea*_thick)/_density)*(UwSkin-_velocity).norm()*(UwSkin-_velocity);  //Mass = Area * thick / density = Vol / Dens
	        //fluidForcehw = rhowater*Chw*((_mass*0.91/1000000)/_density)*(UwSkin-_velocity).norm()*(UwSkin-_velocity);
		}
		//If a floe is in more than one cell, things will change (Multiple Cell Contribution)
		else
		{
			
			//cout <<"Floe Larger than Cell" << endl;
			//4. Use pointIngrain function to determine if point of grid is inside grain
	        //5. Add Forces (use cell area if full, otherwise propArea)
	        //6. Add Moments for each force in Uwg(0) and Uwg(1)
    		Vector2d acting_point;
        	double area_factor = 1;
        	// vector<Vector2d> pointListzero(_pointList.size());
        	// for (size_t i = 0; i < _pointList.size(); i++ )
        	// {
        	// 	pointListzero[i] = _pointList[i]-_position;
        	// }
        	// double total_area = PointsArea(pointListzero);
	        for (size_t i = 0; i < n_cellsy; i++) //At least 2
	        {
	        	for (size_t j = 0; j < n_cellsx; j++) //At least 2
	        	{
	        		//Change x1, x2, y1, y2 to current cell
                    x1 = fluid_coord[ (idx1 + j) + 0 * x_cells](0);
					x2 = fluid_coord[ (idx1 + j + 1) + 0 * x_cells](0);
					y1 = fluid_coord[0 + (idy1 + i) * x_cells](1);
					y2 = fluid_coord[0 + (idy1 + i +  1) * x_cells](1); 
 
                    //Get Cell Points
					Vector2d p11 = Vector2d(x1,y1); 
					Vector2d p12 = Vector2d(x1,y2); 
					Vector2d p21 = Vector2d(x2,y1); 
					Vector2d p22 = Vector2d(x2,y2); 

					bool In11 = pointIngrain(p11); bool In12 = pointIngrain(p12); bool In21 = pointIngrain(p21); bool In22 = pointIngrain(p22);

                    //Cell Completely Inside Grain
					if (In11*In12*In21*In22 == true)
					{
						acting_point = Vector2d(x1 + 0.5*(x2-x1), y1 + 0.5*(y2-y1));
						//area_factor = ((x2-x1)*(y2-y1))/((_mass*(0.91/1000000))/_density); //Assuming mass/density and grid length dimensions coincide
						area_factor = ((x2-x1)*(y2-y1))/(floeArea*_thick);
						//Find function values based on location indices
						fQ11 = Uwg[ (idx1 + j) + (idy1 + i) * x_cells];
						fQ12 = Uwg[ (idx1 + j) + (idy1 + 1 + i) * x_cells];
						fQ21 = Uwg[ (idx1 + 1 + j) + (idy1 + i) * x_cells];
						fQ22 = Uwg[ (idx1 + 1 + j) + (idy1 + 1 + i) * x_cells];

						//Find Value using bilnear interpolation
						Vxx = bilinear(fQ11(0), fQ12(0), fQ21(0), fQ22(0), x1, x2, y1, y2, acting_point(0), acting_point(1));
						Vyy = bilinear(fQ11(1), fQ12(1), fQ21(1), fQ22(1), x1, x2, y1, y2, acting_point(0), acting_point(1));

						UwSkin << Vxx, Vyy;
						//cout << "UwSkin: " << UwSkin << endl;
						//cout << "Vxx: " << Vxx << " Vyy: " << Vyy << endl;
						// cout << "XDir: fQ11: " << fQ11(0) << " fQ12: " << fQ12(0) << " fQ21: " << fQ21(0) << " fQ22: " << fQ22(0) << endl;
						// cout << "YDir: fQ11: " << fQ11(1) << " fQ12: " << fQ12(1) << " fQ21: " << fQ21(1) << " fQ22: " << fQ22(1) << endl;
						// cout << "x1: " << x1 << " x2: " << x2 << " y1: " << y1 << " y2: " << y2 << endl;
						// cout << "xp: " << acting_point(0) << " yp: " << acting_point(1) << endl;
						// cout << "Vxx: " << Vxx << " Vyy: " << Vyy << endl;
					}
					//Cell Not completely in grain (more complex interpolation)
					//First find points inside the cell
					else
					{
						vector<Vector2d> inPoints; //Points in a specific incomplete area cell
						for (size_t idx = 0; idx < _pointList.size(); idx++)
						{
							if (pointInCell(_pointList[idx],x1,x2,y1,y2) == true )
							{
								inPoints.push_back(_pointList[idx]);
							}
						}
						if (inPoints.size() > 3)
						{ 
							//cout <<"ALERT::::Number of points in cell: " << inPoints.size() << endl;
						}
						//cout <<"Number of points in cell: " << inPoints.size() << endl;
						//Add x1, x2, y1, y2 points to list to get centroid and area

						//Add intersection points to get centroid and area

						//Calculate centroid inPoints
						//acting_point = findCentroid(inPoints);
						
						//Calculate Area Contribution
						//area_factor = (PointsArea(inPoints))/(floeArea*_thick);
						
						//Find function values based on location indices
						//fQ11 = Uwg[ (idx1 + j) + (idy1 + i) * x_cells];
						//fQ12 = Uwg[ (idx1 + j) + (idy1 + 1 + i) * x_cells];
						//fQ21 = Uwg[ (idx1 + 1 + j) + (idy1 + i) * x_cells];
						//fQ22 = Uwg[ (idx1 + 1 + j) + (idy1 + 1 + i) * x_cells];

						//Find Value using bilnear interpolation
						//Vxx = bilinear(fQ11(0), fQ12(0), fQ21(0), fQ22(0), x1, x2, y1, y2, acting_point(0), acting_point(1));
						//Vyy = bilinear(fQ11(1), fQ12(1), fQ21(1), fQ22(1), x1, x2, y1, y2, acting_point(0), acting_point(1));

						//UwSkin << Vxx, Vyy;
					}

		        	//CellVelContrib(); //Find UwSkin, mass or area proportion factor and acting point (where UwSkin acts)   
		            temporal_h = rhowater*Chw*((area_factor*(floeArea*_thick))/_density)*(UwSkin-_velocity).norm()*(UwSkin-_velocity);
		            //temporal_h = rhowater*Chw*((area_factor*_mass*0.91/1000000)/_density)*(UwSkin-_velocity).norm()*(UwSkin-_velocity);
		        	fluidForcehw += temporal_h;
		        	cout << "fluidForcehw: " << fluidForcehw << endl;
		        	//fluidMoment += fluidForceva(0)*(midp(0)-_cmLset(0))+fluidForceva(1)*(midp(1)-_cmLset(1)); 
		        	fluidMoment += - temporal_h(0)*(acting_point(0)-floeCentroid(0)) + temporal_h(1)*(acting_point(1)-floeCentroid(1)); //?? + + or  + -
	        	}
	        }

		}






		//SKIN DRAG SIMPLE
		fluidForceha = rhoair*Cha*((floeArea*_thick)/_density)*Ua.norm()*Ua; //use mass0/density to get volumne or in this case to get 2D area
		//fluidForcehw = rhowater*Chw*((floeArea*_thick)/_density)*(Uw-_velocity).norm()*(Uw-_velocity);  //???????????
		//fluidForceha = rhoair*Cha*((_mass*0.91/1000000)/_density)*Ua.norm()*Ua; //use mass0/density to get volumne or in this case to get 2D area
		//fluidForcehw = rhowater*Chw*((_mass*0.91/1000000)/_density)*(Uw-_velocity).norm()*(Uw-_velocity);  //???????????
		fluidForceh = fluidForceha + fluidForcehw;


        //cout << "Fluid Moment after Skin: " << fluidMoment << endl;

		//A constant wind and ocean velocity will not induce moment due to skin drag if applied at center. But using cells that can change.



		//FORM DRAG
		//Input all points of grain to get each normal vector and compare each one to the Ua or Uw and get total Perimeter vertically pressured

		//Adjusted height based on area proportionality
		//hice=function of ( (_mass/_density) ) * hice; But applied on hair and hwater

		hair = _thick*(rhowater-rhoice)/rhowater;   //Let _thick replace hice as a grain property
		hwater = _thick-hair;

		//hair=hice*(rhowater-rhoice)/rhowater;   //Let _thick replace hice as a grain property
		//hwater=hice-hair;

		for (size_t ptidx = 0; ptidx < _pointList.size(); ptidx++) {
		  
			//          for (size_t ii = 0; ii <_pointList.size(); ii++) {
			//               ppos=_pointList[ii]+(_cmLset-_position)  //Update the new position of the center of mass to shift the rest of the points???????????
			//          }
			//          if (ii=_pointList.size()-1)
			//          {
			//             ppv=ppos[ii]-ppos[0];
			//             midp=ppv*0.5+ppos[0];  //Define midpoint for better calculation of moments
			//          }
			//          else
			//          {
			//             ppv=ppos[ii]-ppos[ii+1];
			//             midp=ppv*0.5+ppos[ii+1];
			//          }

			if (ptidx == _pointList.size()-1)
			{
			  ppv=_pointList[ptidx]-_pointList[0];
			  midp=ppv*0.5+_pointList[0];  //Define midpoint for better calculation of moments
			}
			else
			{
			  ppv=_pointList[ptidx]-_pointList[ptidx+1];  //Dist between 2 points
			  midp=ppv*0.5+_pointList[ptidx+1];
			}
			   		 
			ppvn << ppv(1), -1*ppv(0);  //Rotation of 90 degrees for a ccw ppvector s.t. normal faces outwardm, not normal, has a norm ~= 1
			Vector2d temporal_v = Vector2d(0.0 , 0.0);
			
			if(ppvn.dot(Ua) < 0)       //if(ppvn(0)*Ua(0)+ppvn(1)*Ua(1) < 0)
			{
			  temporal_v = rhoair*Cva*hair*abs(ppvn.dot(Ua))*Ua;
			  fluidForceva += temporal_v;   //Accumulation of forces if Ua crashes against face, otherwise, no body drag occurs, it occurs on the projection
			  //fluidForceva+=rhoair*Cva*hair*ppv.norm()*Ua.norm()*Ua;   //Accumulation of forces if Ua crashes against face, otherwise, no body drag occurs
			  //fluidMoment += fluidForceva(0)*(midp(0)-_cmLset(0))+fluidForceva(1)*(midp(1)-_cmLset(1));  //Is this CM of Mass Shifted as we go
			  fluidMoment += - temporal_v(0)*(midp(0)-floeCentroid(0)) + temporal_v(1)*(midp(1)-floeCentroid(1));  //?? +
			}
			   
			//cout << "Form Drag CHECK" << endl;
			//cout << "ppv: " << ppv(0) << " , " << ppv(1) << endl;
			//cout << "ppvn: " << ppvn(0) << " , " << ppvn(1) << endl;
			//cout << "Rel Vel: " << Uw(0)-_velocity(0) << " , " << Uw(1)-_velocity(1) << endl;


	        //Update to defined grid

            //cout << "MIDP: " << midp(0) << " , " << midp(1) << endl; 
	        //Find Uw at midp
	        //cout << "Form Drag Interpolation" << endl;
			//cout << "Midp: " << midp(0) << "  " << midp(1) << endl;	        
	        Vector2d midUwg = bilinear_Interp(midp, fluid_coord, Uwg, x_cells, y_cells, offset); //Interpolate at point using grid info  //Check for performance
	        //cout << "MIDUWG: " << midUwg(0) << " , " << midUwg(1) << endl;
	        //Vector2d midUwg = Uw;

	        if(ppvn.dot(midUwg - _velocity) < 0) 
	        { 
	        	//cout << "Water Form Drag" << endl;
	        	//cout << "Form Drag Water" << endl;
	        	temporal_v = rhowater * Cvw * hwater * abs(ppvn.dot(midUwg-_velocity))*(midUwg-_velocity);  //ppvn norm scales perimeter of floe going against current
	        	fluidForcevw += temporal_v;
	        	//fluidMoment += fluidForcevw(0)*(midp(0)-_cmLset(0)) + fluidForcevw(1)*(midp(1)-_cmLset(1));  //More works
	        	fluidMoment += - temporal_v(0)*(midp(0)-floeCentroid(0)) + temporal_v(1)*(midp(1)-floeCentroid(1));  //?? +
	        }	

			// if(ppvn.dot(Uw-_velocity) < 0) // if(ppvn(0)*(Uw(0)-_velocity(0))+ppvn(1)*(Uw(1)-_velocity(1)) < 0)
			// {
			//   cout << "Water Form Drag" << endl;
			//   fluidForcevw+=rhowater*Cvw*hwater*abs(ppvn.dot(Uw-_velocity))*(Uw-_velocity);
			//   //fluidForcevw+=rhowater*Cvw*hwater*ppv.norm()*(Uw-_velocity).norm()*(Uw-_velocity);
			//   fluidMoment+=fluidForcevw(0)*(midp(0)-_cmLset(0)) + fluidForcevw(1)*(midp(1)-_cmLset(1)); //Should be + for ccw positive moment
			//   //cout << "X: " << fluidForcevw(0) << " Y: " <<  fluidForcevw(1) << " Moment: " << fluidForcevw(0)*(midp(0)-_cmLset(0))+fluidForcevw(1)*(midp(1)-_cmLset(1)) << endl;

			//   //fluidMoment*=0.0; //Inspect Moment
			// }
		}

		//fluidForceva=rhoair*Cva*hair*Perair*Ua.norm()*Ua;     //use Perimeter*exposed_height to get 2D area
		//fluidForcevw= rhowater*Cvw*hwater*Perwater*(Uw-_velocity).norm()*(Uw-_velocity);
		fluidForcev = fluidForceva + fluidForcevw;
		fluidForce = fluidForcev + fluidForceh;

		cout << "Skin Drag Force"  << endl;
		cout << "X: " << fluidForceh(0) << " Y: " <<  fluidForceh(1) << endl;
		cout << "Form Drag Force"  << endl;
		cout << "X: " << fluidForcev(0) << " Y: " <<  fluidForcev(1) << endl;
		cout << "Fluid Moment after Form: " << fluidMoment << endl;

		fluidMoment *= (_mass/_density/floeArea); //TEMPORARY!!! //Calibration Factor to Avoid Huge Moments TODO IMPROVE  (_momentInertia / floeInertia)
        fluidForce *= 10000*(_mass/_density/floeArea);

        //cout << "Mass Factor: " << (_mass/_density/floeArea) << endl;

		// cout << "Position X: " << _position(0) << " Position Y: " << _position(1) << endl;
		// cout << "CM X: " << _cmLset(0) << " CM Y: " << _cmLset(1) << endl;
		// cout << "Point 0 X: " << _pointList[0](0) << " Point 0 Y: " << _pointList[0](1) << endl;
		// cout << "Point 70 X: " << _pointList[70](0) << " Point 70 Y: " << _pointList[70](1) << endl;

		//_momentInertia = 0; _omega = 0; _cmLset;

		//Get grain center of mass and moment of inertia to calculate moment
		//Moment=fluidForcev(0)*(_position(0)-_cmLset(0))+fluidForcev(1)*(_position(1)-_cmLset(1))
		//fluidMoment= SUM  fluidForcev*ptCM
        
    }

	double round_Ocean(Vector2d & pointxy, vector<double> & oceanTemp, const size_t & x_cells, const size_t & y_cells, const Vector2d & offset)
	{
		size_t idx1, idy1;  // j are x cols, i are y rows
		int cell_sizex = int(offset(0))/int(x_cells-1);
		int cell_sizey = int(offset(1))/int(y_cells-1);
		//cout << "cell size x: " << cell_sizex << " cell_sizey: " << cell_sizey << endl;

		double tol_pos = 2.00;
		if ( isnan(pointxy(0)) ||  isnan(pointxy(1)) || pointxy(0) > offset(0) + tol_pos  || pointxy(0) < 0.0 - tol_pos || pointxy(1) > offset(1) + tol_pos || pointxy(1) < 0.0 - tol_pos )
		{
		  //cout <<"WARNING NaN or BAD POSITION VALUES PROVIDED FOR INTERPOLATION, OUTPUT = 0" << endl;
		  return 0.0;
		}

		Vector2d input_p = pointxy; //Input point that accounts for periodic bcs when needed, so that way you just shift the calculation of this point when needed
		//Update periodic BCs for pointxy coordinates
		if (pointxy(0) > offset(0))
		{
		  input_p(0) = pointxy(0) - offset(0);
		}
		if (pointxy(0) < 0.0)
		{
		  input_p(0) = pointxy(0) + offset(0);
		}
		if (pointxy(1) > offset(1))
		{
		  input_p(1) = pointxy(1) - offset(1);
		}
		if (pointxy(1) < 0.0)
		{
		  input_p(1) = pointxy(1) + offset(1);
		}

		//Now we use the periodic coordinates to find the index of col and row for space coordinates
		idx1 = size_t(round(input_p(0)) / cell_sizex);

		//Adjust for ycells good indexing
		idy1 = size_t(round(input_p(1)) / cell_sizey);
		if (idy1 == y_cells - 1)
		{
		  idy1 = y_cells - 2;
		}

		//Use this exact grid point rounded from the input point at local scale
		return oceanTemp[idx1 + idy1*x_cells];
	}
	
	Vector2d round_OceanV(Vector2d & pointxy, vector<Vector2d> & oceanV, const size_t & x_cells, const size_t & y_cells, const Vector2d & offset)
	{
		size_t idx1, idy1;  // j are x cols, i are y rows
		int cell_sizex = int(offset(0))/int(x_cells-1);
		int cell_sizey = int(offset(1))/int(y_cells-1);
		//cout << "cell size x: " << cell_sizex << " cell_sizey: " << cell_sizey << endl;

		double tol_pos = 2.00;
		if ( isnan(pointxy(0)) ||  isnan(pointxy(1)) || pointxy(0) > offset(0) + tol_pos  || pointxy(0) < 0.0 - tol_pos || pointxy(1) > offset(1) + tol_pos || pointxy(1) < 0.0 - tol_pos )
		{
		  //cout <<"WARNING NaN or BAD POSITION VALUES PROVIDED FOR INTERPOLATION, OUTPUT = 0" << endl;
		  Vector2d error_v;
		  error_v << 0.0 , 0.0;
		  return error_v;
		}

		Vector2d input_p = pointxy; //Input point that accounts for periodic bcs when needed, so that way you just shift the calculation of this point when needed
		//Update periodic BCs for pointxy coordinates
		if (pointxy(0) > offset(0))
		{
		  input_p(0) = pointxy(0) - offset(0);
		}
		if (pointxy(0) < 0.0)
		{
		  input_p(0) = pointxy(0) + offset(0);
		}
		if (pointxy(1) > offset(1))
		{
		  input_p(1) = pointxy(1) - offset(1);
		}
		if (pointxy(1) < 0.0)
		{
		  input_p(1) = pointxy(1) + offset(1);
		}

		//Now we use the periodic coordinates to find the index of col and row for space coordinates
		idx1 = size_t(round(input_p(0)) / cell_sizex);

		//Adjust for ycells good indexing
		idy1 = size_t(round(input_p(1)) / cell_sizey);
		if (idy1 == y_cells - 1)
		{
		  idy1 = y_cells - 2;
		}

		//Use this exact grid point rounded from the input point at local scale
		return oceanV[idx1 + idy1*x_cells];
	}


    //This function interpolates the value of Uwg at pointxy using Uwg values at fluid_coord	
    Vector2d bilinear_Interp(const Vector2d & pointxy, const vector<Vector2d> & fluid_coord, const vector<Vector2d> & Uwg, const size_t & x_cells, const size_t & y_cells, const Vector2d & offset)
    {
    	Vector2d Output;
    	size_t idx1, idx2, idy1, idy2;  // j are x cols, i are y rows
    	double Vxx, Vyy;
    	Vector2d fQ11; Vector2d fQ12;  Vector2d fQ21; Vector2d fQ22;
    	double x1, x2, y1, y2;
    	//double fXQ11, fXQ12, fXQ21, fXQ22, fYQ11, fYQ12, fYQ21, fYQ22, x1, x2, y1, y2;
        int cell_sizex = int(offset(0))/int(x_cells-1);
        int cell_sizey = int(offset(1))/int(y_cells-1);
        //cout << "cell size x: " << cell_sizex << " cell_sizey: " << cell_sizey << endl;

        Vector2d input_p = pointxy; //Input point that accounts for periodic bcs when needed
        //Update periodic BCs for pointxy coordinates
        if (pointxy(0) > offset(0))
        {
        	input_p(0) = pointxy(0) - offset(0);
        }
    	if (pointxy(0) < 0.0)
        {
        	input_p(0) = pointxy(0) + offset(0);
        }
        if (pointxy(1) > offset(1))
        {
        	input_p(1) = pointxy(1) - offset(1);
        }
    	if (pointxy(1) < 0.0)
        {
        	input_p(1) = pointxy(1) + offset(1);
        }


    	idx1 = size_t(input_p(0)) / size_t(cell_sizex);
    	idx2 = idx1 + 1;
    	
    	idy1 = size_t(input_p(1)) / size_t(cell_sizey);
    	idy2 = idy1 + 1;
  //   	//Find 4 Points using pointxy as Input
		// // IN Y
		// for (size_t i = 0; i < y_cells; i++) {
		// 	if ( pointxy(1) < fluid_coord[0+i*x_cells](1) )
		// 	{
		// 		idy2 = i;
		// 		idy1 = i-1;
		// 		break;
		// 	}
		// }	
		// //IN X
		// for (size_t i = 0; i < x_cells; i++) {
		// 	if ( pointxy(0) < fluid_coord[i+0*x_cells](0) )
		// 	{
		// 		idx2 = i;
		// 		idx1 = i-1;
		// 		break;
		// 	}

		// }	

        //cout << "inputpx: " << input_p(0) << " inputpy: " << input_p(1) << endl;	
		//Find point values value based on location indices	
		size_t NNx_cells = size_t(offset(0))/(cell_sizex) + 1;
        size_t NNy_cells = size_t(offset(1))/(cell_sizey) + 1;
		// cout << "x_cells: " << NNx_cells << endl;
  //       cout << "y_cells: " << NNy_cells << endl;

		// cout << "idx1: " << idx1 << endl;
		// cout << "idy1: " << idy1 << endl;
		// cout << "idx2: " << idx2 << endl;
		// cout << "idy2: " << idy2 << endl;
		// cout << "fluid_coord size: " << fluid_coord.size() << endl;
		x1 = fluid_coord[idx1 + 0*x_cells](0);
		x2 = fluid_coord[idx2 + 0*x_cells](0);
		y1 = fluid_coord[0 + idy1*x_cells](1);
		y2 = fluid_coord[0 + idy2*x_cells](1);

		//Find function values based on location indices
		// cout << "idx1 + idy1*x_cells: " << idx1 + idy1*x_cells << endl;
		// cout << "idx1 + idy2*x_cells: " << idx1 + idy2*x_cells << endl;
		// cout << "idx2 + idy1*x_cells: " << idx2 + idy1*x_cells << endl;
		// cout << "idx2 + idy2*x_cells: " << idx2 + idy2*x_cells << endl;
		// cout << "UWG size: " << Uwg.size() << endl;
		fQ11 = Uwg[idx1];
		fQ12 = Uwg[idy1];
		fQ21 = Uwg[idx2];
		fQ22 = Uwg[idy2];
		//fQ11 = Uwg[idx1 + idy1*x_cells];
		//fQ12 = Uwg[idx1 + idy2*x_cells];
		//fQ21 = Uwg[idx2 + idy1*x_cells];
		//fQ22 = Uwg[idx2 + idy2*x_cells];

		// cout << "XDir: fQ11: " << fQ11(0) << " fQ12: " << fQ12(0) << " fQ21: " << fQ21(0) << " fQ22: " << fQ22(0) << endl;
		// cout << "YDir: fQ11: " << fQ11(1) << " fQ12: " << fQ12(1) << " fQ21: " << fQ21(1) << " fQ22: " << fQ22(1) << endl;
		// cout << "x1: " << x1 << " x2: " << x2 << " y1: " << y1 << " y2: " << y2 << endl;
		// cout << "xp: " << pointxy(0) << " yp: " << pointxy(1) << endl;
		// cout << "inputpx: " << input_p(0) << " inputpy: " << input_p(1) << endl;				

		//Find Value using bilnear interpolation
        Vxx = bilinear(fQ11(0), fQ12(0), fQ21(0), fQ22(0), x1, x2, y1, y2, input_p(0), input_p(1));
        Vyy = bilinear(fQ11(1), fQ12(1), fQ21(1), fQ22(1), x1, x2, y1, y2, input_p(0), input_p(1));

        //cout << "Vxx: " << Vxx << " Vyy: " << Vyy << endl;

        Output << Vxx, Vyy;
    	return Output;
    }  

    //This function interpolates the value of OceanTemp at pointxy using values at fluid_coord	
    double bilinear_Interp_Ocean(Vector2d & pointxy, const vector<Vector2d> & fluid_coord, vector<double> & oceanTemp, const size_t & x_cells, const size_t & y_cells, const Vector2d & offset)
    {
    	double Output;
    	size_t idx1, idx2, idy1, idy2;  // j are x cols, i are y rows
    	double Vxx; //Value for calc.
    	double fQ11, fQ12, fQ21, fQ22; //Ocean temp values
    	double x1, x2, y1, y2;  //Coordinate values
        int cell_sizex = int(offset(0))/int(x_cells-1);
        int cell_sizey = int(offset(1))/int(y_cells-1);
        //cout << "cell size x: " << cell_sizex << " cell_sizey: " << cell_sizey << endl;

	    if ( isnan(pointxy(0)) ||  isnan(pointxy(1)) )
	    {
	      cout <<"WARNING NaN VALUES PROVIDED FOR BILNEAR INTERPOLATION, OUTPUT = 0" << endl;
	      return 0.0;
	    }

        Vector2d input_p = pointxy; //Input point that accounts for periodic bcs when needed, so that way you just shift the calculation of this point when needed
        //Update periodic BCs for pointxy coordinates
        if (pointxy(0) > offset(0))
        {
        	input_p(0) = pointxy(0) - offset(0);
        }
    	if (pointxy(0) < 0.0)
        {
        	input_p(0) = pointxy(0) + offset(0);
        }
        if (pointxy(1) > offset(1))
        {
        	input_p(1) = pointxy(1) - offset(1);
        }
    	if (pointxy(1) < 0.0)
        {
        	input_p(1) = pointxy(1) + offset(1);
        }


        //Now we use the periodic coordinates to find the index of col and row for space coordinates
    	idx1 = size_t(input_p(0)) / size_t(cell_sizex);
    	idx2 = idx1 + 1;
    	
    	//Adjust for ycells good indexing
    	idy1 = size_t(input_p(1)) / size_t(cell_sizey);
    	if (idy1 == y_cells - 1)
    	{
    		idy1 = y_cells - 2;
    	}
    	idy2 = idy1 + 1;
  //   	//Find 4 Points using pointxy as Input
		// // IN Y
		// for (size_t i = 0; i < y_cells; i++) {
		// 	if ( pointxy(1) < fluid_coord[0+i*x_cells](1) )
		// 	{
		// 		idy2 = i;
		// 		idy1 = i-1;
		// 		break;
		// 	}
		// }	
		// //IN X
		// for (size_t i = 0; i < x_cells; i++) {
		// 	if ( pointxy(0) < fluid_coord[i+0*x_cells](0) )
		// 	{
		// 		idx2 = i;
		// 		idx1 = i-1;
		// 		break;
		// 	}

		// }	

	    //Verify if it's an exact grid point
	    if ( ( fmod(input_p(0), int(cell_sizex)) == 0) && ( fmod(input_p(1), int(cell_sizey))  == 0  )  )
	    {
	       return oceanTemp[idx1 + idy1*x_cells];
	    }

        //cout << "inputpx: " << input_p(0) << " inputpy: " << input_p(1) << endl;	
		//Find point values value based on location indices	
		size_t NNx_cells = size_t(offset(0))/(cell_sizex) + 1;
        size_t NNy_cells = size_t(offset(1))/(cell_sizey) + 1;
		// cout << "x_cells: " << NNx_cells << endl;
  //       cout << "y_cells: " << NNy_cells << endl;

		// cout << "idx1: " << idx1 << endl;
		// cout << "idy1: " << idy1 << endl;
		// cout << "idx2: " << idx2 << endl;
		// cout << "idy2: " << idy2 << endl;
		// cout << "fluid_coord size: " << fluid_coord.size() << endl;

		//These row and column index give us the global coordinates in x and y
		x1 = fluid_coord[idx1 + 0*x_cells](0);
		x2 = fluid_coord[idx2 + 0*x_cells](0);
		y1 = fluid_coord[0 + idy1*x_cells](1);
		y2 = fluid_coord[0 + idy2*x_cells](1);

		//Find function values based on location indices (IFF oceanTemp and fluid_coord have the same size, which they must have)
		// cout << "idx1 + idy1*x_cells: " << idx1 + idy1*x_cells << endl;
		// cout << "idx1 + idy2*x_cells: " << idx1 + idy2*x_cells << endl;
		// cout << "idx2 + idy1*x_cells: " << idx2 + idy1*x_cells << endl;
		// cout << "idx2 + idy2*x_cells: " << idx2 + idy2*x_cells << endl;
		// cout << "oceanTemp size: " << oceanTemp.size() << endl;
		//fQ11 = oceanTemp[idx1];
		//fQ12 = oceanTemp[idy1];
		//fQ21 = oceanTemp[idx2];
		//fQ22 = oceanTemp[idy2];
		fQ11 = oceanTemp[idx1 + idy1*x_cells];
		fQ12 = oceanTemp[idx1 + idy2*x_cells];
		fQ21 = oceanTemp[idx2 + idy1*x_cells];
		fQ22 = oceanTemp[idx2 + idy2*x_cells];

		// cout << "XDir: fQ11: " << fQ11 << " fQ12: " << fQ12 << " fQ21: " << fQ21 << " fQ22: " << fQ22 << endl;
		// cout << "x1: " << x1 << " x2: " << x2 << " y1: " << y1 << " y2: " << y2 << endl;
		// cout << "xp: " << pointxy(0) << " yp: " << pointxy(1) << endl;
		// cout << "inputpx: " << input_p(0) << " inputpy: " << input_p(1) << endl;				

		//Find Value using bilnear interpolation
        Vxx = bilinear(fQ11, fQ12, fQ21, fQ22, x1, x2, y1, y2, input_p(0), input_p(1));

        //cout << "Oxx: " << Vxx << endl;

        Output = Vxx;
    	return Output;
    }  


    double bilinear(const double & fQ11, const double & fQ12, const double & fQ21, const double & fQ22,  const double & x1, const double & x2, const double & y1, const double & y2, const double & x, const double & y)
    {
    	if (x1 == x2 || y1 == y2)
    	{
    		cout <<"WARNING WRONG VALUES PROVIDED FOR BILNEAR INTERPOLATION, OUTPUT = 0" << endl;
    		return 0.0;
    	}

    	if (abs(x) > 1000000.0 ||  abs(y) > 1000000.0 )
	    {
	      cout <<"WARNING TOO HIGH VALUES PROVIDED FOR BILNEAR INTERPOLATION, OUTPUT = 0" << endl;
	      return 0.0;
	    }

	    if (isnan(x) ||  isnan(y) )
	    {
	      cout <<"WARNING NaN VALUES PROVIDED FOR BILNEAR INTERPOLATION, OUTPUT = 0" << endl;
	      return 0.0;
	    }

    	double a0, a1, a2, a3;
    	a0 = (( fQ11 * x2 * y2 )/( (x1 - x2) * (y1 - y2) )) +  (( fQ12 * x2 * y1 )/( (x1 - x2) * (y2 - y1) )) +  (( fQ21 * x1 * y2 )/( (x1 - x2) * (y2 - y1) )) +  (( fQ22 * x1 * y1 )/( (x1 - x2) * (y1 - y2) )) ;
    	a1 = (( fQ11 * y2 )/( (x1 - x2) * (y2 - y1) )) +  (( fQ12 * y1 )/( (x1 - x2) * (y1 - y2) )) +  (( fQ21 * y2 )/( (x1 - x2) * (y1 - y2) )) +  (( fQ22 *  y1 )/( (x1 - x2) * (y2 - y1) )) ;
    	a2 = (( fQ11 * x2 )/( (x1 - x2) * (y2 - y1) )) +  (( fQ12 * x2 )/( (x1 - x2) * (y1 - y2) )) +  (( fQ21 * x1 )/( (x1 - x2) * (y1 - y2) )) +  (( fQ22 *  x1 )/( (x1 - x2) * (y2 - y1) )) ;
    	a3 = (( fQ11 )/( (x1 - x2) * (y1 - y2) )) +  (( fQ12 )/( (x1 - x2) * (y2 - y1) )) +  (( fQ21 )/( (x1 - x2) * (y2 - y1) )) +  (( fQ22 )/( (x1 - x2) * (y1 - y2) )) ;
        // cout << "a0: " << a0 << endl;
        // cout << "a1: " << a1 << endl;
        // cout << "a2: " << a2 << endl;
        // cout << "a3: " << a3 << endl;
    	return a0 + a1*x + a2*y + a3*x*y;
    }

    //Decide if an arbitrary points is inside a grain PointList
    bool pointIngrain(const Vector2d & pointxy)
    {
    	
        if ((pointxy - _position).norm() > _radius )
    	{
    		return false;
    	} 
    	//Find Max and Min Distances
    	//double maxDist = 0.0;  //Should be radius
    	double minDist = _radius*10000000.0;
    	Vector2d cpoint;
    	for (size_t ptidx = 0; ptidx < _pointList.size(); ptidx++) {
    		cpoint = _pointList[ptidx];
    		// if ( (cpoint-_cmLset).norm() > maxDist )
    		// {
    		// 	maxDist = (cpoint-_cmLset).norm();
    		// }
    		if ( (cpoint-_position).norm() < minDist )
    		{
    			minDist = (cpoint-_position).norm();
    		}

    	}

    	// if ((pointxy - _cmLset).norm() > maxDist )
    	// {
    	// 	return false;
    	// }
    	if ((pointxy - _position).norm() <= minDist )
        {
        	return true;
        }

        //Otherwise Need to interpolate by approximation
        if ( (pointxy - _position).norm() <= _radius &&  (pointxy - _position).norm() > minDist)
        {	
        	//Find two closest Grain Points to point of interest
        	double closeD1 = _radius*100.0; double closeD2 = _radius*100.0;
        	size_t closeidx1, closeidx2;
        	//Closest Point
        	for (size_t ptidx = 0; ptidx < _pointList.size(); ptidx++) {
        		cpoint = _pointList[ptidx];
        		if ( (cpoint-pointxy).norm() < closeD1 )
	    		{
	    			closeD1 = (cpoint-_position).norm();
	    			closeidx1 = ptidx;
	    		}
        	}	
            //Second Closest Point
        	for (size_t ptidx = 0; ptidx < _pointList.size(); ptidx++) {
        		if (ptidx != closeidx1)
	        	{
	        		cpoint = _pointList[ptidx];
	        		if ( (cpoint-pointxy).norm() < closeD2 )
		    		{
		    			closeD2 = (cpoint-_position).norm();
		    			closeidx2 = ptidx;
		    		}
		    	}	
        	}
            //Line 1
        	Vector2d pG1 = _pointList[closeidx1];
        	Vector2d pG2 = _pointList[closeidx2];	
            //Line 2 is _centroid and pointxy

            //Solve for intersection if slopes are NOT the same for Lines 1 and 2
            double m1  = (pG2(1) - pG1(1))/(pG2(0) - pG1(1));  double m2  = (_position(1) - pointxy(1))/(_position(0) - pointxy(0));

            if (m1 == m2)
            {
            	cout <<"ERROR SAME SLOPES, ABORT AND FALSE"<< endl;
            	return false;
            }

            double b1 = pG1(1) - m1*pG1(0);  double b2 = pointxy(1) - m2*pointxy(0);

            //Find Intersection Point for Comparison
            Vector2d Inter_Point;

            //Solve system for x and y
            double xsol, ysol;
            xsol = (b2-b1)/(m2-m1);
            ysol = m1*xsol + b1;

            Inter_Point << xsol, ysol;

            //Compare Distance, if more not inside, if equal or less inside
            if ((pointxy - _position).norm() > (Inter_Point - _position).norm())
            {
            	return false;
            }
            if ((pointxy - _position).norm() <= (Inter_Point - _position).norm())
            {
				return true;
            }

        }	
    }

    bool pointInCell(const Vector2d & pointxy, const double & x1, const double & x2, const double & y1, const double & y2)
    {
    	if ( pointxy(0) >= x1 && pointxy(0) <= x2 && pointxy(1) >= y1 && pointxy(1) <= y2 )
    	{
    		return true;
    	}
    	return false;
    }

    //Find MaxMin Positions for all Points of a Grain
    void PointsBounds(double & xmin, double & xmax, double & ymin, double & ymax, size_t & idx1, size_t & idx2, size_t & idy1, size_t & idy2,  const size_t & x_cells, const size_t & y_cells, const Vector2d & offset, bool & per_ind, const vector<Vector2d> & fluid_coord)
    {
    	Vector2d cpoint;
    	for (size_t ptidx = 0; ptidx < _pointList.size(); ptidx++) {
    		cpoint = _pointList[ptidx];
    		
    		//First Process for periodic PBC to keep well constrained problem
	        if (cpoint(0) > offset(0))
	        {
	        	cpoint(0) = cpoint(0) - offset(0);
	        }
	    	if (cpoint(0) < 0.0)
	        {
	        	cpoint(0) = cpoint(0) + offset(0);
	        }
	        if (cpoint(1) > offset(1))
	        {
	        	cpoint(1) = cpoint(1) - offset(1);
	        }
	    	if (cpoint(1) < 0.0)
	        {
	        	cpoint(1) = cpoint(1) + offset(1);
	        }

            //Find Max and Min
    		if ( cpoint(0) > xmax  )
    		{
    			xmax = cpoint(0);
    		}
    		if ( cpoint(1) > ymax )
    		{
    			ymax = cpoint(1);
    		}
    		if ( cpoint(0) < xmin )
    		{
    			xmin = cpoint(0);
    		}
    		if ( cpoint(1) < ymin )
    		{
    			ymin = cpoint(1);
    		}
    	}

        int cell_sizex = int(offset(0))/int(x_cells-1);
        int cell_sizey = int(offset(1))/int(y_cells-1);



    	idx1 = size_t( int(xmin) / cell_sizex );
    	idx2 = size_t( (int(xmax) / cell_sizex) + 1 ); //CHECK JUST IN CASE
    	idy1 = size_t( int(ymin) / cell_sizey );
    	idy2 = size_t( (int(ymax) / cell_sizey) + 1 ); //CHECK JUST IN CASE

    	double x1, x2, y1, y2;
		x1 = fluid_coord[idx1 + 0*x_cells](0);
		x2 = fluid_coord[idx2 + 0*x_cells](0);
		y1 = fluid_coord[0 + idy1*x_cells](1);
		y2 = fluid_coord[0 + idy2*x_cells](1);

    	//Update for Periodic BCS
    	if ( (x2 - x1) == offset(0) ||  (y2 - y1) == offset(1) )
    	{
    		per_ind = true;
    	}
        else
        {
        	per_ind = false; 
        }
    }

    //Find x1 x2 y1 y2 grid Positions for centroid of grain
    void CentroidBounds(double & xmin, double & xmax, double & ymin, double & ymax, size_t & idx1, size_t & idx2, size_t & idy1, size_t & idy2,  const size_t & x_cells, const size_t & y_cells, const Vector2d & offset, bool & per_ind, const vector<Vector2d> & fluid_coord, double & x1, double & x2, double & y1, double & y2)
    {
    	Vector2d cpoint;
    	//for (size_t ptidx = 0; ptidx < _pointList.size(); ptidx++) {
    		cpoint = findCentroid(_pointList);  //_cmLset;
    		
    		//First Process for periodic PBC to keep well constrained problem
	        if (cpoint(0) > offset(0))
	        {
	        	cpoint(0) = cpoint(0) - offset(0);
	        }
	    	if (cpoint(0) < 0.0)
	        {
	        	cpoint(0) = cpoint(0) + offset(0);
	        }
	        if (cpoint(1) > offset(1))
	        {
	        	cpoint(1) = cpoint(1) - offset(1);
	        }
	    	if (cpoint(1) < 0.0)
	        {
	        	cpoint(1) = cpoint(1) + offset(1);
	        }

            //Find Max and Min
    		if ( cpoint(0) > xmax  )
    		{
    			xmax = cpoint(0);
    		}
    		if ( cpoint(1) > ymax )
    		{
    			ymax = cpoint(1);
    		}
    		if ( cpoint(0) < xmin )
    		{
    			xmin = cpoint(0);
    		}
    		if ( cpoint(1) < ymin )
    		{
    			ymin = cpoint(1);
    		}
    	//}

        int cell_sizex = int(offset(0))/int(x_cells-1);
        int cell_sizey = int(offset(1))/int(y_cells-1);



    	idx1 = size_t( int(xmin) / cell_sizex );
    	idx2 = size_t( (int(xmax) / cell_sizex) + 1 ); //CHECK JUST IN CASE
    	idy1 = size_t( int(ymin) / cell_sizey );
    	idy2 = size_t( (int(ymax) / cell_sizey) + 1 ); //CHECK JUST IN CASE

		x1 = fluid_coord[idx1 + 0*x_cells](0);
		x2 = fluid_coord[idx2 + 0*x_cells](0);
		y1 = fluid_coord[0 + idy1*x_cells](1);
		y2 = fluid_coord[0 + idy2*x_cells](1);

        //Aditional bound validations
		if (x1 < cell_sizex || x1 >= offset(0))
		{
			x1 = 0.0;
			x2 = cell_sizex;
		}


		if (y1 < cell_sizey || y1 >= offset(1))
		{
			y1 = 0.0;
			y2 = cell_sizey;
		}

    	//Update for Periodic BCS
    	if ( (x2 - x1) == offset(0) ||  (y2 - y1) == offset(1) )
    	{
    		per_ind = true;
    	}
        else
        {
        	per_ind = false; 
        }
    }


    //Find area Using Points (closed Polygon)
    double PointsArea ( vector<Vector2d> & VecOrigin ) 
    {
    	double Area;
    	size_t n = VecOrigin.size();

    	for (size_t i = 0; i < n-1; i++)
    	{
    		Area += ( VecOrigin[i](0) * VecOrigin[i+1](1) -  VecOrigin[i](1) * VecOrigin[i+1](0) ); 
    	}
    	Area += (VecOrigin[n-1](0) * VecOrigin[0](1) -  VecOrigin[n-1](1) * VecOrigin[0](0) ); 

    	return 0.5*abs(Area);
    }

    
    //Find Centroid Using Points
    Vector2d findCentroid  ( vector<Vector2d> & VecOrigin ) 
    {
    	Vector2d Output;
    	size_t n = VecOrigin.size();
    	double SumX, SumY;

    	for (size_t i = 0; i < n; i++)
    	{
    		SumX += VecOrigin[i](0);
    		SumY += VecOrigin[i](1);
    	}
       
 		Output << (1.0/double(n))*SumX , (1.0/double(n))*SumY;
    	return Output;
    }

    
    // -------------------------------------------------------------------------------------------------------------------------------------------------------------
    //End Fluid Modification


    //Temperature over 1-D Thickness and/or 2D Area
    
    //Change Aug 22, 2022
    void TemperatureModif (double Tair, double Twater, const double & dt, double dh, const double alphaice, double meltTemp, double meltid, double & meltVin, size_t & tstep, const double & KIc,  const double afactor,  const double & fdim, const size_t & START_TEMP, const double & Khor, const double & meltVSun, const size_t & melt_flag, const size_t & dstep, const vector<Vector2d> & fluid_coord, const size_t & x_cells, const size_t & y_cells, vector<double> & oceanTemp, const Vector2d & offset,
                           double & loss_mcv_temp, double & loss_mcv_solar_temp, double & loss_mcv_ocean_temp, double & hlat)
    //void TemperatureModif (double Tair, double Twater, const double & dt, double dh, const double alphaice, double meltTemp, double meltid, double & meltVin, size_t & tstep, const double & KIc,  const double afactor,  const double & fdim, const size_t & START_TEMP, const double & Khor, const double & meltVSun, const size_t & melt_flag, const size_t & dstep, const vector<Vector2d> & fluid_coord, const size_t & x_cells, const size_t & y_cells, vector<double> & oceanTemp, const Vector2d & offset) //double Utemper[80000], double Utemper0[80000] 
    {      
         
          vector<double> Tvec=_grainTemp2D.getLevelset();
          vector<double> Tvec0=_grainTemp2D0.getLevelset();

          vector<double> LsetvecCheck=_lset.getLevelset();

          vector<double> Lsetvec=_lset.getLevelset();
          vector<double> Lsetvec0=_lset0.getLevelset();
          vector<double> InOutvec=_lset.getLevelset();

          //Thickness of ice
          vector<double> Thickvec=_Mthick.getLevelset();
          vector<double> Thickvec0=_Mthick.getLevelset();

       //    //Adjust Boundary Conditions by Quadrant
	      // if (_position(1)<500 && _position(0)<500)
	      // {
	      //    Tair = Tair ;
	      //    Twater = 1.5*Twater ;
	      // } 
	      // else
	      // {    if(_position(1)>500 && _position(0)<500)
	      //       {
	      //          Tair = Tair ;
	      //          Twater = 1.5*Twater ;
	      //       }
	      //       else
	      //       {   if(_position(1)>500 && _position(0)>500)
	      //           {
	      //             Tair = Tair ;
	      //             Twater = 1.5*Twater ;
	      //           }
	      //           else
	      //           {
	      //             Tair = Tair ;
	      //             Twater = 1.5*Twater ;
	      //           }
	      //       }
	      // }

          Tair = Tair ;
          //Twater = 1.5*Twater ; //RECENT CHANGE
          Twater  = 1.0*Twater;

          double meltV = meltVin;
	      //meltV = (-0.000001/20)*(Twater-meltTemp)*alphaice;  //-0.00001 good for melt low speed //Melt
	      //meltV = (-0.000008/20)*(Tair-meltTemp)*alphaice;    //Freeze


	     // //Update Boundary Conditions (INITIAL THERMAL CONDITIONS OF THE ICE FLOE COME FROM readInputFile.h)
	     // _Utemper[0]=Twater;    //Bottom
      //    _Utemper[800-1]=Tair; //Top


        //2D Update Constant Fluid Temperature Outside of the Floe using only a constant Twater (inside update, outside keep constant for now) (Based on LSET Geometry)

        Tvec=_grainTemp2D.getLevelset();



        //Set up geometry with level set geo or with ice temp
        if (tstep == START_TEMP){
        	Lsetvec=_lset.getLevelset();			
        }
        else{
       		Lsetvec=_grainTemp2D.getLevelset(); 	
        }
        

        //Get x y dimensions of level set for convenience

        size_t xlength = 0;  //X dim over which we have ice
        size_t ylength = 0;  //X dim over which we have ice
        //Floe Extents (as a box)
        size_t minicej = 1000000; //Use high value for finding minimum
        size_t maxicej = 0;
        size_t minicei = 1000000; //Use high value for finding minimum
        size_t maxicei = 0;


        //Rectangularize floating ice to find hice and hwater
        for (size_t i = 0; i<_lset.getXdim(); i++) 
	    {   
	        size_t xind = 0; //Just to track if we should add this column       
		    for (size_t j = 0; j<_lset.getYdim(); j++)
		    {
		     	if (LsetvecCheck[(j*_lset.getXdim())+i] <= 0){     //Only if meltTemp is 0
		     	   if (j<minicej) {
		     	   	   minicej = j;
		     	   }
                   if (j>maxicej) {
		     	   	   maxicej = j;
		     	   }
                   if (i<minicei) {
		     	   	   minicei = i;
		     	   }
                   if (i>maxicei) {
		     	   	   maxicei = i;
		     	   }
		     	}
		    
		    }
	    } 
        xlength = maxicei-minicei+1;
        ylength = maxicej-minicej+1;
        //cout <<" Normal   xlength: "<< xlength <<endl;
		//cout <<" Normal   ylength: "<< ylength <<endl;
        
        if (xlength < 1) {
        	xlength = 1;
        }
        if (ylength < 1) {
        	ylength = 1;
        }
        if (xlength ==1 || ylength == 1){
	     	//_thick0=0;  //Let's use thick0 temporarily as a filter
	     	//_thick = 0.0;
	     	cout << "A Grain will melt by too small Level Set" << endl;
	        //world.RemoveGrain;
	     	//_mass=0;
	     	_remove = true;
	     	return;
        }
        else{    //VITAL VALIDATION OF GEOMETRY

			if (xlength <=7 || ylength <= 7){
		        //cout <<" Small   xlength: "<< xlength <<endl;
		        //cout <<" Small   ylength: "<< ylength <<endl;
	        }

	        //tstep=1; //-->Initialize melting or freezing (is the same)
	        //size_t tstarttemp = START_TEMP; 
	        size_t tstarttemp = dstep;

	        if (tstep == tstarttemp) {    //Only for the first step lets use the geometric level set to separate water from ICE  -->  TRY APPLY SINCE READ INPUT FIELD
		         
		         cout << "TEMPERATURE INITIALIZATION ON ICE" << endl;
		         //Save to preserve orig thickness as taper
	        	 _thick0 = _thick;

		         for (size_t ii = 0; ii<_lset.getXdim(); ii++) 
		         {            
		            for (size_t jj = 0; jj<_lset.getYdim(); jj++)
		            {

		               if (_lset.getGridValue2(ii,jj)>0.0) //If outside of the grain let's set Twater, nonsusceptible to ice heat transfer   -->  //Set thickness to neg. value outside
		               {
		                  
		                  Tvec[(jj*_lset.getXdim())+ii] = Twater;  //Melt  
		                  
		                  //cout << "Temp0: " << Tvec[(jj*_lset.getXdim())+ii] << endl;
		                  
		                  InOutvec[(jj*_lset.getXdim())+ii] = 1;     //OUTSIDE

		                  Thickvec[(jj*_lset.getXdim())+ii] = -abs(_thick); //Set non-existent ice to -0.0000001 for contrast for now, it will keep dropping and delimit grain


		                  //Tvec[(jj*_lset.getXdim())+ii] = Tair;  //Freeze
		               }
		               else{
		               	  
		               	  //cout << "Temp0: " << Tvec[(jj*_lset.getXdim())+ii] << endl;
		               	  
		               	  //Small vrs Large floes, large floes do have a cold underlying region, small do not.
		               	  if (melt_flag == 1)
		               	  {
		               	  	Tvec[(jj*_lset.getXdim())+ii] = 0.0;  //Set water under ice just to melting point  //Here locally, meltTemp is 0, but its really -1.8, just shifted
		               	  }
		               	  else
		               	  {
		               	  	Tvec[(jj*_lset.getXdim())+ii] = Twater;  //Floe too small to affect water temperature
		               	  }

		               	  InOutvec[(jj*_lset.getXdim())+ii] = -1;  //INSIDE
		               	  //Thickness inside the grain will stay the same
		               }	
		            }	
		         }
		    }

		    else if (tstep > tstarttemp && tstep <4000000000001){  //for succesive steps we use the Temperature or Heat Level Set to Keep the Ice Temp Evolving and anything above meltTemp becoming water as a constant heat bath  -->
		        //Melting
		        cout << "TEMPERATURE MODIFICATION ON ICE" << endl;
		        if (meltV > 0.0)  //For melting direction //WARNING!!!
		        {	
			        for (size_t ii = 0; ii<_lset.getXdim(); ii++) 
			        {            
			            for (size_t jj = 0; jj<_lset.getYdim(); jj++)
			            {

 			              //Temperature susceptibility for small floes	
 			              if (melt_flag == 2)
		               	  { 		
		               	  	Tvec[(jj*_lset.getXdim())+ii] = Twater; 
		               	  }	
			              //if (_lset.getGridValue2(ii,jj)>0.0) //If outside of the grain let's set Twater, nonsusceptible to ice heat transfer   -->
			              if(_grainTemp2D.getGridValue2(ii,jj)>meltTemp)  //MELT
			              //if(_grainTemp2D.getGridValue2(ii,jj)<meltTemp)  //FREEZE
			              {
			                  
			                  //Tvec[(jj*_lset.getXdim())+ii] = Twater;  //Do not keep updating for now, might be updated by external incoming temperatures
			                  InOutvec[(jj*_lset.getXdim())+ii] = 1;     //OUTSIDE
			                  //Tvec[(jj*_lset.getXdim())+ii] = Tair;  //Freeze
			              }	
			              else{
			               	  InOutvec[(jj*_lset.getXdim())+ii] = -1;  //INSIDE
			              }

				          if ( LsetvecCheck[(jj*_lset.getXdim())+ii]>0.0 && Thickvec[(jj*_lset.getXdim())+ii] > 0.0 ) //You need to redefine thickness in case LS is changed by breakage or other reasons, to avoid bugs
			              {
			               	  //cout << "Thickness in Out!!! LSval: " << LsetvecCheck[(jj*_lset.getXdim())+ii] << " Thicval: "  << Thickvec[(jj*_lset.getXdim())+ii] << endl;
			               	  //cout << "Correcting Thickness after break Temp!!!" << endl;
			               	  Thickvec[(jj*_lset.getXdim())+ii] = -abs(_thick); //Inspect Effect //or Lsetvec //Did it help?? A// Yes
				          }	
				          // if ( LsetvecCheck[(jj*_lset.getXdim())+ii]>0.0 && Thickvec[(jj*_lset.getXdim())+ii] > 0.0 ) //You need to redefine thickness in case LS is changed by breakage or other reasons, to avoid bugs
			           //    {
			           //     	  Thickvec[(jj*_lset.getXdim())+ii] = -_thick; //Inspect Effect //or Lsetvec
				          // }	


			            }	
			        }
			    }
			    //Freezing
			    else
			    {
			    	for (size_t ii = 0; ii<_lset.getXdim(); ii++) 
			        {            
			            for (size_t jj = 0; jj<_lset.getYdim(); jj++)
			            {

			               //if (_lset.getGridValue2(ii,jj)>0.0) //If outside of the grain let's set Twater, susceptible to ice heat transfer   -->
			               if(_grainTemp2D.getGridValue2(ii,jj)>meltTemp)  
			               {
			                  
			                  //Tvec[(jj*_lset.getXdim())+ii] = Twater;  /Don't update water temp since cold isotherm has to keep expanding
			                  InOutvec[(jj*_lset.getXdim())+ii] = 1;     //OUTSIDE
			                  //Tvec[(jj*_lset.getXdim())+ii] = Tair;  //Freeze
			               }	
			               else{
			               	  Tvec[(jj*_lset.getXdim())+ii] = -20.0; //Internal Ice Temp to FREEZE
			               	  //cout << "Freezing Ongoing!!!" << endl;
			               	  InOutvec[(jj*_lset.getXdim())+ii] = -1;  //INSIDE
			               }
			            }	
			        }

			    }    

		    }     
		    
		    //Still not applies
		    else if (tstep==4000000000001) {  //for succesive steps we use the Temperature or Heat Level Set to Keep the Ice Temp Evolving and anything above meltTemp becoming water as a constant heat bath  -->
		          //cout << "Cold Pulse Start step" << endl;	
		          for (size_t ii = 0; ii<_lset.getXdim(); ii++) {   //Reset loop         
		            for (size_t jj = 0; jj<_lset.getYdim(); jj++){

		               if (_grainTemp2D.getGridValue2(ii,jj)>meltTemp) //If outside of the grain let's set Twater, nonsusceptible to ice heat transfer   -->
		                {
		                  
		                   //Tvec[(jj*_lset.getXdim())+ii] = Tair;     //Separate water from Ice
		                   //Tvec[(jj*_lset.getXdim())+ii] = Twater;  //Freeze
		                }
		                else{
		                   //Tvec[(jj*_lset.getXdim())+ii] = -20.0; 	//Convert to Cold Core 
		                }

		              }	
		          }
		          //meltTemp=-20;

		         // for (size_t ii = 0; ii<_lset.getXdim(); ii++) {   //First  Cold Pulse Loop        
		         //    for (size_t jj = 0; jj<_lset.getYdim(); jj++){

		         //       if (_grainTemp2D.getGridValue2(ii,jj)<=meltTemp) //If outside of the grain let's set Tair, susceptible to ice heat transfer, If inside let's fix the temp.   -->
		         //        {
		                  
		         //           //Tvec[(jj*_lset.getXdim())+ii] = Twater;
		         //           Tvec[(jj*_lset.getXdim())+ii] = -20.0;  //Freeze
		         //        }
		         //        else{
		         //           Tvec[(jj*_lset.getXdim())+ii] = Tair; 
		         //           //cout << "Cold Pulse Start" << endl;	
		         //        }

		         //    }	
		        //}
		        //meltV = (0.00008/20)*(Tair-meltTemp)*alphaice; 
		      
		    }      
	          
			else{   //Growth Loop
				meltTemp=-20;
				//cout << "Expand Start" << endl;	
			    for (size_t ii = 0; ii<_lset.getXdim(); ii++) 
			     {            
			        for (size_t jj = 0; jj<_lset.getYdim(); jj++)
			         {

			           if (_grainTemp2D.getGridValue2(ii,jj)>=-20.1) //If outside of the grain let's set Twater, nonsusceptible to ice heat transfer   -->
			            {
			              
			               //Tvec[(jj*_lset.getXdim())+ii] = Twater;
			               //Tvec[(jj*_lset.getXdim())+ii] = -20.0;  //Freeze
			            }	
			         }	
			     }
			     //meltV = (0.00008/20)*(Tair-meltTemp)*alphaice;     
			 }
		         
	  

	         _grainTemp2D.changeLevelset(Tvec);
	         _grainTemp2D0.changeLevelset(Tvec);  //For prev. time step storage of Heat Equation 2D  -->
	         //vector<double> Tvec00=_grainTemp2D0.getLevelset();  //For prior heat comparison in order to know if we should change the level set at all -->
             vector<double> Tvec00=_lset.getLevelset(); //BEFORE BEING MODIFIED


	      //    //FDM 1-D Thermal Loop for Space (time loop uses dt)
		     // dh =_thick0/800.;  //For 80,000 we need thick0 to be at most 0.4 at the input
	      //    _Utemper0=_Utemper; 

	      //    for (int idh = 1; idh<800-1; idh++) 
	      //    {
	      //       _Utemper[idh] = _Utemper0[idh] + alphaice*(dt/(dh*dh))*(_Utemper0[idh+1]-2*_Utemper0[idh]+_Utemper0[idh-1]); //FDM scheme for 1-D diffusion for U_{t}-aU_{xx}=0 FWD t, CNT x
	      //       //_Utemper[idh] = _Utemper0[idh] + 0.005*alphaice*(dt/(dh*dh))*(_Utemper0[idh+1]-2*_Utemper0[idh]+_Utemper0[idh-1]); //FDM scheme for 1-D diffusion for U_{t}-aU_{xx}=0 FWD t, CNT x
	      //       //_Utemper[idh] = _Utemper0[idh] + ((dh*dh)/(2*alphaice))/(dh*dh)*(_Utemper0[idh+1]-2*_Utemper0[idh]+_Utemper0[idh-1]);
	      //       if (_Utemper[idh]>meltTemp)
	      //       {
	      //            meltid=idh;
	      //       }	

	      //    }


     		//cout << "Run thermal PDE" << endl;

	        //2D FDM Thermal Loop for Space (of local ocean heat grids)
	        dh = 1.;  //Space step given level set extension  //TODO, is dt also equal to one for thermal processes???
	        Tvec=_grainTemp2D.getLevelset();
	        // for (size_t ii = 1; ii<_lset.getXdim()-1; ii++) 
	        // {            
	        //     for (size_t jj = 1; jj<_lset.getYdim()-1; jj++)
	        //     {
	        //           //Newly Modified PDE in 2D
	            	  
	        //     	  Tvec[(jj*_lset.getXdim())+ii] = _grainTemp2D0.getGridValue2(ii,jj)+dstep*1.0*(Khor/100000000)*(_grainTemp2D0.getGridValue2(ii+1,jj)+_grainTemp2D0.getGridValue2(ii,jj+1)+_grainTemp2D0.getGridValue2(ii-1,jj)+_grainTemp2D0.getGridValue2(ii,jj-1)-4*_grainTemp2D0.getGridValue2(ii,jj));  //Try later:  - meltV*(_grainTemp2D0.getGridValue2(ii,jj)-meltTemp);
	                  
	        //           //KHor/1000 for 1000 steps, KHor for 100000 for 10000, what for steps 1000000 use Khor/100000000

	        //           //Tvec[(jj*_lset.getXdim())+ii] = _grainTemp2D0.getGridValue2(ii,jj)+(-meltV)*(dh/(alphaice))*(_grainTemp2D0.getGridValue2(ii+1,jj)+_grainTemp2D0.getGridValue2(ii,jj+1)+_grainTemp2D0.getGridValue2(ii-1,jj)+_grainTemp2D0.getGridValue2(ii,jj-1)-4*_grainTemp2D0.getGridValue2(ii,jj));
	        //           //Tvec[(jj*_lset.getXdim())+ii] = _grainTemp2D0.getGridValue2(ii,jj)+(-meltV*50)*(dh/(2*alphaice))*(_grainTemp2D0.getGridValue2(ii+1,jj)+_grainTemp2D0.getGridValue2(ii,jj+1)+_grainTemp2D0.getGridValue2(ii-1,jj)+_grainTemp2D0.getGridValue2(ii,jj-1)-4*_grainTemp2D0.getGridValue2(ii,jj));
         //    		  //cout << "Temp: " << Tvec[(jj*_lset.getXdim())+ii] << endl;
	        //     }	
	        // }


	        //Update Local Temp Grid using Global Ocean Temp Grid, then this Tvec and its LS interact with thickness
	        double lowLx, lowLy, lowLx0, lowLy0, lslocx, lslocy, lslocx0, lslocy0; //Define lower left corner of level set local grid and current lslocx, lslocy points (April 12, 2023)
	        double szx = _lset.getXdim(); //LS size in x direction
	        double szy = _lset.getYdim(); //LS size in y direction 
	        Vector2d pointOcxy; //Shifting point for each cell starting from lower left corner of LS
	        double Tempval; //Value to use in Temp Local Grid interpolated from Global Grid
	        //Shift to global coordinates and move to left bottom corner of level set (April 12, 2023)
	        size_t rotLS = 1; //0 not rotate, 1 or anything else rotate matrix (keep 0 as control when needed) (April 12, 2023)
	        
	        
	        //NO ROTATION (April 12, 2023)
	        if (rotLS == 0){
	            lowLx = _position(0) - 0.5*szx +1.0; //IFF LS has unit of 1 for grid cells  //ADD ROTATION HERE????????
	            lowLy = _position(1) - 0.5*szy +1.0; //IFF LS has unit of 1 for grid cells  //ADD ROTATION HERE????????
	        }
	        //WITH ROTATION (April 12, 2023)
	        else{
	            lowLx0 = - 0.5*szx +1.0; //IFF LS has unit of 1 for grid cells  //ADD ROTATION HERE???????? Equivalent cmx
	            lowLy0 = - 0.5*szy +1.0; //IFF LS has unit of 1 for grid cells  //ADD ROTATION HERE???????? Equivalent cmy
	            //lowLx = ( lowLx0 * cos(_theta) - lowLy0 * sin(_theta) ) + _position(0); 
	            //lowLy = ( lowLx0 * sin(_theta) + lowLy0 * cos(_theta) ) + _position(1); 
	        }
	        
	        //Is this method correct??? From LS coords or ref to global, vice versa uses rotMatrix << cos(_theta), sin(_theta), -sin(_theta), cos(_theta);
// 			rotMatrix << cos(_theta), -sin(_theta), sin(_theta), cos(_theta);

// 			for (size_t i = 0; i < _pointList.size(); i++) {
// 				_pointList[i] = rotMatrix*(_pointList[i] - _cmLset) + _position;
// 			}
	        //cout << "Rounding Ocean Temperature based on proximity to global grid cell" << endl;
	        
	        for (size_t ii = 0; ii<_lset.getXdim(); ii++)  //for (size_t ii = 1; ii<_lset.getXdim()-1; ii++) 
	        {            
	            for (size_t jj = 0; jj<_lset.getYdim(); jj++) //for (size_t jj = 1; jj<_lset.getYdim()-1; jj++)
	            {
	            	
	            	//NO ROTATION (April 12, 2023)
	            	if (rotLS == 0){
	                   lslocx = lowLx + 0.5 + ii; //IFF LS has unit of 1 for grid cells  
	            	   lslocy = lowLy + 0.5 + jj; //IFF LS has unit of 1 for grid cells  
	            	}
	            	//WITH ROTATION (April 12, 2023)
	            	else{
	            	    lslocx0 = lowLx0 + 0.5 + ii; //IFF LS has unit of 1 for grid cells  //ADD ROTATION HERE????????
    	            	lslocy0 = lowLy0 + 0.5 + jj; //IFF LS has unit of 1 for grid cells  //ADD ROTATION HERE????????
    	                lslocx = ( lslocx0 * cos(_theta) - lslocy0 * sin(_theta) ) + _position(0); 
    	                lslocy = ( lslocx0 * sin(_theta) + lslocy0 * cos(_theta) ) + _position(1); 
	            	}
	            	
	            	pointOcxy << lslocx , lslocy; //POint of LS over Global Grid   //ADD ROTATION HERE????????    
	            	//cout << "pointOcxy: " << lslocx << " , " << lslocy << endl;

	               	//For Debugging Purposes, remove grain later on and abort 
				    if ( isnan(lslocx) || isnan(lslocy) )
				    {
				      cout << "Ill-positioned Grain, remove" << endl;
				      _mass = 0.0;
				      return;
				    }

	            	//Bilinear interpolation of this point on Ocean Grid and get temp value (how fast will this be?)
	            	//Tempval = bilinear_Interp_Ocean(pointOcxy, fluid_coord, oceanTemp, x_cells, y_cells, offset);
	            	Tempval = round_Ocean(pointOcxy, oceanTemp, x_cells, y_cells, offset);


	            	//Update each grid cell with this temperature, then you just call this grid when updating thickness
	            	Tvec[(jj*_lset.getXdim())+ii] = Tempval+1.8; //Don't forget using sea ice melt shift for convenience
	            }
	        }



             //Given Dirichlet BCs, borders Tvec have to be updated in other ways

	        _grainTemp2D.changeLevelset(Tvec); //How can we raise temperature more?????
	        //_lset.changeLevelset(Tvec); //Use Temp not +1 -1
	        //_grainTemp2D0.changeLevelset(Tvec0);  //Store for next. time step   -->

	        Thickvec0 = Thickvec; //Get prior step to get change

	        //Find number of out of count elements
	        //cout << "Count Thickness > 0 cells" << endl;
	        size_t orig_outcount = 0;
	        size_t new_outcount = 0;
	        for (size_t ii = 1; ii<_lset.getXdim()-1; ii++) 
	        {            
	            for (size_t jj = 1; jj<_lset.getYdim()-1; jj++)
	            {
	            	if (Thickvec0[(jj*_lset.getXdim())+ii] < 0.0 )
	            	{
	            		orig_outcount++;
	            	}
            	}
        	}

        	//cout << "THICKNESS UPDATE ON ICE" << endl;
        	double ct_hlat = 0.0;
        	double ave_thick = 0;
        	size_t ct_ave_thick = 0;  
        	double ave_solar = 0;
        	double ave_ocean = 0;
        	double ave_solarT = 0;
        	double ave_oceanT = 0;
	        //Use Tvec temperature to process thickness changes
	        for (size_t ii = 1; ii<_lset.getXdim()-1; ii++) 
	        {            
	            for (size_t jj = 1; jj<_lset.getYdim()-1; jj++)
	            {
	            	//Update with 1D diffusion equation in vertical direction //Newly Modified PDE in 1D
	            	
	            	//Older expression
	            	//Thickvec[(jj*_lset.getXdim())+ii] = Thickvec0[(jj*_lset.getXdim())+ii] + dstep*(meltV)*(meltTemp-Tvec[(jj*_lset.getXdim())+ii]) + dstep*meltVSun ; //1.0 okay for 1000 and 10000 what for 1000000
	                
	            	//With thickness average for Fines
	            	Thickvec[(jj*_lset.getXdim())+ii] = Thickvec0[(jj*_lset.getXdim())+ii] + dstep*(meltV)*(meltTemp-Tvec[(jj*_lset.getXdim())+ii]) + dstep*meltVSun  ;
                    //Update Thickness                    //Prior Thickness                 //Denom dt qv / pi*Lf   (Tf - T)                        //Denom (dt / pi*Lf) * (-Q(1-a_i)+(A+B*Tice))

	                //cout << "Temperature Diference: " << (meltTemp-Tvec[(jj*_lset.getXdim())+ii]) << endl;
	                //cout << "Ice Temperature: " << (meltTemp) << endl;
	                //cout << "Ocean Temperature: " << (Tvec[(jj*_lset.getXdim())+ii]) << endl;

	                if (Thickvec[(jj*_lset.getXdim())+ii] <= 0.0)  //Too low or high thick is not ok
	                //if (Thickvec[(jj*_lset.getXdim())+ii] <= 0.0 || Thickvec[(jj*_lset.getXdim())+ii] > 2 )  //Too low or high thick is not ok
	            	{
	            		new_outcount++;
	            		if (Thickvec0[(jj*_lset.getXdim())+ii] > 0.0){
	            		    ct_hlat += 1.0;
	            		    hlat += Thickvec0[(jj*_lset.getXdim())+ii];
	            		}
	            	}
	            	else
	            	{
	            		ave_thick += Thickvec[(jj*_lset.getXdim())+ii];
	            		ct_ave_thick++;
	            		ave_solar += abs(meltVSun);
                        ave_ocean += abs((meltV)*(meltTemp-Tvec[(jj*_lset.getXdim())+ii]));
	            	}

	            }

            }
            _Mthick.changeLevelset(Thickvec); //Update LS vector in Mthick to use later in Reinit2D
            
            //Get average hlat thickness
            if (ct_hlat > 0.0){
                hlat /= ct_hlat;
            }
            else{
                hlat = 0.0;
            }
            
            //Change Aug 22, 2022
            //BEGIN
            //Obtain area for mass loss calculation using points
            vector<Vector2d> VecOrigin = _pointList;
            double AreaMass = 0.0;
            size_t nmass = VecOrigin.size();
            
            for (size_t ii = 0; ii < nmass-1; ii++)
            {
                AreaMass += ( VecOrigin[ii](0) * VecOrigin[ii+1](1) -  VecOrigin[ii](1) * VecOrigin[ii+1](0) ); 
            }
            AreaMass += (VecOrigin[nmass-1](0) * VecOrigin[0](1) -  VecOrigin[nmass-1](1) * VecOrigin[0](0) ); 
            
            if (isnan(AreaMass) == 1)
            {
                cout << "WARNING: Nan in TempState algorithm for inner grain!!" << endl;
                AreaMass = 1.0; //Just to debug
            }
            
            double AreaMeltMass = 0.5*abs(AreaMass); //Area from points to calculate floe melt (using original pointlist)
            double prev_thick = _thick;
            //END
            

            //Update average thickness
            double mult_loss = 1.0;
            double frac_solar = 0.0;
            double frac_ocean = 0.0;
            
            
            if (ct_ave_thick > 0)
            {
            	_thick = ave_thick / double(ct_ave_thick);
            	ave_solarT = ave_solar / double(ct_ave_thick);
            	ave_oceanT = ave_ocean / double(ct_ave_thick);
            	//cout << "Average Thickness Count: " << ct_ave_thick << endl;
            	//cout << "Average Thickness Sum: " << ave_thick << endl;
            	//cout << "Thickness: " << _thick << endl;
            	
            	 if ( (ave_solarT + ave_oceanT) > 0){
                    frac_solar = (ave_solarT)/(ave_solarT + ave_oceanT);
                    frac_ocean = (ave_oceanT)/(ave_solarT + ave_oceanT);
                }
                else{
                    frac_solar = 0.0;
                    frac_ocean = 0.0;
                }


            	//Update mass relative to thick0  !!!NEW CHANGE!! TODO: CHECK!!!!
            	cout << "New thick: " << _thick << " old thickness: " << _thick0 << endl;
            	mult_loss = _thick/_thick0;
            	//mult_loss = 1.0; 
            	_mass *= mult_loss;
            	
            	//Change Aug 22, 2022
                //BEGIN
                loss_mcv_temp += AreaMeltMass * max( (prev_thick - _thick) , 0.0) * 0.001 * _density; 
                loss_mcv_solar_temp += frac_solar * AreaMeltMass * max( (prev_thick - _thick) , 0.0) * 0.001 * _density; 
                loss_mcv_ocean_temp += frac_ocean * AreaMeltMass * max( (prev_thick - _thick) , 0.0) * 0.001 * _density; 
                //END
                //cout << "melt component: " << loss_mcv_temp << endl;
                //cout << "S melt component: " << loss_mcv_solar_temp << endl;
                //cout << "O melt component: " << loss_mcv_ocean_temp << endl;
            }
            else
            {
            	//Change Aug 22, 2022
                //BEGIN
                loss_mcv_temp += AreaMeltMass * max(_thick, 0.0) * 0.001 * _density; 
                loss_mcv_solar_temp += 0.5 * AreaMeltMass * max(_thick, 0.0) * 0.001 * _density; 
                loss_mcv_ocean_temp += 0.5 * AreaMeltMass * max(_thick, 0.0) * 0.001 * _density; 
                //END
            	
            	cout << "WARNING: IN Cells: " << ct_ave_thick << endl;
            	_thick = 0;
            	_mass = 0.0; //Here if we get NO internal cells or cells with positive thickness it's understood thickness is zero and so is mass
            	ave_solar = 0;
            	ave_ocean = 0;
            	return; //We need to leave now, further action is pointless
            }

            //Now we can use Thickvec to update geometry if it happens at every step or at certain steps


	        //Can we use the level set as temperature instead of lset from Now on???????   -->
	        //_lset.changeLevelset(Tvec);


	        //2D Update if a part of the ice will melt
	        int changeshape=0;

            //Define if shape actually changes and more cells are out than before due to melting
            if (new_outcount > orig_outcount)
            {
            	changeshape = 1;
            	//_mass *= mult_loss; //Will be down
            }
            else
            {
                //_mass *= mult_loss;
            }


            //TODO: Inspect use of this section
	        Tvec=_grainTemp2D.getLevelset();    //_lset.changeLevelset(Tvec);   ?????? -->
	        for (size_t ii = 0; ii<_lset.getXdim(); ii++) 
	        {            
	            for (size_t jj = 0; jj<_lset.getYdim(); jj++)
	            {

	               //cout << "Check change shape" << endl;
	               //If outside of the grain let's set Twater, nonsusceptible to ice heat transfer  //////????????????????? -->   Bound in case we keep Twater stable to only find susceptible value to melt ice ???
	               

	               //if (_grainTemp2D.getGridValue2(ii,jj)>=meltTemp-0.01 &&  _grainTemp2D.getGridValue2(ii,jj)<=meltTemp+0.01 )    //MELT

	               if (changeshape = 1) //THICKNESS CHANGE 

	               //if (_grainTemp2D.getGridValue2(ii,jj)<=meltTemp+0.1 &&  _grainTemp2D.getGridValue2(ii,jj)>=meltTemp-0.1 )    //FREEZE
	               {
	                  //cout << "Shape change" << endl;
	                  changeshape=1;  //If this happens we must modify LS, mass, radius, cm, Inertia and Points for each Grain
	                  //Tvec[(jj*_lset.getXdim())+ii] = Twater; //How to alter this Temp matrix in our favor???
	               }	
	               if ( InOutvec[(jj*_lset.getXdim())+ii] < 0 &&  Tvec[(jj*_lset.getXdim())+ii] >=meltTemp){
	                    

	                    InOutvec[(jj*_lset.getXdim())+ii] = 1;
	                    
	                }
	            }	
	        }
	         //_grainTemp2D.changeLevelset(Tvec);
	         //_grainTemp2D0.changeLevelset(Tvec);

            //THIS IS WRONG. LEVEL SETS CANNOT BE PURE BINARY, BUT DISTANCE-BASED FUNCTIONS, USE RE-INIT based on temp or thickness
	        //Update Level Set with Help of Temporal InOut Vector
	        vector<double> Lvec00=_lset.getLevelset();  //For prior geometric comparison in order to know if we should change the level set at all -->
	        for (size_t ii = 0; ii<_lset.getXdim(); ii++) 
	        {            
	            for (size_t jj = 0; jj<_lset.getYdim(); jj++)
	            {
	                //Lsetvec[(jj*_lset.getXdim())+ii] = InOutvec[(jj*_lset.getXdim())+ii];
	            }
	        }
	        //_lset.changeLevelset(Lsetvec); //Now our new Level Set (altered because of the heat equation becomes our geometric regent and then we use to update Dirichlet BC)    	
	 

	        //UPDATE _fracLoc here !!!!!!!!!!
	        //cout <<"Get Frac Loc"<<endl;
	        //Let's avoid breaking grains that are too small
	        if (xlength<=5 || ylength <=5){
               _fracLoc = 0;
	        }
	        else{
	           _fracLoc = minicei+xlength*0.5;
	           //_fracLoc = max(minicei+xlength*0.5, minicej+ylength*0.5);  //IN LONGER AXIS (FLEXURE)
	        }

            //melt_flag = 0; //Temporary for debugging!!! 
	        if (changeshape==1 && melt_flag == 1){ //Change shape just for big grains, which will be a very small change, otherwise only thickness is reduced
	        //if (changeshape==1){ 
	        	//cout << "Exploration of shape change" << endl;

		        //2D Changes for LS *********************** *********************** ****************************** ***********************************
                //Replaced using Reinit Function
		  //       //Modify Geometry if we melt any ice or if any cell of grainTemp2D>=Tmelt, otherwise do nothing  if(changeshape=1){}
		  //       //Assume Rows are Y and Columns are X, we will use Matrix[j][i] or y,x or rows, columns

		  //       //Generate Matrix for Level Set Temperature for Convenience Sake
				// Matrixx MatrixTLS;
				// const size_t N = _lset.getYdim(); //Rows j is y
				// const size_t M = _lset.getXdim(); //Columns i is x
		  //       Tvec=_grainTemp2D.getLevelset();
				// for(size_t j = 0; j < N; ++j) //y
				// {
				//     Row row(M);
				//     for(size_t i = 0; i < M; ++i) //x
				//     {
				//         row[i] = Tvec[(j*_lset.getXdim())+i];
				//     }
				//     MatrixTLS.push_back(row); // push each row after you fill it
				// }
				
			 //    //Generate Matrix for Level Set for Convenience Sake
				// Matrixx MatrixLS;
				// //const size_t N = _lset.getYdim(); //Rows j is y
				// //const size_t M = _lset.getXdim(); //Columns i is x

		  //       Tvec=_lset.getLevelset();
				// for(size_t j = 0; j < N; ++j) //y
				// {
				//     Row row(M);

				//     for(size_t i = 0; i < M; ++i) //x
				//     {
				//         row[i] = Tvec[(j*_lset.getXdim())+i];
				//     }

				//     MatrixLS.push_back(row); // push each row after you fill it
				// }


		  //       //Space Grid Resolution (use half to average left and right)
		  //       double dxinv= _lset.getXdim()*0.5;
		  //       double dyinv= _lset.getYdim()*0.5;

		  //       //Obtain Ghost Matrices from LS in order to do Finite Difference Scheme (1st Fwd)
		  //       Matrixx GhostLSX;
		  //       Matrixx GhostLSY;
				
				// const size_t Nx = _lset.getYdim(); const size_t Ny = _lset.getYdim()+2; //Rows j is y
				// const size_t Mx = _lset.getXdim()+2; const size_t My = _lset.getXdim(); //Columns i is x
				
				// for(size_t j = 0; j < Nx; ++j) //y
				// {
				//     Row row(Mx);
				//     for(size_t i = 0; i < Mx; ++i) //x
				//     {
				//     	if (i==0) 
				//     	{
				//     	     row[i] = MatrixLS[j][i+2];	
				//     	}  else if (i==Mx-1) {
		  //                    row[i] = MatrixLS[j][i-2];	
				//     	} else {
				//              row[i] = MatrixLS[j][i-1];
				//         }
				//     }		    
				//     GhostLSX.push_back(row); // push each row after you fill it
				// }

			 //    for(size_t j = 0; j < Ny; ++j) //y
				// {
				//     Row row(My);
				//     for(size_t i = 0; i < My; ++i) //x
				//     {
		  //               if (j==0) 
				//     	{
				//     	     row[i] = MatrixLS[j+2][i];	
				//     	}  else if (j==Ny-1) {
		  //                    row[i] = MatrixLS[j-2][i];	
				//     	} else {
				//              row[i] = MatrixLS[j-1][i];
				//         }
				//     }
				//     GhostLSY.push_back(row); // push each row after you fill it
				// }

				// //Obtain Central Derivative Matrixes
		  //       Matrixx DerivX;
		  //       Matrixx DerivY;
				
				// const size_t Nxx = _lset.getYdim(); const size_t Nyy = _lset.getYdim()+1; //Rows j is y
				// const size_t Mxx = _lset.getXdim()+1; const size_t Myy = _lset.getXdim(); //Columns i is x
				
				// for(size_t j = 0; j < Nxx; ++j) //y
				// {
				//     Row row(Mxx);
				//     for(size_t i = 0; i < Mxx; ++i) //x
				//     {
				//     	row[i] = dxinv*(GhostLSX[j][i+1]-GhostLSX[j][i]);	
				//     }		    
				//     DerivX.push_back(row); // push each row after you fill it
				// }

				// for(size_t j = 0; j < Nyy; ++j) //y
				// {
				//     Row row(Myy);
				//     for(size_t i = 0; i < Myy; ++i) //x
				//     {
			 //            row[i] = dxinv*(GhostLSY[j+1][i]-GhostLSY[j][i]);	   //Check how to control this magnitude!!!!?????
				//     }
				//     DerivY.push_back(row); // push each row after you fill it
				// }


		  //       //Obtain Left and Right Derivative Matrixes
		  //       Matrixx DerivLX; Matrixx DerivRX;
		  //       Matrixx DerivLY; Matrixx DerivRY;
				// //const size_t N = _lset.getYdim(); //Rows j is y
				// //const size_t M = _lset.getXdim(); //Columns i is x

		  //       for(size_t j = 0; j < _lset.getYdim(); ++j) //y
				// {
				//     Row row(M);
				//     Row row2(M);
				//     for(size_t ii = 0; ii < _lset.getXdim(); ++ii) //x
				//     {
			 //    	     row[ii] = DerivX[j][ii];	
			 //    	     row2[ii] = DerivX[j][ii+1];	
				//     }
				//     DerivLX.push_back(row); // push each row after you fill it
				//     DerivRX.push_back(row2); // push each row after you fill it
				// }

				// for(size_t j = 0; j < N; ++j) //y
				// {
				//     Row row(M);
				//     Row row2(M);
				//     for(size_t i = 0; i < M; ++i) //x
				//     {
				//     {
			 //    	     row[i] = DerivY[j][i];	
			 //    	     row2[i] = DerivY[j+1][i];	
				//     }
				//     }
				//     DerivLY.push_back(row); // push each row after you fill it
				//     DerivRY.push_back(row2); // push each row after you fill it
				// }


		  //       //Obtain Left and Right Product Matrixes
		  //       Matrixx ProdLX; Matrixx ProdRX;
		  //       Matrixx ProdLY; Matrixx ProdRY;
				// //const size_t N = _lset.getYdim(); //Rows j is y
				// //const size_t M = _lset.getXdim(); //Columns i is x
		        
		  //       for(size_t j = 0; j < N; ++j) //y
				// {
				//     Row row(M);
				//     Row row2(M);
				//     for(size_t i = 0; i < M; ++i) //x
				//     {
			 //    	     row[i] = meltV*DerivLX[j][i];	
			 //    	     row2[i] = meltV*DerivRX[j][i];	
				//     }
				//     ProdLX.push_back(row); // push each row after you fill it
				//     ProdRX.push_back(row2); // push each row after you fill it
				// }

				// for(size_t j = 0; j < N; ++j) //y
				// {
				//     Row row(M);
				//     Row row2(M);
				//     for(size_t i = 0; i < M; ++i) //x
				//     {
			 //    	     row[i] = meltV*DerivLY[j][i];	
			 //    	     row2[i] = meltV*DerivRY[j][i];	
				//     }
				//     ProdLY.push_back(row); // push each row after you fill it
				//     ProdRY.push_back(row2); // push each row after you fill it
				// }



		  //       //Obtain Left and Right Magnitude Matrixes
		  //       Matrixx MagLX; Matrixx MagRX;
		  //       Matrixx MagLY; Matrixx MagRY;
				// //const size_t N = _lset.getYdim(); //Rows j is y
				// //const size_t M = _lset.getXdim(); //Columns i is x

				// for(size_t j = 0; j < N; ++j) //y
				// {
				//     Row row(M);
				//     Row row2(M);
				//     for(size_t i = 0; i < M; ++i) //x
				//     {
			 //    	     row[i] = abs(ProdLX[j][i]);	
			 //    	     row2[i] = abs(ProdRX[j][i]);
				//     }
				//     MagLX.push_back(row); // push each row after you fill it
				//     MagRX.push_back(row2); // push each row after you fill it
				// }

				// for(size_t j = 0; j < N; ++j) //y
				// {
				//     Row row(M);
				//     Row row2(M);
				//     for(size_t i = 0; i < M; ++i) //x
				//     {
			 //    	     row[i] = abs(ProdLY[j][i]);	
			 //    	     row2[i] = abs(ProdRY[j][i]);
				//     }
				//     MagLY.push_back(row); // push each row after you fill it
				//     MagRY.push_back(row2); // push each row after you fill it
				// }


		  //       //Obtain Left and Right Flow Bolean Matrixes
		  //       Matrixx FlowLX; Matrixx FlowRX;
		  //       Matrixx FlowLY; Matrixx FlowRY;
				// //const size_t N = _lset.getYdim(); //Rows j is y
				// //const size_t M = _lset.getXdim(); //Columns i is x

		  //       for(size_t j = 0; j < N; ++j) //y
				// {
				//     Row row(M);
				//     Row row2(M);
				//     for(size_t i = 0; i < M; ++i) //x
				//     {
		  //               if ( ( ProdLX[j][i]>=0 && ProdRX[j][i]>=0 ) || ( ProdLX[j][i]>=0 && ProdRX[j][i]<=0 && (MagLX[j][i]>=MagRX[j][i]) )  )
				//     	{
				//     	     row[i] = 1;	
				//     	}  else {
				//              row[i] = 0;
				//         }  
				//         if ( ( ProdLX[j][i]<=0 && ProdRX[j][i]<=0 ) || ( ProdLX[j][i]>=0 && ProdRX[j][i]<=0 && (MagLX[j][i]<MagRX[j][i]) )  )
				//     	{
				//     	     row2[i] = 1;	
				//     	}  else {
				//              row2[i] = 0;
				//         } 
		  //           }
				//     FlowLX.push_back(row); // push each row after you fill it
				//     FlowRX.push_back(row2); // push each row after you fill it
				// }

			 //    for(size_t j = 0; j < N; ++j) //y
				// {
				//     Row row(M);
				//     Row row2(M);
				//     for(size_t i = 0; i < M; ++i) //x
				//     {
			 //            if ( ( ProdLY[j][i]>=0 && ProdRY[j][i]>=0 ) || ( ProdLY[j][i]>=0 && ProdRY[j][i]<=0 && (MagLY[j][i]>=MagRY[j][i]) )  )
				//     	{
				//     	     row[i] = 1;	
				//     	}  else {
				//              row[i] = 0;
				//         }  
				//         if ( ( ProdLY[j][i]<=0 && ProdRY[j][i]<=0 ) || ( ProdLY[j][i]>=0 && ProdRY[j][i]<=0 && (MagLY[j][i]<MagRY[j][i]) )  )
				//     	{
				//     	     row2[i] = 1;	
				//     	}  else {
				//              row2[i] = 0;
				//         } 
			 //        }
				//     FlowLY.push_back(row); // push each row after you fill it
				//     FlowRY.push_back(row2); // push each row after you fill it
				// }



		  //       //Accumulate Magnitude
		  //       Matrixx MagSum;
		  //       //const size_t N = _lset.getYdim(); //Rows j is y
				// //const size_t M = _lset.getXdim(); //Columns i is x

				// for(size_t j = 0; j < N; ++j) //y
				// {
				//     Row row(M);

				//     for(size_t i = 0; i < M; ++i) //x
				//     {
				//         row[i] = ((DerivLX[j][i]*DerivLX[j][i])*FlowLX[j][i]+(DerivRX[j][i]*DerivRX[j][i])*FlowRX[j][i])+((DerivLY[j][i]*DerivLY[j][i])*FlowLY[j][i]+(DerivRY[j][i]*DerivRY[j][i])*FlowRY[j][i]);
				//     }

				//     MagSum.push_back(row); // push each row after you fill it
				// }
				
				// //Sqrt. Magnitude
		  //       Matrixx MagLS;
		  //       //const size_t N = _lset.getYdim(); //Rows j is y
				// //const size_t M = _lset.getXdim(); //Columns i is x
				// for(size_t j = 0; j < N; ++j) //y
				// {
				//     Row row(M);

				//     for(size_t i = 0; i < M; ++i) //x
				//     {
				//         row[i] = sqrt(MagSum[j][i]);
				//     }

				//     MagLS.push_back(row); // push each row after you fill it
				// }

				// //Obtain Change for LS
				// Matrixx DeltaLS;
				// //const size_t N = _lset.getYdim(); //Rows j is y
				// //const size_t M = _lset.getXdim(); //Columns i is x

				// for(size_t j = 0; j < N; ++j) //y
				// {
				//     Row row(M);
				//     for(size_t i = 0; i < M; ++i) //x
				//     {
				//         row[i] = -meltV*MagLS[j][i];
				//     }
				//     DeltaLS.push_back(row); // push each row after you fill it
				// }

				// //Perhaps we must modify meltV in terms of matrix that depends on the 2D heat propagation and then resize all the level set, based on the reduction of melt edges to become >0 but spread to all.  
				// //Implement a row[i] = -meltV*TempM[j][i]*MagLS[j][i] for DeltaLS; We would have direct Temperature Dependence and Geometric Effect

		  //       //Implement change over Real Level Set
		  //       Tvec0=_lset.getLevelset();   //Original
		  //       Tvec=_lset.getLevelset();    //To modify
		  //       for (size_t jj = 0; jj<_lset.getYdim(); jj++) 
		  //        {            
		  //           for (size_t ii = 0; ii<_lset.getXdim(); ii++)
		  //           {
		  //                 Tvec[(jj*_lset.getXdim())+ii]=Tvec0[(jj*_lset.getXdim())+ii]+DeltaLS[jj][ii];   //Melt the grain once a certain value below zero dissapears like -0.5 or -1 if function FIND CONTOUR
		  //                 //Tvec[(jj*_lset.getXdim())+ii]=FlowRY[jj][ii]*Tvec0[(jj*_lset.getXdim())+ii]*0;
		  //                 //Tvec[(jj*_lset.getXdim())+ii]=MatrixLS[jj][ii]+0.001;  //Summing is GOOD
		 
		  //           }	

		  //       }
		  //       _lset.changeLevelset(Tvec);   ///Update New Level Set Based on Normal Gradient Parameters





		        //BEGIN REINIT ATTEMPT
		        //cout << "Perform Reinit" << endl;
                //vector<double> Tvecsm = _lset.getLevelset();  //Use advanced new geometry and reinit
                //vector<double> Tvecsm = _grainTemp2D.getLevelset();  //Or Use Temp as zero contour and the Reinit geometry with it
                vector<double> Tvecsm = _Mthick.getLevelset();  //Or Thickvec

                ArrayXd Tvecsmooth(_lset.getYdim()*_lset.getXdim());

                for (size_t j = 0; j<_lset.getYdim(); j++) 
				{	                 			                 
		            for (size_t ii = 0; ii<_lset.getXdim(); ii++)
		            {
		            	Tvecsmooth((j*_lset.getXdim())+ii) = -Tvecsm[(j*_lset.getXdim())+ii];   //Need to invert if you are using positive Thickness for this Reinit
		            }
		        }    
                        
	                //WARNING-REINCLUDE LSM-LIB	
		        //Essential for lateral melt. Needs LSM-LIB, otherwise use commented block above 
		        reinit2D(Tvecsmooth, _lset.getXdim(), _lset.getYdim());

		        for (size_t jj = 0; jj<_lset.getYdim(); jj++) 
				{	                 			                
		            for (size_t ii = 0; ii<_lset.getXdim(); ii++)
		            {
		            	Tvecsm[(jj*_lset.getXdim())+ii] = Tvecsmooth((jj*_lset.getXdim())+ii);
		            }
		        } 
		       	_lset.changeLevelset(Tvecsm);

		        //END REINIT ATTEMPT
		        //cout << "End Reinit" << endl;




		        //_temper=DeltaLS[20][20];


		        //Link Temp Gradient to Melt Normal Gradient (Second to Last Part)
		        //Full Melt  (Last Part)

		       

		        //Modification of Properties Based on Changes for LS *********************** *********************** ****************************** ***********************************
		        //Create change functions if necessary (URGENT)

		        int changeind=0;  //Validate a significant shape in the level set
		        //double tolerance=0.01; //Include a tolerance to apply only relevant changes to save resources
	            
	            Tvec=_lset.getLevelset();     // declare Tvec=_lset.getLevelset(); if you are going to work with geometric level set instead of heat level set -->
		        //Tvec = _grainTemp2D.getLevelset();  //CANT USE TEMP AS A LEVEL SET FOR GEOMETRY
		        for (size_t jj = 0; jj<_lset.getYdim(); jj++) 
		         {            
		            for (size_t ii = 0; ii<_lset.getXdim(); ii++)
		            {
		               if (Tvec[(jj*_lset.getXdim())+ii] == Tvec00[(jj*_lset.getXdim())+ii]) {       // --> Tvec0 is previous level set from above, Tvec00 is previous geo/heat level set  
		               }
		               else {
		               	changeind=1;  //If at least one cell is different the points should evolve
		               	break;
		               }
		            }	
		         }

	            int zerocounter=0;
	            //Skip point update if no zeros level is detected
		         for (size_t jj = 0; jj<_lset.getYdim(); jj++) 
		         {            
		            for (size_t ii = 0; ii<_lset.getXdim(); ii++)
		            {

	                    //if (meltTemp>=0) {

			               if (Tvec[(jj*_lset.getXdim())+ii]<0.0) //GEO //Make sure meltTemp is cero for geometric level sets
			               //if (Tvec[(jj*_lset.getXdim())+ii]<meltTemp) //MELT //Make sure meltTemp is cero for geometric level sets
		            	   //if (Tvec[(jj*_lset.getXdim())+ii]>meltTemp) //FREEZE
			               {  //Use 0.00000001 for Level Set, Use MeltTemp for Heat Level Set -->
		                      zerocounter++;
			               }

			            // }
			            // else if (meltTemp>=-10){    ///FOR FREEZING CASE ONLY MUST MODIFY IF ICE IS MELTING -->

	              //             if (Tvec[(jj*_lset.getXdim())+ii]<=meltTemp) {  //Use 0.00000001 for Level Set, Use MeltTemp for Heat Level Set -->
		             //          zerocounter++;
			            //    }

			            // } 
	              //       else {    ///FOR FREEZING CASE ONLY MUST MODIFY IF ICE IS MELTING -->
	                         
	              //             if (Tvec[(jj*_lset.getXdim())+ii]=meltTemp) {  //Use 0.00000001 for Level Set, Use MeltTemp for Heat Level Set -->
		             //          zerocounter++;
			               
			            //    }

			            // } 


		            }	
		        }


		        size_t zeroc;
		        Levelset2d glset = _lset;
		    	vector<double> gvec = glset.getLevelset();		
		        double damping= 0.3;

				double LSLimit = 0; //RECENT CHANGE

				vector<Vector2d> zeroCrossings;
				for (size_t j = 0; j < glset.getYdim()-1; j++) {
					for (size_t i = 0; i < glset.getXdim()-1; i++) {
						double val = gvec[(j*glset.getXdim())+i];
						if (val*gvec[((j+1)*glset.getXdim())+i] < (LSLimit)   ) { 
							zeroCrossings.push_back(Vector2d( double(i), double(j)+fabs(val)/(fabs(val) + fabs(gvec[((j+1)*glset.getXdim())+i]))  ));
						}
						if (val*gvec[(j*glset.getXdim())+(i+1)] < (LSLimit)   ) {	
							zeroCrossings.push_back(Vector2d(  double(i)+fabs(val)/(fabs(val) + fabs(gvec[(j*glset.getXdim())+(i+1)])), double(j)  ));
						}
					}
				}	

			    zeroc = zeroCrossings.size();

                size_t minsize = 9; //Minimum Size Before Melting is Preferrable because too small is not convenient for matrix advancement TODO: Check if only do vertical melting
		        if (zerocounter < minsize  || _mass==0.0 || zeroc < 2){    //  || _mass==0.0){    //Is Mass 0 necessary if already erased   //Generate condition to delete grain if our level set has no more value below 0 or below MeltTemp as there is no more ice
		         	 changeind=0;
		         	 //_thick0=0;  //Let's use thick0 temporarily as a filter
		         	 //_thick = 0.0;
		         	 cout << "A Grain has melted due small zero crossings, will go to fines" << endl;
		         	 cout << "zerocounter: " << zerocounter << " mass: " << _mass << " zeroc: " << zeroc << endl;
	                 //world.RemoveGrain;
		         	 //_mass = 0.000;
		         	 _remove = true;
		         	 changeind = 0; //Avoid generating new geometry since we will erase anyway.
		         	 //_radius=0;
		         	 return;
		         }



		        //if (tstep % START_TEMP == 0 && changeind == 1){  //TIME STEP Control of changing shape (for efficiency)
		        if (changeind==1){ //Implement the changes for each single grain, skip if the grain becomes too small to avoid Core Segmentation Faults
		        //cout <<"Really Changing Shape!!!"<<endl;
			        //Change in Mass due to Change of Level Set
			        
			        // compute mass and center of mass (assuming a density of 1)
					double vol = 0;
					//double eps = 1.5+meltTemp;    //tolerance epsilon 1.5 for Geometric Level Set, change using meltTemp for Heat Level Set -->
					double eps = 1.5;  //MODIFY WITH MELT TEMP IF NEEDED
					Vector2d cm = Vector2d (0.,0.);    //Vector3d cm(0., 0., 0.);
					Tvec=_lset.getLevelset();                     //ArrayXd heavi = sHeavi(-lset.getGrid()); Geometric Level Set-->
					//Tvec=_lset.getLevelset();                                                           //Heat Level Set-->
					//Tvec = _grainTemp2D.getLevelset();
					//ArrayXd heavi = sHeavi(-_lset.getLevelset()); //ArrayXd heavi = sHeavi(-lset.getGrid());
					//if (meltTemp>=-10){
						
						for (size_t j = 0; j < _lset.getYdim(); j++) {
								for (size_t i = 0; i < _lset.getXdim(); i++) {   //Substitute of the Heaviside Function to Identify Inside and Outside of Grain
									if (-Tvec[(j*_lset.getXdim())+i] > eps) {  
				                       Tvec[(j*_lset.getXdim())+i] = 1; 
				                    }
									else if (-Tvec[(j*_lset.getXdim())+i] < -eps) {
				                       Tvec[(j*_lset.getXdim())+i] = 0;
									}	
									else{
				                       Tvec[(j*_lset.getXdim())+i] = 0.5*(1. + Tvec[(j*_lset.getXdim())+i]/eps + sin(M_PI*Tvec[(j*_lset.getXdim())+i]/eps)/M_PI); //0.5*(1. + vec(i)/eps + sin(M_PI*vec(i)/eps)/M_PI);
									}	
								}
						}
					//}	
					// else{   //For Freezing
					//    eps=eps-1.7;
	    //                for (size_t j = 0; j < _lset.getYdim(); j++) {
					// 			for (size_t i = 0; i < _lset.getXdim(); i++) {   //Substitute of the Heaviside Function to Identify Inside and Outside of Grain
					// 				if (Tvec[(j*_lset.getXdim())+i] > eps) {  
				 //                       Tvec[(j*_lset.getXdim())+i] = 1; 
				 //                    }
					// 				// else if (Tvec[(j*_lset.getXdim())+i] < eps+3) {
				 //     //                   Tvec[(j*_lset.getXdim())+i] = 0;
					// 				// }	
					// 				else{
					// 				    Tvec[(j*_lset.getXdim())+i] = 0; 	 
				 //                       //Tvec[(j*_lset.getXdim())+i] = 0.5*(1. + Tvec[(j*_lset.getXdim())+i]/eps + sin(M_PI*Tvec[(j*_lset.getXdim())+i]/eps)/M_PI); //0.5*(1. + vec(i)/eps + sin(M_PI*vec(i)/eps)/M_PI);
					// 				}	
					// 			}
					// 	}
					// 	eps = 1.5+meltTemp;
					// }
					//Tvec changes to an almost binary vector to find centroid and moment inertia.

					//size_t iter = 0;  //DO I NEED IT?
					for (size_t j = 0; j < _lset.getYdim(); j++) {
						for (size_t i = 0; i < _lset.getXdim(); i++) {
							vol+=Tvec[(j*_lset.getXdim())+i];     //vol += heavi(iter);
							cm(0) += Tvec[(j*_lset.getXdim())+i]*(double)i;    //cm(0) += heavi(iter)*(double)i;
							cm(1) += Tvec[(j*_lset.getXdim())+i]*(double)j;    //cm(1) += heavi(iter)*(double)j;
							//iter++;
						}
					}

					cm /= vol;

					// cout << "Old Mass: " << _mass << endl;
					// cout << "Old CM: " << _cmLset(0) << " , " << _cmLset(1) <<endl;
					// cout << "Old Density: " << _density << endl;
					// cout << "Old Radius: " << _radius << endl;

					//Preservation of Density (is this really needed)
					//if (_density == 0.91e-6)
					// {}
					// else{
					// 	_density = 0.91e-6;
					// }

					_mass = vol*_density;  //*0.89;  //info._mass = vol;   //Effect of mass?? FACTOR 1/0.9 Increase Consistently 

					//Update for thickness
					mult_loss = _thick/_thick0;
					//mult_loss = 1.0; 
                   	_mass *= mult_loss;
                   	
                   	//Aug-18-2022 Change
                   	//BEGIN
                   	//Rescale mass with thickness since it's not unitary
                    
                    _mass *= _thick0 * 0.001; //CHECK EFFECT!!!!
                   	
                   	//END
                   	

                    //_mass = vol*_density;
					// cout << "New Mass: " << _mass << endl;

					//_mass = vol*dens*0.89;
					//_mass = _mass*0.91e-6;
					// cout << "New cm for grain inner is : " << cm(0) << " , " << cm(1) <<endl;
					// cout << "New density for grain inner is : " << _density <<endl;
					// cout << "New volume for grain inner is : " << vol <<endl;
     //                cout << "New mass for grain inner is : " << _mass <<endl;
	                //_temper=_mass;

			        //Change in Center of Mass due to Change of Level Set
			        //_cmLset //Will it affect possition or what will happen depending on symmetry of melting

			        Vector2d old_cm = _cmLset;
			        Vector2d old_position = _position;

			        //New CM from updated Level Set
			        _cmLset = cm; //info._cm = cm;   //Pretty much the same with slight alterations

			        //Modification of position using old and new cm ?????
			        _position = old_position + (_cmLset - old_cm);

			        //Change in Moment of Inertia: _momentInertia
			        

				    // iter = 0;
					// for (size_t k = 0; k < lset.getZdim(); k++) {
					double Inertia=0;
					double rx=0;
					double ry=0;

					//old cm
					//cout << "Old Moment of Inertia: " << _momentInertia << endl;

					for (size_t j = 0; j < _lset.getYdim(); j++) {
						for (size_t i = 0; i < _lset.getXdim(); i++) {
							rx = (double)i - _cmLset(0);
							ry = (double)j - _cmLset(1);
							//rz = (double)k - cm(2);
							//Inertia += Tvec[(j*_lset.getXdim())+i]*(rx*rx + ry*ry); //Effect of mass???? FACTOR 1/0.9 Increase Consistently
							//Inertia += _density*Tvec[(j*_lset.getXdim())+i]*(rx*rx + ry*ry); //Effect of mass???? FACTOR 1/0.9 Increase Consistently  #Original MoI
							//Inertia += _density*Tvec[(j*_lset.getXdim())+i]*(rx*rx + ry*ry) * mult_loss ; //Effect of mass???? FACTOR 1/0.9 Increase Consistently # MODIF DECEMBER 5, 2022 (MoI Adjustment) v1   (Using relative thick change)
							Inertia += _density*Tvec[(j*_lset.getXdim())+i]*(rx*rx + ry*ry) * Tvecsm[(j*_lset.getXdim())+i] ; //Effect of mass???? FACTOR 1/0.9 Increase Consistently # MODIF DECEMBER 5, 2022 (MoI Adjustment) v2 (Using all the thickness vector cell by cell)
							//I(0,1) += -heavi2(iter)*rx*ry;
						}
					}
					//Inertia /= vol; //???????
			        
			        //_momentInertia=Inertia*0.89;
			        _momentInertia=Inertia;
			        //cout << "New Moment of Inertia: " << _momentInertia << endl;

			        //cout << "New Point Generation" << endl;
			        vector <Vector2d> pneworder(_npoints);
			        bool fail_ind = false; //Assume ok
			        PointsGen(_lset, pneworder, _npoints, _cmLset, fail_ind ); //Assuming MeltTemp is equal to zero, otherwise must be modified!!!!
		            //PointsGen(_grainTemp2D, pneworder, _npoints, _cmLset ); //Assuming MeltTemp is equal to zero, otherwise must be modified!!!!
 
			        if (fail_ind)
			        {
			        	cout <<"WARNING: Wrong Point INTERPOLATION contact problems may occur" << endl;
			        	return; //Leave points the same to avoid mess
			        }
     //                cout << "Print Old Points" << endl;
					// for (size_t i = 0; i < _npoints; i++) {
			  //          cout << _pointList[i](0) << " " << _pointList[i](1) <<endl;
					// }

			        //cout << "Points for Audit" << endl;
	                //Update New Points 
					for (size_t i = 0; i < _npoints; i++) {
			           //_pointList[i]=pneworder[i]+_position-_cmLset;   ///Actually Adding Position is Necessary to ReLocate them
			           _pointList[i] = pneworder[i];   ///Actually Adding Position is Necessary to ReLocate them
			           //cout << pneworder[i](0) << pneworder[i](1) << endl;
					}

                    //Set Points with the Right Rotation as well and Move to Centroid Position
					Matrix2d rotMatrix;
					rotMatrix << cos(_theta), -sin(_theta), sin(_theta), cos(_theta);

					for (size_t i = 0; i < _pointList.size(); i++) {
						_pointList[i] = rotMatrix*(_pointList[i] - _cmLset) + _position;
					}

     //                cout << "Print New Points" << endl;
					// for (size_t i = 0; i < _npoints; i++) {
			  //          cout << _pointList[i](0) << " " << _pointList[i](1) <<endl;
					// }

			        //End


			        //Change in Bbox Radius of Grain: _radius
				
					double maxradius = 0;
					// double minradius = (_pointList[0]-_position).norm();
					for (size_t i = 0; i < _npoints; i++) {
					     if ((_pointList[i]-_position).norm() > maxradius) {     //Should we subtract position to find real radius?
					         maxradius = (_pointList[i]-_position).norm();
					        //cout << maxradius << endl;
					     }
			       //       if((_pointList[i]-_position).norm() < minradius) {
			       //          minradius = (_pointList[i]-_position).norm();
			       // //          //cout << maxradius << endl;
			       //       } 

					}
					_radius=maxradius;
					//cout << "New Radius: " << _radius << endl;
					
        			//Change Aug 22, 2022
                    //BEGIN
                    //Obtain new area for lateral mass loss calculation using points
                    vector<Vector2d> VecOriginN = _pointList;
                    double AreaMassN = 0.0;
                    size_t nmassN = VecOriginN.size();
                    
                    for (size_t i = 0; i < nmassN-1; i++)
                    {
                        AreaMassN += ( VecOriginN[i](0) * VecOriginN[i+1](1) -  VecOriginN[i](1) * VecOriginN[i+1](0) ); 
                    }
                    AreaMassN += (VecOriginN[nmassN-1](0) * VecOriginN[0](1) -  VecOriginN[nmassN-1](1) * VecOriginN[0](0) ); 
                    
                    if (isnan(AreaMassN) == 1)
                    {
                        cout << "WARNING: Nan in TempState algorithm for inner grain!!" << endl;
                        AreaMassN = 1.0; //Just to debug
                    }
                    double AreaMeltMassN = 0.5*abs(AreaMassN); //Area from points to calculate floe melt (using original pointlist)

                    loss_mcv_temp += max((AreaMeltMass - AreaMeltMassN), 0.0) * max( _thick , 0.0) * 0.025 * 0.001 * _density;  //WARNING: Since _thick is the average thickness. We assume a 5% thickness at basal melt zones due to shape change
                    
                    if (ct_ave_thick > 0)
                    {
                    	ave_solarT = ave_solar / double(ct_ave_thick);
                    	ave_oceanT = ave_ocean / double(ct_ave_thick);
                    	//cout << "Average Thickness Count: " << ct_ave_thick << endl;
                    	//cout << "Average Thickness Sum: " << ave_thick << endl;
                    	//cout << "Thickness: " << _thick << endl;
                    	
                    	 if ( (ave_solarT + ave_oceanT) > 0){
                            frac_solar = (ave_solarT)/(ave_solarT + ave_oceanT);
                            frac_ocean = (ave_oceanT)/(ave_solarT + ave_oceanT);
                        }
                        else{
                            frac_solar = 0.0;
                            frac_ocean = 0.0;
                        }
                        loss_mcv_solar_temp += frac_solar * max((AreaMeltMass - AreaMeltMassN), 0.0) * max( _thick , 0.0) * 0.025 * 0.001 * _density;
                        loss_mcv_ocean_temp += frac_ocean * max((AreaMeltMass - AreaMeltMassN), 0.0) * max( _thick , 0.0) * 0.025 * 0.001 * _density;
                        //cout << "melt component: " << max((AreaMeltMass - AreaMeltMassN), 0.0) * max( _thick , 0.0) * 0.025 * 0.001 * _density << endl;
                        //cout << "S melt component: " << loss_mcv_solar_temp << endl;
                        //cout << "O melt component: " << loss_mcv_ocean_temp << endl;
                    }
                    else{
                        loss_mcv_solar_temp += 0.5 * max((AreaMeltMass - AreaMeltMassN), 0.0) * max( _thick , 0.0) * 0.025 * 0.001 * _density;
                        loss_mcv_ocean_temp += 0.5 * max((AreaMeltMass - AreaMeltMassN), 0.0) * max( _thick , 0.0) * 0.025 * 0.001 * _density;
                    }
                    
                    
                    //If necessary this loss_mcl term will be removed and lateral melt will be assumed to be averaged by thickness. It is not exact to assume _thick to be the thickness of the lateral melt loss, must be less because the edges are thinner than the center. //TRY AND SEE
                    //END 
                    
	        	} 
                
        	}

	        // //1D Changes
	        // //Update Thickness if we have melt after temperature change 
	        // if(_thick0-(dh*meltid)>0.000097)  //
	        // {
	        //   _thick=_thick0-(dh*meltid);
	        //   //_thick=_thick0-(dh*meltid);   //-(dh*meltid);  //Constant melting
	        // }
	        // else
	        // {
	        //   _thick=0.000; //Eliminate Grain  
	        // }
	        // _thick0 = _thick; //Update thickness for next step

	        //Update Mass if Thickness is reduced
	        //_mass=(_mass0*0.91/1000000)*(_thick/_thick0); //*(_thick/_thick0);    //????????
	        //_temper=_mass-(_mass0*0.91/1000000.);
	        //Tvec=_grainTemp2D.getLevelset();
	        //_temper=Tvec[(    (_lset.getXdim()*0.25*_lset.getXdim())   +  0.25*_lset.getYdim()   )];               //-->
	        //_temper=Tvec[(    (1*_lset.getXdim())   +  1   )];
	        //_temper=_pointList[_npoints-1](0);
	        //_temper=_grainTemp2D.getGridValue2(_lset.getXdim()*0.5,_lset.getYdim()*0.5); 
	        //_temper=_mass;
        }   //Only do it if the grain is big enough otherwise grain must be melted

        //cout << "End TemperatureModif Okay" << endl;
    	//return;
    }

 //    //Separate function (OPTIONAL)
	// void TemperatureModifBasal (double Tair, double Twater, const double & dt, double dh, const double alphaice, double meltTemp, double meltid, size_t & tstep, const size_t & START_TEMP) {  

	// 	//Shift from meter to km scale
	// 	if (tstep == START_TEMP)
	// 	{
	// 		// cout << "Setting to km scale" << endl;
	// 		// _thick *= 0.001;
	// 		// _thick0 *= 0.001;

	// 	}

	// 	double thick_or = _thick;

	// 	//Update Boundary Conditions (INITIAL THERMAL CONDITIONS OF THE ICE FLOE COME FROM readInputFile.h)
	// 	_Utemper[0]=Twater;    //Bottom
	// 	_Utemper[800-1]=Tair; //Top
	// 	//FDM 1-D Thermal Loop for Space (time loop uses dt)
	// 	dh =_thick0/800.;  //For 80,000 we need thick0 to be at most 0.4 at the input
	// 	_Utemper0=_Utemper; 

	// 	double alphaicef = 1.0e-2;

	// 	for (int idh = 1; idh<800-1; idh++) 
	// 	{
	// 		_Utemper[idh] = _Utemper0[idh] + alphaicef*alphaice*(dt/(dh*dh))*(_Utemper0[idh+1]-2*_Utemper0[idh]+_Utemper0[idh-1]); //FDM scheme for 1-D diffusion for U_{t}-aU_{xx}=0 FWD t, CNT x
	// 		//_Utemper[idh] = _Utemper0[idh] + 0.005*alphaice*(dt/(dh*dh))*(_Utemper0[idh+1]-2*_Utemper0[idh]+_Utemper0[idh-1]); //FDM scheme for 1-D diffusion for U_{t}-aU_{xx}=0 FWD t, CNT x
	// 		//_Utemper[idh] = _Utemper0[idh] + ((dh*dh)/(2*alphaice))/(dh*dh)*(_Utemper0[idh+1]-2*_Utemper0[idh]+_Utemper0[idh-1]);
	// 		if (_Utemper[idh]>meltTemp)
	// 		{
	// 			meltid=idh;
	// 		}	

	// 	}

	// 	//1D Changes
	// 	//Update Thickness if we have melt after temperature change 
	// 	if(_thick0-(dh*meltid)>0.000097)  //
	// 	{
	// 		_thick=_thick0-(dh*meltid);
	// 	//_thick=_thick0-(dh*meltid);   //-(dh*meltid);  //Constant melting
	// 		//Resize mass and we lose mass with basal melt
	// 		_mass *=  _thick / thick_or;
	// 	}
	// 	else
	// 	{
	// 		_thick=0.000; //Eliminate Grain  
	// 		_mass = 0.000;  //Eliminate Grain  
	// 	}
	// 	_thick0 = _thick; //Update thickness for next step

	// }




    void StressModifElastic(const double & rhoice,  const double & rhowater, double Tair, double Twater, double Ticeedge, const double & A_rheo, int flow_ind, size_t & tstep)
    {
        //Get original Stress Level Set and Other Level Sets												
		const size_t N = _lset.getYdim(); //Rows j is y 
		const size_t M = _lset.getXdim(); //Columns i is x    	

    }	
    
    void StressModif(const double & rhoice,  const double & rhowater, double Tair, double Twater, double Ticeedge, const double & A_rheo, size_t & ywsealevel, int flow_ind, size_t & tstep)
    {
    
        //Get original Stress Level Set and Other Level Sets
		//Matrixx MatrixS;   													
		const size_t N = _lset.getYdim(); //Rows j is y 
		const size_t M = _lset.getXdim(); //Columns i is x
		//vector <double>  Svec;                   
		//MatrixXd MatrixS(N,M);                                                         //**-----> MatrixXd MatrixS(N,M);   //using Svec
		//Svec = _grainStress2D.getLevelset();
		

		//vector <double> Tvec = _grainTemp2D0.getLevelset();  //CAREFUL Remove Zero When Doing fully
		vector <double> geovec = _lset.getLevelset();  //CAREFUL Remove Zero When Doing fully

		
		//MatrixS = MatrixDefVec(N, M, Svec);
        MatrixXd MatrixGeo = MatrixDefVec(N, M, geovec);              //**-----> MatrixXd MatrixT(N,M);   //using geovec
        MatrixXd SlopeGlobalX = MatrixDefVec(N, M, geovec);
        MatrixXd SlopeGlobalY = MatrixDefVec(N, M, geovec);
        MatrixXd Cls = MatrixGeo;                                   ///**-----> MatrixXd Cls = MatrixT;   
        
        cout << "Shape" << endl;
        cout << Cls<< endl;
        
        //Initial Limits for Stress Detection
        //size_t x0Lb = 0;         //Initial x position //Anything Left is at -1<= index but we use >=0 for ice
    	//size_t xgl = 122;        //<=123 i x Position of grounding line for reference
    	//size_t y0Db = 18;        //<=19 j  Location of Ocean Floor the Specific Level Set >=19 for ice or >18
    	//size_t ywsealevel = 34;  //<=161 j  y location of sealevel  (Actually 35 in index 1)

		//size_t xg0 = 18; // Or more is Ice
	    //size_t xg1 = 132; // Or Less is Ice (varies)
	    //size_t yg0 = 16; // Or More is Ice (varies)
	    //y0Db = 16;
	    //size_t yg1 = 46; // Or Less is Ice 


	    //Delimit Ice extents 
        //size_t x0 = 18; // Or more is Ice
        //x0Lb = 18;
        //size_t x1 = 132; // Or Less is Ice (varies)
        //size_t y0 = 16; // Or More is Ice (varies)
        //size_t y1 = 46; // Or Less is Ice 
        //size_t ywsealevel = 34;

    	//x0Lb and y0Db are the most important to delimit where to clip ocean floor

        //Modify Stress Matrix using Matrixx Shape
        // for(size_t j = 0; j < N; ++j) //y
		// {
		//     for(size_t i = 0; i < M; ++i) //x
		//     {
		//         MatrixS[j][i] = j*i;  //Add element at each column per row  //Extract each LS Element using this Mechanism
		//     }
		// }  

		//Position of Uphill and Ocean (Depends on level set)
		size_t y_uphill =  172; //Or More is Uphill, less is sidemountains  //169
		size_t y_ocean  =  78;  //Or More is sidemountains, less is ocean, below here it can flow!!!

   
        //Material Tag separation using Geometry and Simple Vertical Separation, a Matrix of Material could be given later
        MatrixXd Cls_Tag = MatrixXd::Zero(N,M);       ///**-----> MatrixXd Cls_Tag = MatrixXd::Constant(N,M,0); 

        //Material Tagging Function  (0 Inside, 1 Outside)
        for(size_t j = 0; j < N; ++j) //y
	    {   
	    	for(size_t i = 0; i < M; ++i) //x
		    {
		     	if ( Cls(j,i) == 1 )   
		     	{
		     		Cls_Tag(j,i) = 1;
		     	}
		    } 

        }
        //Material Tagging Using only Geometry Data



        //Softening function to avoid 3 side free cell and be sure of corner at most B.C.s is done before entering this function
        //size_t weird_ind = 1;   //Indicator of something not soft
        //size_t ttries = 0;      //Number of tries to get to objective (as a while loop stuck countermeasure)
        //size_t maxttries = 50;  //Max number of tries to avoid getting stuck in softening

        //Matrixx Cls_Tagr = Cls_Tag;
        //MatrixXd Cls_Tagr = Cls_Tag;          //**-----> MatrixXd Cls_Tagr = MatrixTag; 

        //cout<<"Before Softening"<<endl;
        //cout<< Cls_Tagr << endl;
        //Cls_Tag = LS_Soft(Cls_Tagr, weird_ind, ttries, maxttries, y0Db, ywsealevel); 
        //cout<<"After Softening"<<endl;
        //cout<< Cls_Tagr << endl;

        //Obtain general extent of Ice or Relevant Components for Continuum Mechanics Analysis (Max Dims of Glacier to eliminate outer rim)
        size_t  maxcol = 0;
        size_t  maxrow = 0;
        size_t  mincol = M*2000000; //Large number for comparison
        size_t  minrow = N*2000000;
        for(size_t j = 0 ; j < N-1 ; ++j) //y
	    {   
	    	for(size_t i = 0; i < M-1 ; ++i) //x
		    {
             	if  ( (Cls_Tag(j,i) < Cls_Tag(j,i+1) ) && i>maxcol )  {
               		maxcol = i;
                } 
                if  ( (Cls_Tag(j,i) > Cls_Tag(j,i+1) ) && i<mincol )  {
               		mincol = i;
                } 
           		if  ( (Cls_Tag(j,i) < Cls_Tag(j+1,i) ) && j>maxrow  ) {
               		maxrow = j;
           		}  
           		if  ( (Cls_Tag(j,i) > Cls_Tag(j+1,i) ) && j<minrow )  {
               		minrow = j;
                } 
		    }
		}    
		//cout<<"Bk1"<<endl; 

        //Update sea level position for stress analysis
        //size_t yw0_sl = ywsealevel - y0Db;

        //Update y_uphill and y_ocean for Stress Analysis
        size_t y_uphill_cut, y_ocean_cut; //Given that minrow is index 0 or minrow-minrow

        if (minrow < y_ocean)
        {
         	y_uphill_cut = y_uphill - minrow; //Given that minrow is index 0 or minrow-minrow
         	y_ocean_cut =  y_ocean - minrow;
        }
        else if (minrow >= y_ocean && minrow <= y_uphill)
        {
        	y_uphill_cut = y_uphill - minrow;
        	y_ocean_cut =  0;
        }	
        else
        {
        	cout << "ERROR: We are going out of glaciar bounds due to cuts. Exit." << endl;
        	exit(1);
        }

        cout << "TIME: " << tstep << endl;
        cout << "minrow cutting" << minrow << endl;

        size_t minrow_abs = minrow; //For absolute control of geometry.


        //Initialize Boundary Condition Matrix
        
        //size_t lsy = maxrow + 1 - ( y0Db + 1 )  ;
        //size_t lsy = maxrow + 1 - ( y0Db + 1 );
        //size_t lsy = maxrow - ( y0Db - 1 ) - 1 ;  //CAREFUL!!! These dimensions are critical for Stress Calculation and for going back to full Level Set
        //size_t lsx = maxcol - (xg0 - 1) ;
        //size_t lsx = maxcol + 1;
        size_t lsx = maxcol - mincol + 1;
        size_t lsy = maxrow - minrow + 1;
        cout << "LSY: " << lsy <<  " LSX: " << lsx << endl; //COMPARE BELOW
        cout << "MaxCol: " << maxcol <<  " MaxRow: " << maxrow << endl; 
        cout << "MinCol: " << mincol <<  " MinRow: " << minrow << endl; 


        //Dimensions in General (let's keep one-to-one for upper view)
    	double DimY = 1;
    	double DimX = (lsy/lsx) * DimY; //To avoid effects of left BC and have a better terminus behavior

        //Matrixx Cls_BC = MatrixDefZero(lsy, lsx);
        MatrixXd Cls_BC = MatrixXd::Constant(lsy,lsx,0);  ///**-----> MatrixXd Cls_BC = MatrixXd::Constant(lsy,lsx,0); 

  //       //Global Slope to Clip versus Glacier Area build using
  //       //MatrixXd SlopeGlobalFullX = _globalTerrain; 
  //       //MatrixXd SlopeGlobalFullY = _globalTerrain;
  //       MatrixXd SlopeGlobalFullX = MatrixXd::Constant(N,M,0);
  //       MatrixXd SlopeGlobalFullY = MatrixXd::Constant(N,M,0);

  //       //Extract for Stresses
  //       for (size_t j = 0; j < N - 1; j++)
  //       {
  //       	for(size_t i = 0; i < M - 1; i++)
  //       	{
  //       		//SlopeGlobalFullX(j,i) = (_globalTerrain(j,i+1) - _globalTerrain(j,i))/(DimX / lsx) ;
  //       	    SlopeGlobalFullY(j,i) = (_globalTerrain(j+1,i) - _globalTerrain(j,i))/(DimY / lsy) ;
  //       	    SlopeGlobalFullX(j,i) = (_globalTerrain(j,i+1) - _globalTerrain(j,i))/(DimY / lsy) ;  //FOR NOW
  //       	} 
  //       }	
  //       MatrixXd SlopeGlobalX = MatrixXd::Constant(lsy,lsx,0);
  //       MatrixXd SlopeGlobalY = MatrixXd::Constant(lsy,lsx,0);

  //       //Pre-process Cls_BC by Clipping from Cls_Tag   //VERY CAREFUL!!!!!!!!!   //Also for GLobal Slope
  //       size_t ii = 0;
  //       size_t jj = 0;
  //       for(size_t j = minrow; j < maxrow + 1; ++j) //y  // or only minrow   maxrow + 1  //This results in lower empty
	 //    {   
	 //    	ii = 0;
	 //    	for(size_t i = mincol; i < maxcol + 1; ++i) //x //or only mincol  maxcol + 1  //This results in left empty
		//     {   
  //               Cls_BC(jj,ii) = Cls_Tag(j,i); 
  //               SlopeGlobalX(jj,ii) = SlopeGlobalFullX(j,i);
  //               SlopeGlobalY(jj,ii) = SlopeGlobalFullY(j,i);
  //               ++ii;
		//     }
		//     ++jj;
		// }    	
  //       cout << "jj: " << jj <<  " ii: " << ii << endl; //COMPARE ABOVE


        //Set BC Tags (Not the same as Material Tags)
		//Accentuate Outside Objects using -1
        
        for(size_t j = 0; j < lsy; ++j) //y
	    {   
	    	for(size_t i = 0; i < lsx; ++i) //x
		    {   
		       if ( int(Cls_BC(j,i)) != 0 )
		       {
		          Cls_BC(j,i) = -1;
		       }
		    }
		}

		// BC TAGGING
		// Now All inner objects have to be classified from 0 up to 28 depending on boundary conditions (see reference Notebook for symbology)
		// BY Default 0 is inner middle cells and -1 are outer cells. But some inner cells are actually boundary cells for BCs. For example 2 is Upper Right Corner uphill (12 for side-mountains, 22 for Ocean and so on),
		// Several ifs are used to sort all types of cells and w.r.t to y_uphill and y_ocean.

        //size_t gr_ln;
        cout << "BC TAGGING BEGIN" << endl;
        for(size_t j = 0; j < lsy; ++j) //y
	    {   
	    	for(size_t i = 0; i < lsx; ++i) //x
		    { 
             	//By Default everything is inner or zero, but must be polished   
				//Outer Cells Are Only Used as Reference for B.C., Skip outer cells that are not used
				if ( int(Cls_BC(j,i)) == 0 ) 
				{   
                    
                    //Up BC (1 represents Up)
                    if ( j >= 0  && j <lsy - 1 && i > 0 && i <lsx -1 ) // All except border
                    {
                    	if ( Cls_BC(j+1,i) == -1  &&  Cls_BC(j, i+1) != -1 && Cls_BC(j, i-1) != -1  )
                    	{
                    		Cls_BC(j,i) = basicClass(1, j, y_uphill_cut, y_ocean_cut);
                    	}
                    }
                    else if  ( j == lsy - 1 && i > 0 && i <lsx -1 ) //Only upper border
                    {
                    	if ( Cls_BC(j, i+1) != -1 && Cls_BC(j, i-1) != -1  )
                    	{
                    		Cls_BC(j,i) = basicClass(1, j, y_uphill_cut, y_ocean_cut);
                    	}
                    }

                    //Up Right BC (2 represents Upper Right)
                    if ( j > 0  && j <lsy - 1 && i > 0 && i <lsx -1 )
                    {
                    	if ( Cls_BC(j+1,i) == -1  &&  Cls_BC(j, i+1) == -1 && Cls_BC(j, i-1) != -1 && Cls_BC(j-1, i) != -1  )
                    	{
                    		Cls_BC(j,i) = basicClass(2, j, y_uphill_cut, y_ocean_cut);
                    	}
                    } 
                    else if ( j > 0  && j <lsy - 1 && i == lsx -1 )
                    {
                    	if ( Cls_BC(j+1,i) == -1  &&  Cls_BC(j-1, i) != -1  )
                    	{
                    		Cls_BC(j,i) = basicClass(2, j, y_uphill_cut, y_ocean_cut);
                    	}
                    }                      
                    else if ( j == lsy - 1 && i > 0 && i <lsx -1 )
                    {
                    	if ( Cls_BC(j, i+1) == -1 && Cls_BC(j, i-1) != -1 )
                    	{
                    		Cls_BC(j,i) = basicClass(2, j, y_uphill_cut, y_ocean_cut);
                    	}
                    }                     
                    else if ( j == lsy - 1 && i == lsx - 1 )
                    {
                    	Cls_BC(j,i) = basicClass(2, j, y_uphill_cut, y_ocean_cut);
                    } 

                    //Up Left BC (3 represents Upper Left)
                    if ( j > 0  && j <lsy - 1 && i > 0 && i <lsx -1 )
                    {
                    	if ( Cls_BC(j+1,i) == -1  &&  Cls_BC(j, i-1) == -1 && Cls_BC(j, i+1) != -1 && Cls_BC(j-1, i) != -1  )
                    	{
                    		Cls_BC(j,i) = basicClass(3, j, y_uphill_cut, y_ocean_cut);
                    	}
                    } 
                    else if ( j > 0  && j <lsy - 1 && i == 0 )
                    {
                    	if ( Cls_BC(j+1,i) == -1  &&  Cls_BC(j-1, i) != -1  )
                    	{
                    		Cls_BC(j,i) = basicClass(3, j, y_uphill_cut, y_ocean_cut);
                    	}
                    }                      
                    else if ( j == lsy - 1 && i > 0 && i <lsx -1 )
                    {
                    	if ( Cls_BC(j, i-1) == -1 && Cls_BC(j, i+1) != -1 )
                    	{
                    		Cls_BC(j,i) = basicClass(3, j, y_uphill_cut, y_ocean_cut);
                    	}
                    }                     
                    else if ( j == lsy - 1 && i == 0 )
                    {
                    	Cls_BC(j,i) = basicClass(3, j, y_uphill_cut, y_ocean_cut);
                    } 

                    //Right BC (4 represents Right)
                    if ( j > 0  && j <lsy - 1 && i >= 0 && i <lsx - 1 )
                    {
                    	if ( Cls_BC(j,i+1) == -1  &&  Cls_BC(j+1, i) != -1 && Cls_BC(j-1, i) != -1  )
                    	{
                    		Cls_BC(j,i) = basicClass(4, j, y_uphill_cut, y_ocean_cut);
                    	}
                    } 
                    else if ( j > 0  && j <lsy - 1 &&  i == lsx - 1 )
                    {
                    	if ( Cls_BC(j+1, i) != -1 && Cls_BC(j-1, i) != -1  )
                    	{
                    		Cls_BC(j,i) = basicClass(4, j, y_uphill_cut, y_ocean_cut);
                    	}
                    }

                    //Left BC (5 represents Left)
                    if ( j > 0  && j <lsy - 1 && i > 0 && i <lsx )
                    {
                    	if ( Cls_BC(j,i-1) == -1  &&  Cls_BC(j+1, i) != -1 && Cls_BC(j-1, i) != -1  )
                    	{
                    		Cls_BC(j,i) = basicClass(5, j, y_uphill_cut, y_ocean_cut);
                    	}
                    } 
                    else if ( j > 0  && j <lsy - 1 &&  i == 0 )
                    {
                    	if ( Cls_BC(j+1, i) != -1 && Cls_BC(j-1, i) != -1  )
                    	{
                    		Cls_BC(j,i) = basicClass(5, j, y_uphill_cut, y_ocean_cut);
                    	}
                    }

                    //Lower Right BC (6 represents Lower Right)
                    if ( j > 0  && j <lsy - 1 && i > 0 && i <lsx -1 )
                    {
                    	if ( Cls_BC(j-1,i) == -1  &&  Cls_BC(j, i+1) == -1 && Cls_BC(j, i-1) != -1 && Cls_BC(j+1, i) != -1  )
                    	{
                    		Cls_BC(j,i) = basicClass(6, j, y_uphill_cut, y_ocean_cut);
                    	}
                    } 
                    else if ( j > 0  && j <lsy - 1 && i == lsx -1 )
                    {
                    	if ( Cls_BC(j-1,i) == -1  &&  Cls_BC(j+1, i) != -1  )
                    	{
                    		Cls_BC(j,i) = basicClass(6, j, y_uphill_cut, y_ocean_cut);
                    	}
                    }                      
                    else if ( j == 0 && i > 0 && i <lsx -1 )
                    {
                    	if ( Cls_BC(j, i+1) == -1 && Cls_BC(j, i-1) != -1 )
                    	{
                    		Cls_BC(j,i) = basicClass(6, j, y_uphill_cut, y_ocean_cut);
                    	}
                    }                     
                    else if ( j == 0 && i == lsx - 1 )
                    {
                    	Cls_BC(j,i) = basicClass(6, j, y_uphill_cut, y_ocean_cut);
                    } 

                    //Lower Left BC (7 represents Lower Left)
                    if ( j > 0  && j <lsy - 1 && i > 0 && i <lsx -1 )
                    {
                    	if ( Cls_BC(j-1,i) == -1  &&  Cls_BC(j, i-1) == -1 && Cls_BC(j, i+1) != -1 && Cls_BC(j+1, i) != -1  )
                    	{
                    		Cls_BC(j,i) = basicClass(7, j, y_uphill_cut, y_ocean_cut);
                    	}
                    } 
                    else if ( j > 0  && j <lsy - 1 && i == 0 )
                    {
                    	if ( Cls_BC(j-1,i) == -1  &&  Cls_BC(j+1, i) != -1  )
                    	{
                    		Cls_BC(j,i) = basicClass(7, j, y_uphill_cut, y_ocean_cut);
                    	}
                    }                      
                    else if ( j == 0 && i > 0 && i <lsx -1 )
                    {
                    	if ( Cls_BC(j, i-1) == -1 && Cls_BC(j, i+1) != -1 )
                    	{
                    		Cls_BC(j,i) = basicClass(7, j, y_uphill_cut, y_ocean_cut);
                    	}
                    }                     
                    else if ( j == 0 && i == 0 )
                    {
                    	Cls_BC(j,i) = basicClass(7, j, y_uphill_cut, y_ocean_cut);
                    } 

                    //Down BC (8 represents Down)
                    if ( j > 0  && j <lsy && i > 0 && i <lsx -1 ) // All except border
                    {
                    	if ( Cls_BC(j-1,i) == -1  &&  Cls_BC(j, i+1) != -1 && Cls_BC(j, i-1) != -1  )
                    	{
                    		Cls_BC(j,i) = basicClass(8, j, y_uphill_cut, y_ocean_cut);
                    	}
                    }
                    else if  ( j == 0 && i > 0 && i <lsx -1 ) //Only upper border
                    {
                    	if ( Cls_BC(j, i+1) != -1 && Cls_BC(j, i-1) != -1  )
                    	{
                    		Cls_BC(j,i) = basicClass(8, j, y_uphill_cut, y_ocean_cut);
                    	}
                    }

				}   

		    }
		} 
		cout << "BC TAGGING END" << endl;

		//DOES BC Setting generates single edges
	    //cout<<"Before Softening ???"<<endl;
        //cout<< Cls_BC << endl;
        // Cls_BC= LS_Soft(Cls_BC, weird_ind, ttries, maxttries, y0Db, ywsealevel); 
        // cout<<"After Softening"<<endl;
        // cout<< Cls_Tagr << endl;
       
		//MatrixPrint(Cls_BC, lsy, lsx);  //To check     
        //cout << Cls_BC << endl;
        //cout << Cls_Tag << endl;
        //Finite Difference Function 
        //Initialize Variables to Reference within the function
        //Matrixx Stress_Out;
        MatrixXd Stress_Out(lsy,lsx); ///**-----> MatrixXd Stress_Out(lsy,lsx);

        //To avoid extreme slopes at BCs
        for(size_t j = 0; j < lsy; ++j) //y
	    {   
	    	for(size_t i = 0; i < lsx; ++i) //x
		    {   
		       //if ( int(Cls_BC(j,i)) == 0 || int(Cls_BC(j,i)) == -1 )
		       if (fabs(SlopeGlobalX(j,i)) > 100  )
		       {
                    //Control Slope at BCs to avoid extreme values
		       	    SlopeGlobalX(j,i) = 0;
		       }
		       else
		       {
		       		
		       }
		       if (fabs(SlopeGlobalY(j,i)) > 100  )
		       {
                    //Control Slope at BCs to avoid extreme values
		       	    SlopeGlobalY(j,i) = 0;
		       }
		       else
		       {
		       		
		       }
		    }
		}

        Stress_Out = FD_Stokes_Glen(rhoice, rhowater, Cls_BC, lsy, lsx, y_uphill_cut, y_ocean_cut, DimX, DimY, A_rheo, 
                                    SlopeGlobalX, SlopeGlobalY, flow_ind, minrow, mincol, MatrixGeo, y_uphill, y_ocean, minrow_abs);  //!!!!!!!!!
        //Stress_Out = Cls_Tag;
        //Stress_Out = Cls_BC;
        //MatrixPrint(Stress_Out, lsy, lsx);
        //Output will be our desired stress Matrix, which directly will affect level set

        //Paste this Matrix to the N x M original size level set matrix for taumax at the right place, just as it was extracted before. No INFO will ???? !!!! (TODO)
      
        //Matrixx taumax;
        MatrixXd taumax = MatrixXd::Zero(N,M);  ///**-----> MatrixXd taumax::Constant(N,M,0);
        
        cout << "block" <<endl;
        taumax.block( minrow, mincol, Stress_Out.rows(), Stress_Out.cols() ) = Stress_Out;  //Paste Just at the Point of insertion If x0Lb were not zero we might have to use -1 with it
        //taumax = MatrixDefZero(N,M);
        //cout << taumax << endl;
        //taumax = Stress_Out_Full; from FINITE DIFFERENCE SCHEME
        
         

        
        //CHANGE!!!!! (TODO)
        cout << "Very Useful Material Tag" << endl;
        cout << Cls_BC<< endl;
        //cout << "Grain Stress taumax" << endl;
        //cout << taumax << endl;

        //To Level Set
        vector<double> Svec;
		Svec = MatrixtoVec(taumax);
        
        //Modify Level Set internal vector
        _grainStress2D.changeLevelset(Svec); 

    }




    /*
    * FINITE DIFFERENCE SCHEME FOR STOKES EQUATION
    *See Long-Chen Reference, Notes and Presentation of MAC Scheme Solving
    *The scheme here is oriented for Neumann B.C. Up and Right (with changes in x due to water level) , Dirichlet B.C. Down and Left (all zero except for bcLu for a constant glacier inflow)
    *We get a spatially implicit scheme that must be solved by conjugate gradient or other method for sparse matrices for elliptic equations
    *Picard iteration convergence is proposed for a while loop but must be reviewed to prevent extremely long time (due to ice Glen Law behavior viscosity that depends on stress), so it is NOT active right now
    *This scheme can be modified to solve Hookes Law, which would be a much simpler stencil
    */

    //**----->UNCOMMENT ALL FX.
    MatrixXd FD_Stokes_Glen(const double & rhoice,  const double & rhowater, MatrixXd & Cls_BC, const size_t & lsy, const size_t & lsx, const size_t & y_uphill_cut, const size_t & y_ocean_cut, const double & DimX, const double & DimY, const double & A_rheo, 
    	                    MatrixXd & SlopeGlobalX, MatrixXd & SlopeGlobalY, int flow_ind, size_t & minrow, size_t & mincol, MatrixXd & MatrixGeo, const size_t & y_uphill, const size_t & y_ocean, const size_t & minrow_abs) 
    {
    	
    	//MatrixXd Stress_Out;
    	//Stress_Out = MatrixDefZero(lsy,lsx);
    	MatrixXd Stress_Out = MatrixXd::Zero(lsy, lsx);   ///**-----> MatrixXd Stress_Out(lsy,lsx);
        MatrixXd Cls_taumax = MatrixXd::Zero(lsy, lsx);
        MatrixXd Cls_taumax0 = MatrixXd::Zero(lsy, lsx);

        //Get Terrain Info for Slope Body Forces //*--> Only in Y right now
        //Used fromn MatrixXd SlopeGlobalFull = _globalTerrain; 
       
        //GRID RESOLUTION (insert exception just in case)
        //double hx = (DimX / lsx);
		double hy = (DimY / lsy); 
        double hx = hy; //FOR NOW
         
		//Define Boundary Conditions (Left and Down Dirichlet, Up and Right Neumann, see above)
		double ice_x_speed = -20; //Lateral Glacier Constant Flow  //Orig 300 Left Direction
		double ice_y_speed = -800; //Lateral Glacier Constant Flow  //Orig 300 Down Direction
		//Ice Weight Foce
		double grav = -9.81;
		double F2 = -rhoice * grav;
		//3 break 100   6 with 600, good for better propagation
		double w_adjust = 600; //Calibrate Weight if needed 20 orig

		//Dirichlet U and V 
		//double bcLu = 0,                bcRu = ice_x_speed,       bcUu = 0,                   bcDu = 0;
        //double bcLv = 0,                bcRv = 0,                 bcUv = ice_y_speed,         bcDv = 0;

        double bcLu = ice_x_speed,                bcRu = ice_x_speed,           bcUu = ice_x_speed,                   bcDu = ice_x_speed;
        double bcLv = ice_y_speed,                bcRv = ice_y_speed,                 bcUv = ice_y_speed,         bcDv = ice_y_speed;
     
        //Neumann U and V (placeholders if function is used instead) (free Stress condition)
        double bcLun = 0, bcRun = 0, bcUun = 0, bcDun = 0;
        double bcLvn = 0, bcRvn = 0, bcUvn = 0, bcDvn = 0;

        //Error for Picard Loop
        double MaxError = 2;
        double MaxError0 = 2;
        double tol = 0.04; //To avoid loop for the meantime !!!!!!! (TODO) 0.0000001
        size_t kt = 0;  //Number of Picard Loops
        size_t Maxtries = 0;  //Limiting Number of tries to avoid getting stuck
        double eta_val0 = 1;  //Initial Trial Value for Viscosity
        
        //Matrixx Cls_eta;      //Current Viscosity
        //Matrixx Cls_eta0;     //Prior Iteration Viscosity
        MatrixXd Cls_eta(lsy, lsx); ///**-----> MatrixXd Cls_eta  = MatrixXd::Constant(lsy,lsx,0); 
        MatrixXd Cls_eta0(lsy, lsx); ///**-----> MatrixXd Cls_eta0 = MatrixXd::Constant(lsy,lsx,0); 
 
        //Error Matrixes
        MatrixXd U_ls_error = MatrixXd::Zero(lsy, lsx);
        MatrixXd V_ls_error = MatrixXd::Zero(lsy, lsx);
        MatrixXd P_ls_error = MatrixXd::Zero(lsy, lsx);
        MatrixXd U_ls_error0 = MatrixXd::Zero(lsy, lsx);
        MatrixXd V_ls_error0 = MatrixXd::Zero(lsy, lsx);
        MatrixXd P_ls_error0 = MatrixXd::Zero(lsy, lsx);
 
        //So far all matrixes have been lsx*lsy and relatively Dense. But A, BT and B for Finite Diff. Stokes
        //result in VERY SPARSE AND LARGE MATRICES

        while (MaxError > tol)   //START OF WHILE LOOP
        {
        	//Matrixx Cls_U = MatrixDefZero(lsy,lsx);  //Assume BC points are already implicit for left
            //Matrixx Cls_V = MatrixDefZero(lsy,lsx);  //Assume BC points are already implicit for bottom
            //Matrixx Cls_P = MatrixDefZero(lsy,lsx);  //P(1,1) and P(2,1) go to U(1,1) and so on. Use Quadratic P Interp. at Border or 0, check performance.
            MatrixXd Cls_U(lsy,lsx);   ///**-----> MatrixXd Cls_U(lsy,lsx);
            MatrixXd Cls_V(lsy,lsx);   ///**-----> MatrixXd Cls_V(lsy,lsx);
            MatrixXd Cls_P(lsy,lsx);   ///**-----> MatrixXd Cls_P(lsy,lsx);


            //SPARSE	*******************************************	
            //To Assemble together with uF and vF
			//Matrixx Cls_A = MatrixDefZero(2*lsy*lsx,2*lsy*lsx);
			//Matrixx Cls_B = MatrixDefZero(2*lsy*lsx,lsy*lsx);
            MatrixXd Cls_A  = MatrixXd::Zero(2*lsy*lsx, 2*lsy*lsx); 
            MatrixXd Cls_B  = MatrixXd::Zero(2*lsy*lsx, lsy*lsx); 


			//To Assemble with pF
			//Transpose of Cls_B
			//Matrixx Cls_BT = MatrixDefZero(lsy*lsx,2*lsy*lsx);
			//Matrixx Cls_Z = MatrixDefZero(lsy*lsx,lsy*lsx);

            //Defined using tranpose
            //MatrixXd Cls_BT  = Cls_B.transpose(); 
            //MatrixXd Cls_Z  = MatrixXd::Zero(lsy*lsx, lsy*lsx); 
			//SPARSE	*******************************************


            //Row Vector to Assemble Answer
			//Matrixx Cls_uF = MatrixDefZero(lsy*lsx,1);
			//Matrixx Cls_vF = MatrixDefZero(lsy*lsx,1);
			//Matrixx Cls_pF = MatrixDefZero(lsy*lsx,1); 

			VectorXd Cls_uF = VectorXd::Zero((lsy*lsx));///**-----> VectorXd Cls_uF(lsy*lsx);
            VectorXd Cls_vF = VectorXd::Zero((lsy*lsx));///**-----> VectorXd Cls_vF(lsy*lsx);
            VectorXd Cls_pF = VectorXd::Zero((lsy*lsx)); ///**-----> VectorXd Cls_pF(lsy*lsx); 

            //MatrixXd Cls_uF(lsy*lsx,1);///**-----> VectorXd Cls_uF(lsy*lsx);
            //MatrixXd Cls_vF(lsy*lsx,1);///**-----> VectorXd Cls_vF(lsy*lsx);
            //MatrixXd Cls_pF(lsy*lsx,1); ///**-----> VectorXd Cls_pF(lsy*lsx); 



			//Define Viscosity Matrix (Zeros Represent Outer Values)
			if (kt == 0)
			{
			   //Cls_eta = MatrixDefConst(lsy, lsx, eta_val0); 
			    Cls_eta  = MatrixXd::Constant( lsy, lsx, eta_val0 );///**-----> MatrixXd Cls_eta  = MatrixXd::Constant(lsy,lsx, eta_val0); 
			}
			else
			{
			    Cls_eta = Cls_eta0;
			}
			//Set Outer Values to Zero (or NaN)
            for(size_t j = 0; j < lsy; ++j) //y
	    	{   
	    		for(size_t i = 0; i < lsx; ++i) //x
		    	{ 
		    		if ( int (Cls_BC(j,i) ) < 0 )
		    		{
						Cls_eta(j,i) = 0; 
			            //cout << "Initial Viscosity" << endl; 
			            //cout << Cls_eta << endl;  						
		    		}
		    	}	
            }
            //cout << "Initial Viscosity for Picard Step: " << kt << endl; 
			//cout << Cls_eta << endl;  	

            //No Horizontal Internal Forces to add to uF
            //Find Vertical Internal Forces to add to vF (due to weight)  //Find icetop to get Thickness of Weight

            //External BODY FORCES in X direction
   //          size_t icetop;

			// for(size_t j = 0; j < lsy; ++j) //y
			// {	
			//     for(size_t i = 0; i < lsx; ++i) //x
			//     {	
			//         if  (int (Cls_BC(j,i)) >=0)
			// 		{	
			// 			icetop = j;
			// 			if ( int(Cls_BC(lsy-1,i)) >=0 ) 
			// 			{	
			// 			  icetop = lsy-1;
			// 			}
			// 			else
			// 			{	
			// 			    for (size_t jj = j; jj < (lsy-1); ++jj)
			// 			    {	
			// 			        if ( int(Cls_BC(jj,i)) >= 0)
			// 			        { 
			// 			            icetop = jj;
			// 			        }    
			// 			        else
			// 			        {	
			// 			            break; 
			// 			        }
			// 			    } 
			// 			}    
			// 		    Cls_vF(j*lsx + i) = F2 * (icetop-j) * hy * w_adjust;   //CAREFUL
			// 		}    
			//     }   
			// } 

			//Separator for convenience sake
			size_t sepxls = lsx; 
            size_t sepyls = lsx;

            //Other multiples
            size_t st =lsx*lsy;   //Partial A Indexing, Half Diagonal
            size_t st2 =2*lsx*lsy;  //For B and BT

            /*
            * Definition of Stiffness Matrix A and B, and values for uF and vF
            */

            //Make Cls_MM into a Sparse Matrix for Efficiency in the Solver by using Triplets
            //Cls_MM_Sparse
            typedef Eigen::SparseMatrix<double> SpMat;
            typedef Eigen::Triplet<double> T;
            // Triplet List for Rigidity Matrix
            std::vector<T> tripletList;
            //tripletList.reserve(9*Cls_MM.rows());
            //tripletList.reserve(13*9*lsx*lsx*lsy*lsy);
            tripletList.reserve(13*3*lsx*lsy);  //VERIFY IF CORRECT

            // for (size_t j = 0; j < Cls_MM.rows(); ++j)
            // {
            //  	for (size_t i = 0; i < Cls_MM.cols(); ++i)
            //  	{
            //  		//Assemble Triplet
            //  		if (int(Cls_MM(j,i)) == 0)
            //  		{

            //  		}
            //  		else
            //  		{
            //  		    //Assemble Triplet 
            //  		    tripletList.push_back(T(j,i,Cls_MM(j,i))); 
            //  		}

            //  	}
            // }

            //For convenient value
            int Cell;
            double Slope;

            //For U

            for (size_t i = 0; i < lsx*lsy ; ++i)  //Read Left to Right, Down to Up
            {
	            double fi = double(i);                              //CAREFUL!!!!
	            double flsx = double (lsx); 
	            size_t xval =   i - (int(floor(fi/flsx))*lsx);
	            size_t yval =   int(floor(fi/flsx)) + 0;

	         //       //Average Viscosity (OPTIONAL) Effects still unknown (TODO) 
	         //       if round(Cls_BC (yval,xval),1) == (0 || 1 || 2 || 3 || 4 || 5 || 15 || 11 || 12 || 22)
	         //   		etaC = meanf(Cls_eta(yval,xval+1) , Cls_eta(yval,xval) );
	         //   		etaC= Cls_eta(  floor((i-1)/lsx)+1, i - (floor((i-1)/lsx)*lsx) ); 
	   			   // else if round(Cls_BC (yval,xval),1) == (6 || 16 || 7 || 8 || 9 || 19 || 10)
	         //   		etaCR = meanf(Cls_eta(yval,xval) , Cls_eta(yval,xval-1) ); 
	         //   		etaCR= Cls_eta(  floor((i-1)/lsx)+1, i - (floor((i-1)/lsx)*lsx) ); 
	   			   // else if round(Cls_BC (yval,xval),1) == 20
	         //   		etaCR = meanf(Cls_eta(yval,xval) , Cls_eta(yval,xval-1) );  
	         //   		etaCR= Cls_eta(  floor((i-1)/lsx)+1, i - (floor((i-1)/lsx)*lsx) ); 
	   			   // end      
	            
	            //Start with AT CELL values
	            double etaC =  Cls_eta(yval,xval);  //For middle cell
	            double etaCR = Cls_eta(yval,xval); //For Right BC 

	            Cell = int(Cls_BC(yval,xval));
                Slope = SlopeGlobalX(yval,xval)*F2;

                //Stencil Definition with multiple cases
                //Middle Cell
                if (Cell == 0)  
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*4 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    tripletList.push_back( T( i+st2,   i,  hy ) ); 

      			    Cls_uF(i)     += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)
 
                }

                //UP CELL (Uphill) (Dirichlet Stag Up because X)   1s
                else if (Cell == 1)  
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    tripletList.push_back( T( i+st2,   i,  hy ) ); 
                    
                    //On uphill we have a fixed velocity from the Top               
      				Cls_uF(i) += etaC*(2/(hx*hy))*bcUu;	
      				Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)
                }
                //UP CELL (Mountain-side) (Dirichlet Stag Up because X)
                else if (Cell == 11)  
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    tripletList.push_back( T( i+st2,   i,  hy ) ); 
                    
                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcUu = 0 here!!!!)    				
      				Cls_uF(i)     += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)
                }
                //UP CELL (Ocean)  (Neumann Stag Up because X)
                else if (Cell == 21)  
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*3 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    tripletList.push_back( T( i+st2,   i,  hy ) ); 
                    
      				Cls_uF(i) += etaC*(1/hx)*bcUun; 		//Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->   				
      				Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)
                }


                //UP RIGHT CELL (Uphill) (Dirichlet Stag Up and Normal Dirichlet on Right) 2s
                else if (Cell == 2)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    //tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    //tripletList.push_back( T( i+st2,   i,  hy ) ); 

                    //On uphill we have a fixed velocity from the Up - Right
      				Cls_uF(i) += etaC*(1/(hx*hy))*bcRu + etaC*(2/(hx*hy))*bcUu;                  	
      				Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)
                }                
                //UP RIGHT CELL (Mountain-side) (Dirichlet Stag Up and Normal Dirichlet on Right)
                else if (Cell == 12)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    //tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    //tripletList.push_back( T( i+st2,   i,  hy ) ); 

                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcRu = 0 , bcUu = 0 here!!!!)    				    			    
      				Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now) 
                } 
                //UP RIGHT CELL (Ocean) (Neumann Stag Up and Normal Neumann on Right)
                else if (Cell == 22)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaCR*1.5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaCR*1.0 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaCR*0.5 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+st2, 0.5*hy ) ); 
      			    //tripletList.push_back( T( i , i-1+st2,  -0.5*hy ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+st2, i, -0.5*hy ) ); 
      			    //tripletList.push_back( T( i-1+st2,   i,  0.5*hy ) ); 

      				Cls_uF(i)  =  0.5 * (Cls_uF(i)+ etaC*(1/hx)*bcUun ) + etaC*(1/hy)*bcRun;   //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *--> 
      				Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now) 
      			} 

                //UP LEFT CELL (Uphill)  (Dirichlet Stag Up and Normal Dirichlet on Left)  3s 
                else if (Cell == 3)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    tripletList.push_back( T( i+st2,   i,  hy ) ); 

                    //On uphill we have a fixed velocity from the Up - left
      				Cls_uF(i) += etaC*(1/(hx*hy))*bcLu + etaC*(2/(hx*hy))*bcUu;                  	
      				Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)
                }   
                //UP LEFT CELL (Mountain-side)  (Dirichlet Stag Up and Normal Dirichlet on Left)              
                else if (Cell == 13)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    tripletList.push_back( T( i+st2,   i,  hy ) ); 

                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcLu = 0 , bcUu = 0 here!!!!)    				    			    
      				Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)  
                }  
                //UP LEFT CELL (Ocean)  (Neumann Stag Up and Normal Neumann on Left)              
                else if (Cell == 23)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaCR*1.5 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaCR*1.0 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaCR*0.5 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+st2, 0.5*hy ) ); 
      			    tripletList.push_back( T( i , i-1+st2,  -0.5*hy ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+st2, i, -0.5*hy ) ); 
      			    tripletList.push_back( T( i-1+st2,   i,  0.5*hy ) ); 

      				Cls_uF(i)  =  0.5 * (Cls_uF(i)+ etaC*(1/hx)*bcUun ) + etaC*(1/hy)*bcLun;   //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *--> 
      				Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)
                }  

                //RIGHT CELL (Uphill) (Dirichlet Normal because X direction and uphill)  4s 
                else if (Cell == 4)
                {
                	//FOR A
      			    tripletList.push_back( T( i , i,      etaC*4 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    //tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    //tripletList.push_back( T( i+st2,   i,  hy ) ); 

                    //On uphill we have a fixed velocity from the right
      			    Cls_uF(i) += etaC*(1/(hx*hy))*bcRu;
      			    Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)
                }
                //RIGHT CELL (Mountain-side) (Dirichlet Normal because X direction and mountainside) 
                else if (Cell == 14)
                {
                	//FOR A
      			    tripletList.push_back( T( i , i,      etaC*4 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    //tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    //tripletList.push_back( T( i+st2,   i,  hy ) ); 

                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcRu = 0 here!!!!)    				
      			    Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)
                }
                //RIGHT CELL (Ocean) (Neumann Normal Right because X direction and ocean) 
                else if (Cell == 24)
                {
                	//FOR A
      			    tripletList.push_back( T( i , i,      etaC*2 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*0.5 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*0.5 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+1+st2, 0.5*hy ) ); 
      			    //tripletList.push_back( T( i , i+st2,  -0.5*hy ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+1+st2, i, -0.5*hy ) ); 
      			    //tripletList.push_back( T( i+st2,   i,  0.5*hy ) ); 

      			    Cls_uF(i)  = 0.5 * ( Cls_uF(i)) + etaC*(1/hy)*bcRun;  //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *--> 
      			    Cls_uF(i) += 0.5 * etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)
                }

                //LEFT CELL (Uphill) (Dirichlet Normal because X direction and uphill) 5s
                else if (Cell == 5)
                {
                	//FOR A
      			    tripletList.push_back( T( i , i,      etaC*4 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    tripletList.push_back( T( i+st2,   i,  hy ) ); 

                    //On uphill we have a fixed velocity from the left
      			    Cls_uF(i) += etaC*(1/(hx*hy))*bcLu;
      			    Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)
                }
                //LEFT CELL (Mountain-side) (Dirichlet Normal because X direction and mountainside) 
                else if (Cell == 15)
                {
                	//FOR A
      			    tripletList.push_back( T( i , i,      etaC*4 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    tripletList.push_back( T( i+st2,   i,  hy ) ); 

                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcLu = 0 here!!!!)    				
      			    Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)
                }	
                //LEFT CELL (Ocean) (Neumann Normal Left because X direction and ocean) 
                else if (Cell == 25)
                {
                	//FOR A
      			    tripletList.push_back( T( i , i,      etaC*2 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*0.5 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*0.5 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+1+st2, 0.5*hy ) ); 
      			    tripletList.push_back( T( i , i+st2,  -0.5*hy ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+1+st2, i, -0.5*hy ) ); 
      			    tripletList.push_back( T( i+st2,   i,  0.5*hy ) ); 
		  				
      				Cls_uF(i)  = 0.5 * (Cls_uF(i)) + etaC*(1/hx)*bcLun;  //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *--> 
      				Cls_uF(i) += 0.5 * etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)
                }


                //DOWN RIGHT CELL (Uphill) (Dirichlet Stag Down and Normal Dirichlet on Right)   6s
                else if (Cell == 6)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    //tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    //tripletList.push_back( T( i+st2,   i,  hy ) ); 

                    //On uphill we have a fixed velocity from the down - right
      				Cls_uF(i) += etaC*(2/(hx*hy))*bcRu + etaC*(1/(hx*hy))*bcDu;                  	
      				Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)
                }
                //DOWN RIGHT CELL (Mountain-side) (Dirichlet Stag Down and Normal Dirichlet on Right)   
                else if (Cell == 16)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    //tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    //tripletList.push_back( T( i+st2,   i,  hy ) ); 

                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcRu = 0 , bcDu = 0 here!!!!)    				    			    
      				Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)                 	
                }
                //DOWN RIGHT CELL (Ocean) (Neumann Stag Down and Normal Neumann on Right)   
                else if (Cell == 26)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaCR*1.5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaCR*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaCR*0.5 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+st2, 0.5*hy ) ); 
      			    //tripletList.push_back( T( i , i-1+st2,  -0.5*hy ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+st2, i, -0.5*hy ) ); 
      			    //tripletList.push_back( T( i-1+st2,   i,  0.5*hy ) ); 

      				Cls_uF(i)  = 0.5 * (Cls_uF(i)+ etaCR*(1/hx)*bcDun) + etaCR*(1/hy)*bcRun;   //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->                	
                    Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)  
                }                


                //DOWN LEFT CELL (Uphill)  (Dirichlet Stag Down and Normal Dirichlet on Left) 7s
                else if (Cell == 7)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    tripletList.push_back( T( i+st2,   i,  hy ) ); 

                    //On uphill we have a fixed velocity from the down - left
      				Cls_uF(i) += etaC*(1/(hx*hy))*bcLu + etaC*(2/(hx*hy))*bcDu;                  	
      				Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)
                }
                //DOWN LEFT CELL (Mountain-side)  (Dirichlet Stag Down and Normal Dirichlet on Left) 
                else if (Cell == 17)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    tripletList.push_back( T( i+st2,   i,  hy ) ); 

                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcLu = 0 , bcDu = 0 here!!!!)    				    			    
      				Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)      			    
                }
                //DOWN LEFT CELL (Ocean) (Neumann Stag Down and Normal Neumann on Right)   
                else if (Cell == 27)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaCR*1.5 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaCR*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaCR*0.5 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+st2, 0.5*hy ) ); 
      			    tripletList.push_back( T( i , i-1+st2,  -0.5*hy ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+st2, i, -0.5*hy ) ); 
      			    tripletList.push_back( T( i-1+st2,   i,  0.5*hy ) ); 

      				Cls_uF(i)  = 0.5 * (Cls_uF(i)+ etaCR*(1/hx)*bcDun) + etaCR*(1/hy)*bcLun;   //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->                	
                    Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)                 	
                }                   

                //DOWN CELL (Uphill)  (Dirichlet Stag Down because X and uphill)   8s
                else if (Cell == 8)  
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    tripletList.push_back( T( i+st2,   i,  hy ) ); 
                    
                    //On uphill we have a fixed velocity from the bottom
               
      				Cls_uF(i) += etaC*(2/(hx*hy))*bcDu;	
      				Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)
                }
                //DOWN CELL (Mountain-side)  (Dirichlet Stag Down because X and mountain-side)   8s
                else if (Cell == 18)  
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    tripletList.push_back( T( i+st2,   i,  hy ) ); 
                    
                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcDu = 0 here!!!!)    				
      				Cls_uF(i)  += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)

                }
                //DOWN CELL (Ocean)  (Neumann Stag Down because X and opening to ocean)   8s
                else if (Cell == 28)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*3 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+1+st2, hy ) ); 
      			    tripletList.push_back( T( i , i+st2,  -hy ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+1+st2, i, -hy ) ); 
      			    tripletList.push_back( T( i+st2,   i,  hy ) ); 
                    
                    //On uphill we have a fixed velocity from the bottom
               
      				Cls_uF(i) += etaC*(1/hx)*bcDun;	   //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *--> 
      				Cls_uF(i) += etaC*(1/(hx*hy)) * Slope * 0;  //BODY FORCE only in Y *--> (for now)                	
                }	

                else if (Cell == -1)
                {
                    //Do nothing.
                }               
                else
                {
                	cout << "Wrong Cell Selection at U for A,B,BT,uF" << endl;
                	cout << Cls_BC << endl;
                	exit(1);
                }

            }

            //For V

            //Useful variables to initialize
            size_t top, rowpos, colpos; 

            for (size_t i = st; i < 2*lsx*lsy ; ++i)  //Read Left to Right, Down to Up
            {
	            double fi = double(i-st);                              //CAREFUL!!!!
	            double flsx = double (lsx); 
	            size_t xval =   i-st - (int(floor(fi/flsx))*lsx);
	            size_t yval =   int(floor(fi/flsx))+0;

	          ////Average Viscosity (OPTIONAL) Effects still unknown (TODO) 
			  //  if round(Cls_BC (yval,xval),1) == (0 || 1 || 2 || 5.5 || 7 || 8 || 9 || 19 || 10 || 20 || 11 || 12 || 22)
			  //       etaC = meanf(Cls_eta(yval+1,xval) , Cls_eta(yval,xval) );
			  //       etaC= Cls_eta( floor((i-st-1)/lsx)+1, i-st - (floor((i-st-1)/lsx)*lsx) ); 
			  //  elseif round(Cls_BC (yval,xval),0) == (4||5||6)     
			  //       etaCU = meanf(Cls_eta(yval,xval) , Cls_eta(yval-1,xval) );
			  //       etaCU= Cls_eta( floor((i-st-1)/lsx)+1, i-st - (floor((i-st-1)/lsx)*lsx) ); 
			  // elseif round(Cls_BC (yval,xval),1) == (3)     
			  //       etaCU = meanf(Cls_eta(yval,xval) , Cls_eta(yval-1,xval) );
			  //       etaCU= Cls_eta( floor((i-st-1)/lsx)+1, i-st - (floor((i-st-1)/lsx)*lsx) ); 
			  //  elseif round(Cls_BC (yval,xval),1) == (16)     
			  //       etaCU = meanf(Cls_eta(yval,xval) , Cls_eta(yval-1,xval) );  
			  //       etaCU= Cls_eta( floor((i-st-1)/lsx)+1, i-st - (floor((i-st-1)/lsx)*lsx) ); 
			  //  end    
	            
	            //Start with AT CELL values
	            double etaC =  Cls_eta(yval,xval);  //For middle cell
	            double etaCU = Cls_eta(yval,xval); //For Up BC 

	            Cell = int(Cls_BC(yval,xval));
	            Slope = SlopeGlobalY(yval,xval)*F2;

                //Stencil Definition with multiple cases
                if (Cell == 0)  //Middle Cell
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*4 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    tripletList.push_back( T( i+st,     i,  -hx ) ); 

      			    Cls_vF(i-st) += etaC * (1/(hx*hy)) * Slope;  //Body Forces only in Y direction or for Velocity V (for now) *-->
                }

                //UP CELL (Uphill) (Dirichlet Normal Up because Y)   1s
                else if (Cell == 1)  
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*4 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    //tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    //tripletList.push_back( T( i+st,     i,  -hx ) ); 
                    
                    //On uphill we have a fixed velocity from the Top               
      				Cls_vF(i-st) += etaC * (1/(hx*hy)) * bcUv;	
      				Cls_vF(i-st) += etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now)
                } 
                //UP CELL (Mountain-side) (Dirichlet Normal Up because Y)
                else if (Cell == 11)  
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*4 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    //tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    //tripletList.push_back( T( i+st,     i,  -hx ) ); 
                    
                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcUv = 0 here!!!!)    				
      				Cls_vF(i-st) += etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now)
                }
                //UP CELL (Ocean)  (Neumann Normal Up because Y)
                else if (Cell == 21)  
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*2 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*0.5 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*0.5 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+lsx+st, -0.5*hx ) ); 
      			    //tripletList.push_back( T( i , i+st,      0.5*hx ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+lsx+st, i,   0.5*hx ) ); 
      			    //tripletList.push_back( T( i+st,     i,  -0.5*hx ) ); 

      			    Cls_vF(i-st) = 0.5 * Cls_vF(i-st) + etaC*(1/hx)*bcUvn;   //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->   				
      				Cls_vF(i-st) += 0.5 * etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now)
                }


                //UP RIGHT CELL (Uphill) (Dirichlet Normal Up and Stag Dirichlet on Right because Y) 2s
                else if (Cell == 2)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    //tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    //tripletList.push_back( T( i+st,     i,  -hx ) ); 

                    //On uphill we have a fixed velocity from the Top Right    
      			    Cls_vF(i-st) += etaC*(2/(hx*hy))*bcRv + etaC*(1/(hx*hy))*bcUv;   
      				Cls_vF(i-st) += etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now)      			     
                }                
                //UP RIGHT CELL (Mountain-side) (Dirichlet Normal Up and Stag Dirichlet on Right because Y)
                else if (Cell == 12)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    //tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    //tripletList.push_back( T( i+st,     i,  -hx ) ); 

                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcRv = 0 , bcUv = 0 here!!!!)    				    			    
      				Cls_vF(i-st) += etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now)   
                } 
                //UP RIGHT CELL (Ocean) (Neumann Normal Up and Stag Neumann on Right because Y)
                else if (Cell == 22)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*1.5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1   ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*0.5 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+lsx+st, -0.5*hx ) ); 
      			    //tripletList.push_back( T( i , i+st,      0.5*hx ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+lsx+st, i,   0.5*hx ) ); 
      			    //tripletList.push_back( T( i+st,     i,  -0.5*hx ) ); 

      			    Cls_vF(i-st)  = 0.5*Cls_vF(i-st) + etaCU*0.5*(1/hy)*bcRvn + etaCU*(1/hx)*bcUvn;  //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->
      			    Cls_vF(i-st) += 0.5 * etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now) 
      			} 

                //UP LEFT CELL (Uphill)  (Dirichlet Normal Up and Normal Neumann on Left)  3s 
                else if (Cell == 3)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    //tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    //tripletList.push_back( T( i+st,     i,  -hx ) ); 

                    //On uphill we have a fixed velocity from the Top Left    
      			    Cls_vF(i-st) += etaC*(2/(hx*hy))*bcLv + etaC*(1/(hx*hy))*bcUv;   
      				Cls_vF(i-st) += etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now)   
                }   
                //UP LEFT CELL (Mountain-side)  (Dirichlet Normal Up and Stag Dirichlet on Left)              
                else if (Cell == 13)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    //tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    //tripletList.push_back( T( i+st,     i,  -hx ) ); 

                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcLv = 0 , bcUv = 0 here!!!!)    				    			    
      				Cls_vF(i-st) += etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now) 
                }  
                //UP LEFT CELL (Ocean)  (Neumann Normal Up and Stag Neumann on Left)              
                else if (Cell == 23)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*1.5 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1   ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*0.5 ) ); 
                    //FOR B
                    //tripletList.push_back( T( i , i+lsx+st, -0.5*hx ) ); 
      			    //tripletList.push_back( T( i , i+st,      0.5*hx ) ); 
                    //FOR BT
                    //tripletList.push_back( T( i+lsx+st, i,   0.5*hx ) ); 
      			    //tripletList.push_back( T( i+st,     i,  -0.5*hx ) ); 

      			    Cls_vF(i-st)  = 0.5*Cls_vF(i-st) + etaCU*0.5*(1/hy)*bcLvn + etaCU*(1/hx)*bcUvn;  //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->
      			    Cls_vF(i-st) += 0.5 * etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now) 
                }  

                //RIGHT CELL (Uphill) (Dirichlet Stag because Y direction and uphill)  4s 
                else if (Cell == 4)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    tripletList.push_back( T( i+st,     i,  -hx ) ); 

                    //On uphill we have a fixed velocity from the right
      			    Cls_vF(i-st) += etaC*(2/(hx*hy))*bcRv; 
      			    Cls_vF(i-st) +=  etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now) 
                }
                //RIGHT CELL (Mountain-side) (Dirichlet Stag because Y direction and mountainside) 
                else if (Cell == 14)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    tripletList.push_back( T( i+st,     i,  -hx ) ); 

                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcRv = 0 here!!!!)    				
      			    Cls_vF(i-st) +=  etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now) 
                }
                //RIGHT CELL (Ocean) (Neumann Stag Right because Y direction and ocean) 
                else if (Cell == 24)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*3 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    tripletList.push_back( T( i+st,     i,  -hx ) ); 

      			    Cls_vF(i-st) += etaC*(1/hy)*bcRvn;   //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->
      			    Cls_vF(i-st) += etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now)       			    
                }

                //LEFT CELL (Uphill) (Dirichlet Stag because Y direction and uphill) 5s
                else if (Cell == 5)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    tripletList.push_back( T( i+st,     i,  -hx ) ); 

                    //On uphill we have a fixed velocity from the left
      			    Cls_vF(i-st) += etaC*(2/(hx*hy))*bcLv; 
      			    Cls_vF(i-st) +=  etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now) 
                }
                //LEFT CELL (Mountain-side) (Dirichlet Stag because Y direction and mountainside) 
                else if (Cell == 15)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    tripletList.push_back( T( i+st,     i,  -hx ) ); 

                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcLv = 0 here!!!!)    				
      			    Cls_vF(i-st) +=  etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now) 
                }	
                //LEFT CELL (Ocean) (Neumann Stag Left because Y direction and ocean) 
                else if (Cell == 25)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*3 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i-lsx, -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    tripletList.push_back( T( i+st,     i,  -hx ) ); 

      			    Cls_vF(i-st) += etaC*(1/hy)*bcLvn;   //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->
      			    Cls_vF(i-st) += etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now)    
                }


                //DOWN RIGHT CELL (Uphill) (Dirichlet Normal Down and Stag Dirichlet on Right)   6s
                else if (Cell == 6)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    tripletList.push_back( T( i+st,     i,  -hx ) ); 

                    //On uphill we have a fixed velocity from the down - right
      			    Cls_vF(i-st)  += etaC*(2/(hx*hy))*bcRv + etaC*(1/(hx*hy))*bcDv;                  	
      			    Cls_vF(i-st) +=  etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now) 
                }
                //DOWN RIGHT CELL (Mountain-side) (Dirichlet Normal Down and Stag Dirichlet on Right)   
                else if (Cell == 16)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    tripletList.push_back( T( i+st,     i,  -hx ) ); 

                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcRv = 0 , bcDv = 0 here!!!!)    				    			    
      			    Cls_vF(i-st) +=  etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now)                  	
                }
                //DOWN RIGHT CELL (Ocean) (Neumann Normal Down and Stag Neumann on Right because Y)   
                else if (Cell == 26)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*1.5 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*0.5 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+lsx+st, -0.5*hx ) ); 
      			    tripletList.push_back( T( i , i+st,      0.5*hx ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+lsx+st, i,   0.5*hx ) ); 
      			    tripletList.push_back( T( i+st,     i,  -0.5*hx ) ); 

      			    Cls_vF(i-st)  =  0.5*Cls_vF(i-st) + etaC*0.5*(1/hy)*bcRvn + etaC*(1/hx)*bcDvn;   //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->                	
      			    Cls_vF(i-st) += 0.5 * etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now)    
                }                


                //DOWN LEFT CELL (Uphill)  (Dirichlet Stag Down and Normal Dirichlet on Left) 7s
                else if (Cell == 7)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    tripletList.push_back( T( i+st,     i,  -hx ) ); 

                    //On uphill we have a fixed velocity from the down - left
      			    Cls_vF(i-st)  += etaC*(2/(hx*hy))*bcLv + etaC*(1/(hx*hy))*bcDv;                  	
      			    Cls_vF(i-st) +=  etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now) 
                }
                //DOWN LEFT CELL (Mountain-side)  (Dirichlet Stag Down and Normal Dirichlet on Left) 
                else if (Cell == 17)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*5 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    tripletList.push_back( T( i+st,     i,  -hx ) ); 

                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcLv = 0 , bcDv = 0 here!!!!)    				    			    
      			    Cls_vF(i-st) +=  etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now)     			    
                }
                //DOWN LEFT CELL (Ocean) (Neumann Normal Down and Stag Neumann on Right)   
                else if (Cell == 27)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*1.5 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*0.5 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+lsx+st, -0.5*hx ) ); 
      			    tripletList.push_back( T( i , i+st,      0.5*hx ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+lsx+st, i,   0.5*hx ) ); 
      			    tripletList.push_back( T( i+st,     i,  -0.5*hx ) ); 

      			    Cls_vF(i-st)  =  0.5*Cls_vF(i-st) + etaC*0.5*(1/hy)*bcLvn + etaC*(1/hx)*bcDvn;   //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->                	
      			    Cls_vF(i-st) += 0.5*etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now)                	
                }                   

                //DOWN CELL (Uphill)  (Dirichlet Stag Down because X and uphill)   8s
                else if (Cell == 8)  
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*4 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    tripletList.push_back( T( i+st,     i,  -hx ) ); 
       
                    //On uphill we have a fixed velocity from the bottom              
      				Cls_vF(i-st) += etaC*(1/(hx*hy))*bcDv; 	
      				Cls_vF(i-st) += etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now)
                }
                //DOWN CELL (Mountain-side)  (Dirichlet Stag Down because X and mountain-side)   8s
                else if (Cell == 18)  
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*4 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*1 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+lsx+st, -hx ) ); 
      			    tripletList.push_back( T( i , i+st,      hx ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+lsx+st, i,   hx ) ); 
      			    tripletList.push_back( T( i+st,     i,  -hx ) ); 
       
                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcDv = 0 here!!!!)    				
      				Cls_vF(i-st) += etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now)

                }
                //DOWN CELL (Ocean)  (Neumann Stag Down because X and opening to ocean)   8s
                else if (Cell == 28)
                {
                    //FOR A
      			    tripletList.push_back( T( i , i,      etaC*2 ) ); 
      			    tripletList.push_back( T( i , i-1,   -etaC*0.5 ) ); 
      			    tripletList.push_back( T( i , i+1,   -etaC*0.5 ) ); 
      			    tripletList.push_back( T( i , i+lsx, -etaC*1 ) ); 
                    //FOR B
                    tripletList.push_back( T( i , i+lsx+st, -0.5*hx ) ); 
      			    tripletList.push_back( T( i , i+st,      0.5*hx ) ); 
                    //FOR BT
                    tripletList.push_back( T( i+lsx+st, i,   0.5*hx ) ); 
      			    tripletList.push_back( T( i+st,     i,  -0.5*hx ) ); 

      			    Cls_vF(i-st) =  0.5 * Cls_vF(i-st) + etaC*(1/hx)*bcDvn;  //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *--> 
     				Cls_vF(i-st) += 0.5 * etaC * (1/(hx*hy)) * Slope;  //BODY FORCE only in Y *--> (for now)              	
                }	

                else if (Cell == -1)
                {
                    //Do nothing 
                }   
                else
                {
                	cout << "Wrong Cell Selection at V for A,B,BT,vF" << endl;
                	cout << Cls_BC << endl;
                	exit(1);
                }
            }


           /*
            * Definition of Stiffness Matrix BT is done above -B.transpose() (if not transpose of B for NonSymmetry), here is only for values of pF based on BCs
            */

            //For U

            for (size_t i = 0; i < lsx*lsy ; ++i)  //Read Left to Right, Down to Up
            {
	            double fi = double(i);                              //CAREFUL!!!!
	            double flsx = double (lsx); 
	            size_t xval =   i - (int(floor(fi/flsx))*lsx);
	            size_t yval =   int(floor(fi/flsx))+0;

				////Average Viscosity (OPTIONAL) Effects still unknown (TODO) 
				// %NORMAL U
				// if Cls_BC (yval,xval) == (0 || 1 || 2 || 3 || 4 || 5 || 15 || 11 || 12 || 22)
				//      etaC = meanf(Cls_eta(yval,xval+1) , Cls_eta(yval,xval) );
				//      etaC= Cls_eta(  floor((i-1)/lsx)+1, i - (floor((i-1)/lsx)*lsx) ); 
				// elseif Cls_BC (yval,xval) == (6 || 16 || 7 || 8 || 9 || 19 || 10 || 20)
				//      etaCR = meanf(Cls_eta(yval,xval) , Cls_eta(yval,xval-1) );
				//      etaCR= Cls_eta(  floor((i-1)/lsx)+1, i - (floor((i-1)/lsx)*lsx) ); 
				// end       
	            
	            //Start with AT CELL values
	            double etaC =  Cls_eta(yval,xval);  //For middle cell
	            double etaCR = Cls_eta(yval,xval); //For Right BC 

                
                Cell = int(Cls_BC(yval,xval));
                Slope = SlopeGlobalX(yval,xval)*F2;


                //Stencil Definition with multiple cases
                //Middle Cell
                if (Cell == 0)  
                {
                	//No BCs
                }

                //UP CELL (Uphill) (Dirichlet Stag Up because X)   1s
                else if (Cell == 1)  
                {
					//X Direction  
                }
                //UP CELL (Mountain-side) (Dirichlet Stag Up because X)
                else if (Cell == 11)  
                {
					//X Direction  
                }
                //UP CELL (Ocean)  (Neumann Stag Up because X)
                else if (Cell == 21)  
                {
					//X Direction  
                }


                //UP RIGHT CELL (Uphill) (Dirichlet Stag Up and Normal Dirichlet on Right) 2s
                else if (Cell == 2)
                {
					Cls_pF(i) += etaC*(1/hx)*bcRu;
                }                
                //UP RIGHT CELL (Mountain-side) (Dirichlet Stag Up and Normal Dirichlet on Right)
                else if (Cell == 12)
                {
                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcRu = 0 , bcUu = 0 here!!!!) 
                } 
                //UP RIGHT CELL (Ocean) (Neumann Stag Up and Normal Neumann on Right)
                else if (Cell == 22)
                {
				    Cls_pF(i) += etaC*2*bcRun;   //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->  
      			} 

                //UP LEFT CELL (Uphill)  (Dirichlet Stag Up and Normal Dirichlet on Left)  3s 
                else if (Cell == 3)
                {
					Cls_pF(i) += etaC*(1/hx)*bcLu;
                }   
                //UP LEFT CELL (Mountain-side)  (Dirichlet Stag Up and Normal Dirichlet on Left)              
                else if (Cell == 13)
                {
                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcLu = 0 , bcUu = 0 here!!!!)  
                }  
                //UP LEFT CELL (Ocean)  (Neumann Stag Up and Normal Neumann on Left)              
                else if (Cell == 23)
                {
				    Cls_pF(i) += 0.5 * etaC*2*bcLun;   //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->  
                }  

                //RIGHT CELL (Uphill) (Dirichlet Normal because X direction and uphill)  4s 
                else if (Cell == 4)
                {
					Cls_pF(i) += etaC*(1/hx)*bcRu;
                }
                //RIGHT CELL (Mountain-side) (Dirichlet Normal because X direction and mountainside) 
                else if (Cell == 14)
                {
                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcRu = 0 here!!!!) 
                }
                //RIGHT CELL (Ocean) (Neumann Normal Right because X direction and ocean) 
                else if (Cell == 24)
                {
				    Cls_pF(i) += etaC*2*bcRun;   //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->  
                }

                //LEFT CELL (Uphill) (Dirichlet Normal because X direction and uphill) 5s
                else if (Cell == 5)
                {
					Cls_pF(i) += etaC*(1/hx)*bcLu;
                }
                //LEFT CELL (Mountain-side) (Dirichlet Normal because X direction and mountainside) 
                else if (Cell == 15)
                {
                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcLu = 0  here!!!!) 
                }	
                //LEFT CELL (Ocean) (Neumann Normal Left because X direction and ocean) 
                else if (Cell == 25)
                {
				    Cls_pF(i) += 0.5 * etaC*2*bcLun;   //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->  
                }


                //DOWN RIGHT CELL (Uphill) (Dirichlet Stag Down and Normal Dirichlet on Right)   6s
                else if (Cell == 6)
                {
					Cls_pF(i) += etaC*(1/hx)*bcRu;
                }
                //DOWN RIGHT CELL (Mountain-side) (Dirichlet Stag Down and Normal Dirichlet on Right)   
                else if (Cell == 16)
                {
                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcRu = 0 , bcDu = 0 here!!!!)              	
                }
                //DOWN RIGHT CELL (Ocean) (Neumann Stag Down and Normal Neumann on Right)   
                else if (Cell == 26)
                {
				    Cls_pF(i) += etaC*2*bcRun;   //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->  
                }                


                //DOWN LEFT CELL (Uphill)  (Dirichlet Stag Down and Normal Dirichlet on Left) 7s
                else if (Cell == 7)
                {
				    Cls_pF(i) += etaC*(1/hx)*bcLu;
                }
                //DOWN LEFT CELL (Mountain-side)  (Dirichlet Stag Down and Normal Dirichlet on Left) 
                else if (Cell == 17)
                {
                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcLu = 0 , bcDu = 0 here!!!!)    				    			        			    
                }
                //DOWN LEFT CELL (Ocean) (Neumann Stag Down and Normal Neumann on Right)   
                else if (Cell == 27)
                {
				    Cls_pF(i) += 0.5 * etaC*2*bcLun;   //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->                	
                }                   

                //DOWN CELL (Uphill)  (Dirichlet Stag Down because X and uphill)   8s
                else if (Cell == 8)  
                {
					//X Direction  
                }
                //DOWN CELL (Mountain-side)  (Dirichlet Stag Down because X and mountain-side)   8s
                else if (Cell == 18)  
                {
					//X Direction  

                }
                //DOWN CELL (Ocean)  (Neumann Stag Down because X and opening to ocean)   8s
                else if (Cell == 28)
                {
					//X Direction             	
                }


                else if (Cell == -1)
                {
                    //Do nothing 
                }   
                else
                {
                	cout << "Wrong Cell Selection at U for pF" << endl;
                	cout << Cls_BC << endl;
                	exit(1);
                }
            }

            //For V

            for (size_t i = 0; i < lsx*lsy ; ++i)  //Read Left to Right, Down to Up
            {
	            double fi = double(i);                              //CAREFUL!!!!
	            double flsx = double (lsx); 
	            
	            size_t xval =   i - (int(floor(fi/flsx))*lsx);
	            size_t yval =   int(floor(fi/flsx))+0;


	          	////Average Viscosity (OPTIONAL) Effects still unknown (TODO) 
				// if Cls_BC (yval,xval) == (0 || 1 || 2 || 15 || 7 || 8 || 9 || 19 || 10 || 20 || 11 || 12 || 22)
				//     etaC = meanf(Cls_eta(yval+1,xval) , Cls_eta(yval,xval) );
				//     etaC= Cls_eta(  floor((i-1)/lsx)+1, i - (floor((i-1)/lsx)*lsx) );
				// elseif Cls_BC (yval,xval) == (3 || 4 || 5 || 6 || 16)     
				//     etaCU = meanf(Cls_eta(yval,xval) , Cls_eta(yval-1,xval) ); 
				//     etaCU= Cls_eta(  floor((i-1)/lsx)+1, i - (floor((i-1)/lsx)*lsx) );
				// end  
	            
	            //Start with AT CELL values
	            double etaC =  Cls_eta(yval,xval);  //For middle cell
	            double etaCU = Cls_eta(yval,xval); //For Right BC 

                
                Cell = int(Cls_BC(yval,xval));
                Slope = SlopeGlobalY(yval,xval)*F2;

                //Stencil Definition with multiple cases
                if (Cell == 0)  //Middle Cell
                {
                	//No BC Contribution to pF	
                }

                //UP CELL (Uphill) (Dirichlet Normal Up because Y)   1s
                else if (Cell == 1)  
                {
					Cls_pF(i) += etaC*(1/hy)*bcUv;
                } 
                //UP CELL (Mountain-side) (Dirichlet Normal Up because Y)
                else if (Cell == 11)  
                {
                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcUv = 0 here!!!!)    				
                }
                //UP CELL (Ocean)  (Neumann Normal Up because Y)
                else if (Cell == 21)  
                {
					Cls_pF(i) += etaC*2*bcUvn; //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *--> 
                }


                //UP RIGHT CELL (Uphill) (Dirichlet Normal Up and Stag Dirichlet on Right because Y) 2s
                else if (Cell == 2)
                {
					Cls_pF(i) += etaC*(1/hy)*bcUv;    			     
                }                
                //UP RIGHT CELL (Mountain-side) (Dirichlet Normal Up and Stag Dirichlet on Right because Y)
                else if (Cell == 12)
                {
                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcRv = 0 , bcUv = 0 here!!!!)    				    			    
                } 
                //UP RIGHT CELL (Ocean) (Neumann Normal Up and Stag Neumann on Right because Y)
                else if (Cell == 22)
                {
					Cls_pF(i) += etaC*2*bcUvn; //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *--> 
      			} 

                //UP LEFT CELL (Uphill)  (Dirichlet Normal Up and Normal Neumann on Left)  3s 
                else if (Cell == 3)
                {
					Cls_pF(i) += etaC*(1/hy)*bcUv;
                }   
                //UP LEFT CELL (Mountain-side)  (Dirichlet Normal Up and Stag Dirichlet on Left)              
                else if (Cell == 13)
                {
                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcLv = 0 , bcUv = 0 here!!!!)    				    			    
                }  
                //UP LEFT CELL (Ocean)  (Neumann Normal Up and Stag Neumann on Left)              
                else if (Cell == 23)
                {
					Cls_pF(i) += etaC*2*bcUvn; //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *--> 
                }  

                //RIGHT CELL (Uphill) (Dirichlet Stag because Y direction and uphill)  4s 
                else if (Cell == 4)
                {
					//Y Direction
                }
                //RIGHT CELL (Mountain-side) (Dirichlet Stag because Y direction and mountainside) 
                else if (Cell == 14)
                {
					//Y Direction 
                }
                //RIGHT CELL (Ocean) (Neumann Stag Right because Y direction and ocean) 
                else if (Cell == 24)
                {
					//Y Direction   			    
                }

                //LEFT CELL (Uphill) (Dirichlet Stag because Y direction and uphill) 5s
                else if (Cell == 5)
                {
					//Y Direction
                }
                //LEFT CELL (Mountain-side) (Dirichlet Stag because Y direction and mountainside) 
                else if (Cell == 15)
                {
                    //FOR A
					//Y Direction
                }	
                //LEFT CELL (Ocean) (Neumann Stag Left because Y direction and ocean) 
                else if (Cell == 25)
                {
					//Y Direction   
                }

                //DOWN RIGHT CELL (Uphill) (Dirichlet Normal Down and Stag Dirichlet on Right)   6s
                else if (Cell == 6)
                {
					Cls_pF(i) += etaC*(1/hy)*bcDv;
                }
                //DOWN RIGHT CELL (Mountain-side) (Dirichlet Normal Down and Stag Dirichlet on Right)   
                else if (Cell == 16)
                {
                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcRv = 0 , bcDv = 0 here!!!!)    				    			                  	
                }
                //DOWN RIGHT CELL (Ocean) (Neumann Normal Down and Stag Neumann on Right because Y)   
                else if (Cell == 26)
                {
					Cls_pF(i) += 0.5 * etaC*2*bcDvn; //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->  
                }                

                //DOWN LEFT CELL (Uphill)  (Dirichlet Stag Down and Normal Dirichlet on Left) 7s
                else if (Cell == 7)
                {
					Cls_pF(i) += etaC*(1/hy)*bcDv;
                }
                //DOWN LEFT CELL (Mountain-side)  (Dirichlet Stag Down and Normal Dirichlet on Left) 
                else if (Cell == 17)
                {
                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcLv = 0 , bcDv = 0 here!!!!)    				    			    			    
                }
                //DOWN LEFT CELL (Ocean) (Neumann Normal Down and Stag Neumann on Right)   
                else if (Cell == 27)
                {
					Cls_pF(i) += 0.5 * etaC*2*bcDvn; //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->              	
                }                   

                //DOWN CELL (Uphill)  (Dirichlet Stag Down because X and uphill)   8s
                else if (Cell == 8)  
                {
					Cls_pF(i) += etaC*(1/hy)*bcDv;
                }
                //DOWN CELL (Mountain-side)  (Dirichlet Stag Down because X and mountain-side)   8s
                else if (Cell == 18)  
                {
                    //We have 0 velocity on border of mountain-side so all BCs are zero for this part of the grid, so no need to add their contribution on Force Vector (i. e.  bcDv = 0 here!!!!)    				
                }
                //DOWN CELL (Ocean)  (Neumann Stag Down because X and opening to ocean)   8s
                else if (Cell == 28)
                {
					Cls_pF(i) += 0.5 * etaC*2*bcDvn; //Neumann BC. can be free stress (now) or we can add a hydrostatic pressure if needed *-->           	
                }	
                
                else if (Cell == -1)
                {
                    //Do nothing 
                }   
                else
                {
                	cout << "Wrong Cell Selection at V for pF" << endl;
                	cout << Cls_BC << endl;
                	exit(1);
                }

            }    
            //Assemble Matrices to Prepare Outer Cell REMOVAL in Order to Solve Only For Internal Cells
            //Start with Matrix A

            // //Could Use Matrix Cls_BC for elements int(Cls_BC) < 0
            vector<double> Cls_BC_ind(lsx*lsy);  //0 is In, 1 is Out
            size_t count_hit = 0;
            size_t two_count_hit = 0;
            size_t count_hit_full = 0;
            vector<size_t> hit_list;  //For lsx*lsxy
            vector<size_t> hit_list_asc;  //For lsx*lsxy
            vector<size_t> hit_list2; //For 2*lsx*lsy
            vector<size_t> hit_list2_asc; //For 2*lsx*lsy
            vector<size_t> hit_list_full; //For 3*lsx*lsy
            vector<size_t> hit_list_asc_full; //For 3*lsx*lsy
            //VectorXi kill_list;

            // FOR MATRIX B IN COLS
            for (size_t i = 0; i < lsx*lsy; ++i)
            { 
            	double fi = double(i);                              //CAREFUL!!!!
	            double flsx = double (lsx); 
	            size_t xval =   i - (int(floor(fi/flsx))*lsx);
	            size_t yval =   int(floor(fi/flsx))+0;

	            if (  int(Cls_BC(yval,xval)) >= 0 )
	            {
                	Cls_BC_ind[i] = 0;
	            }
	            else
	            {
                    Cls_BC_ind[i] = 1;  //Positive for Outer Cells to Kill
                    count_hit++; //Amount of Rows and/Or Columns to Drop
                    two_count_hit+=2; //Amount of Rows and/Or Columns to Drop for A
                    count_hit_full+=3; //Amount of Rows and/Or Columns to Drop for all Rigidity Matrix
                   
                    hit_list_asc.push_back(i);  //ASCENDING ORDER
                    hit_list2_asc.push_back(i);
                    hit_list_asc_full.push_back(i);
                   
                    hit_list.insert(hit_list.begin(),i);  //REVERSE ORDER
                    hit_list2.insert(hit_list2.begin(),i);
                    hit_list_full.insert(hit_list_full.begin(),i);
	            }
            }

            //FOR MATRIX A IN BOTH DIRECTIONS AND B IN ROWS
            for (size_t i = 0; i < lsx*lsy; ++i)
            { 
            	double fi = double(i);                              //CAREFUL!!!!
	            double flsx = double (lsx); 
	            size_t xval =   i - (int(floor(fi/flsx))*lsx);
	            size_t yval =   int(floor(fi/flsx))+0;

	            if (  int(Cls_BC(yval,xval)) >= 0 )
	            {

	            }
	            else
	            {
                    //REVERSE ORDER
                    hit_list2_asc.push_back(i+(lsx*lsy));
                    hit_list_asc_full.push_back(i+(lsx*lsy));
                    hit_list2.insert(hit_list2.begin(),i+(lsx*lsy));
                    hit_list_full.insert(hit_list_full.begin(),i+(lsx*lsy));
	            }
            }

            //ONE LAST TIME FOR FULL RIGIDITY
            for (size_t i = 0; i < lsx*lsy; ++i)
            { 
            	double fi = double(i);                              //CAREFUL!!!!
	            double flsx = double (lsx); 
	            size_t xval =   i - (int(floor(fi/flsx))*lsx);
	            size_t yval =   int(floor(fi/flsx))+0;

	            if (  int(Cls_BC(yval,xval)) >= 0 )
	            {

	            }
	            else
	            {
                    //REVERSE ORDER
                    hit_list_asc_full.push_back(i+(2*lsx*lsy));
                    hit_list_full.insert(hit_list_full.begin(),i+(2*lsx*lsy));

	            }
            }

            //Print For Revision (Shape of LS-SET)
            MatrixXd Test_BC_ind = MatrixDefVec(lsy, lsx, Cls_BC_ind);
            //MatrixPrint(Test_BC_ind, lsy, lsx);
            cout<< Test_BC_ind <<endl;
            //VecPrint(hit_list);
            //VecPrint(hit_list2);

            double start_delete = omp_get_wtime();
            double end_delete; //MIGHT NOT BE USED
            cout<< "lsx: " << lsx << " lsy: " << lsy  << " count_hit: " << count_hit << endl; 

            //MatrixXd Cls_BT;
            //MatrixXd Cls_Z;

            
            // Do not ERASE anything if count_hit is zero(perfect rectangle) 
            if (count_hit > 0)
            {
	            //For A 

	  

	            //METHOD 1
	            //cout<<"A FIRST Deletion"<<endl;
	            //Cls_A = Row_Deletion(Cls_A, hit_list2);
	            //cout<<"A Second Deletion"<<endl;
	            //MatrixXd Cls_A2 = Cls_A.transpose();
	            //Cls_A2 = Row_Deletion(Cls_A2, hit_list2); 
	            //Cls_A = Cls_A2.transpose();
	            //cout<<"A Deletion Complete"<<endl;

	            //METHOD2

	            // MatrixXd Cls_AA( Cls_A.rows() - two_count_hit, Cls_A.cols() - two_count_hit );
	            // Cls_AA = Row_Del(Cls_A, two_count_hit, two_count_hit,  hit_list2_asc,  hit_list2_asc); 
	            // Cls_A = Cls_AA; 

	            
	            // cout << "Rows Cls_A: " << Cls_A.rows() << " Cols Cls_A: " <<Cls_A.cols() << endl;

	            //For B

	            //METHOD 1
	            //cout<<"B FIRST Deletion"<<endl;
	            //Cls_B = Row_Deletion(Cls_B, hit_list2);
	            //cout<<"B Second Deletion"<<endl;
	            //MatrixXd Cls_B2 = Cls_B.transpose();
	            //Cls_B2 = Row_Deletion(Cls_B2, hit_list); 
	            //Cls_B = Cls_B2.transpose();
	            //cout<<"B Deletion Complete"<<endl;

	            //METHOD 2

	            // MatrixXd Cls_BB( Cls_A.rows() - two_count_hit, Cls_B.cols() - count_hit );
	            // Cls_BB = Row_Del(Cls_B, two_count_hit, count_hit,  hit_list2_asc,  hit_list_asc); 
	            // Cls_B = Cls_BB; 
	 

	            cout << "Rows Cls_B: " << Cls_B.rows() << " Cols Cls_B: " <<Cls_B.cols() << endl;

	            //For BBTT (Symmetric Negative)
	            // Cls_BT = -(Cls_B.transpose());
	            // cout<<"BT Complete"<<endl;
	            // cout << "Rows Cls_BT: " << Cls_BT.rows() << " Cols Cls_B: " <<Cls_BT.cols() << endl;

	            //For ZZ
	            //Cls_Z  = MatrixXd::Zero(Cls_BT.rows(), Cls_B.cols()); 



	            //For u_F, v_F and p_F
	            cout << "Vector Row Deletion" << endl;

	            //METHOD 1
	            //Cls_uF = Row_Deletion_Vec(Cls_uF, hit_list); 
	            //Cls_vF = Row_Deletion_Vec(Cls_vF, hit_list); 
	            //Cls_pF = Row_Deletion_Vec(Cls_pF, hit_list);
	            
	            //METHOD 2

	            // VectorXd Cls_uFa( Cls_uF.size()-count_hit );
	            // VectorXd Cls_vFa( Cls_vF.size()-count_hit );
	            // VectorXd Cls_pFa( Cls_pF.size()-count_hit );

	            //cout << "U: " << endl; 
	            VectorXd Cls_uFa = Row_Del_Vec( Cls_uF, count_hit, hit_list_asc );
	            //cout << "V: " << endl; 
	            VectorXd Cls_vFa = Row_Del_Vec( Cls_vF, count_hit, hit_list_asc );
	            //cout << "P: " << endl; 
	            VectorXd Cls_pFa = Row_Del_Vec( Cls_pF, count_hit, hit_list_asc );

	            Cls_uF = Cls_uFa;
	            Cls_vF = Cls_vFa;
	            Cls_pF = Cls_pFa;

	            cout << "Rows Cls_uF: " << Cls_uF.size() << endl;  
	            cout << "Rows Cls_vF: " << Cls_vF.size() << endl;  
	            cout << "Rows Cls_pF: " << Cls_pF.size() << endl;  


            }
            else
            {
            	// MatrixXd Cls_AA;
            	// MatrixXd Cls_BB;  
            	VectorXd Cls_uFa, Cls_vFa, Cls_pFa;
            	//Cls_BT = -(Cls_B.transpose());
	            //cout<<"BT Complete"<<endl;
	            //Cls_Z  = MatrixXd::Zero(Cls_BT.rows(), Cls_B.cols()); 
            }

            //Assemble System of Equations
            // cout << "Matrix Assembly Part 1" << endl;
            // MatrixXd Cls_M_UP( Cls_A.rows() , Cls_A.cols() + Cls_B.cols() ); 
            // cout << "Matrix Assembly Part 2" << endl;
            // MatrixXd Cls_M_DOWN( Cls_BT.rows() , Cls_A.cols() + Cls_B.cols() );  
            // cout << "Matrix Assembly Part 3" << endl;
            // Cls_M_UP << Cls_A , Cls_B;
            // Cls_M_DOWN << Cls_BT ,  Cls_Z;            

            // cout << "Matrix Assembly Part 4" << endl;
            // MatrixXd Cls_M0(  Cls_A.rows() + Cls_BT.rows() , Cls_A.cols() + Cls_B.cols() );  
            // Cls_M0 << Cls_M_UP , Cls_M_DOWN;
            // MatrixXd Cls_M = (1/(hx*hy))*Cls_M0;  //Apply Scaling Factors to A, B, BT and Z
            // cout << "Rows Cls_M: " << Cls_M.rows() << " Cols Cls_M: " <<Cls_M.cols() << endl;

            cout << "Matrix Assembly Part 5" << endl;
            VectorXd Cls_FF(3*(lsx*lsy - count_hit));
            Cls_FF << Cls_uF, Cls_vF, Cls_pF;

            // cout << "Matrix Assembly Part 6" << endl;            
            
            // //MatrixXd Cls_MM = Cls_M.transpose()*Cls_M; //No need for lscg only for normal cg
            // MatrixXd Cls_MM = Cls_M;  //Faster
            
            //cout << "Rows Cls_MM: " << Cls_MM.rows() << " Cols Cls_MM: " <<Cls_MM.cols() << " Size Cls_F: " << Cls_F.size()<<endl;

            //cout << "Matrix Assembly Part 7" << endl;
            //VectorXd Cls_FF; 
            
            //Cls_FF = (Cls_M.transpose())*Cls_F; //No need for lscg only for normal cg 
            //Cls_FF = Cls_F; //Faster
            
            cout << "Matrix Assembly Part 8" << endl;
            VectorXd Cls_x = VectorXd::Zero( 3*(lsx*lsy - count_hit) );
            cout<< "Size X: " << Cls_x.size() << " Size F: " << Cls_FF.size() << endl; 


            
            // //Make Cls_MM into a Sparse Matrix for Efficiency in the Solver by using Triplets
            // //Cls_MM_Sparse
            // typedef Eigen::SparseMatrix<double> SpMat;
            // typedef Eigen::Triplet<double> T;
            // std::vector<T> tripletList;
            // tripletList.reserve(9*Cls_MM.rows());

            // for (size_t j = 0; j < Cls_MM.rows(); ++j)
            // {
            //  	for (size_t i = 0; i < Cls_MM.cols(); ++i)
            //  	{
            //  		//Assemble Triplet
            //  		if (int(Cls_MM(j,i)) == 0)
            //  		{

            //  		}
            //  		else
            //  		{
            //  		    //Assemble Triplet 
            //  		    tripletList.push_back(T(j,i,Cls_MM(j,i))); 
            //  		}

            //  	}
            // }
            
            cout<<"Sparse Creation Start"<<endl;
            //Prepare Sparse Creation for Solver, all in a single matrix no concatenation for Rigidity Matrix

            SpMat Cls_MM_Sparse(  3*(lsx*lsy) , 3*(lsx*lsy)  );            
            Cls_MM_Sparse.setFromTriplets(tripletList.begin(), tripletList.end());

            //Apply scaling Factor for grid
            Cls_MM_Sparse = (1/(hx*hy))*Cls_MM_Sparse;

            cout<<"Sparse Creation Done"<<endl;

            //Now Create Sparse Matrix to Eliminate elements:

            
            SpMat Cls_MM_Sparse_trim (  3*(lsx*lsy - count_hit) , 3*(lsx*lsy - count_hit)  ); 
            
            //IF WE HAVE IRREGULAR GEOMETRY
            if (count_hit > 0)
            {
            	Eigen::SparseMatrix<double> Cls_clear; 
            	Cls_clear = Del_Sparse(3*lsx*lsy, count_hit_full, hit_list_asc_full);
            	//Modify Cls_MM_Sparse to Get New Trimmed Matrix:
                Cls_MM_Sparse_trim = Eigen::SparseMatrix<double>(Cls_clear.transpose()) * Cls_MM_Sparse * Cls_clear;
            }
            //If we have rectangular geometry and we don't need to delete.
            else
            {
            	Eigen::SparseMatrix<double> Cls_clear;
            	Cls_MM_Sparse_trim = Cls_MM_Sparse;
            }	

            //cout<<Cls_MM_Sparse<<endl;

   //          //Alternatively
   //          SparseMatrix<double> Cls_MM_Sparse(Cls_MM.rows(),Cls_MM.rows());        // default is column major
			// mat.reserve(VectorXi::Constant(cols,9));  //9 Due to Stencil
			// for (size_t j = 0; j < Cls_MM.rows(); ++j)
   //          {
   //           	for (size_t i = 0; i < Cls_MM.cols(); ++i)
   //           	{
   //           		//Assemble Triplet
   //           		if (int(Cls_MM(j,i)) == 0)
   //           		{
   //           			Cls_MM_Sparse.insert(j,i) = Cls_MM(j,i); 	  // alternative: mat.coeffRef(i,j) += v_ij;		
   //           		}
   //           	}
   //          } 		                 
			// Cls_MM_Sparse.makeCompressed();                        // optional


            //Solve for x

   //          //METHOD 1 Conjugate Gradient Using Sparse (ONLY SYMM)
   //          //SpMat<double> A(n,n);
			// // fill A and b
			// Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper > cg;
			// cg.compute(Cls_MM_Sparse);
			// Cls_x = cg.solve(Cls_FF);
			// cout << Cls_x << endl;
			// std::cout << "cg #iterations:     " << cg.iterations() << std::endl;
			// std::cout << "cg estimated error: " << cg.error()      << std::endl;
			// // update b, and solve again
			// //Cls_x = cg.solve(Cls_FF);
			// //cout << Cls_x << endl;


			//METHOD 2 LEASTSQ CG (WORK FOR NON SYMM) (FASTER, NO A'*A Needed.) (SAME SOLUTION AS CG Method 1, but direct)
			// fill A and b
			Eigen::LeastSquaresConjugateGradient<SpMat> lscg;
			lscg.compute(Cls_MM_Sparse_trim);
			Cls_x = lscg.solve(Cls_FF);
			std::cout << "lscg #iterations:     " << lscg.iterations() << std::endl;
			std::cout << "lscg estimated error: " << lscg.error()      << std::endl;
			// update b, and solve again
			//Cls_x = lscg.solve(Cls_FF);
			cout << "Answer Cls_x: " <<endl;
			//cout << Cls_x << endl;



			end_delete = omp_get_wtime();
	        cout << "TIME FOR DELETION AND/OR CONJ. GRAD.: " << (end_delete - start_delete) << endl;
			// cout << "Cls_FF Vector: " <<endl;
			// cout << Cls_FF << endl;

			// //METHOD 3 BiCGSTAB  (Different Answer and not faster ummmm)
			// /* ... fill A and b ... */  
		    // Eigen::BiCGSTAB<SpMat> solver;
		    // solver.compute(Cls_MM_Sparse);
		    // Cls_x = solver.solve(Cls_FF);
		    // std::cout << "BiCGSTAB #iterations:     " << solver.iterations() << std::endl;
		    // std::cout << "BiCGSTAB estimated error: " << solver.error()      << std::endl;
		    // /* ... update b ... */
		    // //Cls_x = solver.solve(Cls_FF);; // solve again
      //       cout << Cls_x << endl;


            //Extract U, V and P
            VectorXd U_part(Cls_uF.size()); 
            U_part	= Cls_x.head(Cls_uF.size());
            VectorXd V_part(Cls_vF.size()); 
            V_part = Cls_x.segment(Cls_uF.size(), Cls_vF.size()); 
            VectorXd P_part(Cls_pF.size());
            P_part = Cls_x.tail(Cls_pF.size());

            MatrixXd U_ls = MatrixXd::Zero(lsy,lsx);
            MatrixXd V_ls = MatrixXd::Zero(lsy,lsx);
            MatrixXd P_ls = MatrixXd::Zero(lsy,lsx);

            size_t shrinker = 0;
            VectorXd U_temp;
            VectorXd V_temp;
            VectorXd P_temp;

            cout << "Reassembling Answer as U, V and P matrices within Cls_BC Domain" << endl;
            //  for (size_t j = 0 ; j < lsy; ++j)
            //  {
            //   	for (size_t i = 0 ; i < lsx; ++i)
            //   	{

            //          U_ls(j,i) = U_part(j*lsx+i);
			        // V_ls(j,i) = V_part(j*lsx+i);
			        // P_ls(j,i) = U_part(j*lsx+i);

            //   	}
            //  }

            for (size_t j = 0 ; j < lsy; ++j)
            {
             	for (size_t i = 0 ; i < lsx; ++i)
             	{

             	   U_temp =  U_part; 
             	   V_temp =  V_part; 
             	   P_temp =  P_part;  
                   if (int(Cls_BC(j,i)) < 0)
                   {
                      //Do Nothing, Out of Bounds.
                   }
                   else
                   {
                   	    shrinker++;
	                   	U_ls(j,i) = U_part(0);
				        V_ls(j,i) = V_part(0);
				        P_ls(j,i) = U_part(0);
                        
                        if (shrinker <= Cls_uF.size() - 1 )  //WARNING THIS IS ASSUMING SAME SIZE FOR ALL!!!!
                        {
					        U_part = U_temp.tail(U_temp.size() - 1);
					        V_part = V_temp.tail(V_temp.size() - 1);
					        P_part = P_temp.tail(P_temp.size() - 1);
					    }    	 
                   }  
             	}
            }
            // cout << "U_ls" << endl;
            // cout << U_ls << endl;
            // //cout << "C_ls" << endl;
            // //cout << Cls_BC << endl;
            // cout << "Shrinker" << endl;
            // cout << shrinker << endl;
            // cout << "V_ls" << endl;
            // cout << V_ls << endl;
            // cout << "P_ls" << endl;
            // cout << P_ls << endl;            

            //Generate Deviatoric Stress Components
            cout << "Deviatoric Stress" << endl;

            //U_ls, V_ls and P_ls are the level set answers. Combine with Cls_BC to get
			//stress and failure surface among other data.
			MatrixXd Cls_tau11 = MatrixXd::Zero(lsy,lsx);
            MatrixXd Cls_tau13 = MatrixXd::Zero(lsy,lsx);
            MatrixXd Cls_tau33 = MatrixXd::Zero(lsy,lsx);

            int CellV;
            cout << "Data" << endl;
            cout << Cls_BC.rows() << " , " << Cls_BC.cols() << endl;
			//Generate deviatoric stresses

			//OLDER
			for (size_t j = 0 ; j < lsy; ++j)
			{	
			    for (size_t i = 0 ; i < lsx; ++i)
			    {    
			       CellV = int(Cls_BC(j,i));

                   if (CellV == 0) //Middle
			       {	
						Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i+1)-U_ls(j,i));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j+1,i)-U_ls(j,i)) + (1/hx)*(V_ls(j,i+1)-V_ls(j,i)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j+1,i)-V_ls(j,i));
			       }
			       else if (CellV == 1) // Up Uphill
                   {
						Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i+1)-U_ls(j,i));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j,i)-U_ls(j-1,i)) + (1/hx)*(V_ls(j,i+1)-V_ls(j,i)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j,i)-V_ls(j-1,i));
                   }
			       else if (CellV == 11) // Up Mountainside
                   {
						Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i+1)-U_ls(j,i));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j,i)-U_ls(j-1,i)) + (1/hx)*(V_ls(j,i+1)-V_ls(j,i)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j,i)-V_ls(j-1,i));
                   }
			       else if (CellV == 21) // Up Ocean
                   {
						Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i+1)-U_ls(j,i));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j,i)-U_ls(j-1,i)) + (1/hx)*(V_ls(j,i+1)-V_ls(j,i)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j,i)-V_ls(j-1,i));
                   }
			       else if (CellV == 2) // Up Right Uphill
                   {
					    Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i)-U_ls(j,i-1));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j,i)-U_ls(j-1,i)) + (1/hx)*(V_ls(j,i)-V_ls(j,i-1)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j,i)-V_ls(j-1,i));
                   }
			       else if (CellV == 12) //Up Right Mountainside 
                   {
					    Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i)-U_ls(j,i-1));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j,i)-U_ls(j-1,i)) + (1/hx)*(V_ls(j,i)-V_ls(j,i-1)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j,i)-V_ls(j-1,i));
                   }
			       else if (CellV == 22) //Up Right Ocean
                   {
					    Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i)-U_ls(j,i-1));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j,i)-U_ls(j-1,i)) + (1/hx)*(V_ls(j,i)-V_ls(j,i-1)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j,i)-V_ls(j-1,i));
                   }
			       else if (CellV == 3) // Up Left Uphill
                   {
                    	Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i+1)-U_ls(j,i));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j,i)-U_ls(j-1,i)) + (1/hx)*(V_ls(j,i+1)-V_ls(j,i)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j,i)-V_ls(j-1,i));
                   }
			       else if (CellV == 13) // Up Left Mountainside
                   {
                    	Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i+1)-U_ls(j,i));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j,i)-U_ls(j-1,i)) + (1/hx)*(V_ls(j,i+1)-V_ls(j,i)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j,i)-V_ls(j-1,i));
                   }
			       else if (CellV == 23) //Up Left Ocean
                   {
                    	Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i+1)-U_ls(j,i));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j,i)-U_ls(j-1,i)) + (1/hx)*(V_ls(j,i+1)-V_ls(j,i)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j,i)-V_ls(j-1,i));
                   }
			       else if (CellV == 4) //Right Uphill
                   {
					    Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i)-U_ls(j,i-1)); 
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j+1,i)-U_ls(j,i)) + (1/hx)*(V_ls(j,i)-V_ls(j,i-1)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j+1,i)-V_ls(j,i));
                   }
			       else if (CellV == 14) // Right Mountainside
                   {
					    Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i)-U_ls(j,i-1)); 
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j+1,i)-U_ls(j,i)) + (1/hx)*(V_ls(j,i)-V_ls(j,i-1)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j+1,i)-V_ls(j,i));
                   }
			       else if (CellV == 24) //Right Ocean
                   {
					    Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i)-U_ls(j,i-1)); 
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j+1,i)-U_ls(j,i)) + (1/hx)*(V_ls(j,i)-V_ls(j,i-1)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j+1,i)-V_ls(j,i));
                   }
			       else if (CellV == 5) //Left Uphill
                   {
					    Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i+1)-U_ls(j,i));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j+1,i)-U_ls(j,i)) + (1/hx)*(V_ls(j,i+1)-V_ls(j,i)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j+1,i)-V_ls(j,i));
                   }
                   else if (CellV == 15) //Left Mountainside
                   {
					    Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i+1)-U_ls(j,i));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j+1,i)-U_ls(j,i)) + (1/hx)*(V_ls(j,i+1)-V_ls(j,i)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j+1,i)-V_ls(j,i));
                   }
			       else if (CellV == 25) //Left Ocean
                   {
					    Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i+1)-U_ls(j,i));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j+1,i)-U_ls(j,i)) + (1/hx)*(V_ls(j,i+1)-V_ls(j,i)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j+1,i)-V_ls(j,i));
                   }
			       else if (CellV == 6) //Down Right Uphill
                   {
					    Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i)-U_ls(j,i-1));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j+1,i)-U_ls(j,i)) + (1/hx)*(V_ls(j,i)-V_ls(j,i-1)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j+1,i)-V_ls(j,i));
                   }
                   else if (CellV == 16) //Down Right Mountainside
                   {
					    Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i)-U_ls(j,i-1));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j+1,i)-U_ls(j,i)) + (1/hx)*(V_ls(j,i)-V_ls(j,i-1)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j+1,i)-V_ls(j,i));
                   }
                   else if (CellV == 26) //Down Right Ocean
                   {
					    Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i)-U_ls(j,i-1));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j+1,i)-U_ls(j,i)) + (1/hx)*(V_ls(j,i)-V_ls(j,i-1)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j+1,i)-V_ls(j,i));
                   }
                   else if (CellV == 7) //Down Left Uphill
                   {
						Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i+1)-U_ls(j,i));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j+1,i)-U_ls(j,i)) + (1/hx)*(V_ls(j,i+1)-V_ls(j,i)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j+1,i)-V_ls(j,i));
                   }
                   else if (CellV == 17) //Down Left Mountainside
                   {
						Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i+1)-U_ls(j,i));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j+1,i)-U_ls(j,i)) + (1/hx)*(V_ls(j,i+1)-V_ls(j,i)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j+1,i)-V_ls(j,i));
                   }
                   else if (CellV == 27) //Down Left Ocean
                   {
						Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i+1)-U_ls(j,i));
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j+1,i)-U_ls(j,i)) + (1/hx)*(V_ls(j,i+1)-V_ls(j,i)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j+1,i)-V_ls(j,i));
                   }
                   else if (CellV == 8) //Down Uphill
                   {
					    Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i+1)-U_ls(j,i)); 
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j+1,i)-U_ls(j,i)) + (1/hx)*(V_ls(j,i+1)-V_ls(j,i)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j+1,i)-V_ls(j,i));
                   }
                   else if (CellV == 18) //Down Mountainside
                   {
					    Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i+1)-U_ls(j,i)); 
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j+1,i)-U_ls(j,i)) + (1/hx)*(V_ls(j,i+1)-V_ls(j,i)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j+1,i)-V_ls(j,i));
                   }
                   else if (CellV == 28) //Down Ocean
                   {
					    Cls_tau11(j,i) = Cls_eta(j,i)*2*(1/hx)*(U_ls(j,i+1)-U_ls(j,i)); 
					    Cls_tau13(j,i) = Cls_eta(j,i)*( (1/hy)*(U_ls(j+1,i)-U_ls(j,i)) + (1/hx)*(V_ls(j,i+1)-V_ls(j,i)) );
					    Cls_tau33(j,i) = Cls_eta(j,i)*2*(1/hy)*(V_ls(j+1,i)-V_ls(j,i));
                   }           
                   else if(CellV == -1) //OUT
                   {
                   		//Do nothing
                   }      
                   else
                   {
                   		cout << "Something fishy about Shear Stress Calc." << endl;
                        exit(1);
                   }   

			    }
			}

			//cout << "Assume Symmetry"  << endl;
			MatrixXd Cls_tau31 = Cls_tau13;	

            //Generate Von Mises and Tau Max Stress (Tau Max will be exported)

            cout << "Von Mises And TauMax"  << endl;
            MatrixXd Cls_vmises0;
            Cls_vmises0 = 1.5 * (Cls_tau11.cwiseProduct(Cls_tau11) + Cls_tau33.cwiseProduct(Cls_tau33) + 2 * Cls_tau13.cwiseProduct(Cls_tau13)) ; 
            MatrixXd Cls_vmises = (Cls_vmises0.array().sqrt()).matrix();
            Cls_taumax = (1/sqrt(3))*Cls_vmises;


            //Generate New Viscosity Using Stresses from U and V

            //Revised Version
			cout << "New Viscosity Calc"  << endl;
			double coefv = (double)-1/3;
			double visc_ct = 0.5*pow(A_rheo, coefv);
			//double visc_ct = 0.5*1;
             
            double hhyy = pow(hy,2.0);
            double hhxx = pow(hx,2.0); 
			
			for (size_t j = 0 ; j < lsy; ++j)
			{	
			    for (size_t i = 0 ; i < lsx; ++i)
			    {    
			       CellV = int(Cls_BC(j,i));
			       

			       if (CellV == 0) //Middle
			       {	
			            Cls_eta(j,i) = visc_ct * pow(  ( (2/(hhxx))* pow( (U_ls(j,i+1)-U_ls(j,i)), 2) + (2/(hhyy))*pow( (V_ls(j+1,i)-V_ls(j,i)) , 2) + (1/(hhyy))*pow( (U_ls(j+1,i)-U_ls(j,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i+1)-V_ls(j,i)) , 2 ) + (2/(hx*hy))*(U_ls(j+1,i)-U_ls(j,i))*(V_ls(j,i+1)-V_ls(j,i)) ) , coefv );
			       }
			       else if (CellV == 1) // Up Uphill
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i+1)-U_ls(j,i)), 2) + (2/(hhyy))*pow( (V_ls(j,i)-V_ls(j-1,i)) , 2) + (1/(hhyy))*pow( (U_ls(j,i)-U_ls(j-1,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i+1)-V_ls(j,i)) , 2) + (2/(hx*hy))*(U_ls(j,i)-U_ls(j-1,i))*(V_ls(j,i+1)-V_ls(j,i)) ) , coefv );
                   }
			       else if (CellV == 11) // Up Mountainside
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i+1)-U_ls(j,i)), 2) + (2/(hhyy))*pow( (V_ls(j,i)-V_ls(j-1,i)) , 2) + (1/(hhyy))*pow( (U_ls(j,i)-U_ls(j-1,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i+1)-V_ls(j,i)) , 2) + (2/(hx*hy))*(U_ls(j,i)-U_ls(j-1,i))*(V_ls(j,i+1)-V_ls(j,i)) ) , coefv );
                   }
			       else if (CellV == 21) // Up Ocean
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i+1)-U_ls(j,i)), 2) + (2/(hhyy))*pow( (V_ls(j,i)-V_ls(j-1,i)) , 2) + (1/(hhyy))*pow( (U_ls(j,i)-U_ls(j-1,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i+1)-V_ls(j,i)) , 2) + (2/(hx*hy))*(U_ls(j,i)-U_ls(j-1,i))*(V_ls(j,i+1)-V_ls(j,i)) ) , coefv );
                   }
			       else if (CellV == 2) // Up Right Uphill
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i)-U_ls(j,i-1)) , 2) + (2/(hhyy))*pow( (V_ls(j,i)-V_ls(j-1,i)) , 2) + (1/(hhyy))*pow( (U_ls(j,i)-U_ls(j-1,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i)-V_ls(j,i-1)) , 2) + (2/(hx*hy))*(U_ls(j,i)-U_ls(j-1,i))*(V_ls(j,i)-V_ls(j,i-1)) ) , coefv );
                   }
			       else if (CellV == 12) //Up Right Mountainside 
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i)-U_ls(j,i-1)) , 2) + (2/(hhyy))*pow( (V_ls(j,i)-V_ls(j-1,i)) , 2) + (1/(hhyy))*pow( (U_ls(j,i)-U_ls(j-1,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i)-V_ls(j,i-1)) , 2) + (2/(hx*hy))*(U_ls(j,i)-U_ls(j-1,i))*(V_ls(j,i)-V_ls(j,i-1)) ) , coefv );
                   }
			       else if (CellV == 22) //Up Right Ocean
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i)-U_ls(j,i-1)) , 2) + (2/(hhyy))*pow( (V_ls(j,i)-V_ls(j-1,i)) , 2) + (1/(hhyy))*pow( (U_ls(j,i)-U_ls(j-1,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i)-V_ls(j,i-1)) , 2) + (2/(hx*hy))*(U_ls(j,i)-U_ls(j-1,i))*(V_ls(j,i)-V_ls(j,i-1)) ) , coefv );
                   }
			       else if (CellV == 3) // Up Left Uphill
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i+1)-U_ls(j,i)), 2) + (2/(hhyy))*pow( (V_ls(j,i)-V_ls(j-1,i)) , 2) + (1/(hhyy))*pow( (U_ls(j,i)-U_ls(j-1,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i+1)-V_ls(j,i)) , 2) + (2/(hx*hy))*(U_ls(j,i)-U_ls(j-1,i))*(V_ls(j,i+1)-V_ls(j,i)) ) , coefv );
                   }
			       else if (CellV == 13) // Up Left Mountainside
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i+1)-U_ls(j,i)), 2) + (2/(hhyy))*pow( (V_ls(j,i)-V_ls(j-1,i)) , 2) + (1/(hhyy))*pow( (U_ls(j,i)-U_ls(j-1,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i+1)-V_ls(j,i)) , 2) + (2/(hx*hy))*(U_ls(j,i)-U_ls(j-1,i))*(V_ls(j,i+1)-V_ls(j,i)) ) , coefv );
                   }
			       else if (CellV == 23) //Up Left Ocean
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i+1)-U_ls(j,i)), 2) + (2/(hhyy))*pow( (V_ls(j,i)-V_ls(j-1,i)) , 2) + (1/(hhyy))*pow( (U_ls(j,i)-U_ls(j-1,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i+1)-V_ls(j,i)) , 2) + (2/(hx*hy))*(U_ls(j,i)-U_ls(j-1,i))*(V_ls(j,i+1)-V_ls(j,i)) ) , coefv );
                   }
			       else if (CellV == 4) //Right Uphill
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i)-U_ls(j,i-1)) , 2) + (2/(hhyy))*pow( (V_ls(j+1,i)-V_ls(j,i)) , 2) + (1/(hhyy))*pow( (U_ls(j+1,i)-U_ls(j,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i)-V_ls(j,i-1)) , 2) + (2/(hx*hy))*(U_ls(j+1,i)-U_ls(j,i))*(V_ls(j,i)-V_ls(j,i-1)) ) , coefv );
                   }
			       else if (CellV == 14) // Right Mountainside
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i)-U_ls(j,i-1)) , 2) + (2/(hhyy))*pow( (V_ls(j+1,i)-V_ls(j,i)) , 2) + (1/(hhyy))*pow( (U_ls(j+1,i)-U_ls(j,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i)-V_ls(j,i-1)) , 2) + (2/(hx*hy))*(U_ls(j+1,i)-U_ls(j,i))*(V_ls(j,i)-V_ls(j,i-1)) ) , coefv );
                   }
			       else if (CellV == 24) //Right Ocean
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i)-U_ls(j,i-1)) , 2) + (2/(hhyy))*pow( (V_ls(j+1,i)-V_ls(j,i)) , 2) + (1/(hhyy))*pow( (U_ls(j+1,i)-U_ls(j,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i)-V_ls(j,i-1)) , 2) + (2/(hx*hy))*(U_ls(j+1,i)-U_ls(j,i))*(V_ls(j,i)-V_ls(j,i-1)) ) , coefv );
                   }
			       else if (CellV == 5) //Left Uphill
                   {
						Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i+1)-U_ls(j,i)), 2) + (2/(hhyy))*pow( (V_ls(j+1,i)-V_ls(j,i)) , 2) + (1/(hhyy))*pow( (U_ls(j+1,i)-U_ls(j,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i+1)-V_ls(j,i)) , 2) + (2/(hx*hy))*(U_ls(j+1,i)-U_ls(j,i))*(V_ls(j,i+1)-V_ls(j,i)) ) , coefv );
                   }
                   else if (CellV == 15) //Left Mountainside
                   {
						Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i+1)-U_ls(j,i)), 2) + (2/(hhyy))*pow( (V_ls(j+1,i)-V_ls(j,i)) , 2) + (1/(hhyy))*pow( (U_ls(j+1,i)-U_ls(j,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i+1)-V_ls(j,i)) , 2) + (2/(hx*hy))*(U_ls(j+1,i)-U_ls(j,i))*(V_ls(j,i+1)-V_ls(j,i)) ) , coefv );
                   }
			       else if (CellV == 25) //Left Ocean
                   {
						Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i+1)-U_ls(j,i)), 2) + (2/(hhyy))*pow( (V_ls(j+1,i)-V_ls(j,i)) , 2) + (1/(hhyy))*pow( (U_ls(j+1,i)-U_ls(j,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i+1)-V_ls(j,i)) , 2) + (2/(hx*hy))*(U_ls(j+1,i)-U_ls(j,i))*(V_ls(j,i+1)-V_ls(j,i)) ) , coefv );
                   }
			       else if (CellV == 6) //Down Right Uphill
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i)-U_ls(j,i-1)) , 2) + (2/(hhyy))*pow( (V_ls(j+1,i)-V_ls(j,i)) , 2) + (1/(hhyy))*pow( (U_ls(j+1,i)-U_ls(j,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i)-V_ls(j,i-1)) , 2) + (2/(hx*hy))*(U_ls(j+1,i)-U_ls(j,i))*(V_ls(j,i)-V_ls(j,i-1)) ) , coefv );
                   }
                   else if (CellV == 16) //Down Right Mountainside
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i)-U_ls(j,i-1)) , 2) + (2/(hhyy))*pow( (V_ls(j+1,i)-V_ls(j,i)) , 2) + (1/(hhyy))*pow( (U_ls(j+1,i)-U_ls(j,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i)-V_ls(j,i-1)) , 2) + (2/(hx*hy))*(U_ls(j+1,i)-U_ls(j,i))*(V_ls(j,i)-V_ls(j,i-1)) ) , coefv );
                   }
                   else if (CellV == 26) //Down Right Ocean
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i)-U_ls(j,i-1)) , 2) + (2/(hhyy))*pow( (V_ls(j+1,i)-V_ls(j,i)) , 2) + (1/(hhyy))*pow( (U_ls(j+1,i)-U_ls(j,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i)-V_ls(j,i-1)) , 2) + (2/(hx*hy))*(U_ls(j+1,i)-U_ls(j,i))*(V_ls(j,i)-V_ls(j,i-1)) ) , coefv );
                   }
                   else if (CellV == 7) //Down Left Uphill
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i+1)-U_ls(j,i)), 2) + (2/(hhyy))*pow( (V_ls(j+1,i)-V_ls(j,i)) , 2) + (1/(hhyy))*pow( (U_ls(j+1,i)-U_ls(j,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i+1)-V_ls(j,i)) , 2) + (2/(hx*hy))*(U_ls(j+1,i)-U_ls(j,i))*(V_ls(j,i+1)-V_ls(j,i)) ) , coefv );
                   }
                   else if (CellV == 17) //Down Left Mountainside
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i+1)-U_ls(j,i)), 2) + (2/(hhyy))*pow( (V_ls(j+1,i)-V_ls(j,i)) , 2) + (1/(hhyy))*pow( (U_ls(j+1,i)-U_ls(j,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i+1)-V_ls(j,i)) , 2) + (2/(hx*hy))*(U_ls(j+1,i)-U_ls(j,i))*(V_ls(j,i+1)-V_ls(j,i)) ) , coefv );
                   }
                   else if (CellV == 27) //Down Left Ocean
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i+1)-U_ls(j,i)), 2) + (2/(hhyy))*pow( (V_ls(j+1,i)-V_ls(j,i)) , 2) + (1/(hhyy))*pow( (U_ls(j+1,i)-U_ls(j,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i+1)-V_ls(j,i)) , 2) + (2/(hx*hy))*(U_ls(j+1,i)-U_ls(j,i))*(V_ls(j,i+1)-V_ls(j,i)) ) , coefv );
                   }
                   else if (CellV == 8) //Down Uphill
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i+1)-U_ls(j,i)), 2) + (2/(hhyy))*pow( (V_ls(j+1,i)-V_ls(j,i)) , 2) + (1/(hhyy))*pow( (U_ls(j+1,i)-U_ls(j,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i+1)-V_ls(j,i)) , 2) + (2/(hx*hy))*(U_ls(j+1,i)-U_ls(j,i))*(V_ls(j,i+1)-V_ls(j,i)) ) , coefv );
                   }
                   else if (CellV == 18) //Down Mountainside
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i+1)-U_ls(j,i)), 2) + (2/(hhyy))*pow( (V_ls(j+1,i)-V_ls(j,i)) , 2) + (1/(hhyy))*pow( (U_ls(j+1,i)-U_ls(j,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i+1)-V_ls(j,i)) , 2) + (2/(hx*hy))*(U_ls(j+1,i)-U_ls(j,i))*(V_ls(j,i+1)-V_ls(j,i)) ) , coefv );
                   }
                   else if (CellV == 28) //Down Ocean
                   {
			            Cls_eta(j,i) = visc_ct*pow( ( (2/(hhxx))*pow( (U_ls(j,i+1)-U_ls(j,i)), 2) + (2/(hhyy))*pow( (V_ls(j+1,i)-V_ls(j,i)) , 2) + (1/(hhyy))*pow( (U_ls(j+1,i)-U_ls(j,i)) , 2 ) + (1/(hhxx))*pow( (V_ls(j,i+1)-V_ls(j,i)) , 2) + (2/(hx*hy))*(U_ls(j+1,i)-U_ls(j,i))*(V_ls(j,i+1)-V_ls(j,i)) ) , coefv );
                   }           
                   else if(CellV == -1) //OUT
                   {
                   		//Do nothing
                   }      
                   else
                   {
                   		cout << "Something fishy about visc_calc" << endl;
                        exit(1);
                   }         

			    }
			}

			Cls_eta0 = Cls_eta;  
            
            cout << "Final Viscosity for Picard Step: " << kt << endl; 
			cout << Cls_eta << endl;  	
            
            //cout << "Cls_BC" << endl; 
			//cout << Cls_BC << endl; 
			//cout << "Cls_taumax" << endl; 
			//cout << Cls_taumax << endl;   


			if (flow_ind == 1)
			{
                //Use velocities in X and Y to modify Thickness (z) and Geometry of Glacier Level Set
                flowAdjust(U_ls, V_ls, minrow, mincol, hy, hx, MatrixGeo, Cls_BC, y_uphill, y_ocean, minrow_abs);  //Modify _lset using this function and _Thickness
			}
			else
			{
				//Do nothing regarding flow (conserve LS Shape) in order to coincide flow shape change and stress so we can evolve damage and fracture
			}			      


            //Find New MaxError to Analyze Glen Convergence
            cout << "Max Error calculation"  << endl;

            U_ls_error =  U_ls;
            V_ls_error =  V_ls;
            P_ls_error =  P_ls;


            if (kt == 0)
            {
			   U_ls_error0 = 0*U_ls_error; 
			   V_ls_error0 = 0*V_ls_error; 
			   P_ls_error0 = 0*P_ls_error; 
			}   

			double Erroru, Errorv, Errorp;

            Erroru = ((U_ls_error-U_ls_error0).array().abs()).matrix().norm() / ((U_ls_error).array().abs()).matrix().norm();
            Errorv = ((V_ls_error-V_ls_error0).array().abs()).matrix().norm() / ((V_ls_error).array().abs()).matrix().norm();
            Errorp = ((P_ls_error-P_ls_error0).array().abs()).matrix().norm() / ((P_ls_error).array().abs()).matrix().norm();

            //Should be maximum, but we will take mean just for now:
            MaxError = (0.33333333333333) * (Erroru + Errorv + Errorp);

            cout<<"Errors for step: "<< kt <<endl;
            cout<<"Error U: "<< Erroru << " Error V: "<< Errorv << " Error P: "<< Errorp << " Error Max: "<< MaxError << endl;

            //Update for later
            U_ls_error0 = U_ls_error;
            V_ls_error0 = V_ls_error;
            P_ls_error0 = P_ls_error;

            //Increase time step
            kt++;

        	//Emergency Exit for While Loop (should report it did not work)
        	if (kt >= Maxtries)
        	{
                cout << "WARNING: Picard Loop was Stopped on Purpose. 'Not reaching tolerance'. \n"; 
                break;
        	} 
        	if (MaxError<tol && MaxError!= 0) 
        	{
        		cout << "Picard Loop CONVERGED after " << kt << " iterations for tolerance: " << tol << endl; 
        	}
        	if (MaxError >= MaxError0)
        	{
        		cout << "WARNING:'Error going up'. \n"; 
        		// cout << "WARNING: Picard iteration was Stopped on Purpose. 'Error going up'. \n"; 
        		// Cls_taumax = Cls_taumax0; 
                // break;
        	}
        	MaxError0 = MaxError;
        	Cls_taumax0 = Cls_taumax; 

        }  // END OF WHILE LOOP

        Stress_Out = Cls_taumax;

        return Stress_Out;
    }

    void flowAdjust(MatrixXd & U_ls, MatrixXd & V_ls, size_t & minrow, size_t & mincol, double & hy, double & hx, MatrixXd & MatrixGeo, MatrixXd & Cls_BC, const size_t y_uphill, const size_t y_ocean, const size_t minrow_abs)
    {
  //   	//Convenient dimension information
  //   	size_t lsy = U_ls.rows();
  //   	size_t lsx = U_ls.cols();
  //   	size_t N = MatrixGeo.rows();
  //   	size_t M = MatrixGeo.cols();

  //   	int CellV;
  //   	double mult  = 0.000000001;  //Multiple to calibrate ice thinning or thickening
  //   	double multf = 10; //Multiple to calibrate ice flow
  //   	double delta = 1;   //Assume 1 for now
  //   	double a_s = 0.0;     //Rate of thickness surface and basal loss and gain  (melt vs snow combination)
  //   	double a_xf = 0.0;     //Rate of x direction front loss ior gain (melt, snow, etc.)
  //   	double a_yf = 0.0;     //Rate of y direction front loss ior gain (melt, snow, etc.)
    
  //       //Store as prior time step state
  //   	MatrixXd Prior_Thick = _grainThickness;
  //   	MatrixXd Prior_Geo = MatrixGeo;

  //   	//MatrixXd Surf_Prior = _globalTerrain + Prior_Thick;   //Surface position =  Terrain elevation + Ice Thickness, if surface changes due to flow, so does thickness, assuming terrain is constant.
  //   	MatrixXd Surf_Prior = MatrixGeo;
  //   	MatrixXd Surf_New = Surf_Prior;    	

  //   	//Matrix of Advance or Loss of distance for cell in x and y
  //   	MatrixXd Moving_x = MatrixXd::Zero(lsy,lsx);
  //   	MatrixXd Moving_y = MatrixXd::Zero(lsy,lsx);


  //       //FLOW MODIFIES SHAPE AND ADVECTS DAMAGE AND THICKNESS

  //       //SHAPE CAN BE CHANGED LAST USING U and V

  //   	//Update Thickness of Glacier using Prior Thickness and Velocity. Use data before flow shape alteration
  //   	//Take into consideration BCs
  //   	//Define areas where thickness will be constant or would all thickness change

  //   	//Update Damage too using U and V and sorrounding Damage, all in all damaged material flows downhill so damage does not stay constant either, but it moves in 2D not in straight line

  //       cout <<"Get new thickness and moving amount for form change: " << endl;
  //       size_t jj, ii; 
  //   	for (size_t j = 0; j < lsy; j++)
  //   	{
  //   		for (size_t i = 0; i < lsx; i++)
  //   		{
		// 	   //Indices for Global Matrices
		// 	   jj = j + minrow;
		// 	   ii = i + mincol;

		// 	   CellV = int(Cls_BC(j,i));

  //              if (CellV == 0) //Middle
		//        {	
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii+1) - Surf_Prior(jj,ii)) / hx) - V_ls(j,i)*((Surf_Prior(jj+1,ii) - Surf_Prior(jj,ii)) / hy ) + a_s);    
		//        }
		//        else if (CellV == 1) // Up Uphill
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii+1) - Surf_Prior(jj,ii)) / hx) - V_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj-1,ii)) / hy ) + a_s);    
  //              }
		//        else if (CellV == 11) // Up Mountainside
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii+1) - Surf_Prior(jj,ii)) / hx) - V_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj-1,ii)) / hy ) + a_s);   
  //              }
		//        else if (CellV == 21) // Up Ocean
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii+1) - Surf_Prior(jj,ii)) / hx) - V_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj-1,ii)) / hy ) + a_s);   
  //                   Moving_y(j,i) =  multf * delta * ( V_ls(j,i) +  a_yf );
  //              }
		//        else if (CellV == 2) // Up Right Uphill
  //              {
		//        		Surf_New(j,i) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj,ii-1)) / hx) - V_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj-1,ii)) / hy ) + a_s);   
  //              }
		//        else if (CellV == 12) //Up Right Mountainside 
  //              {
		//        		Surf_New(j,i) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj,ii-1)) / hx) - V_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj-1,ii)) / hy ) + a_s);   
  //              }
		//        else if (CellV == 22) //Up Right Ocean
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj,ii-1)) / hx) - V_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj-1,ii)) / hy ) + a_s);   
  //                   Moving_x(j,i) =  multf * delta * ( U_ls(j,i) +  a_xf );                        
  //                   Moving_y(j,i) =  multf * delta * ( V_ls(j,i) +  a_yf );
  //              }
		//        else if (CellV == 3) // Up Left Uphill
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii+1) - Surf_Prior(jj,ii)) / hx) - V_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj-1,ii)) / hy ) + a_s);   
  //              }
		//        else if (CellV == 13) // Up Left Mountainside
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii+1) - Surf_Prior(jj,ii)) / hx) - V_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj-1,ii)) / hy ) + a_s);   
  //              }
		//        else if (CellV == 23) //Up Left Ocean
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii+1) - Surf_Prior(jj,ii)) / hx) - V_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj-1,ii)) / hy ) + a_s);   
  //                   Moving_x(j,i) =  multf * delta * ( U_ls(j,i) +  a_xf );    
  //                   Moving_y(j,i) =  multf * delta * ( V_ls(j,i) +  a_yf );
  //              }
		//        else if (CellV == 4) //Right Uphill
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj,ii-1)) / hx) - V_ls(j,i)*((Surf_Prior(jj+1,ii) - Surf_Prior(jj,ii)) / hy ) + a_s);    

  //              }
		//        else if (CellV == 14) // Right Mountainside
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj,ii-1)) / hx) - V_ls(j,i)*((Surf_Prior(jj+1,ii) - Surf_Prior(jj,ii)) / hy ) + a_s);    
  //              }
		//        else if (CellV == 24) //Right Ocean
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj,ii-1)) / hx) - V_ls(j,i)*((Surf_Prior(jj+1,ii) - Surf_Prior(jj,ii)) / hy ) + a_s);    
  //                   Moving_x(j,i) =  multf * delta * ( U_ls(j,i) +  a_xf );    
  //              }
		//        else if (CellV == 5) //Left Uphill
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii+1) - Surf_Prior(jj,ii)) / hx) - V_ls(j,i)*((Surf_Prior(jj+1,ii) - Surf_Prior(jj,ii)) / hy ) + a_s);    
  //              }
  //              else if (CellV == 15) //Left Mountainside
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii+1) - Surf_Prior(jj,ii)) / hx) - V_ls(j,i)*((Surf_Prior(jj+1,ii) - Surf_Prior(jj,ii)) / hy ) + a_s);    
  //              }
		//        else if (CellV == 25) //Left Ocean
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii+1) - Surf_Prior(jj,ii)) / hx) - V_ls(j,i)*((Surf_Prior(jj+1,ii) - Surf_Prior(jj,ii)) / hy ) + a_s);    
  //                   Moving_x(j,i) =  multf * delta * ( U_ls(j,i) +  a_xf );    
  //              }
		//        else if (CellV == 6) //Down Right Uphill
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj,ii-1)) / hx) - V_ls(j,i)*((Surf_Prior(jj+1,ii) - Surf_Prior(jj,ii)) / hy ) + a_s);    
  //              }
  //              else if (CellV == 16) //Down Right Mountainside
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj,ii-1)) / hx) - V_ls(j,i)*((Surf_Prior(jj+1,ii) - Surf_Prior(jj,ii)) / hy ) + a_s);    
  //              }
  //              else if (CellV == 26) //Down Right Ocean
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii) - Surf_Prior(jj,ii-1)) / hx) - V_ls(j,i)*((Surf_Prior(jj+1,ii) - Surf_Prior(jj,ii)) / hy ) + a_s);    
  //                   Moving_x(j,i) =  multf * delta * ( U_ls(j,i) +  a_xf );    
  //                   Moving_y(j,i) =  multf * delta * ( V_ls(j,i) +  a_yf );                   
  //              }
  //              else if (CellV == 7) //Down Left Uphill
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii+1) - Surf_Prior(jj,ii)) / hx) - V_ls(j,i)*((Surf_Prior(jj+1,ii) - Surf_Prior(jj,ii)) / hy ) + a_s);    
  //              }
  //              else if (CellV == 17) //Down Left Mountainside
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii+1) - Surf_Prior(jj,ii)) / hx) - V_ls(j,i)*((Surf_Prior(jj+1,ii) - Surf_Prior(jj,ii)) / hy ) + a_s);    
  //              }
  //              else if (CellV == 27) //Down Left Ocean
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii+1) - Surf_Prior(jj,ii)) / hx) - V_ls(j,i)*((Surf_Prior(jj+1,ii) - Surf_Prior(jj,ii)) / hy ) + a_s);    
  //                   Moving_x(j,i) =  multf * delta * ( U_ls(j,i) +  a_xf );    
  //                   Moving_y(j,i) =  multf * delta * ( V_ls(j,i) +  a_yf );                   
  //              }
  //              else if (CellV == 8) //Down Uphill
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii+1) - Surf_Prior(jj,ii)) / hx) - V_ls(j,i)*((Surf_Prior(jj+1,ii) - Surf_Prior(jj,ii)) / hy ) + a_s);    
  //              }
  //              else if (CellV == 18) //Down Mountainside
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii+1) - Surf_Prior(jj,ii)) / hx) - V_ls(j,i)*((Surf_Prior(jj+1,ii) - Surf_Prior(jj,ii)) / hy ) + a_s);    
  //              }
  //              else if (CellV == 28) //Down Ocean
  //              {
		//        		Surf_New(jj,ii) = Surf_Prior(jj,ii) + mult * delta * ( -U_ls(j,i)*((Surf_Prior(jj,ii+1) - Surf_Prior(jj,ii)) / hx) - V_ls(j,i)*((Surf_Prior(jj+1,ii) - Surf_Prior(jj,ii)) / hy ) + a_s);    
  //                   Moving_y(j,i) =  multf * delta * ( V_ls(j,i) +  a_yf );
  //              }           
  //              else if(CellV == -1) //OUT
  //              {
  //              		//Do nothing
  //              }      
  //              else
  //              {
  //              		cout << "Something fishy about Thick/Flow Calc." << endl;
  //                   exit(1);
  //              }  

  //   		}
  //   	}
  //   	cout << "Just for reference" << endl;
  //   	cout << "U Velocity" << endl;
  //   	cout << U_ls << endl;
  //   	cout << "V Velocity" << endl;
  //   	cout << V_ls << endl;

  //   	_grainThickness = Surf_New; // - _globalTerrain;
  //   	//_grainThickness =  _globalTerrain;


  //   	//Eliminate Negative Thickness
  //       for (size_t j = 0; j < N; j++)  ///CORRECT CORRECT!!!!
  //   	{
  //   		for (size_t i = 0; i < M; i++)
  //   		{
  //   			if (_grainThickness(j,i)<=0)
  //   			{
  //   				_grainThickness(j,i) = 0;	
  //   			}
  //   	    }
  //   	}  

  //       //Get geometry as Matrix
  //   	vector <double> geovec = _lset.getLevelset();  
  //       MatrixXd NewGeo = MatrixDefVec(N, M, geovec);    

  //       //Get temp to Modify
  //       vector<double> Tvec = _grainTemp2D.getLevelset();
  //       MatrixXd NewTemp = MatrixDefVec(N, M, Tvec);          

  //   	//Limit for ice flow advance
  //   	double limit_v = 0.0001;
  //   	double limit_v2 = 0.00000000000000000000000000000000000000001;

  //   	//Indicator of flow
  //   	size_t ind_flow = 0;

  //   	//Original Level Set as Matrix to Bound flow
  //   	vector <double> Geovec0 = _lset0.getLevelset();
  //   	MatrixXd Orig_geo = MatrixDefVec(_lset0.getYdim(), _lset0.getXdim(), Geovec0);
  //       cout << "Orig GEO" << endl;
  //       cout << Orig_geo << endl;

  //   	//Initial Damage at flow induction
  //   	double init_damage = 0.05;


  //       //Update Level Set Using Speed of Ice
  //   	for (size_t j = 0; j < lsy; j++) 
  //   	{
  //   		//Indices for Global Matrices
  //   		jj = j + minrow;
            
  //   		for (size_t i = 0; i < lsx; i++)
  //   		{
  //   	        //Indices for Global Matrices
		// 	    ii = i + mincol;

  //               CellV = int(Cls_BC(j,i));
                
  //               //Prevent invalid cells
  //               if ( (jj>0 && jj <N-1) && (ii>0 &&  ii<M-1) )
  //               {
	 //    	        if (CellV == 21 && Moving_y(j,i) > limit_v )  //Find flow extreme Up Ocean
		//     		{
		//     			MatrixGeo(jj+1,ii) = -1;
		//     			_grainThickness(jj+1,ii) = 1;   //Or inherit???
		//     			_grainDamage(jj+1,ii) = init_damage;      //Or inherit???
		//     			NewTemp(jj+1,ii) = NewTemp(j,i); //Inherit
		//     		}
		// 	       	else if (CellV == 22) //Up Right Ocean
	 //               	{
	 //               		if(Moving_y(j,i) > limit_v)
	 //               		{
		// 	    			MatrixGeo(jj+1,ii) = -1;
		// 	    			_grainThickness(jj+1,ii) = 1;
		// 	    			_grainDamage(jj+1,ii) = init_damage;
		// 	    			NewTemp(jj+1,ii) = NewTemp(j,i);
	 //               		}
	 //               		if(Moving_x(j,i) > limit_v)
	 //               		{
		// 	    			MatrixGeo(jj,ii+1) = -1;
		// 	    			_grainThickness(jj,ii+1) = 1; 
		// 	    			_grainDamage(jj,ii+1) = init_damage; 
		// 	    			NewTemp(jj,ii+1) = NewTemp(j,i);
	 //               		}
	 //                }	
		// 			else if (CellV == 23) //Up Left Ocean
		// 			{
	 //               		if(Moving_y(j,i) > limit_v)
	 //               		{
		// 	    			MatrixGeo(jj+1,ii) = -1;
		// 	    			_grainThickness(jj+1,ii) = 1;
		// 	    			_grainDamage(jj+1,ii) = init_damage; 
		// 	    			NewTemp(jj+1,ii) = NewTemp(j,i);
	 //               		}
	 //               		if(Moving_x(j,i) < -limit_v)
	 //               		{
		// 	    			MatrixGeo(jj,ii-1) = -1;
		// 	    			_grainThickness(jj,ii-1) = 1;  
		// 	    			_grainDamage(jj,ii-1) = init_damage; 
		// 	    			NewTemp(jj,ii-1) = NewTemp(j,i);
	 //               		}
		// 			}
		// 			else if (CellV == 24 && Moving_x(j,i) > limit_v) //Right Ocean
		// 			{
		// 				if (minrow_abs > y_ocean)
	 //               		{
	 //               			if ( Orig_geo(jj,ii+1) < 0)
	 //               			{
		// 		    			MatrixGeo(jj,ii+1) = -1;
		// 		    			_grainThickness(jj,ii+1) = 1; 
		// 		    			_grainDamage(jj,ii+1) = init_damage;   
		// 		    			NewTemp(jj,ii+1) = NewTemp(j,i);
	 //               			}
	 //               		}
	 //               		else
	 //               		{
		// 	    			MatrixGeo(jj,ii+1) = -1;
		// 	    			_grainThickness(jj,ii+1) = 1; 
		// 	    			_grainDamage(jj,ii+1) = init_damage;   
		// 	    			NewTemp(jj,ii+1) = NewTemp(j,i);
	 //               		}	

		// 			}
		// 			else if (CellV == 25 && Moving_x(j,i) < -limit_v) //Left Ocean
		// 			{
		// 				if (minrow_abs > y_ocean)
	 //               		{
	 //               			if ( Orig_geo(jj,ii-1) < 0)
	 //               			{
		// 	               		MatrixGeo(jj,ii-1) = -1;
		// 		    			_grainThickness(jj,ii-1) = 1;
		// 		    			_grainDamage(jj,ii-1) = init_damage;   
		// 		    			NewTemp(jj,ii-1) = NewTemp(j,i); 
	 //               			}
	 //               		}
	 //               		else
	 //               		{
		//                		MatrixGeo(jj,ii-1) = -1;
		// 	    			_grainThickness(jj,ii-1) = 1;
		// 	    			_grainDamage(jj,ii-1) = init_damage;   
		// 	    			NewTemp(jj,ii-1) = NewTemp(j,i); 
	 //               		}	
		// 			}
		// 			else if (CellV == 26) //Down Right Ocean
		// 			{
	 //               		if (minrow_abs > y_ocean)
	 //               		{
	 //               			if ( Orig_geo(jj-1,ii) < 0)
	 //               			{
		// 	            		//if(Moving_y(j,i) < -limit_v)
		// 	               		//{
		// 			    			MatrixGeo(jj-1,ii) = -1;
		// 			    			_grainThickness(jj-1,ii) = 1;
		// 			    			_grainDamage(jj-1,ii) = init_damage; 
		// 			    			NewTemp(jj-1,ii) = NewTemp(j,i);
		// 	               		//}
		// 		    		}	
		// 		    		if ( Orig_geo(jj,ii+1) < 0)
	 //               			{
		// 	               		//if(Moving_x(j,i) > limit_v)
		// 	               		//{
		// 			    			MatrixGeo(jj,ii+1) = -1;
		// 			    			_grainThickness(jj,ii+1) = 1;  
		// 			    			_grainDamage(jj,ii+1) = init_damage; 
		// 			    			NewTemp(jj,ii+1) = NewTemp(j,i);
		// 	               		//}
		// 		    		}	
		// 	    		}	
		// 	    		else
		// 	    		{	
		// 		    		//if(Moving_y(j,i) < -limit_v)
		//                		//{
		// 		    			MatrixGeo(jj-1,ii) = -1;
		// 		    			_grainThickness(jj-1,ii) = 1;
		// 		    			_grainDamage(jj-1,ii) = init_damage; 
		// 		    			NewTemp(jj-1,ii) = NewTemp(j,i);
		//                		//}
		//                		//if(Moving_x(j,i) > limit_v)
		//                		//{
		// 		    			MatrixGeo(jj,ii+1) = -1;
		// 		    			_grainThickness(jj,ii+1) = 1;  
		// 		    			_grainDamage(jj,ii+1) = init_damage; 
		// 		    			NewTemp(jj,ii+1) = NewTemp(j,i);
		//                		//}
		// 	    		}
		// 			}
		// 			else if (CellV == 27) //Down Left Ocean
		// 			{
		//                	if (minrow_abs > y_ocean)
	 //               		{
	 //               			if ( Orig_geo(jj-1,ii) < 0)
	 //               			{
		//                			//if(Moving_y(j,i) < -limit_v)
		// 	               		//{
		// 			    			MatrixGeo(jj-1,ii) = -1;
		// 			    			_grainThickness(jj-1,ii) = 1;
		// 			    			_grainDamage(jj-1,ii) = init_damage; 
		// 			    			NewTemp(jj-1,ii) = NewTemp(j,i);
		// 	               		//}
  //                           }	
  //                           if ( Orig_geo(jj,ii-1) < 0)
	 //               			{				    			
		// 	               		//if(Moving_x(j,i) < -limit_v)
		// 	               		//{
		// 			    			MatrixGeo(jj,ii-1) = -1;
		// 			    			_grainThickness(jj,ii-1) = 1;  
		// 			    			_grainDamage(jj,ii-1) = init_damage; 
		// 			    			NewTemp(jj,ii-1) = NewTemp(j,i);
		// 	               		//} 
		// 	               	}	 
	 //               		}	
	 //               		else
	 //               		{
	 //               			//if(Moving_y(j,i) < -limit_v)
		//                		//{
		// 		    			MatrixGeo(jj-1,ii) = -1;
		// 		    			_grainThickness(jj-1,ii) = 1;
		// 		    			_grainDamage(jj-1,ii) = init_damage; 
		// 		    			NewTemp(jj-1,ii) = NewTemp(j,i);
		//                		//}
		//                		//if(Moving_x(j,i) < -limit_v)
		//                		//{
		// 		    			MatrixGeo(jj,ii-1) = -1;
		// 		    			_grainThickness(jj,ii-1) = 1;  
		// 		    			_grainDamage(jj,ii-1) = init_damage; 
		// 		    			NewTemp(jj,ii-1) = NewTemp(j,i);
		//                		//}  
	 //               		}              
		// 			}
		// 			//else if ( CellV == 28   &&   Moving_y(j,i) < -limit_v2) //Down Ocean  but also other down flow could exist
		// 			else if ( CellV == 28 ) //Down Ocean  but also other down flow could exist, regardless of speed on frontier cells
		// 			{
		// 			    if (minrow_abs > y_ocean)
	 //               		{
	 //               			if ( Orig_geo(jj-1,ii) < 0)
	 //               			{
		// 	               		MatrixGeo(jj-1,ii) = -1;
		// 		    			_grainThickness(jj-1,ii) = 1;
		// 		    			_grainDamage(jj-1,ii) = init_damage; 
		// 		    			NewTemp(jj-1,ii) = NewTemp(j,i);
		// 		    		}	
	 //               		}
	 //               		else
	 //               		{
		//                		MatrixGeo(jj-1,ii) = -1;
		// 	    			_grainThickness(jj-1,ii) = 1;
		// 	    			_grainDamage(jj-1,ii) = init_damage; 
		// 	    			NewTemp(jj-1,ii) = NewTemp(j,i);
	 //               		}
		// 			} 
		// 			else
		// 			{
		// 			    //Do nothing.
		// 			}
		// 		}
		// 		else
		// 		{
		// 			cout << "WARNING: Limit of LS Grid Reached, reduce flow speed or increase grid size" << endl;
		// 		}	
  //   	    }
  //   	}    

  //   	//Soften before returning
  //   	size_t maxt = 50;
  //   	MatrixXd MatrixGeoS = LS_Soft_Simp(MatrixGeo, maxt);

		// //SOFTENED THICKNESS, DAMAGE AND TEMPERATURE TOO
		// for (size_t jj = 0; jj < N; jj++)
		// {
		// 	for (size_t ii = 0; ii < M; ii++)
		// 	{ 
		// 		if (MatrixGeoS(jj,ii) > 0)
		// 		{  
		// 			 NewTemp(jj,ii) = NewTemp(0,0); //OuterTemp Equiv,
		// 			 _grainDamage(jj,ii) = -1;
		// 			 _grainThickness(jj,ii) = 0;
		// 		}
		// 		if (MatrixGeoS(jj,ii) <= 0 && MatrixGeo(jj,ii)>0)
		// 		{
		// 			 NewTemp(jj,ii) = NewTemp(jj+1,ii);
		// 			 _grainDamage(jj,ii) = 0;
		// 			 _grainThickness(jj,ii) = 1;
		// 		}	
		// 	}     
		// }

  //   	//Return geometry to update Geometry Level Set and Temp. too  -->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //       vector<double> newgeometry = MatrixtoVec(MatrixGeoS);
  //       vector<double> newtemperature = MatrixtoVec(NewTemp);
  //       change_Full_lset(newgeometry, M, N);
  //       change_Full_Temp_lset(newtemperature, M, N);

  //       //Update all other Geometric Properties due to flow deformation
        
  //       // compute mass and center of mass (assuming a density of 1)
		// double vol = 0;
		// //double eps = 1.5+meltTemp;   
		// double eps = 0.015;            
		// Vector2d cm = Vector2d (0.,0.);   
		// vector<double> TTvec=_lset.getLevelset();                   

	 //    //Binarization and Smoothing of Level Set (1 in and 0 out)
		// for (size_t j = 0; j < _lset.getYdim(); j++) {
		// 		for (size_t i = 0; i < _lset.getXdim(); i++) {   //Substitute of the Heaviside Function to Identify Inside and Outside of Grain
		// 			if (-TTvec[(j*_lset.getXdim())+i] > eps) {  
  //                      TTvec[(j*_lset.getXdim())+i] = 1; 
  //                   }
		// 			else if (-TTvec[(j*_lset.getXdim())+i] < -eps) {
  //                      TTvec[(j*_lset.getXdim())+i] = 0;
		// 			}	
		// 			else{
  //                      TTvec[(j*_lset.getXdim())+i] = 0.5*(1. + TTvec[(j*_lset.getXdim())+i]/eps + sin(M_PI*TTvec[(j*_lset.getXdim())+i]/eps)/M_PI); //0.5*(1. + vec(i)/eps + sin(M_PI*vec(i)/eps)/M_PI);
		// 			}	
		// 		}
		// }

		// //For the Source Grain Our Center of Mass has to be shifted so the grain does not shift, this is because the Source Grain DOES NOT MOVE only whatever breaks from it
		// //double cmxorig = _cmLset(0);
		// //double cmyorig = _cmLset(1);
		
  //       Vector2d g1cm = _cmLset; //Save older cm for position, before updating it

		// for (size_t j = 0; j < _lset.getYdim(); j++) {
		// 	for (size_t i = 0; i < _lset.getXdim(); i++) {
		// 		vol+=TTvec[(j*_lset.getXdim())+i];     //vol += heavi(iter);
		// 		cm(0) += TTvec[(j*_lset.getXdim())+i]*(double)i;    //cm(0) += heavi(iter)*(double)i;
		// 		cm(1) += TTvec[(j*_lset.getXdim())+i]*(double)j;    //cm(1) += heavi(iter)*(double)j;
		// 		//iter++;
		// 	}
		// }

		// cm /= vol;

		// //Preservation of Density
		// if (_density == 0.91e-6)
		// {}
		// else{
		// 	//_density = 0.91e-6;
		// }

		// _mass = vol*_density*0.89;  //info._mass = vol;   //Effect of mass?? FACTOR 1/0.9 Increase Consistently
  //       _temper=_mass;

  //       //Change in Center of Mass due to Change of Level Set
  //       //_cmLset //Will it affect possition or what will happen depending on symmetry of melting

  //       //Fixed Center of Mass
  //       //cm(0) = cmxorig;
  //       //cm(1) = cmyorig;
  //       _cmLset = cm; //info._cm = cm;   //Pretty much the same with slight alterations Update LSET
        

  //       //Change in Moment of Inertia
	 //    // iter = 0;
		// // for (size_t k = 0; k < lset.getZdim(); k++) {
		// double Inertia=0;
		// double rx=0;
		// double ry=0;

		// for (size_t j = 0; j < _lset.getYdim(); j++) {
		// 	for (size_t i = 0; i < _lset.getXdim(); i++) {
		// 		rx = (double)i - _cmLset(0);
		// 		ry = (double)j - _cmLset(1);
		// 		//rz = (double)k - cm(2);
		// 		Inertia += TTvec[(j*_lset.getXdim())+i]*(rx*rx + ry*ry); //Effect of mass???? FACTOR 1/0.9 Increase Consistently
		// 		//I(0,1) += -heavi2(iter)*rx*ry;
		// 	}
		// }
  //       _momentInertia=Inertia*0.89;

		// vector <Vector2d> pneworder(_npoints);
  //       PointsGen(_lset, pneworder, _npoints, _cmLset ); 	

  //       //cout<<"Point Generation Adjust"<<endl;
  //       //_fracProps.PointsGen(_lset , pneworder, _npoints, _cmLset);

  //       const double cos1 = cos(_theta);
  //       const double sin1 = sin(_theta);
  //       double cosine1 = cos1;
  //       double sine1 = sin1;
  //       //Update Position
  //       Vector2d positionnew = _position-Vector2d((g1cm(0)-_cmLset(0))*cosine1 - (g1cm(1)-_cmLset(1))*sine1,
  //                                                                (g1cm(0)-_cmLset(0))*sine1 + (g1cm(1)-_cmLset(1))*cosine1);

        
  //       //Update New Points and Move to Centroid Position (THE SLOWEST THE CHANGE THE SMOOTHER THE TRANSITION)
		// for (size_t i = 0; i < _npoints; i++) {
  //          //_pointList[i] = pneworder[i]+_position-_cmLset;   ///Actually Adding Position is Necessary to ReLocate them
		//    //pneworder[i] = pneworder[i]+positionnew-_cmLset-g1cm;  //Works with no change in Pos of grain 0
		//    //pneworder[i] = pneworder[i] + positionnew - _cmLset - _position + g1cm ;// - g1cm;

  //          //pneworder[i] = pneworder[i] - g1cm + _cmLset; //SECOND
  //          pneworder[i] = pneworder[i] + (  - (_position  -   g1cm)   +   (positionnew  -     _cmLset)   );    //WINNING SO FAR               
  //          //pneworder[i] = pneworder[i]  +   _cmLset;
  //          //pneworder[i] = pneworder[i]; //THIRD
		//    //pneworder[i] = pneworder[i]+_position-_cmLset;
		//    //_pointList[i]=pneworder[i]-_cmLset;
		// }
		// _pointList = pneworder;

		// //????
		// _position = positionnew;  //WE DO NOT NEED TO UPDATE POSITION OF GRAIN ZERO BECAUSE IT SHOULD STAY FIXED AS IT REPRESENT DETACHABLE CLIFF

  //       //End

  //       //Change in Bbox Radius of Grain
  //       //_radius
	
		// double maxradius = 0;
		// // double minradius = (_pointList[0]-_position).norm();
		// for (size_t i = 0; i < _npoints; i++) {
		//      if ((_pointList[i]-_position).norm() > maxradius) {     //Should we subtract position to find real radius?
		//          maxradius = (_pointList[i]-_position).norm();
		//         //cout << maxradius << endl;
		//      }
  //      //       if((_pointList[i]-_position).norm() < minradius) {
  //      //          minradius = (_pointList[i]-_position).norm();
  //      // //          //cout << maxradius << endl;
  //      //       } 

		// }
		// _radius=maxradius;  // -->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  //       //Just for view
  //       //_grainThickness = Moving_x;
  //       //_grainThickness = Moving_y;
  //       //_grainThickness.block( minrow, mincol, lsy, lsx ) = Moving_x;
  //       //_grainThickness.block( minrow, mincol, lsy, lsx ) = Moving_y;
  //       //_grainThickness.block( minrow, mincol, lsy, lsx ) = V_ls;
  //       _grainThickness.block( minrow, mincol, lsy, lsx ) = Cls_BC;
    }

    //Basic Locator Classifier
    double basicClass (int value, size_t & j, const size_t & y_uphill_cut, const size_t & y_ocean_cut)
    {
    	if (j >= y_uphill_cut)
    	{
    		//Value stays the same.

    	}
    	else if (j >= y_ocean_cut && j < y_uphill_cut)
    	{
    		//Increases by 10 for mountain side
    		value += 10;

    	}
    	else if (j < y_ocean_cut)
    	{
    		//Increases by 20 for ocean
            value += 20;
    	}
    	else
    	{
    		//ERROR MESSAGE
    		cout << "ERROR in basicClass." << endl;
    		exit(1);
    	}

    	if ( value == 14 &&  y_ocean_cut == 0) //Convert other intermediate regimes to OPEN Neumann BC due to cut
    	{
    		value = 24;
    	}
    	if ( value == 15 &&  y_ocean_cut == 0) //Convert other intermediate regimes to OPEN Neumann BC due to cut
    	{
    		value = 25;
    	}
    	if ( value == 16 &&  y_ocean_cut == 0) //Convert other intermediate regimes to OPEN Neumann BC due to cut
    	{
    		value = 26;
    	}
    	if ( value == 17 &&  y_ocean_cut == 0) //Convert other intermediate regimes to OPEN Neumann BC due to cut
    	{
    		value = 27;
    	}
    	if ( value == 18 &&  y_ocean_cut == 0) //Convert other intermediate regimes to OPEN Neumann BC due to cut
    	{
    		value = 28;
    	}
    	
    	return (double) value;

    }

    //Eigen-Matrix created using a lset vector
    MatrixXd  MatrixDefVec(const size_t & Rows, const size_t & Cols, vector <double> & Svec) {
        size_t N = Rows;
    	size_t M = Cols;
    	MatrixXd MatrixS(N,M);

        for(size_t j = 0; j < N; ++j) //y
		{
		    //Row row(M);

		    for(size_t i = 0; i < M; ++i) //x
		    {
		        //row[i] = Svec[(j*M)+i];  //Add element at each column per row  //Extract each LS Element using this Mechanism
		        MatrixS(j,i) = Svec[(j*M)+i]; 
		    }

		    //MatrixS.push_back(row); // push each row after you fill it
		}
  		return MatrixS;
    }

    //Eigen-Matrix to Vector for Lset
    vector <double> MatrixtoVec(MatrixXd & Mat) //const size_t & Rows, const size_t & Cols)
    {
        size_t N = Mat.rows();
    	size_t M = Mat.cols();
        vector<double> VEC(N*M);


        for(size_t j = 0; j < N; ++j) //y
		{
		    for(size_t i = 0; i < M; ++i) //x
		    {
		       VEC[(j*M)+i] = Mat(j,i);  //Add element at each column per row  //Extract each LS Element using this Mechanism
		    }
		}  
		return VEC;    
    }

    //Eigen-Matrix softening function
    MatrixXd LS_Soft(MatrixXd & Cls_Tag0, size_t & weird_ind, size_t & ttries, const size_t & maxttries, const size_t & y0Db, const size_t & ywsealevel) //const size_t & SZY, const size_t & SZX)
    {
	    
	    size_t SZY = Cls_Tag0.rows();
    	size_t SZX = Cls_Tag0.cols();
	    MatrixXd Cls_Tag = Cls_Tag0; 

	    while (weird_ind == 1) 
	    {
	       weird_ind = 0; 
	       for (size_t j = 0 ; j < SZY ;  ++j)  //Priorly ii
	       {
	           for (size_t i = 0 ; i < SZX ;  ++i)   //Priorly jj
	           {

		            if ( int(Cls_Tag0(j,i)) == 0 &&  i>0  &&  j>y0Db  &&  j<(SZY-1)  &&  i<(SZX-1)  ) {   //CAREFUL!!!!
		                 
		                //Define value to replace (Air or Water)
		                size_t switchvalue;
		                if (j <= ywsealevel) { 
  							switchvalue = 1;
		                }  
		                else{
		                    switchvalue = 2;
		                }   
		                 
		                //Remove Weird Right
		                if ( int(Cls_Tag0(j+1,i)) != 0  &&  int(Cls_Tag0(j-1,i)) != 0  &&  int(Cls_Tag0(j,i+1)) != 0) {
		                    Cls_Tag(j,i) = switchvalue;
		                    weird_ind = 1; 
		                }    
		                 
		                //Remove Weird Left
		                if ( int(Cls_Tag0(j+1,i)) != 0 && int(Cls_Tag0(j-1,i)) != 0 && int(Cls_Tag0(j,i-1)) != 0) {
		                    Cls_Tag(j,i) = switchvalue;
		                    weird_ind = 1; 
		                }    
		                 
		                //Remove Weird Up
		                if (int(Cls_Tag0(j+1,i)) != 0 && int(Cls_Tag0(j,i+1)) != 0 && int(Cls_Tag0(j,i-1)) != 0) {
		                    Cls_Tag(j,i) = switchvalue;
		                    weird_ind = 1; 
		                }   
		                //Remove Weird Down 
		                if (int(Cls_Tag0(j-1,i)) != 0 && int(Cls_Tag0(j,i+1)) != 0 && int(Cls_Tag0(j,i-1)) != 0) {
		                    Cls_Tag(j,i) = switchvalue;
		                    weird_ind = 1; 
		                }  
		                //Remove Floating Cells
		                if (int(Cls_Tag0(j+1,i)) != 0 && int(Cls_Tag0(j-1,i)) != 0 && int(Cls_Tag0(j,i+1)) != 0 && int(Cls_Tag0(j,i-1)) != 0) {
		                    Cls_Tag(j,i) = switchvalue;
		                    weird_ind = 1; 
		                }  
		            }    
	            }
	        } 

	       //Saving Grace to Avoid Getting Stuck in While Loop, but it should not be used at all
	       ttries = ttries + 1;
	       if (ttries > maxttries) {
	          weird_ind = 0; 
	       } 

	    }
	    return Cls_Tag;      

    }

	//Eigen-Matrix softening function
	MatrixXd LS_Soft_Simp(MatrixXd & Cls_Tag0, const size_t & maxttries) //const size_t & SZY, const size_t & SZX)
	{

		size_t SZY = Cls_Tag0.rows();
		size_t SZX = Cls_Tag0.cols();
		size_t ttries = 0;
		size_t weird_ind = 1;
		MatrixXd Cls_Tag = Cls_Tag0; 

		//Define value to replace (Out)
		size_t switchvalue;
		switchvalue = 1;

		while (weird_ind == 1) 
		{
		   weird_ind = 0; 
		   for (size_t j = 0 ; j < SZY ;  ++j)  //Priorly ii
		   {
		       for (size_t i = 0 ; i < SZX ;  ++i)   //Priorly jj
		       {

		          if ( int(Cls_Tag(j,i)) == -1 ) {   //CAREFUL!!!!
		               
		               
		              //Remove Weird Right
		              if ( (int)Cls_Tag(j+1,i) == 1  &&  (int)Cls_Tag(j-1,i) == 1  &&  (int)Cls_Tag(j,i+1) == 1) {
		                  Cls_Tag(j,i) = switchvalue;
		                  weird_ind = 1; 
		                }    
		               
		              //Remove Weird Left
		              if ( int(Cls_Tag(j+1,i)) == 1 && int(Cls_Tag(j-1,i)) == 1 && int(Cls_Tag(j,i-1)) == 1) {
		                  Cls_Tag(j,i) = switchvalue;
		                  weird_ind = 1; 
		                }    
		               
		              //Remove Weird Up
		              if (int(Cls_Tag(j+1,i)) == 1 && int(Cls_Tag(j,i+1)) == 1 && int(Cls_Tag(j,i-1)) == 1) {
		                  Cls_Tag(j,i) = switchvalue;
		                  weird_ind = 1; 
		                }   
		              //Remove Weird Down 
		              if (int(Cls_Tag(j-1,i)) == 1 && int(Cls_Tag(j,i+1)) == 1 && int(Cls_Tag(j,i-1)) == 1) {
		                  Cls_Tag(j,i) = switchvalue;
		                  weird_ind = 1; 
		                }  
		              //Remove Floating Cells
		              if (int(Cls_Tag(j+1,i)) == 1 && int(Cls_Tag(j-1,i)) == 1 && int(Cls_Tag(j,i+1)) == 1 && int(Cls_Tag(j,i-1)) == 1) {
		                  Cls_Tag(j,i) = switchvalue;
		                  weird_ind = 1; 
		                }  
		            }    
		        }
		    } 

		   //Saving Grace to Avoid Getting Stuck in While Loop, but it should not be used at all
		   ttries = ttries + 1;
		   if (ttries > maxttries) {
		      weird_ind = 0; 
		   } 

		}

		//Also eliminate single holes
		for (size_t j = 1 ; j < SZY-1 ;  ++j)  //Priorly ii
		{
	       for (size_t i = 1 ; i < SZX-1 ;  ++i)   //Priorly jj
	       {

	       		if ( (int) Cls_Tag(j,i) == 1  && ( (int) Cls_Tag(j+1,i) == -1 && (int) Cls_Tag(j-1,i) == -1 && (int) Cls_Tag(j,i+1) == -1 && (int) Cls_Tag(j,i-1) == -1 ) )
	       		{
	       			Cls_Tag(j,i) = -1;
	       		}	

	       		//FIND BIG HOLES PART ( TODO )



	       }
	    }   
		       	
		return Cls_Tag;      

	}

    void VecPrint(vector<size_t> & Vec)
    {   
    	size_t SZY = Vec.size(); 
     	for (size_t i = 0 ; i < SZY ;  ++i)  //Priorly ii
     	{
     		cout << Vec[i] << ";\n";
     	}
    }


    //SPARSE MATRIX MULTIPLICATION TO DELETE ROWS AND COLUMNS FOR FULL SQUARE RIGIDITY (size is either rows or cols for square matrix)
    Eigen::SparseMatrix<double> Del_Sparse (size_t size_Matrix, size_t & minus_cells, vector<size_t> & hit_list_asc )
    {
    	Eigen::SparseMatrix<double> Out( size_Matrix, size_Matrix-minus_cells );
        
        typedef Eigen::Triplet<double> T;
    	std::vector<T> tripletList;
        //tripletList.reserve(9*lsx*lsx*lsy*lsy);  //Check reserve effect on performance
        tripletList.reserve(size_Matrix); 

    	size_t cellN = 0;
        size_t hitdex = 0;

        for (size_t i = 0; i < size_Matrix; ++i)
        {
            if ( i != hit_list_asc[hitdex] )
            {
                 //Triplet Modification
                 tripletList.push_back( T( i, cellN, 1) ) ; 
                 //Shift columns as we go
                 cellN++; 
            } 
            else
            {
            	//Move to next skippable cell 
                hitdex++;
            }

        }

        Out.setFromTriplets(tripletList.begin(), tripletList.end());
        return Out;
    }

    //USE MATRIX MULTIPLICATION TO DELETE ROWS (AND COLUMNS) at a matrix (ACTUALLY MATRIX MULT DELETES COLS)
    //Hit list must be in ascending order
    MatrixXd Row_Del(MatrixXd & MatBase, size_t & minus_rows, size_t & minus_cols, vector<size_t> & hit_list_rows, vector<size_t> & hit_list_cols )
    {
     	MatrixXd Terminator1 = MatrixXd::Zero ( MatBase.cols() , MatBase.cols() - minus_cols ); 
        MatrixXd Terminator2 = MatrixXd::Zero ( MatBase.rows() , MatBase.rows() - minus_rows );

        //Terminator1
        size_t colN = 0;
        size_t hitdex = 0;

        for (size_t i = 0; i < Terminator1.rows(); ++i)
        {
            if ( i != hit_list_cols[hitdex] )
            {
                Terminator1(i,colN) = 1;
                colN++; 
            } 
            else
            {
                hitdex++;
            }

        }

        //Terminator2 
        colN = 0;
        hitdex = 0;
        for (size_t i = 0; i < Terminator2.rows(); ++i)
        {
            if ( i != hit_list_rows[hitdex] )
            {
                Terminator2(i,colN) = 1;
                colN++; 
            } 
            else
            {
                hitdex++;
            }

        }

	    MatrixXd Mid1 = MatBase * Terminator1;
	    MatrixXd Mid2 = (Mid1.transpose() * Terminator2).transpose();

        return Mid2;
    }

    VectorXd Row_Del_Vec( VectorXd & VecBase, size_t & minusel, vector<size_t> hit_list  )
    {
    	//Dense
    	//MatrixXd Terminator1 = MatrixXd::Zero ( VecBase.size(), VecBase.size() - minusel );

        //Sparse
        Eigen::SparseMatrix<double> Terminator1(  VecBase.size(), VecBase.size() - minusel );
        
        typedef Eigen::Triplet<double> T;
    	std::vector<T> tripletList;
        tripletList.reserve(VecBase.size()); 

        //Terminator1
        size_t colN = 0;
        size_t hitdex = 0;
        for (size_t i = 0; i < Terminator1.rows(); ++i)
        {
            if ( i != hit_list[hitdex] )
            {
                //Dense
                //Terminator1(i,colN) = 1;
                //Triplet Sparse Modification
                tripletList.push_back( T( i, colN, 1) ) ; 
                colN++; 
            } 
            else
            {
                hitdex++;
            }

        }

        Terminator1.setFromTriplets(tripletList.begin(), tripletList.end()); 
    	VectorXd Mid1 = Terminator1.transpose()*VecBase;

    	return Mid1;

    }


    MatrixXd Row_Deletion(MatrixXd & TEST_ERASE, vector<size_t> hit_list)
    {
        MatrixXd Temporal_Test1;
        MatrixXd Temporal_Test2;
        size_t posdel;
        for (size_t i =0; i < hit_list.size(); ++i)  //Because it is in inverse order
        {

           posdel = hit_list[i];
           if (int(hit_list[i]) == 0)
           {
              Temporal_Test1 = TEST_ERASE.bottomLeftCorner(TEST_ERASE.rows()-1,TEST_ERASE.cols());
              TEST_ERASE = Temporal_Test1;    
           }	
           else if ( int(hit_list[i]) == int(TEST_ERASE.size()-1) )
           {
              Temporal_Test1 = TEST_ERASE.topLeftCorner(TEST_ERASE.rows()-1,TEST_ERASE.cols());
              TEST_ERASE = Temporal_Test1;                  
           }
           else
           {
              Temporal_Test1 = TEST_ERASE.topLeftCorner(posdel,TEST_ERASE.cols());
              Temporal_Test2 = TEST_ERASE.bottomLeftCorner(TEST_ERASE.rows()-posdel-1,TEST_ERASE.cols());
              MatrixXd Temporal_Test3(Temporal_Test1.rows()+Temporal_Test2.rows(), Temporal_Test1.cols());  	
              Temporal_Test3 << Temporal_Test1 , Temporal_Test2;   //----
              TEST_ERASE = Temporal_Test3;
           }
        } 
        return TEST_ERASE;
    }    

    VectorXd Row_Deletion_Vec(VectorXd & TEST_ERASE, vector<size_t> hit_list)
    {
        VectorXd Temporal_Test1;
        VectorXd Temporal_Test2;
        size_t posdel;
        for (size_t i =0; i < hit_list.size(); ++i)  //Because it is in inverse order
        {

           posdel = hit_list[i];
           if (int(hit_list[i]) == 0)
           {
              Temporal_Test1 = TEST_ERASE.tail(TEST_ERASE.rows()-1);
              TEST_ERASE = Temporal_Test1;    
           }	
           else if ( int(hit_list[i]) == int(TEST_ERASE.size()-1) )
           {
              Temporal_Test1 = TEST_ERASE.head(TEST_ERASE.rows()-1);
              TEST_ERASE = Temporal_Test1;                  
           }
           else
           {
           	  Temporal_Test1 = TEST_ERASE.head(posdel);
           	  Temporal_Test2 = TEST_ERASE.tail(TEST_ERASE.rows()-posdel-1);
              VectorXd Temporal_Test3(Temporal_Test1.rows()+Temporal_Test2.rows());  	
              Temporal_Test3 << Temporal_Test1 , Temporal_Test2;   //----
              TEST_ERASE = Temporal_Test3;
           }
        } 
        return TEST_ERASE;
    }   

    VectorXd Row_Insertion_Vec(VectorXd & TEST_ERASE, vector<size_t> hit_list)  //Use list in ascending order of nodes eliminated
    {

        size_t posadd;
        for (size_t i =0; i < hit_list.size(); ++i)  //Because it is in inverse order
        {
           VectorXd Temporal_Test1;
           VectorXd Temporal_Test2;
           posadd = hit_list[i];
           if (int(hit_list[i]) == 0)
           {
              Temporal_Test1 << 0;
              Temporal_Test2 = TEST_ERASE;
              VectorXd Temporal_Test3(Temporal_Test1.rows()+Temporal_Test2.rows());  	
              Temporal_Test3 << Temporal_Test1 , Temporal_Test2;   //----
              TEST_ERASE = Temporal_Test3;  
           }	
           else if ( int(hit_list[i]) == int(TEST_ERASE.size()-1) )
           {
              Temporal_Test1 = TEST_ERASE;
              Temporal_Test2 << 0;
              VectorXd Temporal_Test3(Temporal_Test1.rows()+Temporal_Test2.rows());  	
              Temporal_Test3 << Temporal_Test1 , Temporal_Test2;   //----
              TEST_ERASE = Temporal_Test3;                 
           }
           else
           {
           	  Temporal_Test1 = TEST_ERASE.head(posadd);
           	  Temporal_Test2 << 0;
           	  VectorXd Temporal_Test3 = TEST_ERASE.tail(TEST_ERASE.rows()-posadd);
              VectorXd Temporal_Test4(Temporal_Test1.rows()+Temporal_Test2.rows()+Temporal_Test3.rows());  	
              Temporal_Test4 << Temporal_Test1 , Temporal_Test2, Temporal_Test3;   //----
              TEST_ERASE = Temporal_Test4;
           }
        } 
        return TEST_ERASE;
    }   


    //Quick mean for 2 numbers
    double meanf(double & Val1, double & Val2)
    {
        double out = (Val1 + Val2) * 0.5;
        return out;
    }

    //Lateral Hydrostatic Pressure
    double bcRunf(const double & rhowater, double & hy, const size_t & yw0_sl, size_t & yval)
    {
        double g = 9.81; //Internal gravity definition
        double out = 0;

        if (yval <= yw0_sl)
        {
        	out = -rhowater * g * (yw0_sl-yval) * hy;
        }
    	return out;
    }

    //BUOYANCY PRESSURE
	//TRY SAME AS LATERAL HYDROSTATIC BUT IN UP DIRECTION for idx at limit, depending on x position
	//Some blocks might have right and up pressures as melting progresses. Must provide top position of ice below sea level or is it equivalent to be at sea level for displacement
    double bcDvnf(const double & rhowater, double & hy, const size_t & yw0_sl, size_t & top, size_t & yval)
    {
        double g = 9.81; //Internal gravity definition
        double out = 0;
        double dist;

        if (yval <= yw0_sl)
        {
        	dist = min( yw0_sl - yval , top - yval );
        	out = rhowater * g * dist * hy;
        }
    	return out;
    }     	


	// Finds contact data between *this and other - for output purposes
	CData findContactData(const Grain2d & other, const double & dt, const double & offset) const {

		// Declare temp variables
		CData 			   cData;								 // output
		Vector2d 		   force;								 // total force
		Vector2d		   ptOtherCM; 		                     // point wrt the center of mass of other in real space
		Vector2d		   ptOtherLset; 	                     // point in the reference config of other's level set
		double		       penetration;	                         // penetration amount (is negative by convention of the level set)
		Vector2d		   normal; 			                     // surface normal pointing from other to *this
		const double       cos2 = cos(other.getTheta());
		const double       sin2 = sin(other.getTheta());
		Vector2d		   ptThisCM; 		                     // point wrt the center of mass of this in real space
		Vector2d		   df; 				                     // force increment from a single point
		Vector2d		   tangent; 		                     // surface tangent
		double		       sdot; 			                     // relative velocity of a point of contact
		double 			   ndot;
		double		       ds;									 // relative displacement of a point of contact
		double 			   dn;
		Vector2d		   Fs; 				                     // vector of frictional/shear force
		Vector2d           relativeVel;							 // vector of relative velocity between the two grains
        double Kn = (_kn+other.getKn())/2.;						 // average normal contact stiffness between the particles in contact
		double Ks = (_ks+other.getKs())/2.;						 // average shear contact stiffness between the particles in contact
		double Fsmag;											 // magnitude of shear force
		double cres = .43;										 // auxilliary variable for computing contact damping coeff below
		double GamaN = -2*sqrt(_kn*_mass*other.getMass()/(_mass+other.getMass()))*log(cres)/sqrt(M_PI*M_PI + log(cres)*log(cres));
		Vector2d		   forceNormal;
		double 			   sigMoment;

		Vector2d offsetVec;
		if (other.getPosition()(0) > _position(0)) {
			offsetVec << -offset, 0.;
		}
		else {
			offsetVec << offset, 0.;
		}

		// Initialize!
		force.fill(0.0);
		normal.fill(0.0);

		// iterate through all of the points of *this and check for contact for each one
		for (size_t ptidx = 0; ptidx < _pointList.size(); ptidx++) {

			ptThisCM = _pointList[ptidx] - _position; //OFF Melting, ON Breakage

			ptOtherCM = _pointList[ptidx] - other.getPosition() - offsetVec;
			ptOtherLset(0) =  ptOtherCM(0)*cos2 + ptOtherCM(1)*sin2;
			ptOtherLset(1) = -ptOtherCM(0)*sin2 + ptOtherCM(1)*cos2;
			ptOtherLset += other.getCmLset();

			//	if there is penetration, get contact info into the cData struct
			if ( other.getLset().isPenetration(ptOtherLset, penetration, normal) ) {

				normal << normal(0)*cos2 - normal(1)*sin2, normal(0)*sin2 + normal(1)*cos2;
				tangent << -normal(1), normal(0);
				relativeVel << other.getVelocity()(0) - other.getOmega()*ptOtherCM(1) - (_velocity(0) - _omega*ptThisCM(1)),
				               other.getVelocity()(1) + other.getOmega()*ptOtherCM(0) - (_velocity(1) + _omega*ptThisCM(0));
				df = penetration*normal*Kn - GamaN*normal.dot(relativeVel)*normal;
				force -= df;
				if (_nodeShears[ptidx] > 0) {
					Fsmag = std::min(_nodeShears[ptidx], df.norm()*_mu );
				}
				else {
					Fsmag = std::max(_nodeShears[ptidx], -df.norm()*_mu );
				}
				Fs = -tangent*Fsmag;
				force += Fs;

				cData._cpairs.push_back(Vector2i(_id, other.getId()));
				cData._forces.push_back(force);
				cData._normals.push_back(normal);
				cData._clocs.push_back(_pointList[ptidx]);
			}
		} // end for loop iterating through contact points
		return cData;
	} // findContactData
	
	// Finds contact data between *this and other - for output purposes regular DEM
	CData findContactDataDEM(const Grain2d & other, const double & dt, const double & offset) const {

		CData cData;				// output
		Vector2d force;			// total force
		Vector2d ptThisCM; 		// point wrt the center of mass of *this in real space
		Vector2d ptOtherCM; 		// point wrt the center of mass of other in real space
		Vector2d normal; 			// surface normal pointing out of other in the reference config (initially) and then pointing out of *this in real space (after)
		double Fsmag;
		Vector2d shear;
		Vector2d Fs, Fn;

		Vector2d dr = _position - other.getPosition();
		double drMag = dr.norm();
		double sigma = _radius + other.getRadius();
		double penetration = drMag - sigma ; //negative if contact 
		size_t ptidx = other.getId();

		Vector2d offsetVec;
		if (other.getPosition()(0) > _position(0)) {
			offsetVec << -offset, 0.;
		}
		else {
			offsetVec << offset, 0.;
		}

		// Initialize!
		force.fill(0.0);
		normal.fill(0.0);

		//	if there is penetration, get contact info into the cData struct
		if (penetration<=0.){

			// update force: normal force contribution
			// note: penetration is negative which negates the negative sign in eq (2) of KWL
			normal = -dr/drMag;
			ptThisCM = _radius * normal;
			ptOtherCM = ptThisCM + _position - other.getPosition() - offsetVec;
			Fn = penetration*normal*_kn;
			force += Fn;

            //Shear force correction
			if (_nodeNormals.find(ptidx) != _nodeNormals.end()){
			    std::map<size_t, Vector2d> tempShearv2 = _nodeShearsv2;
				shear = tempShearv2[ptidx];
				Fsmag = min(Fn.norm()*_mu, shear.norm() ); // eq (14)
				if(Fsmag > 0.0){
					Fs = Fsmag*shear/shear.norm(); // eq (13)
					force += Fs;
				}
			}

			cData._cpairs.push_back(Vector2i(_id, other.getId()));
			cData._forces.push_back(force);
			cData._normals.push_back(normal);
			cData._clocs.push_back(ptThisCM + _position);
		}
		return cData;
	} // findContactDataDEM regular DEM


	// Finds contact data between *this and other - for output purposes
	CData findContactDataIDX(const Grain2d & other, const double & dt, const double & offset, vector<size_t> & idx) const {

		// Declare temp variables
		CData 			   cData;								 // output
		Vector2d 		   force;								 // total force
		Vector2d		   ptOtherCM; 		                     // point wrt the center of mass of other in real space
		Vector2d		   ptOtherLset; 	                     // point in the reference config of other's level set
		double		       penetration;	                         // penetration amount (is negative by convention of the level set)
		Vector2d		   normal; 			                     // surface normal pointing from other to *this
		const double       cos2 = cos(other.getTheta());
		const double       sin2 = sin(other.getTheta());
		Vector2d		   ptThisCM; 		                     // point wrt the center of mass of this in real space
		Vector2d		   df; 				                     // force increment from a single point
		Vector2d		   tangent; 		                     // surface tangent
		double		       sdot; 			                     // relative velocity of a point of contact
		double 			   ndot;
		double		       ds;									 // relative displacement of a point of contact
		double 			   dn;
		Vector2d		   Fs; 				                     // vector of frictional/shear force
		Vector2d           relativeVel;							 // vector of relative velocity between the two grains
        double Kn = (_kn+other.getKn())/2.;						 // average normal contact stiffness between the particles in contact
		double Ks = (_ks+other.getKs())/2.;						 // average shear contact stiffness between the particles in contact
		double Fsmag;											 // magnitude of shear force
		double cres = .43;										 // auxilliary variable for computing contact damping coeff below
		double GamaN = -2*sqrt(_kn*_mass*other.getMass()/(_mass+other.getMass()))*log(cres)/sqrt(M_PI*M_PI + log(cres)*log(cres));
		Vector2d		   forceNormal;
		double 			   sigMoment;

		Vector2d offsetVec;
		if (other.getPosition()(0) > _position(0)) {
			offsetVec << -offset, 0.;
		}
		else {
			offsetVec << offset, 0.;
		}

		// Initialize!
		force.fill(0.0);
		normal.fill(0.0);

		// iterate through all of the points of *this and check for contact for each one
		for (size_t ptidx = 0; ptidx < _pointList.size(); ptidx++) {

			ptThisCM = _pointList[ptidx] - _position;
			ptOtherCM = _pointList[ptidx] - other.getPosition() - offsetVec;
			ptOtherLset(0) =  ptOtherCM(0)*cos2 + ptOtherCM(1)*sin2;
			ptOtherLset(1) = -ptOtherCM(0)*sin2 + ptOtherCM(1)*cos2;
			ptOtherLset += other.getCmLset();

			//	if there is penetration, get contact info into the cData struct
			if ( other.getLset().isPenetration(ptOtherLset, penetration, normal) ) {

				normal << normal(0)*cos2 - normal(1)*sin2, normal(0)*sin2 + normal(1)*cos2;
				tangent << -normal(1), normal(0);
				relativeVel << other.getVelocity()(0) - other.getOmega()*ptOtherCM(1) - (_velocity(0) - _omega*ptThisCM(1)),
				               other.getVelocity()(1) + other.getOmega()*ptOtherCM(0) - (_velocity(1) + _omega*ptThisCM(0));
				df = penetration*normal*Kn - GamaN*normal.dot(relativeVel)*normal;
				force -= df;
				if (_nodeShears[ptidx] > 0) {
					Fsmag = std::min(_nodeShears[ptidx], df.norm()*_mu );
				}
				else {
					Fsmag = std::max(_nodeShears[ptidx], -df.norm()*_mu );
				}
				Fs = -tangent*Fsmag;
				force += Fs;

				cData._cpairs.push_back(Vector2i(_id, other.getId()));
				cData._forces.push_back(force);
				cData._normals.push_back(normal);
				cData._clocs.push_back(_pointList[ptidx]);
				idx.push_back(ptidx); //Vector of contacts points in pointList for this grain
			}
		} // end for loop iterating through contact points
		return cData;
	} // findContactData


	// Compute kinetic energy of the grain
	double computeKineticEnergy() const {
		double ke = 0;
		ke += .5*_mass*_velocity.squaredNorm();
		ke += .5*_momentInertia*_omega*_omega;
		return ke;
	}

	// Explicit time integration
	void takeTimestep(const Vector2d & force, const double & moment, const double & gDamping, const double & dt) {
	    //Apply this reset for a time step after breakage and then ignore?
	    //Vector2d forcea = Vector2d(0.,0.);
        //double momenta = 0.0;
		// cout << "Force update " << force.transpose() <<endl;
		// cout << "Moment update " << moment <<endl;
		//_velocity = 1/(1+gDamping*dt/2)*( (1-gDamping*dt/2)*_velocity + dt*forcea/_mass   );
		//_omega = 1/(1+gDamping*dt/2)*( (1-gDamping*dt/2)*_omega + dt*momenta/_momentInertia);
		
        //cout << "TT Force Drag: " << endl;
        //cout << "X: " << force(0) << " Y: " << force(1) << endl;
        //cout << "TT Moment Drag: " << moment << endl;

		_velocity = 1/(1+gDamping*dt/2)*( (1-gDamping*dt/2)*_velocity + dt*force/_mass   );
		_omega = 1/(1+gDamping*dt/2)*( (1-gDamping*dt/2)*_omega + dt*moment/_momentInertia);
		double cosd = cos(_omega*dt);
		double sind = sin(_omega*dt);
		
		//cout << "TT Velocity: " << endl;
        //cout << "X: " << _velocity(0) << " Y: " << _velocity(1) << endl;
        //cout << "TT Omega: " << _omega << endl;
		
		
		// must update the points
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] << (_pointList[ptid](0)-_position(0))*cosd - (_pointList[ptid](1)-_position(1))*sind,
								(_pointList[ptid](0)-_position(0))*sind + (_pointList[ptid](1)-_position(1))*cosd;
			_pointList[ptid] += _position + _velocity*dt;
		}
		_position = _position + dt*_velocity;
		_theta = _theta + dt*_omega;
		
		//cout << "TT Position: " << endl;
        //cout << "X: " << _position(0) << " Y: " << _position(1) << endl;
        //cout << "TT Theta: " << _theta << endl;
	}
	
	//Control Velocity from Data
	void takeTimestepVelocity(double & speed, double & angle, const double & moment, const double & gDamping, const double & dt) {
        double dts = 1.0; //Actually it is 1 second
        double rad_angle = angle*3.141592653589793/180.0; //Convert from degrees to radians
		_velocity(0) = speed * cos(rad_angle);
		_velocity(1) = speed * sin(rad_angle);
		
		_omega = 1/(1+gDamping*dt/2)*( (1-gDamping*dt/2)*_omega + dt*moment/_momentInertia);
		double cosd = cos(_omega*dt);
		double sind = sin(_omega*dt);
		
		// must update the points
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
		    //Rotation
			_pointList[ptid] << (_pointList[ptid](0)-_position(0))*cosd - (_pointList[ptid](1)-_position(1))*sind,
								(_pointList[ptid](0)-_position(0))*sind + (_pointList[ptid](1)-_position(1))*cosd;
			//Displacement
			_pointList[ptid] += _position + _velocity*dts;
		}
		_position = _position + dts*_velocity;
		_theta = _theta + dt*_omega;

	}
	
	
	//For load grains
	void takeTimestepLoad(const Vector2d & force, const double & moment, const double & gDamping, const double & dt, Vector2d & tryvel) {
	    //Apply this reset for a time step after breakage and then ignore?
	    //Vector2d forcea = Vector2d(0.,0.);
        //double momenta = 0.0;
		// cout << "Force update " << force.transpose() <<endl;
		// cout << "Moment update " << moment <<endl;
		//_velocity = 1/(1+gDamping*dt/2)*( (1-gDamping*dt/2)*_velocity + dt*forcea/_mass   );
		//_omega = 1/(1+gDamping*dt/2)*( (1-gDamping*dt/2)*_omega + dt*momenta/_momentInertia);
		
        //cout << "TT Force Drag: " << endl;
        //cout << "X: " << force(0) << " Y: " << force(1) << endl;
        //cout << "TT Moment Drag: " << moment << endl;

		//_velocity = 1/(1+gDamping*dt/2)*( (1-gDamping*dt/2)*_velocity  );
		_velocity = tryvel; //FIX NOW
		_omega = 1/(1+gDamping*dt/2)*( (1-gDamping*dt/2)*_omega   );
		double cosd = cos(_omega*dt);
		double sind = sin(_omega*dt);
		
		//cout << "TT Velocity: " << endl;
        //cout << "X: " << _velocity(0) << " Y: " << _velocity(1) << endl;
        //cout << "TT Omega: " << _omega << endl;
		
		
		// must update the points
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] << (_pointList[ptid](0)-_position(0))*cosd - (_pointList[ptid](1)-_position(1))*sind,
								(_pointList[ptid](0)-_position(0))*sind + (_pointList[ptid](1)-_position(1))*cosd;
			_pointList[ptid] += _position + _velocity*dt;
		}
		_position = _position + dt*_velocity;
		_theta = _theta + dt*_omega;
		
		//cout << "TT Position: " << endl;
        //cout << "X: " << _position(0) << " Y: " << _position(1) << endl;
        //cout << "TT Theta: " << _theta << endl;
	}
	

	void moveGrain(const Vector2d & amount) {
		_position = _position + amount;
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] += amount;
		}
	}

	//Point Generating Function using a Level Set
    void PointsGen(const Levelset2d & glset, vector<Vector2d> & gpts, const int & npts, const Vector2d & newcm, bool & fail_ind){  //const int & npts
    	bool fail_indi = false; //Assume it's okay
    	//cout <<"Point Gen" << endl;
    	vector<double> gvec = glset.getLevelset();		
        double damping= 0.3;
        //size_t pointsmeltdense = 400;

        //_npoints = pointsmeltdense;
		vector<Vector2d> pointsnew(npts); //npts
		
		double MeltTemp = 0;

        cout<<"Zero Cross"<<endl;
    	//  place points interpolating around 0 contour or meltTemp
		vector<Vector2d> zeroCrossings;
		for (size_t j = 0; j < glset.getYdim()-1; j++) {
			for (size_t i = 0; i < glset.getXdim()-1; i++) {
				double val = gvec[(j*glset.getXdim())+i];
				if (val*gvec[((j+1)*glset.getXdim())+i] < (MeltTemp)   ) { 
					zeroCrossings.push_back(Vector2d( double(i), double(j)+fabs(val)/(fabs(val) + fabs(gvec[((j+1)*glset.getXdim())+i]))  ));
				}
				if (val*gvec[(j*glset.getXdim())+(i+1)] < (MeltTemp)   ) {	
					zeroCrossings.push_back(Vector2d(  double(i)+fabs(val)/(fabs(val) + fabs(gvec[(j*glset.getXdim())+(i+1)])), double(j)  ));
				}
			}
		}	

	    double step = double(zeroCrossings.size())/double(npts);

        // cout<<"Temporal pointsnew"<<endl;
        cout<<"Numberpoints: "<< npts <<endl;
        cout<<"ZeroC size: "<< zeroCrossings.size() <<endl;
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
				for (size_t i = 0; i < glset.getXdim(); i++) {   //Substitute of the sDirac Function to Identify Insidre and Outside of Grain
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
					//pointsnew[i] << 5+(double)rand()/RAND_MAX,5+(double)rand()/RAND_MAX; //or << //PSOE I (ok)
				}
				forces[i] << 0., 0.; // or <<
			}
	// 		// iterate through points
			for (size_t i = 0; i < npts; i++) {
 			//drag points closer to surface using gradient/signed distance
	        	//pointsnew[i] -= .8*sdlist[i]*glist[i];                                          //PSOE II
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
		        //pointsnew[i] += ddt*v[i];                                                 //PSOE III
			}
		}
		
		// additional iterations to move points to zero contour
		for (int tt = 0; tt < 200; tt++) {   //200
			for (size_t i = 0; i < npts; i++) {
				if (glset.findValueGradient(pointsnew[i], sd, grad)) {                   //change _lset for Geometric Level Set, use _grainTemp2D for Heat Level Set -->   //Does it work for Binary Lset mmmmm unlikely
					//pointsnew[i] -= .8*sd*grad;                                        //PSOE IV (not ok)
				}
			}
		}

		//Check if pointsnew were done ok
		bool nanp_before = false;
		cout << "Check Nan Points in Point Gen!!!!" << endl;
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
			cout << "CONFIRMED: Nan Points in Point Gen!!!!" << endl;
			return;  //Skip to next healthy grain
		}
        
        //REORDER POINTS CLOCKWISELY

        //FIND CLOCKWISE ANGLE FROM CENTROID TO EACH POINT, THEN ORDER THE POINTS FROM LOWEST TO HIGHEST ANGLE (ALTERNATIVE 2)
        
        cout<<"Sorting CW"<<endl;

        vector<Vector2d> pneworder=pointsnew;

        vector<Vector2d> pangle=pointsnew;
        vector<Vector2d> pangle2=pointsnew;

        vector<int> indicesorder(npts);
        size_t minindex=0;
        const double pii=atan(1)*4;
        Vector2d savior = Vector2d (0.,0.);

        for (size_t i = 0; i < npts; i++) {       //Loop for Setting Up Angles
        	
             pangle[i](0) = (atan (  (  pneworder[i](1)-newcm(1) ) /  (  pneworder[i](0)-newcm(0) ) ) )*(180/pii); //Get CCW angle in degrees
             pangle[i](1) = i; //Point index storage for convenience
             
             //Get CW angle in degrees and in the right quadrant measured as CW from theta=0
             if (   (pneworder[i](0)-newcm(0)) >= 0 && (pneworder[i](1)-newcm(1)) <= 0  ) {
             	   pangle[i](0)=-pangle[i](0); 
             } else if (   (pneworder[i](0)-newcm(0)) < 0 && (pneworder[i](1)-newcm(1)) < 0  ) {
             	   pangle[i](0)=180-pangle[i](0);

             }  else if (   (pneworder[i](0)-newcm(0)) < 0 && (pneworder[i](1)-newcm(1)) > 0  ) {
             	   pangle[i](0)=180-pangle[i](0);

             }  else { //if (   (pneworder[i](0)-_cmLset(0)) > 0 && (pneworder[i](1)-_cmLset(1)) > 0  ) {
                   pangle[i](0)=360-pangle[i](0);
             } //else {
                 //cout << "grain " << _id << " has a problem at its newpointagle" << endl; 
             //}
             pangle2[i](0)=pangle[i](0);
             pangle2[i](1)=pangle[i](1);
        }
        
 		cout<<"Sorting CW Part 2"<<endl;
        //Sorting Loop CW
        double closeness=0.0;
        for (size_t i = 0; i < npts; i++) {          //sort(v.begin(), v.end()); 
        	minindex=0;
        	if (i == npts-1){      //Check Closing Process for Points
                pangle2[i] = pangle[0];   //Should be only pangle2[i] = pangle[1]; //This is only to get the last remaining point
        	}
        	//MODIFICATION 2: 1 or 2 seem to be nice.   Try using 5 at 90 degress to fully have CW
        	else if (i<2){   // We can map (2)3 points clockwise to set up right clockwise direction and then go with nearest distance instead of angle 0 , 3  or 7
            //else if (i < npts-1) {
                for (size_t j = 0; j < pangle.size() ; j++) {	
                	if ( pangle[minindex](0)>pangle[j](0) ) {    //Improve criterion
                   		minindex=j;
                	}
                }	
               savior=pangle[minindex];
               pangle.erase(pangle.begin() + minindex);
               pangle2[i]= savior;
        	}
            else{   //This is to sort all the other points using nearest distance to next instead of angle ordering, given regular spacing this should be fine 
            	 
                //Should a competition for nearest distance be established for improving quality?????? Given preference to closer on Y???


            	closeness=(pneworder[pangle2[i-1](1)]-pneworder[pangle[minindex](1)]).norm();
            	for (size_t j = 1; j < pangle.size() ; j++)   //pangle.size()
            	{     
            		   

            		//MODIFICATION 3: Use CM for your convenience                                              
                    if (   ((pneworder[pangle2[i-1](1)]-pneworder[pangle[j](1)]).norm())<closeness )
                    { 
                    //if (   ((pneworder[pangle2[i-1](1)]-pneworder[pangle[j](1)]).norm())<closeness  &&   pneworder[pangle[j](1)](0) > _cmLset(0)   &&    pneworder[pangle[j](1)](1) > _cmLset(1)  ){                           
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
        for (size_t i = 0; i < npts; i++) {
        	
        	minindex=pangle2[i](1);
        	
        	//minindex=indicesorder[i];
        	//cout << "pangle " << pangle2[i](0) << " grain no." << _id << endl;

        	gpts[i]=pointsnew[minindex];

        }	

    }


	// Change methods
	void changeMu(const double & newmu) {
		_mu = newmu;
	}

	void changePos(const Vector2d & pos) {
		Vector2d disp = pos - _position;
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] += disp;
		}
		_position = pos;
	}
	
	//Use only to initialize initial positions for measuring displacements
	void definePos0(const Vector2d & new_pos) {
		_position0 = new_pos;
	}
	

	void changeRot(const double & rot) {
		double dtheta = rot - _theta;
		_theta = rot;
		double cosd = cos(dtheta);
		double sind = sin(dtheta);
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] << (_pointList[ptid](0)-_position(0))*cosd - (_pointList[ptid](1)-_position(1))*sind,
									  (_pointList[ptid](0)-_position(0))*sind + (_pointList[ptid](1)-_position(1))*cosd;
			_pointList[ptid] += _position;
		}
	}

	void changeKn(const double & kn) {
		_kn = kn;
	}

	void changeKs(const double & ks) {
		_ks = ks;
	}

	void changeDensity(const double & density) {
		_mass *= density/_density;
		_momentInertia *= density/_density;
		_density = density;
	}

	void changeId(const size_t & id) {
		_id = id;
	}

	void scaleRadius(const double & scale){
		_radius += scale;
	}

	// erase friction history
	void clearFriction() {
		for (size_t i = 0; i < _nodeShears.size(); i++) {
			_nodeShears[i] = 0.;
			// _nodeNormals[i] << 0,0,0;
			_nodeContact[i] = 0;
		}
	}

	// For wall grain purposes
	void changeVelocity(const Vector2d & velocity) {
		_velocity = velocity;
	}
	void changeOmega(const double & omega) {
		_omega = omega;
	}
	void changeTheta(const double & theta) {
		_theta = theta;
	}	
	void changeId(const unsigned int newId) {
		_id = newId;
	}


	void changeMass (const double & mnew){ 
	 	_mass= mnew;
	}

	void changeRadius (const double & rnew){ 
	 	_radius= rnew;
	}

	void changeTemper (const double & tempernew){ 
	 	_temper= tempernew;
	}

	void changeThickness (const double & thicknew){ 
	 	_thick= thicknew;
	}
	
	void changeThickness0 (const double & thicknew){ 
	 	_thick0 = thicknew;
	}




    void changeFracFlag(const bool & fracflag) {
		_fracFlag = fracflag;
	}
	void changeFracLoc(const int & fracloc) {
		_fracLoc = fracloc;
	}

	void change_Mthick(vector<double> & lsvector, size_t & xdim, size_t & ydim)
	{
         Levelset2d lsetnew(lsvector, xdim, ydim);
         _Mthick = lsetnew;
    }   

	void change_Full_lset(vector<double> & lsvector, size_t & xdim, size_t & ydim)
	{
         Levelset2d lsetnew(lsvector, xdim, ydim);
         _lset = lsetnew;
    }   

    void change_Full_lset0(vector<double> & lsvector, size_t & xdim, size_t & ydim)
	{
         Levelset2d lsetnew(lsvector, xdim, ydim);
         _lset0 = lsetnew;
    }    

    void change_Full_Temp_lset(vector<double> & lsvector, size_t & xdim, size_t & ydim)
	{
         Levelset2d lsetnew(lsvector, xdim, ydim);
         _grainTemp2D = lsetnew;
    }   

     void changeInertia(double & newinertia)
	{
        _momentInertia = newinertia;
    }     

    void changeCM(Vector2d & newCM)
	{
        _cmLset = newCM;
    }     

    void changePoints(vector<Vector2d> & newpoints)
	{
        _pointList = newpoints;
    }    

    void changeNPoints (int & newnumber)
    {
    	//_npoints = newnumber;  //DEACTIVATED
    }

    void changeDamage (vector<Vector3d> & newDamage)
    {
    	_grainDamage = newDamage;
    }   

    void changeInternalStress (vector<double> & lsvector, size_t & xdim, size_t & ydim)
    {
    	Levelset2d lsetnew(lsvector, xdim, ydim);
    	_grainStress2D = lsetnew;
    }

    void changegrainThickness (MatrixXd & newThickness)
    {
    	_grainThickness = newThickness;
    }

    void changeYield(const double & yield){
		_yield = yield;
	}
	void changeCritStress(const double & VM){
		_critStress = VM;
	}
	
	void changeTfail(const double & new_time_fail){
		_time_fail = new_time_fail;
	}
	void changeTorigin(const double & new_origin_time_fail){
		_origin_time_fail = new_origin_time_fail;
	}



	// Helper methods
	const double & getMass() const {
		return _mass;
        //return _mass; //*(_thick/_thick0);
	}

	const double & getMoI() const {
		return _momentInertia;		
	}

	const Vector2d & getPosition() const {
		return _position;
	}
	const Vector2d & getVelocity() const {
		return _velocity;
	}	
	const Vector2d & getPosition0() const {
		return _position0;
	}


	const double & getTemperature() const {
		return _temper;
	}
	const double & getThickness() const {
		return _thick;
	}
    const double & getThickness0() const {
		return _thick0;
	}

	const Levelset2d & getgrainTemp2D() const {
		return _grainTemp2D;
	}
	const vector<double> & getUtemper() const {
		return _Utemper;
	}

	const Levelset2d & getgrainStress2D() const {
		return _grainStress2D;
	}

	const vector<Vector3d> & getgrainDamage() const {
		return _grainDamage;
	}


	const MatrixXd & getgrainThickness() const {
		return _grainThickness;
	}	

	const Levelset2d & getMthick() const {
		return _Mthick;
	}


	const double & getDensity() const {
		return _density;
	}


	const double & getTheta() const {
		return _theta;
	}
	const double & getOmega() const {
		return _omega;
	}
	const Vector2d & getCmLset() const {
		return _cmLset;
	}
	const double & getRadius() const {
		return _radius;
	}
	const Levelset2d & getLset() const {
		return _lset;
	}
	const Levelset2d & getLset0() const {
		return _lset0;
	}
	const double & getKn() const {
		return _kn;
	}
	const double & getKs() const {
		return _ks;
	}
	const double & getMu() const {
		return _mu;
	}

    const bool & getFracFlag() const{
		return _fracFlag;
	}
	const int & getFracLoc() const{
		return _fracLoc;
	}

	const double & getYield() const{
		return _yield;
	}
	const double & getCritStress() const {
		return _critStress;
	}
	const bool & getRemove() const {
		return _remove;
	}
	const double & getDiameter() const {
		return _d;
	}


	const size_t & getId() const {
		return _id;
	}
    const size_t & getmorphologyID()const {
        return _morphologyID;
    }
    const vector<Vector2d> getPointList() const {
    	return _pointList;
    }
    const vector<Vector2d> & getRefPointList() const {
		return _refpointList;
	}
    const int & getnpoints() const {
    	return _npoints;
    }
	const vector<size_t> & getNodeContact() const {
		return _nodeContact;
	}
	const vector<double> & getNodeShears() const {
		return _nodeShears;
	}


	// const vector<Vector2d> & getNodeNormals() const {
	// 	return _nodeNormals;
	// }

	// non const
	vector<double> & getNodeShearsNonConst() {
		return _nodeShears;
	}
	vector<size_t> & getNodeContactNonConst() {
		return _nodeContact;
	}
	
	const double & getTfail()
	{
	    return _time_fail;
	}
	const double & getTorigin()
	{
	    return _origin_time_fail;
	}
	
	//Bonded Particle Method Properties (May 7, 2023)
    // bond functions //FULL LS-DEM
	Vector3i createInterparticleBond(Grain2d & other) { // 0/1, this id , other Id
		// loop through all the grid point, find each cohesive distance, find the smallest one, determine if creating a bond
		double value = 0.;
		Vector2d normal; normal.fill(0);
		Vector3i bondInfoIJ; bondInfoIJ.fill(0);

        Vector2d ptOtherCM; // point wrt the center of mass of other in real space
		Vector2d ptOtherLset;
		const double cos2 = cos(other.getTheta());
		const double sin2 = sin(other.getTheta());

		for (size_t ptindex = 0; ptindex < _pointList.size(); ptindex++) {
			ptOtherCM = _pointList[ptindex] - other.getPosition();
			ptOtherLset(0) =  ptOtherCM(0)*cos2 + ptOtherCM(1)*sin2;
			ptOtherLset(1) = -ptOtherCM(0)*sin2 + ptOtherCM(1)*cos2;
			ptOtherLset += other.getCmLset();
			if (other.getLset().isCohesion(ptOtherLset, _cohesiveDistance, value, normal)) {
				// bond created
				normal = -Vector2d(normal(0)*cos2 - normal(1)*sin2, normal(0)*sin2 + normal(1)*cos2);
				bondInfoIJ(0) = 1;
				bondInfoIJ(1) = _id;
				bondInfoIJ(2) = other.getId();

				_bondpts.push_back(ptindex);
				_bondInfo.push_back(bondInfoIJ);
				_bondForceShear.push_back(Vector2d::Zero());
				_bondForceNormal.push_back(Vector2d::Zero());
				_bondThisMomentNormal.push_back(0.);
				_fSig.push_back(0.);
				_fTau.push_back(0.);
				_sigrF.push_back(1.0);
				_taurF.push_back(1.0);
				_bondNormals.push_back(normal);
			}
		}

		return bondInfoIJ;
	}
	
    // bond functions //Regular DEM //For centroid bonds
	Vector3i createInterparticleBondDEM(Grain2d & other) { // 0/1, this id , other Id
		// loop through all the grid point, find each cohesive distance, find the smallest one, determine if creating a bond
		Vector2d normal; normal.fill(0);
		Vector3i bondInfoIJ; bondInfoIJ.fill(0);
		size_t bondpt1, bondpt2;
		double distmag;
		Vector2d dist;
		
		distmag = (_position - other.getPosition()).norm();

		if (distmag-_radius - other.getRadius() <= _cohesiveDistance) {
			// bond created
			dist = other.getPosition() - _position;
			normal = dist / distmag; // point outside of this into other
			bondInfoIJ(0) = 1;
			bondInfoIJ(1) = _id;
			bondInfoIJ(2) = other.getId();
			
			bondpt1 = _refpointList.size();
			_refpointList.push_back(Vector2d(_radius*normal));
			bondpt2 = other.getRefPointList().size();
			other.getRefPointListNonConst().push_back(Vector2d(other.getRadius()*-normal));
			
			//Also update pointList for plotting purposes
			Vector2d pref = Vector2d(_radius*normal);
			double cosd = cos(_theta);
    		double sind = sin(_theta);
    		Vector2d plist;
    		plist(0) = pref(0)*cosd - pref(1)*sind;
    		plist(1) = pref(0)*sind + pref(1)*cosd;
    		plist += _position;
    		_pointList.push_back(plist);

			_bondpts.push_back(bondpt1);
			_bondOtherpts.push_back(bondpt2);
			_bondInfo.push_back(bondInfoIJ);
			_bondForceShear.push_back(Vector2d::Zero());
			_bondForceNormal.push_back(Vector2d::Zero());
			_bondThisMomentNormal.push_back(0.);
			_bondThisMomentShear.push_back(0.);
			_fSig.push_back(0.);
			_fTau.push_back(0.);
			_sigrF.push_back(1.0);
			_taurF.push_back(1.0);
			_bondNormals.push_back(normal);
		}

		return bondInfoIJ;
	}
	
	// bond functions //Regular DEM //For multiple bonds on surface
	Vector3i createInterparticleBondDEM2(Grain2d & other) { // 0/1, this id , other Id
		// loop through all the grid point, find each cohesive distance, find the smallest one, determine if creating a bond
		Vector2d normal; normal.fill(0);
		Vector3i bondInfoIJ; bondInfoIJ.fill(0);
		double distmag;
		Vector2d dist;
		double sigma;
		
        // For all points instead of simple bond
		for (size_t ptindex = 0; ptindex < _pointList.size(); ptindex++) {
			dist = _pointList[ptindex] - other.getPosition();
			distmag = dist.norm();
			normal = dist / distmag;
			sigma = distmag - other.getRadius();

			if (sigma <= _cohesiveDistance) {
				// bond created
				bondInfoIJ(0) = 1;
				bondInfoIJ(1) = _id;
				bondInfoIJ(2) = other.getId();

				_bondpts.push_back(ptindex);
				_bondInfo.push_back(bondInfoIJ);
				_bondForceShear.push_back(Vector2d::Zero());
				_bondForceNormal.push_back(Vector2d::Zero());
				_bondThisMomentNormal.push_back(0.);
				_bondThisMomentShear.push_back(0.);
				_fSig.push_back(0.);
				_fTau.push_back(0.);
				_sigrF.push_back(1.0);
				_taurF.push_back(1.0);
				_bondNormals.push_back(normal);
			}
		}
		
		return bondInfoIJ;
	}
	
	// bond functions //Regular DEM //For ONE bond on surface points
	Vector3i createInterparticleBondDEM3(Grain2d & other) { // 0/1, this id , other Id
		// loop through all the grid point, find each cohesive distance, find the smallest one, determine if creating a bond
		Vector2d normal; normal.fill(0);
		Vector3i bondInfoIJ; bondInfoIJ.fill(0);
		double distmag = 100000000;   //Distance to minimize
		Vector2d dist;
		double sigma; 
		Vector2d dist_choice; dist_choice.fill(0);
		size_t minf = 0;
		size_t minl = 0;
		
        // Loop through all leader points to find closest to follower
		for (size_t ptindex = 0; ptindex < _pointList.size(); ptindex++) {
		    // Loop through all follower points to find the nearest point possible (only 1)
			for (size_t ptindexf = 0; ptindexf < other.getPointList().size(); ptindexf++) {
			
    			dist = _pointList[ptindex] - other.getPointList()[ptindexf];
    			
    			if ( dist.norm() < distmag && dist.norm() > 0  ){  //Condition to avoid same point!!!
    			    distmag = dist.norm();
    			    minl = ptindex;
    			    minf = ptindexf;
    			    dist_choice =  dist;
    			}
			      
			}
		}
		
		//Control
// 		cout << "Min distance: " << distmag << endl;
// 		cout << "Min distance v: " << dist_choice(0) << " " << dist_choice(1) << endl;
// 		cout << "Min leader index: " << minl << endl;
// 		cout << "Min follower index: " << minf << endl;
		normal = dist_choice / distmag;
	    sigma = distmag;
		
		//Assume no or small overlap
		if (sigma <= _cohesiveDistance) {
				// bond created
				bondInfoIJ(0) = 1;
				bondInfoIJ(1) = _id;
				bondInfoIJ(2) = other.getId();

				_bondpts.push_back(minl);
				_bondptsf.push_back(minf);
				_bondInfo.push_back(bondInfoIJ);
				_bondForceShear.push_back(Vector2d::Zero());
				_bondForceNormal.push_back(Vector2d::Zero());
				_bondThisMomentNormal.push_back(0.);
				_bondThisMomentShear.push_back(0.);
				_fSig.push_back(0.);
				_fTau.push_back(0.);
				_sigrF.push_back(1.0);
				_taurF.push_back(1.0);
				_bondNormals.push_back(normal);
		}
		if (sigma < 0)
		{
		    cout << "WARNING: BONDS FORMED WITH PENETRATING POINTS" << endl;
		}
		
		return bondInfoIJ;
	}	
	
	//Weaken or delete bonds
	void deleteInterparticleBondDEM(size_t & bondidx) {   // 0/1, this id , other Id
	    //Turn back off
		_bondInfo[bondidx][0] = 0;
		_bondForceNormal[bondidx] = Vector2d(0,0);
		_bondForceShear[bondidx] = Vector2d(0,0);
		_bondThisMomentNormal[bondidx] = 0;
		_bondThisMomentShear[bondidx] = 0;
		_fSig[bondidx] = 0;
		_fTau[bondidx] = 0;
		_bondNormals[bondidx] = Vector2d(0,0);
		_sigrF[bondidx] = 0;
		_taurF[bondidx] = 0;
		
		return;
	}	
	
    //Bond get functions
	const vector<Vector3i> & getBondInfo() const {
		return _bondInfo;
	}
	
    void changeBondInfo0(vector<Vector3i> & bondInfo0) {
		_bondInfo = bondInfo0;
	}
	
	void changeBondInfo(vector<Vector3i> & bondInfo0) {
		_bondInfo = bondInfo0;
	}
	
	vector<size_t> & getbondpts() {
		return _bondpts;
	}
	
	vector<size_t> & getbondptsf() {
		return _bondptsf;
	}
	
    vector<Vector2d> & getbondForceNormal() {
		return _bondForceNormal;
	}
	
	vector<Vector2d> & getbondForceShear() {
		return _bondForceShear;
	}
	
	vector<double> & getbondThisMomentNormal() {
		return _bondThisMomentNormal;
	}
	
	vector<double> & getbondThisMomentShear() {
		return _bondThisMomentNormal;
	}
	
	vector<Vector2d> & getbondNormals() {
		return _bondNormals;
	}
	
	vector<double> & getfSig() {
		return _fSig;
	}
	
	vector<double> & getfTau() {
		return _fTau;
	}
	
	double & getCohesiveDistance() {
		return _cohesiveDistance;
	}
	const bool & getLoadGrain() {
		return _load_grain;
	}
	//Get to change bond reduction factor SIG
	vector<double> & getsigrF() {   //Locate specific bond with id and change value
       return _sigrF;
	}	
	
	//Get to change bond reduction factor TAU
	vector<double> & gettaurF() {   //Locate specific bond with id and change value
       return _taurF;
	}	

	
	//Bond change functions
	void changeCohesiveDistance(const double & cohesivedistance) {
		_cohesiveDistance = cohesivedistance;
	}
	void changeKnBond(const double & knbond) {
		_kn_bond = knbond;
	}
	void changeKsBond(const double & ksbond) {
		_ks_bond = ksbond;
	}
	void changeBondRadius(const double & bondradius) {
		_bondRadius = bondradius;
		_bondAreas = _bondRadius * 2.; // John's paper, t = 1
		_bondMomentsInertia = 2./3.*pow(_bondRadius,3);
	}
	void changeSigC(const double & sigc) {
		_sigC = sigc;
	}
	void changeTauC(const double & tauc) {
		_tauC = tauc;
	}
	void changeLoadGrain(const bool & new_load_grain) {
		_load_grain = new_load_grain;
	}
	
	//Nearest neighbor for regular DEM
	void changeNNList(const vector<size_t> & nnList){
		_nnList = nnList;
	}
	const vector<size_t> & getNNList() const{
		return _nnList;
	}
	const vector<size_t> & getBondOtherPts() const {
		return _bondOtherpts;
	}
	vector<size_t> & getBondOtherptsNonConst() {
		return _bondOtherpts;
	}
	vector<Vector2d> & getRefPointListNonConst() {
		return _refpointList;
	}
	const std::map<size_t, Vector2d>  & getNodeNormals() const {
		return _nodeNormals;
	}
	std::map<size_t, Vector2d>  & getNodeNormalsNonConst()  {
		return _nodeNormals;
	}
	const std::map<size_t, Vector2d>  & getNodeShearsv2() const {
		return _nodeShearsv2;
	}
	std::map<size_t, Vector2d>  & getNodeShearsv2NonConst()  {
		return _nodeShearsv2;
	}
	//Change bond reduction factor SIG
	void changesigrF(size_t & bondidx, double & new_val) {   //Locate specific bond with id and change value
		_sigrF[bondidx] = new_val;
		return;
	}	
	
	//Change bond reduction factor TAU
	void changetaurF(size_t & bondidx, double & new_val) {   //Locate specific bond with id and change value
		_taurF[bondidx] = new_val;
		return;
	}	

private:

	double 			 _mass;				// particle mass
	Vector2d 		 _position; 			// location of center of mass in real space
	Vector2d 		 _velocity;			// velocity of center of mass

    double           _mass0;             //Initial Mass
	double           _temper;           //temperture of grain
	double           _thick;         //thickness of grain
	double           _temper0;           //Init. temperture of grain
	double           _thick0;         //Init. thickness of grain
	
	Levelset2d 		 _Mthick;				//Thickness distribution in grain

	vector<double>   _Utemper;  //Temperature vector
    vector<double>   _Utemper0;  //Init. Temperature vector
    Levelset2d 		 _grainTemp2D;				// Temperature of level set of grain
    Levelset2d 		 _grainTemp2D0;				// Initial Temperature of level set of grain

    Levelset2d 		 _grainStress2D;				// Relevant 2D Stress fo Failure or Deformation
    vector<Vector3d> _grainDamage; 
    //MatrixXd         _grainDamage;                  // Damage Value from 0 to 1 to track deterioration of ice
    MatrixXd         _grainThickness;               // Thickness grid for shelf and cliff

	double 			 _momentInertia;		// moment of inertia in principal frame (purely diagonal terms)
	double 			 _theta;				// particle orientation
	double 			 _omega;				// particle angular velocity
	Vector2d		 _cmLset; 			// center of mass wrt the level set reference configuration (cm: at (0,0), I: diagonal)
	vector<Vector2d> _pointList; 	// list of points comprising the grain in real space (translated and rotated)
	int              _npoints;       //Number of points to represent the particle
	double 			 _radius;			// radius of bounding sphere
	Levelset2d 		 _lset;				// level set of grain
	Levelset2d 		 _lset0;			// Original level set of grain
	double			 _kn;				// normal stiffness
	double			 _ks;				// shear stiffness
	double			 _mu;				// interparticle friction

	bool             _fracFlag;         //Fracture flag for grain
	int              _fracLoc;           //Fracture horizontal location for grain (assuming vertical plane)


	double			 _density;			// particle density (default 1?)
	size_t      	 _morphologyID;		// ID representing the morphology type of the particle
	size_t 			 _ncontacts;			// number of contacts of the grain (wals + other grains)
	size_t 			 _id;				// ID (numbering) of the grain
	vector<double>	 _nodeShears;	// shears at each node
	vector<size_t>	 _nodeContact;  	// index of grain the node is contacting

	double 			 _scalePixToMM;

	vector<Vector2d> _refpointList; 	// list of points in reference config (center of mass is at (0,0,0) and I is diagonal)
	double 		     _yield;
	double		     _critStress;
	double		     _rsq; // squared radius
	bool			 _remove; //flag for removal
	double		     _d;
	
	double           _time_fail;   //Expected time for failure
	double           _origin_time_fail;   //Start of counter for Expected time for failure
	// vector<Vector2d> _nodeNormals;
	
	//Bonded Particle Method Properties (May 7, 2023)
    //bond properties
	vector<Vector3i>	_bondInfo;  // info on bonds (0->broken/1->intact ,  id, slave grain id)
	vector<size_t>		_bondpts;
	vector<size_t>		_bondptsf;
	vector<Vector2d>	_bondForceNormal;
	vector<Vector2d>	_bondForceShear;
	vector<double>		_bondThisMomentNormal;
	vector<Vector2d>	_bondNormals;
	vector<double>		_fSig;              //Stress acting on the bond
	vector<double>		_fTau;
	vector<double>	    _sigrF;              //Vector of bond critical stress reduction in case bonds weaken due to temperature or other factors, multiplies _sigC but depends on individual bonds. By default is 1.0;
	vector<double>	    _taurF;              //Vector of bond critical stress reduction in case bonds weaken due to temperature or other factors, multiplies _sigC but depends on individual bonds. By default is 1.0;
	double				_sigC;				// critical stress for bond fracture in pressure
	double				_tauC;
	double				_kn_bond; 		// normal stiffness for the bond
	double				_ks_bond;
	double				_cohesiveDistance;
	double				_bondRadius;
	double				_bondAreas;
	double				_bondMomentsInertia;
	
	vector<size_t>		_nnList;  //For regular DEM
	double	         	_vr; //verlet radius
	vector<size_t>      _bondOtherpts; //Bond for detection of follower grain
	vector<double>		_bondThisMomentShear;
	std::map<size_t, Vector2d> _nodeShearsv2; // shear force with each other particle
	std::map<size_t, Vector2d> _nodeNormals;// contact normals at each node
	
	bool                _load_grain;  //Apply load only on this tagged grain (sort of Neumann BC, just like fixed_grain is Dirichlet)
 	Vector2d            _position0;   //Initial position to measure displacement
	
};

#endif /* Grain2D_H_ */
