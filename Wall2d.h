/*
 * Wall2d.h
 *
 *  Created on: Jun 7, 2014
 *      Author: Reid
 */

#ifndef WALL2D_H_
#define WALL2D_H_

#include "definitions.h"
#include "Grain2d.h"

class Wall2d {
public:
	Wall2d() {_kn = 0; _ks = 0; _d = 0; _mu = 0; _id = 0;}
	
	Wall2d(const Vector2d & position, const Vector2d & normal, const double & kn, const double & ks, const double & mu, const size_t & id):
	_position(position), _normal(normal), _kn(kn), _ks(ks), _mu(mu), _id(id) {
		_velocity << 0.,0;
		_d = -_normal.dot(_position);
	}
	
	Wall2d(const Vector2d & position, const Vector2d & velocity, const Vector2d & normal, const double & kn, const double & ks, const double & mu, const size_t & id):
	_position(position), _velocity(velocity), _normal(normal), _kn(kn), _ks(ks), _mu(mu), _id(id) {
		_d = -_normal.dot(_position);
	}
	
	bool bcircleWallIntersection(Grain2d & grain) const {
		return _normal.dot(grain.getPosition())+_d < grain.getRadius();
	}
	
	bool findWallForceMoment(Grain2d & grain, Vector2d & force, double & moment, size_t & ncontacts, Vector3d & stress, const double & dt) const {
		// zero the input vars
		force << 0., 0.;
		moment = 0.;
		stress << 0,0,0;
		ncontacts = 0;
		bool 	 checkflag = false;			// changes to true if penetration exists
		double 	 penetration; 				// penetration amount
		Vector2d df;						// force increment on grain
		Vector2d v;							// relative velocity
		Vector2d tangent; 		            // surface tangent
		double   ds;						// relative displacement of a point of contact
		double sdot;						// projection of relative velocity into tangential direction
		Vector2d Fs;						// shear force
		Vector2d ptcm;						// point wrt the center of mass of the grain in real space
		double Fsmag;						// magnitude of shear force
        
        
        
        
       Vector2d relativeVel = Vector2d (0.,0.);
        
        

		for (size_t ptidx = 0; ptidx < grain.getPointList().size(); ptidx++) {
			penetration = _normal.dot(grain.getPointList()[ptidx]) + _d;
			if ( penetration < 0 ) {
				ncontacts++;
				ptcm = grain.getPointList()[ptidx] - grain.getPosition();
				ptcm = ptcm + penetration*ptcm/ptcm.norm();
				checkflag = true;
				

                
                
                
                double cres = .40;                                         // auxilliary variable for computing contact damping coeff below   .3
                double GamaN = -2*sqrt(_kn*grain.getMass())*log(cres)/sqrt(M_PI*M_PI + log(cres)*log(cres));
                
                
//              //Dampen wall contact
                relativeVel << grain.getVelocity()(0) - grain.getOmega()*grain.getCmLset()(1) - (_velocity(0) - grain.getOmega()*grain.getCmLset()(1)),
                               grain.getVelocity()(1) + grain.getOmega()*grain.getCmLset()(0) - (_velocity(1) + grain.getOmega()*grain.getCmLset()(0));
                // viscoelastic contact law
                df = -penetration*_normal*_kn - GamaN*_normal.dot(relativeVel)*_normal;
//                   
                
                
                
//              df = -penetration*_normal*_kn;
				force += df;
				// moment -= ptcm.cross(df);
				moment += ptcm(0)*df(1) - ptcm(1)*df(0);

//				stress(0) += df(0)*ptcm(0);
//				stress(1) += df(1)*ptcm(1);
//				stress(2) += 0.5*(df(1)*ptcm(0) + df(0)*ptcm(1));

				// Friction contribution
				// Find relative tangential kinematics
				v << _velocity(0) - grain.getVelocity()(0) + grain.getOmega()*ptcm(1),
					  _velocity(1) - grain.getVelocity()(1) - grain.getOmega()*ptcm(0);
				tangent << -_normal(1), _normal(0);
				sdot = tangent.dot(-v);
				ds = sdot*dt;
				// TODO: Compare the above with the following to make sure they are equivalent
				// ds = (v - v.dot(_normal)*_normal)*dt; 
				
				grain.getNodeShearsNonConst()[ptidx] += ds*_kn; // technically should be ks but walls don't have ks
				grain.getNodeContactNonConst()[ptidx] = _id;

				if (grain.getNodeShearsNonConst()[ptidx] > 0) {
					Fsmag = std::min(grain.getNodeShearsNonConst()[ptidx], df.norm()*_mu );
				}
				else {
					Fsmag = std::max(grain.getNodeShearsNonConst()[ptidx], -df.norm()*_mu );
				}
				Fs = -tangent*Fsmag;
				grain.getNodeShearsNonConst()[ptidx] = Fsmag;
				force += Fs;
				stress(0) += (df(0)+Fs(0))*ptcm(0);
				stress(1) += (df(1)+Fs(1))*ptcm(1);
				stress(2) += 0.5*((df(1)+Fs(1))*ptcm(0) + (df(0)+Fs(0))*ptcm(1));
			    moment += ptcm(0)*(df(1)+Fs(1)) - ptcm(1)*(df(0)+Fs(0));
			}
			else if (grain.getNodeContactNonConst()[ptidx] == _id ){
				grain.getNodeContactNonConst()[ptidx] = INT_MAX;
				grain.getNodeShearsNonConst()[ptidx] = 0.0;
			}
		}
		return checkflag;
	}

	size_t FindWallForceForFrac(const Grain2d & grain, double & force){
		double 		penetration;
		Vector2d 	df;
		Vector2d 	forcevector;
		Vector2d 	tangent;
		force = 0;
		size_t 		contactpoint =1;

		for (size_t p = 0; p < grain.getPointList().size(); p++){

			penetration = _normal.dot(grain.getPointList()[p]) + _d;
			if (penetration < 0){
				double forcei;
				df = -penetration*_normal*_kn;
				tangent << -_normal(1), _normal(0);
				forcevector = df-tangent*grain.getNodeShears()[p];
				forcei = forcevector.norm();

				if (forcei > force) {
					force = forcei;
					contactpoint = p;
				}

			}
		}
		return contactpoint;
	}
	
	// netForce is in the direction of the normal (outward of wall)
	void takeForceTimestep(const double & netForce, const size_t & ncontacts) {
		const double _alpha = 0.9;
		if (ncontacts < 2)
			moveWall( _normal*_alpha*netForce/_kn/2. );
		else
			moveWall( _normal*_alpha*netForce/_kn/double(ncontacts) );
	}
	
	Vector2d takeForceTimestep2(const double & netForce, const size_t & ncontacts) {
		const double _alpha = 0.9;
		if (ncontacts < 2) {
			moveWall( _normal*_alpha*netForce/_kn/2. );
			return _normal*_alpha*netForce/_kn/2.;
		}
		else {
			moveWall( _normal*_alpha*netForce/_kn/double(ncontacts) );
			return _normal*_alpha*netForce/_kn/double(ncontacts);
		}
	}
	
	void rotateWall(double amount) {
		_normal << cos(amount)*_normal(0)-sin(amount)*_normal(1), sin(amount)*_normal(0)+cos(amount)*_normal(1);
		_d = -(_normal.dot(_position) );
	}
	void moveWall(Vector2d amount) {
		_position += amount;
		_d = -(_normal.dot(_position) );
	}
	void changeKn(const double & kn) {
		_kn = kn;
	}
	void changeKs(const double & ks) {
		_ks = ks;
	}
	void changeMu(const double & mu) {
		_mu = mu;
	}
	void changeVelocity(Vector2d vel) {
		_velocity = vel;
	}
	
	const Vector2d & getPosition() const {
		return _position;
	}
	const Vector2d & getVelocity() const {
		return _velocity;
	}
	const Vector2d & getNormal() const {
		return _normal;
	}
	const double & getKn() const {
		return _kn;
	}
	const size_t & getId() const {
		return _id;
	}

private:
	Vector2d _position;
	Vector2d _velocity;
	Vector2d _normal; // outward normal
	double 	_kn;
	double 	_ks;
	double 	_mu;
	size_t	_id;
	double	_d;
};


#endif /* WALL2D_H_ */
