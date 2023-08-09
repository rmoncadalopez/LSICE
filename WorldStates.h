/*
 * WorldStates.h
 *
 *  Created on: October 25, 2016
 *  Author: Reid Kawamoto
 */

#ifndef WORLDSTATES_H_
#define WORLDSTATES_H_

struct GrainState2d {

	GrainState2d() {}
	GrainState2d(const size_t & ngrains, const size_t & nwalls) {
		_grainStress.resize(ngrains);
		_grainCoord.resize(ngrains);
		_grainForces.resize(ngrains);
		_grainMoments.resize(ngrains);
		_wallForces.resize(nwalls);
		_wallContacts.resize(nwalls);
	}
	GrainState2d(const Vector3d & svtemp, const vector<Vector2d> & gftemp, const vector<double> & gmtemp, const vector<Vector2d> & wtemp):
		_stressVoigt(svtemp), _grainForces(gftemp), _grainMoments(gmtemp), _wallForces(wtemp) {}

	void resize(const int & ngrains, const int & nwalls) {    //or size_t for nwalls???
		
		_grainCoord.resize(ngrains);
		_grainStress.resize(ngrains);
		_grainForces.resize(ngrains);
		_grainMoments.resize(ngrains);
		_grainContacts.resize(ngrains);
        _wallForces.resize(nwalls);
        _wallContacts.resize(nwalls);
	}

	void reset() {
		_stressVoigt << 0., 0., 0.;
		for (size_t i = 0; i < _grainForces.size(); i++) {
			_grainCoord[i].clear();
			_grainCoord[i].resize(0);
			_grainStress[i] << 0, 0, 0;
			_grainForces[i] << 0., 0.;
			_grainMoments[i] = 0.;
			_grainContacts[i] = 0;
		}
        for (size_t i = 0; i < _wallForces.size(); i++) {
            _wallForces[i] << 0, 0;
            _wallContacts[i] = 0;
        }
	}

	void operator+=(const GrainState2d & w) {
		_stressVoigt += w._stressVoigt; 
		for (size_t i = 0; i < _grainForces.size(); i++) {
			_grainCoord[i].insert(_grainCoord[i].end(),w._grainCoord[i].begin(),w._grainCoord[i].end());
			_grainStress[i] += w._grainStress[i];
			_grainForces[i] += w._grainForces[i];
			_grainMoments[i] += w._grainMoments[i];
			_grainContacts[i] += w._grainContacts[i];
		}
        for (size_t i = 0; i < _wallForces.size(); i++) {
            _wallForces[i] += w._wallForces[i];
            _wallContacts[i] += w._wallContacts[i];
        }
	}

	vector<vector<size_t>> _grainCoord; //coordination number of grains
	vector<Vector3d> _grainStress; //average particle stress
	Vector3d 			_stressVoigt;			// macroscopic stress of assembly
	vector<Vector2d> 	_grainForces;			// forces on grain
	vector<double>   	_grainMoments;			// moments on grain
	vector<size_t> 		_grainContacts;			// number of contacts
    vector<Vector2d>    _wallForces; // forces on wall
    vector<size_t>      _wallContacts;
    
};

struct CData {
	
	void operator+=(const CData & c) {
		_cpairs.insert( _cpairs.end(),	c._cpairs.begin(),	c._cpairs.end());
		_forces.insert( _forces.end(),	c._forces.begin(),	c._forces.end());
		_normals.insert(_normals.end(),	c._normals.begin(),	c._normals.end());
		_clocs.insert(_clocs.end(),		c._clocs.begin(),		c._clocs.end());
		_nodes.insert(_nodes.end(),		c._nodes.begin(),		c._nodes.end());
	}

	// maybe add more variables like branch vectors in the future
	vector<Vector2i> _cpairs;				// contact pairs
	vector<Vector2d> _forces;				// contact forces
	vector<Vector2d> _normals;				// contact normals
	vector<Vector2d> _clocs;				// location of contact points
	vector<size_t> _nodes; // nodes of _cpairs[i](0)
};

#endif /* WORLDSTATES_H_ */
