/*
 * readInputFile.h
 *
 *  Created on: Jun 3, 2014
 *      Modified by: Liuchi
 */

#ifndef READINPUTFILE_H_
#define READINPUTFILE_H_

#include "definitions.h"
#include "Levelset2d.h"
#include "Grain2d.h"

// creates a vector of grain objects from input files
vector<Grain2d> generateGrainsFromFiles(string morphologyMaterialIDs,
										string morphologyDir,
	                                   	string Positions,
	                                   	string Velocities,
	                                   	string Temperatures,
	                                   	string Thickness) {

	string line;
	string line_property;
	string line_velocity;
	string line_position;

    string line_temper;
	string line_thick;

	string partial;
	string propertyfile;
	istringstream iss;


	ifstream file_information(morphologyMaterialIDs.c_str());		// construct an ifstream and open the grain morphology file
	ifstream file_position(Positions.c_str());						// construct another ifstream and open the grain position file
	ifstream file_velocity(Velocities.c_str());						// construct another ifstream and open the grain velocity file

    ifstream file_temper(Temperatures.c_str());						// construct another ifstream and open the grain temperature file
	ifstream file_thick(Thickness.c_str());					    	// construct another ifstream and open the grain thickness file
	

	getline(file_information, line);								// read a line (first) of the morphology file (number of particles)
    iss.str(line);													// copy line string to into iss (we basically bind iss to the line we just read)
    getline(iss, partial, ' ');										// extracts characters until delimiter ' ' is found; the latter is extracted and discarded 
	size_t numberOfGrains = atoi(line.c_str());						// converts the string to an integer
	//cout<<"Number of Grains INPUT: "<< numberOfGrains << endl;
	iss.clear();													// clear the error state of the stream (e.g. end-of-file -> no error)
    char tempfilename[100];

    // Initialize the vector of grain objects
	vector<Grain2d> grainList(numberOfGrains);

	// temp stuff
	Vector2d point;
	Vector2d position;
	double theta;
	Vector2d velocity;
	double omega;

	// Go through each grain 
	for (size_t grainidx = 0; grainidx < numberOfGrains; grainidx++) {

        // Read morphology index for each particle - each index has each own property .dat file
        getline(file_information, line);
        iss.str(line);
        getline(iss, partial, ' ');
		int morphologyID = atoi(partial.c_str());
        iss.clear();
        sprintf(tempfilename, "grainproperty%d.dat", morphologyID);
        
        propertyfile = morphologyDir + tempfilename;
        ifstream file_gp(propertyfile.c_str());

        // mass
        getline(file_gp, line_property);
        double mass = atof(line_property.c_str());
        double mass0 = atof(line_property.c_str());


        // moment of inertia
        getline(file_gp, line_property);
		double momentOfInertia = atof(line_property.c_str());

		// cmLset (center of mass)
		getline(file_gp, line_property);
		Vector2d cmLset;
		iss.str(line_property);
		getline(iss, partial, ' ');
		cmLset(0) = atof(partial.c_str());  //Multiply by zero to reset origin?
		getline(iss, partial, ' ');
		cmLset(1) = atof(partial.c_str());  //Multiply by zero to reset origin?
		iss.clear();

		// number of points on the grain surface (INTEGER)
		getline(file_gp, line_property);
		int npoints = atoi(line_property.c_str());

		// the point coordinates
		getline(file_gp, line_property);
		vector<Vector2d> pointList(npoints);
		iss.str(line_property);
		for (int ptidx = 0; ptidx < npoints; ptidx++) {
			getline(iss, partial, ' ');
			point(0) = atof(partial.c_str());
			getline(iss, partial, ' ');
//            point(1) = atof(partial.c_str());  //Mirror Flip Points
//            point(1)=-point(1);
            point(1) = atof(partial.c_str());
			pointList[ptidx] = point;
		}
		iss.clear();

		vector<Vector2d> refpointList(npoints);
		refpointList = pointList;

		// bounding box radius
		getline(file_gp, line_property);
		double bboxRadius = atof(line_property.c_str());
		//cout<<"Radius INPUT: "<< bboxRadius << endl;

		// level set dimensions (INTEGERS)
		getline(file_gp, line_property);
		iss.str(line_property);
		getline(iss, partial, ' ');
		int xdim = atoi(partial.c_str());
		getline(iss, partial, ' ');
		int ydim = atoi(partial.c_str());
		iss.clear();

		// level set
		getline(file_gp, line_property);
		vector<double> lsetvec(xdim*ydim);
		iss.str(line_property);
		for (int i = 0; i < xdim*ydim; i++) {
			getline(iss, partial, ' ');
			lsetvec[i] = atof(partial.c_str());
		}
		iss.clear();

		// kn
		getline(file_gp, line_property);
		double kn = atof(line_property.c_str());

        // ks
		getline(file_gp, line_property);
		double ks = atof(line_property.c_str());

        // mu
		getline(file_gp, line_property);
		double mu = atof(line_property.c_str());

		// cohesive distance
		getline(file_gp, line_property);
		double cohesiveDistance = atof(line_property.c_str());

		// bond modulus
		getline(file_gp, line_property);
		double bondModulus = atof(line_property.c_str());

		// bond area
		getline(file_gp, line_property);
		double bondArea = atof(line_property.c_str());

		// bond moment of inertia
		getline(file_gp, line_property);
		double bondMomentInertia = atof(line_property.c_str());

		// critical bond pressure
		getline(file_gp, line_property);
		double sigC = atof(line_property.c_str());

		// critical bond shear
		getline(file_gp, line_property);
		double tauC = atof(line_property.c_str());

        // Clear string for property file ready for next grain - do we need this?
        propertyfile.clear();
        line_property.clear();

	    // Read position and theta from position file
	    getline(file_position, line_position);
		iss.str(line_position);
		getline(iss, partial, ' ');
		position(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		position(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		theta = atof(partial.c_str());
		iss.clear();

		// Read velocity and omega from velocity file
		getline(file_velocity, line_velocity);
		iss.str(line_velocity);
		getline(iss, partial, ' ');
		velocity(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		velocity(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		omega = atof(partial.c_str());;
		iss.clear();


		//Read Temperature from File
        getline(file_temper, line_temper);
		double temper = atof(line_temper.c_str());

		//Read Thickness from File
        getline(file_thick, line_thick);
		double thick = atof(line_thick.c_str());

        //double mass0=mass;
        double temper0=temper;
        double thick0=thick;

        // 2D Thickness of Level Set (Initialize all Level Set with an initial thickness. This will be changed later)
			vector<double> MthickV(xdim*ydim);
			for (size_t i = 0; i < xdim*ydim; i++) {
				MthickV[i] = thick;
			}


        vector<double> Utemper(800);
        //double Utemper[80000];
        for (int iii=0; iii<(800); iii++)
        {
           Utemper[iii]=temper0;
        }

        vector<double> Utemper0(800);
        Utemper0=Utemper;

        // 2D Temp of Level Set (Initialize all Level Set with an initial temperature. This will be differentiated later)
		vector<double> grainTemp2Dvec(xdim*ydim);
		for (size_t i = 0; i < xdim*ydim; i++) {
			grainTemp2Dvec[i] = temper;
		}
		vector<double> grainTemp2Dvec0(xdim*ydim);
		grainTemp2Dvec0=grainTemp2Dvec;

		// 2D Internal Continuum Stress Level Set (for Ice Rheology and Getting Rheology/Elastic Stresses)
        vector<double> grainInnerStress2D(xdim*ydim);
        for (size_t i = 0; i < xdim*ydim; i++) {
			grainInnerStress2D[i] = 0.0 ;
		}
	   	
		// Create level set objects
		Levelset2d lset(lsetvec, xdim, ydim);
		Levelset2d lset0(lsetvec, xdim, ydim);
	   Levelset2d grainTemp2D(grainTemp2Dvec, xdim, ydim);
		Levelset2d grainTemp2D0(grainTemp2Dvec0, xdim, ydim);

		Levelset2d grainStress2D(grainInnerStress2D, xdim, ydim);
		Levelset2d Mthick(MthickV, xdim, ydim);

		//Create Damage Matrix for tracking evolution of damaged material
		//MatrixXd grainDamage = MatrixXd::Zero(ydim, xdim);
		vector<Vector3d> grainDamage;


		//Create Thickness Matrix for tracking evolution of damage and thickness RIGHT NOW INITIALLY CONSTANT FOR ALL (UNITARY) *-->
		double Thicc = 1.0;
		MatrixXd grainThickness = MatrixXd::Constant(ydim, xdim, Thicc);		


		bool fracFlag = false;
        int fracLoc = 0; //Should this be initialized here????
        
        //Temporary for time failure
        const double time_fail = 10.000;
        const double origin_time_fail = 0.000;


		// Update grain object in the vector that was created at the beginning of this function
		grainList[grainidx] = Grain2d(mass, position, velocity, mass0, temper, thick, Mthick, temper0, thick0, Utemper, Utemper0, grainTemp2D, grainTemp2D0, momentOfInertia, theta, omega, cmLset, pointList, npoints, bboxRadius, lset, lset0, grainidx, morphologyID, kn, ks, mu, fracFlag, fracLoc, grainStress2D, grainDamage, grainThickness, refpointList, time_fail, origin_time_fail);
	}

	return grainList;
} // end generateGrainsFromFiles

// creates a vector of grain objects from input files
Grain2d generateGrainFromFiles(string morphologyMaterialIDs,
										string morphologyDir,
	                                   	string Positions,
	                                   	string Velocities,
	                                   	string Temperatures,
	                                   	string Thickness) {

	string line;
	string line_property;
	string line_velocity;
	string line_position;

    string line_temper;
	string line_thick;


	string partial;
	string propertyfile;
	istringstream iss;

	ifstream file_information(morphologyMaterialIDs.c_str());		// construct an ifstream and open the grain morphology file
	ifstream file_position(Positions.c_str());						// construct another ifstream and open the grain position file
	ifstream file_velocity(Velocities.c_str());						// construct another ifstream and open the grain velocity file

    ifstream file_temper(Temperatures.c_str());						// construct another ifstream and open the grain temperature file
	ifstream file_thick(Thickness.c_str());					    	// construct another ifstream and open the grain thickness file


    char tempfilename[100];

	// temp stuff
	Vector2d point;
	Vector2d position;
	double theta;
	Vector2d velocity;
	double omega;

	// Read morphology index for each particle - each index has each own property .dat file
	getline(file_information, line);
	getline(file_information, line);
	iss.str(line);
	getline(iss, partial, ' ');
	int morphologyID = atoi(partial.c_str());
	iss.clear();
	sprintf(tempfilename, "grainproperty%d.dat", morphologyID);

	propertyfile = morphologyDir + tempfilename;
	ifstream file_gp(propertyfile.c_str());

	// mass
	getline(file_gp, line_property);
	double mass = atof(line_property.c_str());
	double mass0 = atof(line_property.c_str());

	// moment of inertia
	getline(file_gp, line_property);
	double momentOfInertia = atof(line_property.c_str());

	// cmLset (center of mass)
	getline(file_gp, line_property);
	Vector2d cmLset;
	iss.str(line_property);
	getline(iss, partial, ' ');
	cmLset(0) = atof(partial.c_str());
	getline(iss, partial, ' ');
	cmLset(1) = atof(partial.c_str());
	iss.clear();

	// number of points on the grain surface (INTEGER)
	getline(file_gp, line_property);
	int npoints = atoi(line_property.c_str());

	// the point coordinates
	getline(file_gp, line_property);
	vector<Vector2d> pointList(npoints);
	iss.str(line_property);
	for (int ptidx = 0; ptidx < npoints; ptidx++) {
		getline(iss, partial, ' ');
		point(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		point(1) = atof(partial.c_str());
		pointList[ptidx] = point;
	}
	iss.clear();

	vector<Vector2d> refpointList(npoints);
	refpointList = pointList;

	// bounding box radius
	getline(file_gp, line_property);
	double bboxRadius = atof(line_property.c_str());

	// level set dimensions (INTEGERS)
	getline(file_gp, line_property);
	iss.str(line_property);
	getline(iss, partial, ' ');
	int xdim = atoi(partial.c_str());
	getline(iss, partial, ' ');
	int ydim = atoi(partial.c_str());
	iss.clear();

	// level set
	getline(file_gp, line_property);
	vector<double> lsetvec(xdim*ydim);
	iss.str(line_property);
	for (int i = 0; i < xdim*ydim; i++) {
		getline(iss, partial, ' ');
		lsetvec[i] = atof(partial.c_str());
	}
	iss.clear();

	// kn
	getline(file_gp, line_property);
	double kn = atof(line_property.c_str());

	// ks
	getline(file_gp, line_property);
	double ks = atof(line_property.c_str());

	// mu
	getline(file_gp, line_property);
	double mu = atof(line_property.c_str());

	// cohesive distance
	getline(file_gp, line_property);
	double cohesiveDistance = atof(line_property.c_str());

	// bond modulus
	getline(file_gp, line_property);
	double bondModulus = atof(line_property.c_str());

	// bond area
	getline(file_gp, line_property);
	double bondArea = atof(line_property.c_str());

	// bond moment of inertia
	getline(file_gp, line_property);
	double bondMomentInertia = atof(line_property.c_str());

	// critical bond pressure
	getline(file_gp, line_property);
	double sigC = atof(line_property.c_str());

	// critical bond shear
	getline(file_gp, line_property);
	double tauC = atof(line_property.c_str());

	// Clear string for property file ready for next grain - do we need this?
	propertyfile.clear();
	line_property.clear();

	// Read position and theta from position file
	getline(file_position, line_position);
	iss.str(line_position);
	getline(iss, partial, ' ');
	position(0) = atof(partial.c_str());
	getline(iss, partial, ' ');
	position(1) = atof(partial.c_str());
	getline(iss, partial, ' ');
	theta = atof(partial.c_str());
	iss.clear();

	// Read velocity and omega from velocity file
	getline(file_velocity, line_velocity);
	iss.str(line_velocity);
	getline(iss, partial, ' ');
	velocity(0) = atof(partial.c_str());
	getline(iss, partial, ' ');
	velocity(1) = atof(partial.c_str());
	getline(iss, partial, ' ');
	omega = atof(partial.c_str());;
	iss.clear();

	//Read Temperature from File
    getline(file_temper, line_temper);
	double temper = atof(line_temper.c_str());

	//Read Thickness from File
    getline(file_thick, line_thick);
	double thick = atof(line_thick.c_str());

    //double mass0=mass;
    double temper0=temper;
    double thick0=thick;


    // 2D Thickness of Level Set (Initialize all Level Set with an initial thickness. This will be changed later)
 	 vector<double> MthickV(xdim*ydim);
    for (size_t i = 0; i < xdim*ydim; i++) {
		MthickV[i] = thick;
	 }

    vector<double> Utemper(800);
    //double Utemper[80000];
    for (int iii=0; iii<(800); iii++)
    {
        Utemper[iii]=temper0;
    }

    vector<double> Utemper0(800);
    Utemper0=Utemper;

    // 2D Temp of Level Set (Initialize all Level Set with an initial temperature. This will be differentiated later)
	vector<double> grainTemp2Dvec(xdim*ydim);
	for (size_t i = 0; i < xdim*ydim; i++) {
		grainTemp2Dvec[i] = temper;
	}
	vector<double> grainTemp2Dvec0(xdim*ydim);
	grainTemp2Dvec0=grainTemp2Dvec;

	// 2D Internal Continuum Stress Level Set (for Ice Rheology and Getting Rheology/Elastic Stresses)
    vector<double> grainInnerStress2D(xdim*ydim);
    for (size_t i = 0; i < xdim*ydim; i++) {
		grainInnerStress2D[i] = 0.0 ;
	}
   	
	// Create level set objects
	Levelset2d lset(lsetvec, xdim, ydim);
	Levelset2d lset0(lsetvec, xdim, ydim);
    Levelset2d grainTemp2D(grainTemp2Dvec, xdim, ydim);
	Levelset2d grainTemp2D0(grainTemp2Dvec0, xdim, ydim);

	Levelset2d grainStress2D(grainInnerStress2D, xdim, ydim);
	Levelset2d Mthick(MthickV, xdim, ydim);


	//Create Damage Matrix for tracking evolution of damaged material
	//MatrixXd grainDamage = MatrixXd::Zero(ydim, xdim);
	vector<Vector3d> grainDamage;


	//Create Thickness Matrix for tracking evolution of damage and thickness RIGHT NOW INITIALLY CONSTANT FOR ALL (UNITARY) *-->
	double Thicc = 1.0;
	MatrixXd grainThickness = MatrixXd::Constant(ydim, xdim, Thicc);		


	bool fracFlag = false;
    int fracLoc = 0; //Should this be initialized here????
    
    //Temporary for time failure
    const double time_fail = 10.000;
    const double origin_time_fail = 0.000;


	// Update grain object in the vector that was created at the beginning of this function
	Grain2d grain = Grain2d(mass, position, velocity, mass0, temper, thick, Mthick, temper0, thick0, Utemper, Utemper0, grainTemp2D, grainTemp2D0, momentOfInertia, theta, omega, cmLset, pointList, npoints, bboxRadius, lset, lset0, 0, morphologyID, kn, ks, mu, fracFlag, fracLoc, grainStress2D, grainDamage, grainThickness, refpointList, time_fail, origin_time_fail);


	return grain;
} // end generateGrainsFromFiles


// creates a vector of grain objects from a single input file (e.g. from level set imaging)
vector<Grain2d> generateGrainsFromFile(string filename,
										string Temperatures,
										string Thickness) {

	string line;

    string line_temper;
	string line_thick;

	string partial;
	istringstream iss;
	ifstream file(filename.c_str());								


    ifstream file_temper(Temperatures.c_str());						// construct another ifstream and open the grain temperature file
	ifstream file_thick(Thickness.c_str());					    	// construct another ifstream and open the grain thickness file


	getline(file, line);
	size_t numberOfGrains = atoi(line.c_str());

    // Initialize the vector of grain objects
	vector<Grain2d> grainList(numberOfGrains);

	// temp stuff
	Vector2d point;
	Vector2d position;
	Vector2d velocity;

	// morphologyID is of no importance here since each grain is special; just pass any ID to the constructor
	size_t morphologyID = 0;

	// Go through each grain 
	for (size_t grainidx = 0; grainidx < numberOfGrains; grainidx++) {

        // mass
        getline(file, line);
        double mass = atof(line.c_str());
        double mass0 = atof(line.c_str());

		// position
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		position(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		position(1) = atof(partial.c_str());
		iss.clear();

		// velocity
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		velocity(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		velocity(1) = atof(partial.c_str());
		iss.clear();

        // moment of inertia
        getline(file, line);
		double momentOfInertia = atof(line.c_str());

		// theta
		getline(file, line);
		double theta = atof(line.c_str());

		// omega
		getline(file, line);
		double omega = atof(line.c_str());

		// cmLset (center of mass)
		getline(file, line);
		Vector2d cmLset;
		iss.str(line);
		getline(iss, partial, ' ');
		cmLset(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		cmLset(1) = atof(partial.c_str());
		iss.clear();

		// number of points on the grain surface (INTEGER)
		getline(file, line);
		int npoints = atoi(line.c_str());

		// the point coordinates
		getline(file, line);
		vector<Vector2d> pointList(npoints);
		iss.str(line);
		for (int ptidx = 0; ptidx < npoints; ptidx++) {
			getline(iss, partial, ' ');
			point(0) = atof(partial.c_str());
			getline(iss, partial, ' ');
			point(1) = atof(partial.c_str());
			pointList[ptidx] = point;
		}
		iss.clear();

		vector<Vector2d> refpointList(npoints);
		refpointList = pointList;

		// bounding box radius
		getline(file, line);
		double bboxRadius = atof(line.c_str());

		// level set dimensions (INTEGERS)
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		int xdim = atoi(partial.c_str());
		getline(iss, partial, ' ');
		int ydim = atoi(partial.c_str());
		iss.clear();

		// level set
		getline(file, line);
		vector<double> lsetvec(xdim*ydim);
		iss.str(line);
		for (int i = 0; i < xdim*ydim; i++) {
			getline(iss, partial, ' ');
			lsetvec[i] = atof(partial.c_str());
		}
		iss.clear();

		// kn
		getline(file, line);
		double kn = atof(line.c_str());

        // ks
		getline(file, line);
		double ks = atof(line.c_str());

        // mu
		getline(file, line);
		double mu = atof(line.c_str());

        // cohesive distance
		getline(file, line);
		double cohesiveDistance = atof(line.c_str());

		// bond modulus
		getline(file, line);
		double bondModulus = atof(line.c_str());

		// bond area
		getline(file, line);
		double bondArea = atof(line.c_str());

		// bond moment of inertia
		getline(file, line);
		double bondMomentInertia = atof(line.c_str());

		// critical bond pressure
		getline(file, line);
		double sigC = atof(line.c_str());

		// critical bond shear
		getline(file, line);
		double tauC = atof(line.c_str());

		//Read Temperature from File
	    getline(file_temper, line_temper);
		double temper = atof(line_temper.c_str());

		//Read Thickness from File
	    getline(file_thick, line_thick);
		double thick = atof(line_thick.c_str());

        //double mass0=mass;
		double temper0=temper;
      double thick0=thick;


       // 2D Thickness of Level Set (Initialize all Level Set with an initial thickness. This will be changed later)
	 	 vector<double> MthickV(xdim*ydim);
	    for (size_t i = 0; i < xdim*ydim; i++) {
			MthickV[i] = thick;
		 }



        vector<double> Utemper(800);
        //double Utemper[80000];
        for (int iii=0; iii<(800); iii++)
        {
           Utemper[iii]=temper0;
        }

        vector<double> Utemper0(800);
        Utemper0=Utemper;

        // 2D Temp of Level Set (Initialize all Level Set with an initial temperature. This will be differentiated later)
		vector<double> grainTemp2Dvec(xdim*ydim);
		for (size_t i = 0; i < xdim*ydim; i++) {
			grainTemp2Dvec[i] = temper;
		}
		vector<double> grainTemp2Dvec0(xdim*ydim);
		grainTemp2Dvec0=grainTemp2Dvec;

		// 2D Internal Continuum Stress Level Set (for Ice Rheology and Getting Rheology/Elastic Stresses)
        vector<double> grainInnerStress2D(xdim*ydim);
        for (size_t i = 0; i < xdim*ydim; i++) {
			grainInnerStress2D[i] = 0.0 ;
		}
	   	
		// Create level set objects
		Levelset2d lset(lsetvec, xdim, ydim);
		Levelset2d lset0(lsetvec, xdim, ydim);
	   Levelset2d grainTemp2D(grainTemp2Dvec, xdim, ydim);
		Levelset2d grainTemp2D0(grainTemp2Dvec0, xdim, ydim);

		Levelset2d grainStress2D(grainInnerStress2D, xdim, ydim);
		Levelset2d Mthick(MthickV, xdim, ydim);


		//Create Damage Matrix for tracking evolution of damaged material
		//MatrixXd grainDamage = MatrixXd::Zero(ydim, xdim);
		vector<Vector3d> grainDamage;


		//Create Thickness Matrix for tracking evolution of damage and thickness RIGHT NOW INITIALLY CONSTANT FOR ALL (UNITARY) *-->
		double Thicc = 1.0;
		MatrixXd grainThickness = MatrixXd::Constant(ydim, xdim, Thicc);		


		bool fracFlag = false;
      int fracLoc = 0; //Should this be initialized here????
      
         //Temporary for time failure
        const double time_fail = 10.000;
        const double origin_time_fail = 0.000;
      

		// Update grain object in the vector that was created at the beginning of this function
		grainList[grainidx] = Grain2d(mass, position, velocity, mass0, temper, thick, Mthick, temper0, thick0, Utemper, Utemper0, grainTemp2D, grainTemp2D0, momentOfInertia, theta, omega, cmLset, pointList, npoints, bboxRadius, lset, lset0, grainidx, morphologyID, kn, ks, mu, fracFlag, fracLoc, grainStress2D, grainDamage, grainThickness, refpointList, time_fail, origin_time_fail);
	}

	return grainList;
} // end generateGrainsFromFile


vector<Vector2d> readPositionFile(string filename, size_t ngrains) {
	ifstream file(filename.c_str());
	string  line;
	string partial;
	istringstream iss;
	vector<Vector2d> positions;
	positions.resize(ngrains);
	for (size_t grainidx = 0; grainidx < ngrains; grainidx++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		positions[grainidx](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		positions[grainidx](1) = atof(partial.c_str());
		iss.clear();
	}
	return positions;
} // end readPositionFile


vector<double> readRotationFile(string filename, size_t ngrains) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	vector<double> rotations;
	rotations.resize(ngrains);
	for (size_t grainidx = 0; grainidx < ngrains; grainidx++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		rotations[grainidx] = atof(partial.c_str());
		iss.clear();
	}
	return rotations;
} // end readRotationFile


// For imposing strain boundary conditions through wall displacements
vector<double> readWallPositionFile(string filename, size_t nwalls) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	vector<double> wallPositions;
	wallPositions.resize(nwalls);
	for (size_t wallidx = 0; wallidx < nwalls; wallidx++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		wallPositions[wallidx] = atof(partial.c_str());
		iss.clear();
	}
	return wallPositions;
} // end readWallPositionFile


vector<Matrix<double, 5, 1> > readContactParameters(string filename, size_t contacttype) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	vector<Matrix<double, 5, 1> > contactParameters;
	contactParameters.resize(contacttype);
	for (size_t grainidx = 0; grainidx < contacttype; grainidx++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		contactParameters[grainidx](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		contactParameters[grainidx](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		contactParameters[grainidx](2) = atof(partial.c_str());
		getline(iss, partial, ' ');
		contactParameters[grainidx](3) = atof(partial.c_str());
		getline(iss, partial, ' ');
		contactParameters[grainidx](4) = atof(partial.c_str());
		iss.clear();
	}
	return contactParameters;	
} // end readContactParameters


vector<double> readMaximumValues(string filename) {
    // return the minimum, maximum positions of all grains: x_min, x_max, y_min, y_max, max_radius;
    ifstream file(filename.c_str());
    string line;
    vector<double> extremevalues;
    for (size_t index = 0; index < 5; index++) {
        getline(file, line);
        extremevalues.push_back(atof(line.c_str()));
    }
    return extremevalues;
} // end readMaximumValues


vector<Vector3d> readCorrectionFile(string filename, size_t ngrains) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	vector<Vector3d> corrections;
	corrections.resize(ngrains);
	for (size_t grainidx = 0; grainidx < ngrains; grainidx++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		corrections[grainidx](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		corrections[grainidx](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		corrections[grainidx](2) = atof(partial.c_str());
		iss.clear();

	}
	return corrections;
} // end readCorrectionFile


vector<double> readScalarFieldFile(string filename, size_t n) {
	ifstream file(filename.c_str());
	string  line;
	vector<double> field;
	for (size_t idx = 0; idx < n; idx++) {
		getline(file, line);
		field.push_back(atof(line.c_str()));
	}
	return field;
}

MatrixXd readStressFile (string filename, string filenumber)
{
	
    ifstream file(filename.c_str());		// construct an ifstream and open the stress file
    ifstream file2(filenumber.c_str());		// construct an ifstream and open the number element file
    string   line;
    string 	partial;
	istringstream iss;

	getline(file2, line);
	iss.str(line);
	getline(iss, partial, ' ');
	int nelem = atoi(partial.c_str());
	iss.clear();

	MatrixXd stress_mat = MatrixXd::Zero(nelem,5);
    

	// the point coordinates
	for (int i = 0; i < nelem; i++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		stress_mat(i,0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		stress_mat(i,1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		stress_mat(i,2) = atof(partial.c_str());
		getline(iss, partial, ' ');
		stress_mat(i,3) = atof(partial.c_str());
		getline(iss, partial, ' ');
		stress_mat(i,4) = atof(partial.c_str());
		iss.clear();
	}
	
	return stress_mat;
}

//Smaller for Efficiency
MatrixXd readStressFile2 (string filename, string filenumber)
{
	
    ifstream file(filename.c_str());		// construct an ifstream and open the stress file
    ifstream file2(filenumber.c_str());		// construct an ifstream and open the number element file
    string   line;
    string 	partial;
	istringstream iss;

	getline(file2, line);
	iss.str(line);
	getline(iss, partial, ' ');
	int nelem = atoi(partial.c_str());
	iss.clear();

	MatrixXd stress_mat = MatrixXd::Zero(nelem,3);
    

	// the point coordinates
	for (int i = 0; i < nelem; i++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		stress_mat(i,0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		stress_mat(i,1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		stress_mat(i,2) = atof(partial.c_str());
		iss.clear();
	}
	
	return stress_mat;
}

vector<Vector2d> readPointFile (string filename, string filenumber)
{
	
    ifstream file(filename.c_str());		// construct an ifstream and open the stress file
    ifstream file2(filenumber.c_str());		// construct an ifstream and open the number element file
    string   line;
    string 	partial;
	istringstream iss;

	getline(file2, line);
	iss.str(line);
	getline(iss, partial, ' ');
	int nelem = atoi(partial.c_str());
	iss.clear();

	vector<Vector2d> points_out(nelem);
    

	// the point coordinates
	for (int i = 0; i < nelem; i++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		points_out[i](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		points_out[i](1) = atof(partial.c_str());
		iss.clear();
	}
	
	return points_out;
}

vector<Vector3d> readElemFile (string filename, string filenumber)
{
	
	ifstream file(filename.c_str());		// construct an ifstream and open the stress file
	ifstream file2(filenumber.c_str());		// construct an ifstream and open the number element file
	string   line;
	string 	partial;
	istringstream iss;

	getline(file2, line);
	iss.str(line);
	getline(iss, partial, ' ');
	int nelem = atoi(partial.c_str());
	iss.clear();

	vector<Vector3d> elems_out(nelem);
    

	// the point coordinates
	for (int i = 0; i < nelem; i++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		elems_out[i](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		elems_out[i](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		elems_out[i](2) = atof(partial.c_str());
		iss.clear();
	}
	
	return elems_out;
}

//Smaller for Efficiency
void readSimInputFile (string filename, size_t caseNo, size_t & melt_step, size_t & break_step, double & IlimDiamMax, double & IlimDiamMin, double & Vert_q, double & Q_atmos)
{
	
   ifstream file(filename.c_str());		// construct an ifstream and open the input file
   string   line;
   string 	partial;
	istringstream iss;

	//Declare inputs or uses addresses
	double input0, input1, input2, input3, input4, input5;

   //Get data from text, move to the right case line   	
	for (size_t i = 0; i < caseNo; i++)
	{
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		iss.clear();
   }

	getline(file, line);
	iss.str(line);
	getline(iss, partial, ' ');
	input0 = atof(partial.c_str());
	melt_step = size_t(input0);
	getline(iss, partial, ' ');
	input1 = atof(partial.c_str());
	break_step = size_t(input1);
	getline(iss, partial, ' ');
	input2 = atof(partial.c_str());
   IlimDiamMax = input2;
	getline(iss, partial, ' ');
	input3 = atof(partial.c_str());
	IlimDiamMin = input3;
	getline(iss, partial, ' ');
	input4 = atof(partial.c_str());
	Vert_q = input4;
	getline(iss, partial, ' ');
	input5 = atof(partial.c_str());
	Q_atmos = input5;

	iss.clear();
	
	//Print Results
	cout << "Case Results for Case #: " << caseNo << endl;
	cout << "Melt Step: " << melt_step << " MAXD: " << IlimDiamMax << " MIND: " << IlimDiamMin << " Break Step: " << break_step << " Qvert: " << Vert_q << " QATM: " << Q_atmos << endl;

}

//Read Conc. File (1 line, 2 columns)
void readConcFile (string filename, double & fine_conc, double & coarse_conc)
{
	ifstream file(filename.c_str());		// construct an ifstream and open the input file
   string   line;
   string 	partial;
	istringstream iss;

	//Declare inputs or uses addresses
	double input0, input1;

	getline(file, line);
	iss.str(line);
	getline(iss, partial, ' ');
	input0 = atof(partial.c_str());
	fine_conc = input0;
	getline(iss, partial, ' ');
	input1 = atof(partial.c_str());
	coarse_conc = input1;

	iss.clear();

}

vector<Vector2d> readOutTempFile (string filename)
{
	
	//First count of number of lines in the file
	int count = 0;
	string linect;
	ifstream filect(filename.c_str());
	while (getline(filect, linect))
	{
		count++;
	}

	ifstream file(filename.c_str());		// construct an ifstream and open the stress file
	string   line;
	string 	partial;
	istringstream iss;

	vector<Vector2d> out_temp(count);
    

	// Reading file lines
	for (int i = 0; i < count; i++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		out_temp[i](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		out_temp[i](1) = atof(partial.c_str());
		iss.clear();

		std::cout << out_temp[i](1) << std::endl;
	}
	
	return out_temp;
}

vector<Vector3d> readWaveFile (string filename)
{
	
	//First count of number of lines in the file
	int count = 0;
	string linect;
	ifstream filect(filename.c_str());
	while (getline(filect, linect))
	{
		count++;
	}

	ifstream file(filename.c_str());		// construct an ifstream and open the stress file
	string   line;
	string 	partial;
	istringstream iss;

	vector<Vector3d> wave(count);
    

	// Reading file lines
	for (int i = 0; i < count; i++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		wave[i](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		wave[i](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		wave[i](2) = atof(partial.c_str());
		iss.clear();

		std::cout << wave[i](2) << std::endl;
	}
	
	return wave;
}

vector<Vector3d> readFloeVelFile (string filename)
{
	
	//First count of number of lines in the file
	int count = 0;
	string linect;
	ifstream filect(filename.c_str());
	while (getline(filect, linect))
	{
		count++;
	}

	ifstream file(filename.c_str());		// construct an ifstream and open the stress file
	string   line;
	string 	partial;
	istringstream iss;

	vector<Vector3d> vel(count);
    

	std::cout << "Importing floe velocities!!" << std::endl;
	// Reading file lines
	for (int i = 0; i < count; i++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		vel[i](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		vel[i](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		getline(iss, partial, ' ');
		getline(iss, partial, ' ');
		vel[i](2) = atof(partial.c_str());
		iss.clear();

		std::cout << vel[i](0) << " " << vel[i](1) << " " << vel[i](2) << std::endl;
	}
	
	return vel;
}

vector<Vector2d> arbSplitterFile (string filename, size_t & xdom)
{
	
	double finalx = 0;
	double finaly = 0;
	
	//First count of number of lines in the file
	int count = 0;
	string linect;
	ifstream filect(filename.c_str());
	while (getline(filect, linect))
	{
		count++;
	}

	ifstream file(filename.c_str());		// construct an ifstream and open the point list file
	string   line;
	string 	partial;
	istringstream iss;

	vector<Vector2d> pointList(count-1);
    

	// Reading file lines
	for (int i = 0; i < count; i++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		if (i < count - 1){
    		pointList[i](0) = atof(partial.c_str());
    		getline(iss, partial, ' ');
    		pointList[i](1) = atof(partial.c_str());
    		iss.clear();
		}
		else{
		    finalx = atof(partial.c_str());
    		getline(iss, partial, ' ');
    		finaly = atof(partial.c_str());
    		iss.clear();
		}
	}
	
	if (finalx + finaly >= 0)
	{
	    xdom = 1;  //X direction of line, horizontal fracture
	}
	else
	{
	    xdom = 0; //Y direction of line, vertical fracture    
	}

	return pointList;
}

//Realistic Ocean Data

//Grid Vector for Square Grid in x and y
vector<double> readRealGridFile(string filename)
{
	
	//First count of number of lines in the file
	int count = 0;
	string linect;
	ifstream filect(filename.c_str());
	while (getline(filect, linect))
	{
		count++;
	}

	ifstream file(filename.c_str());		// construct an ifstream and open the stress file
	string   line;
	string 	partial;
	istringstream iss;

	vector<double> gridxy(count);
    
    //std::cout<<" Importing Grid with length: "<< count << std::endl;
	
	// Reading file lines
	for (int i = 0; i < count; i++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		gridxy[i] = atof(partial.c_str());
		iss.clear();

		//std::cout << gridxy[i] << std::endl;
	}
	
	return gridxy;
}

//Matrix with same and x and y dimensions transformed into a vector of doubles for convenience
vector<double> readRealValueFile (string filename)
{
	
	//First count of number of lines in the file
	int count = 0;
	string linect;
	ifstream filect(filename.c_str());
	while (getline(filect, linect))
	{
		count++;
	}

	ifstream file(filename.c_str());		// construct an ifstream and open the stress file
	string   line;
	string 	partial;
	istringstream iss;

    //Assuming square grid and data received as Ny lines with Nx columns for each line and Ny = Nx = count
	vector<double> gridV(count*count);
    

	// Reading file lines
	size_t ct = 0;
	//Read "count" lines
	for (int i = 0; i < count; i++) {
		//std::cout << "Reading line number: " << i << std::endl;
		getline(file, line);
		iss.str(line);
		
		//Read "count" columns in each line
        for (int j = 0; j < count; j++){
		    getline(iss, partial, ' ');
		    gridV[ct] = atof(partial.c_str());
		    //std::cout << gridV[ct] << std::endl; 
		    ct++;
        }
		iss.clear();
	}
	
	return gridV;
}



#endif // READINPUTFILE_H_ 
