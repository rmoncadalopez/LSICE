/*
 *      main_Periodic_Ocean.cpp
 *      For 2D Modeling of Sea-Ice
 *      A group of grains floating in the ocean and subject to wind currents are either under Periodic B.C. or within walls
 *
 *      Authors: Konstantinos Karapiperis - Liuchi Li /  Modified by: John Harmon - Rigoberto Moncada
 *      FLuid & Thermal Additions: Rigoberto Moncada
 *      Fracture Components: Based on John Harmon and Adapted by Rigoberto Moncada
 */

#include "definitions.h"
#include "Levelset2d.h"
#include "Grain2d.h"
#include "readInputFile.h"
#include "World2d.h"
#include "Fracture.h"

//Main units for the program
/*
 Mass = kg
 Length =  km  #Start with grain Points and LS and units as km
 Time = seconds or hours  1 month = 2592000 s, 1 day = 86400 s
 
 Density = [Mass]/[L]**3 --> kg/km^3 1 kg/m3 = 1e9 kg / km3
 Volume = km^3
 Area = km^2
 Moment of Inertia = kg * m2    1 kg * m2 = 1e-6 kg * km2 3D
 Thermal Diffusivity = W / m*K   1 kg * m/ (K * s3) = 1e-3 kg*km/(K*s^3)
 1 Watt = 1kg * m2/s3.    or  1e-6 kg * km2 / s3
 1 W / m2 = 1 kg/s3 
 Stress = Pa = N / m2 = kg * m*s^-2 * m^-2 = 1 kg /(m*s2) = 1e3 kg / km*s2 = 1 Pa, 1 kPa = 1e6 kg/ (km*s2), 1 Gpa = 1e12 etc.
 */

//Modification of constants
/*
 
 */


//Find area Using Points (closed Polygon)
double PointsAreaM (vector<Vector2d> VecOrigin) 
{
    double Area;
    size_t n = VecOrigin.size();

    for (size_t i = 0; i < n-1; i++)
    {
        Area = Area + ( VecOrigin[i](0) * VecOrigin[i+1](1) -  VecOrigin[i](1) * VecOrigin[i+1](0) ); 
    }
    Area = Area + (VecOrigin[n-1](0) * VecOrigin[0](1) -  VecOrigin[n-1](1) * VecOrigin[0](0) ); 

    double result = 0.5*abs(Area);
    return result;
}

//Find if Ocean is exposed to ice or atmosphere (improve for PBC)
bool under_ice(const Grain2d & grainS, Vector2d & ptDetect)
{
    bool result_cover = false; //We assume no ice cover by default and prove the opposite.
    double      penetration;    // penetration amount (is negative by convention of the level set)
    Vector2d        normal;             // surface normal
    Vector2d        ptOtherCM;      // point wrt the center of mass of other in real space
    Vector2d        ptOtherLset;    // point in the reference config of other's level set

    //Start with a simple BBox radius criterion for simplicity (improve for PBC)
    if (  (ptDetect-grainS.getPosition()).norm() <= grainS.getRadius()  )  //No need for PBC for initialization, but needs for shifting grains
    {
        
        //You need to shift ptDetect to the right coordinate system before checking penetration
        const double cosdet = cos(grainS.getTheta());
        const double sindet = sin(grainS.getTheta());
        ptOtherCM = ptDetect - grainS.getPosition();
        ptOtherLset(0) =  ptOtherCM(0)*cosdet + ptOtherCM(1)*sindet;
        ptOtherLset(1) = -ptOtherCM(0)*sindet + ptOtherCM(1)*cosdet;
        ptOtherLset += grainS.getCmLset();


        //If within BBox radius now check if grid point penetrates ice level set
        if ( grainS.getLset().isPenetration(ptOtherLset, penetration, normal) )  
        {
            result_cover = true;
        }
    }

    return result_cover;
}

//Round Ocean Temperature for Output
double round_Ocean(Vector2d & pointxy, vector<double> & oceanTemp, const size_t & x_cells, const size_t & y_cells, const Vector2d & offset)
{
	size_t idx1, idy1;  // j are x cols, i are y rows
	int cell_sizex = int(offset(0))/int(x_cells-1);
	int cell_sizey = int(offset(1))/int(y_cells-1);
	//cout << "cell size x: " << cell_sizex << " cell_sizey: " << cell_sizey << endl;

	double tol_pos = 2.00;
	if ( isnan(pointxy(0)) ||  isnan(pointxy(1)) || pointxy(0) > offset(0) + tol_pos  || pointxy(0) < 0.0 - tol_pos || pointxy(1) > offset(1) + tol_pos || pointxy(1) < 0.0 - tol_pos )
	{
	  cout <<"WARNING NaN or BAD POSITION VALUES PROVIDED FOR INTERPOLATION, OUTPUT = 0" << endl;
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


//MAIN FUNCTION
int main(int argc, char * argv[]) {
    
    // MPI initialization
    int numprocessors, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    char name[50]; int count;
    MPI_Get_processor_name(name, &count);
    
    //Define Major Criteria
    bool create_folders = true; //true use Case folder, false use default

    double yearImport = atof(argv[4]);
    int yearImport_i = yearImport; //Make into integer
    
    //Choose Case (different case, different run, maybe specified from outside)
    // To read caseImport program must be called as Name 10 10 caseImport
    double caseImport = atof(argv[3]); //164697 30  //164757 40  //164577 10
    size_t caseNo = caseImport;
    //size_t caseNo = 4; //164637 20 same
    
    size_t year_case; //1 is 2020, 2 or else is 2018
    if (yearImport_i == 2020)
    {
        year_case = 1;
    }
    else
    {
        year_case = 2;
    }
    
    //Provide SIM PARAMETERS
    size_t melt_step, break_step;
    double IlimitMaxDiam, IlimitMinDiam, Vert_q, Q_atmos, mixed_layerh_in, kappa_in;

    //Import Sim Parameters Text File
    char tempfname00[300];
    
    //Input folder Location Define
    string inputFolderLoc;
    string inputGFolderLoc;
    
    size_t fluid_mode = 1;  //1 is super simple for concentration analysis, //2 is more arbitrary using Uwg // Anything else is not covered
    
    if (year_case == 1)
    {
        inputFolderLoc = "./Input/SeaIceWGJ2020/"; //2020
        inputGFolderLoc = "./Input/grainsIceWGJ/"; //2020
    }
    else
    {
        if (caseNo == 1033){
            inputFolderLoc = "./Input/SeaIceWGJ2018_MoreUniform1/"; //2018 1033 1.1 alpha
        }
        else if (caseNo == 1034){
            inputFolderLoc = "./Input/SeaIceWGJ2018_MoreUniform2/"; //2018 1034
        }
        else if (caseNo == 1035){
            inputFolderLoc = "./Input/SeaIceWGJ2018_LessUniform1/"; //2018 1035
        }
        else if (caseNo == 1036){
            inputFolderLoc = "./Input/SeaIceWGJ2018_LessUniform2/"; //2018 1036
        }
        else if (caseNo == 1077){
            inputFolderLoc = "./Input/SeaIceWGJ2018_Fracture/"; //2018 1036
            fluid_mode = 2;  //For more complex movement using an Uwg grid
        }
        else{
            inputFolderLoc = "./Input/SeaIceWGJ2018/"; //2018
        }
        inputGFolderLoc = "./Input/grainsIceWGJ2018/"; //2020
    }
    sprintf(tempfname00, (inputFolderLoc + "SimInputParamsv2.dat").c_str());  //DIR  //use 2 #2018 make sure right folder has this file //PROB CHANGE OFF
    //sprintf(tempfname00, (inputFolderLoc + "SimInputParamsv2_prob2.dat").c_str());  //DIR  //use 2 #2018 make sure right folder has this file //PROB CHANGE ON
    //sprintf(tempfname00, (inputFolderLoc + "SimInputNCasesZoom.dat").c_str());  //DIR  //use 2 #2020
    //sprintf(tempfname00, "./Input/SeaIceWGJ2018/SimInputNCasesZoom.dat");  //DIR  //use 2 #2018
    
    //sprintf(tempfname00, "./Input/SeaIceOcean/SimInput.dat");  //DIR
    std::string file_inputs = tempfname00;
    
    //Run Input Function
    //readSimInputFile(file_inputs, caseNo, melt_step, break_step, IlimitMaxDiam, IlimitMinDiam, Vert_q, Q_atmos);   
    readSimInputFileMore(file_inputs, caseNo, melt_step, break_step, IlimitMaxDiam, IlimitMinDiam, Vert_q, Q_atmos, mixed_layerh_in, kappa_in);

    const size_t START_TEMP = melt_step; //Can also be used as Step if melt Step is Constant
    double limitMaxDiam = IlimitMaxDiam;
    const double THERM_COEFF = 0.005; //MeltV could also be a constant //0.200 //0.005 for 100 steps //High speed makes unstable grains
    const double MELT_MASS = 0.06; //0.02 eg for mass //Divided Break and Melt //10000 for point size  //SHOULD IT BE LARGER FOR FASTER SLOPE????
    double limitMinDiam = IlimitMinDiam;
    const size_t BREAK_STEP = break_step;
    int BREAK_PROB = break_step; // Leave constant for now  //PROB CHANGE STAY REF
    
    //December 20, 2022 MODIF (Temporal)
    double VERTQ = Vert_q; //OFF To use qv
    //double VERTQ = 25.0; //ON To use kappa
    if (caseNo == 1077){
        VERTQ = 2.0;
    }
    
    double QATM = Q_atmos;
    
    //Import initial concentration
    char tempfnameconc[300];
    sprintf(tempfnameconc, (inputFolderLoc+"Initial_Concentration.dat").c_str() );  //DIR  //use 2
    //sprintf(tempfnameconc, "./Input/SeaIceWGJ2018/Initial_Concentration.dat");  //DIR  //use 2
    std::string file_conc = tempfnameconc;
    double initial_fine, initial_coarse;
    readConcFile(file_conc, initial_fine, initial_coarse);
    //Imported as %, turn to decimal
    cout << "Initial fine: " << initial_fine << " Initial Coarse: " << initial_coarse << endl;
    initial_fine /= 100.0;
    initial_coarse /= 100.0;
    
    // Get morphology, init. position and init. velocity input files
    char tempfname[300];

    //Grains    
    sprintf(tempfname, (inputFolderLoc + "morphologies.dat").c_str() );
    //sprintf(tempfname, "./Input/SeaIceWGJ2018/morphologies.dat");
    //sprintf(tempfname, "./Input/SeaIceOcean/morphologies.dat");
    //sprintf(tempfname, "./Input/mainCompression3/morphologies.dat");
    string file_morph = tempfname;
    sprintf(tempfname, (inputFolderLoc + "initialpositions.dat").c_str() );
    //sprintf(tempfname, "./Input/SeaIceWGJ2018/initialpositions.dat");
    //sprintf(tempfname, "./Input/SeaIceOcean/initialpositions.dat");
    //sprintf(tempfname, "./Input/mainCompression3/initialpositions.dat");
    string file_pos = tempfname;
    sprintf(tempfname, (inputFolderLoc + "initialvelocities.dat").c_str() );
    //sprintf(tempfname, "./Input/SeaIceWGJ2018/initialvelocities.dat");
    //sprintf(tempfname, "./Input/SeaIceOcean/initialvelocities.dat");
    //sprintf(tempfname, "./Input/mainCompression3/initialvelocities.dat");
    string file_vel = tempfname;
    sprintf(tempfname, (inputGFolderLoc).c_str());
    //sprintf(tempfname, "./Input/grainsIceWGJ2018/");
    //sprintf(tempfname, "./Input/grainsAug/");
    //sprintf(tempfname, "./Input/grains/");
    //sprintf(tempfname, "./Input/grainsComp2/");
    string morph_dir = tempfname;

    //Get temperature-related quantitities
    sprintf(tempfname, (inputFolderLoc + "initialtemper.dat").c_str() );
    //sprintf(tempfname, "./Input/SeaIceWGJ2018/initialtemper.dat");
    //sprintf(tempfname, "./Input/SeaIceOcean/initialtemper.dat");
    //sprintf(tempfname, "./Input/mainCompression3/initialtemper.dat");
    string file_temper = tempfname;
    sprintf(tempfname, (inputFolderLoc + "initialthick.dat").c_str() );
    //sprintf(tempfname, "./Input/SeaIceWGJ2018/initialthick.dat");
    //sprintf(tempfname, "./Input/SeaIceOcean/initialthick.dat");
    //sprintf(tempfname, "./Input/mainCompression3/initialthick.dat");
    string file_thick = tempfname;
    
    //Import environmental data
    sprintf(tempfname, (inputFolderLoc + "Temp_data.dat").c_str() );
    //sprintf(tempfname, "./Input/SeaIceWGJ2018/Temp_data.dat");
    string file_out_temper = tempfname;
    sprintf(tempfname, (inputFolderLoc + "WaveH_data.dat").c_str() );
    //sprintf(tempfname, "./Input/SeaIceWGJ2018/WaveH_data.dat");
    string file_wave = tempfname;
    sprintf(tempfname, (inputFolderLoc + "WaveH_data_Day.dat").c_str() );
    //sprintf(tempfname, "./Input/SeaIceWGJ2018/WaveH_data_Day.dat");
    string file_waveD = tempfname;
    sprintf(tempfname, (inputFolderLoc + "WaveLength_data.dat").c_str() );
    //sprintf(tempfname, "./Input/SeaIceWGJ2018/WaveLength_data.dat");
    string file_waveL = tempfname;
    sprintf(tempfname, (inputFolderLoc + "WaveLength_data_Day.dat").c_str() );
    //sprintf(tempfname, "./Input/SeaIceWGJ2018/WaveLength_data_Day.dat");
    string file_waveLD = tempfname;
    sprintf(tempfname, (inputFolderLoc + "Vel_Change.dat").c_str() );
    string file_floevel = tempfname;

    vector<Vector2d> Out_temperature = readOutTempFile(file_out_temper);
    vector<Vector2d> Wave_HD = readOutTempFile(file_waveD);  //For a two-column file
    vector<Vector3d> WaveH = readWaveFile(file_wave); //For a three-column file
    vector<Vector2d> Wave_LD = readOutTempFile(file_waveLD);
    vector<Vector3d> WaveL = readWaveFile(file_waveL);
    vector<Vector3d> floeVel = readFloeVelFile(file_floevel);
    
    //Calculate average floe speed for scaling
    double ave_floeVel = 0.0;
    for (size_t ifl = 0; ifl < floeVel.size(); ifl++) {
	ave_floeVel += floeVel[ifl](1);
    }	
    ave_floeVel /= floeVel.size();	
	

    // //Read Grains Walls Files
    // sprintf(tempfname, "./Input/SeaIce/morphologiesWallF.dat");
    // string file_morphW = tempfname;
    // sprintf(tempfname, "./Input/SeaIce/initialpositionsWallFrev.dat");
    // string file_posW = tempfname;
    // sprintf(tempfname, "./Input/SeaIce/initialvelocitiesW.dat");
    // string file_velW = tempfname;
    // sprintf(tempfname, "./Input/grainsWallsF/");
    // string morph_dirW = tempfname;

    // //Get temperature-related quantitities
    // sprintf(tempfname, "./Input/SeaIce/initialtemperW.dat");
    // string file_temperW = tempfname;
    // sprintf(tempfname, "./Input/SeaIce/initialthickW.dat");
    // string file_thickW = tempfname;
    

    
    //Read Initial Wall Files  Are they even necessary????
    sprintf(tempfname, "./Input/SeaIce/initialpositionsTopWall.dat");
    string file_wallT = tempfname;
    sprintf(tempfname, "./Input/SeaIce/initialpositionsBottomWall.dat");
    string file_wallB = tempfname;
    sprintf(tempfname, "./Input/SeaIce/initialpositionsLeftWall.dat");
    string file_wallL = tempfname;
    sprintf(tempfname, "./Input/SeaIce/initialpositionsRightWall.dat");
    string file_wallR = tempfname;
    


    // Generate grains
//    vector<Grain2d> grains0 = generateGrainsFromFiles(file_morph, morph_dir, file_pos, file_vel);
//    vector<Grain2d> grains;
//    grains.push_back(grains0[0]);
    vector<Grain2d> grains = generateGrainsFromFiles(file_morph, morph_dir, file_pos, file_vel, file_temper, file_thick);
    size_t ngrains = grains.size();

    
    vector<Grain2d> grainsWall;
    cout << "No grains wall so size is: " << grainsWall.size() << endl;
    // vector<Grain2d> grainsWall = generateGrainsFromFiles(file_morphW, morph_dirW, file_posW, file_velW, file_temperW, file_thickW);
    // size_t ngrainsWall = grainsWall.size();

    // Change grain properties (for ICE)
    //double rho = 0.91e-6;                      // g/pixel^3 Density of Ice (Herman, 2011)
    
    //If pixels are assumed to be km2
    double rho = 910e9;  // Density of Ice = 0.91 kg/m3 //Change Aug 22, 2022 rho = 0.91e9 iff 1 m thick adjust but thickness varies, so changed to rho = 910e9 kg / km^3
    
    double mu = 0.8;  //Internal Friction Coefficient of Ice (Weiss, 2013)
    //double kn = 8e2; //2e2 //6e3  //8e2                         // g/s^2 1e4 stabilizes here  (Original 1e6) (Young Modulus E=6e9 Pa)
    double kn = 6e12;   //(Young Modulus E=6e9 Pa to km)
    
    double ks = 0.8*kn;                           // g/s^2


    double PrConv   = 1e-6;  // MPa per LSDEM pressure

    // grain properties
//  double yield = 80/PrConv;
    double yield = .4/PrConv; // in LSDEM pressure units (MPa / (MPa/LSDEM))
    double d0 = 28.6;
    double m=3;


    // Apply grain properties indicated above
    double randmax = double(RAND_MAX);
    for (size_t i = 0; i < ngrains; i++) {
        grains[i].changeDensity(rho);
        grains[i].changeMu(mu);
        //grains[i].changeKn(kn);
        grains[i].changeKn(ks);
        
        //Aug-18-2022 Change Mass Scaled by thickness
        //BEGIN
        grains[i].changeKn(kn * 1.0); //Change Aug 22, 2022 //Adjust for force changes due to mass changes, divided by 1000 (thickness) but multiplied by 1000 (change in rho) See what happens with 1.0 force field factor. //PHASE II //Leave as this for now. See effect of current
        
        double mass_temp = grains[i].getMass() * grains[i].getThickness0() * 0.001; //Scale mass area * density ( (area * thick = 1 = Volume) * density * proportional to thickness) (Thus, area is to adjust thickness mass, final term is volume mass now not area mass) (If thickness changes, rescale _mass, if area changes, REDO)
        grains[i].changeMass(mass_temp);
        //END


        double Pf = 0.6;

        double Sig_c = yield*pow((-1.*pow((d0/2*grains[i].getRadius()),3)*log(1.-Pf)),1./m); //2*R =Diameter
        grains[i].changeYield(1e-9*Sig_c);


       //Apply a random velocity to get grains moving
       // grains[i].changeVelocity(50.*Vector2d(rand()/randmax,rand()/randmax));  //Random change of grains as to model random external forcing (must be separated for wind and water) Originally 10000
       // grains[i].changeOmega(10*rand()/randmax);  //Function for randomly setting angular velocity (if needed)
    }

    // size_t ID_Diff = 100000; //For ID differentiation grains from grainWalls, depends on size of grains

    // for (size_t i = 0; i < ngrainsWall; i++) {
    //     grainsWall[i].changeDensity(rho);
    //     grainsWall[i].changeMu(mu);
    //     grainsWall[i].changeKn(kn);
    //     grainsWall[i].changeKn(ks);
    //     size_t Id_temp = grainsWall[i].getId();
    //     grainsWall[i].changeId(Id_temp+ID_Diff);
    // }



    string variant("0");


    // frac properties
    double RoMratio = 100;
    double maxRoM   = 100;
    size_t minpts  = 100;
    FracProps2D FracProps = FracProps2D(variant,RoMratio,maxRoM,minpts);


    //Update properties such as Mass to be consistent
    double new_mass;
    double new_I;
    Vector2d new_cm;
    Levelset2d ref_lset;
    for (size_t i = 0; i < ngrains; i++) {
        //ref_lset = grains[i].getLset();
        //FracProps.findGrainProps0(ref_lset, new_mass, new_cm, new_I);
        //grains[i].changeInertia(new_I);
        //grains[i].changeMass(new_mass);
        //grains[i].changeCM(new_cm);

    }
    
    


    //Start Wall Modification
    //--------------------------------------------------------------------------------------------------------------------------------------------------
    
    double height = 8000;  //2645.3
    // double width = 8000;
    
    
    // // create walls   CHANGE WALL POSITIONS???????
    vector<Wall2d> walls; // 0bot 1left 2right
    walls.resize(2);
    // // bottom wall
    walls[0] = Wall2d(Vector2d(-80000,-80000), Vector2d(0,1), kn, 0.9*kn, 0.5, INT_MAX-3);
    
    // double rot_angle = 70.0 * (3.141592653589793/180.00);
    // double rot_angle2 = -70.0 * (3.141592653589793/180.00);
    // double v11 = 0 * cos(rot_angle) - 1 * sin(rot_angle);
    // double v12 = 0 * sin(rot_angle) + 1 * cos(rot_angle);
    // double v21 = 0 * cos(rot_angle2) - 1 * sin(rot_angle2);
    // double v22 = 0 * sin(rot_angle2) + 1 * cos(rot_angle2);
    // //cout << v11 << endl;
    // //cout << v12 << endl;
    // //cout << v21 << endl;
    // //cout << v22 << endl;
    // //walls[0] = Wall2d(Vector2d(250.00,1150.00), Vector2d(v21,v22), kn, 0.9*kn, 0.5, INT_MAX-3);
    // //walls[1] = Wall2d(Vector2d(900.00,1150.00), Vector2d(v11,v12), kn, 0.9*kn, 0.5, INT_MAX-2);
    
    // // top wall
    walls[1] = Wall2d(Vector2d(-400000,height), Vector2d(0,-1), kn, 0.9*kn, 0.5, INT_MAX-2);
    // // left wall
    // //walls[2] = Wall2d(Vector2d(-2000,-2000), Vector2d(1,0), kn, 0.9*kn, 0.5, INT_MAX-1);
    // // right wall
    // //walls[3] = Wall2d(Vector2d(width,-2000), Vector2d(-1,0), kn, 0.9*kn, 0.5, INT_MAX-0);
    
    //--------------------------------------------------------------------------------------------------------------------------------------------------
    //End Wall Modification
    
    
    
    // Set-up I/O
    char posRotFile[300];
    char velocityFile[300];
    char forceFile[300];
    char momentFile[300];

    
    char temperFile[300];
    char thickFile[300];
    char sampleLSFile[300];
    char pointsoutFile[300];
    char pointsoutLSFile[300];
    char normAreaFile[300];
    char lossAreaFile[300];
    char normConcFile[300];
    char fluxAreaFile[300];
    char numbergoutFile[300];
    char DiametersoutFile[300];
    char npointspergrainFile[300];
    char oceanTGridFile[300];
    char oceanTGridDimFile[300];
    char oceanTGridCoordFile[300];
    char testParamsFile[300];
    char posFineFile[300];
    char numberFineFile[300];
    
    //Aug-18-2022 Change
    //BEGIN (reference down later)
    char massLossFile[300];
    //END

    char wallTFile[300];  //But Walls are set to be static why update their positions???
    char wallBFile[300];
    char wallLFile[300];
    char wallRFile[300];
        
    string folderLoc, folderLoc2, folderLoc3, folderLoc4, folderLoc5, year_subfolder;
    
    if (year_case == 1)
    {
        //year_subfolder = "SeaIce_";
        //year_subfolder = "SeaIce2020N_";
        year_subfolder = "SeaIce2020_";
    }
    else
    {
        //year_subfolder = "SeaIce_";  //Actually for 2018
        year_subfolder = "SeaIce2018_"; //PROB CHANGE OFF
        //year_subfolder = "SeaIce2018Prob_"; //PROB CHANGE ON
        //year_subfolder = "SeaIce2018h_";
    }
    
    //Seed for multiple cases!!!
    bool seed_on = true;


    //Extension for ensemble runs //Seed Remove or add 
    double seedImport; size_t seedNo; string seedString;
    if (seed_on == true){
        seedImport = atof(argv[5]);
        seedNo = seedImport;
        seedString = (std::to_string(seedNo));  
    }
    
    if (create_folders)
    {
        //Create folders for each respective case and subfolders as well
        cout << "Creating New Folders" << endl;
        
        string caseString = (std::to_string(caseNo));
        
        //Seed
        if (seed_on == true){
            folderLoc = "./Output/"+year_subfolder+caseString+"_"+seedString+"/";
            folderLoc2 = "./Output/"+year_subfolder+caseString+"_"+seedString+"/ContactInfo/";
            folderLoc3 = "./Output/"+year_subfolder+caseString+"_"+seedString+"/GSD/";
            folderLoc4 = "./Output/"+year_subfolder+caseString+"_"+seedString+"/Damage/";
            folderLoc5 = "./Output/"+year_subfolder+caseString+"_"+seedString+"/Vertical_Thickness/";
        }
        else{
        //No seed
            folderLoc = "./Output/"+year_subfolder+caseString+"/";
            folderLoc2 = "./Output/"+year_subfolder+caseString+"/ContactInfo/";
            folderLoc3 = "./Output/"+year_subfolder+caseString+"/GSD/";
            folderLoc4 = "./Output/"+year_subfolder+caseString+"/Damage/";
            folderLoc5 = "./Output/"+year_subfolder+caseString+"/Vertical_Thickness/";
        }

        //Turn off when necessary
        mkdir(folderLoc.c_str(), 0700);
        mkdir(folderLoc2.c_str(), 0700);
        mkdir(folderLoc3.c_str(), 0700);
        mkdir(folderLoc4.c_str(), 0700);
        mkdir(folderLoc5.c_str(), 0700);
    }
    else
    {
        //Use to not create new folders
        folderLoc =  "./Output/SeaIce/";  //Default for controlling output
        folderLoc2 = "./Output/SeaIce/ContactInfo/";
        folderLoc3 = "./Output/SeaIce/GSD/";
        folderLoc4 = "./Output/SeaIce/Damage/";
        folderLoc5 = "./Output/SeaIce/Vertical_Thickness/";
    }
    

    
    sprintf(posRotFile,(folderLoc+"positions"+variant+".dat").c_str());
    sprintf(velocityFile,(folderLoc+"velocities"+variant+".dat").c_str());
    sprintf(forceFile,(folderLoc+"forces"+variant+".dat").c_str());
    sprintf(momentFile,(folderLoc+"moments"+variant+".dat").c_str());


    sprintf(temperFile,(folderLoc+"temperatures"+variant+".dat").c_str());
    sprintf(thickFile,(folderLoc+"thickness"+variant+".dat").c_str());

    sprintf(sampleLSFile,(folderLoc+"samplels"+variant+".dat").c_str());
    sprintf(pointsoutFile,(folderLoc+"pointsout"+variant+".dat").c_str());
    sprintf(pointsoutLSFile,(folderLoc+"pointsoutLS"+variant+".dat").c_str());
    sprintf(normAreaFile,(folderLoc+"normarea"+variant+".dat").c_str());
    sprintf(lossAreaFile,(folderLoc+"loss_area"+variant+".dat").c_str());
    sprintf(normConcFile,(folderLoc+"normConc"+variant+".dat").c_str());
    sprintf(fluxAreaFile,(folderLoc+"fluxarea"+variant+".dat").c_str());
    sprintf(numbergoutFile,(folderLoc+"numberg"+variant+".dat").c_str());
    sprintf(DiametersoutFile,(folderLoc+"Diameters"+variant+".dat").c_str());
    sprintf(npointspergrainFile,(folderLoc+"npointsperg"+variant+".dat").c_str());
    sprintf(oceanTGridFile,(folderLoc+"oceanTGrid"+variant+".dat").c_str());
    sprintf(oceanTGridDimFile,(folderLoc+"oceanTGridDim"+variant+".dat").c_str());
    sprintf(oceanTGridCoordFile,(folderLoc+"oceanTGridCoord"+variant+".dat").c_str());
    sprintf(testParamsFile,(folderLoc+"testParams"+variant+".dat").c_str());
    sprintf(posFineFile,(folderLoc+"posFines"+variant+".dat").c_str()); 
    sprintf(numberFineFile,(folderLoc+"numberFines"+variant+".dat").c_str());
    
    //Aug-18-2022 Change
    //BEGIN (reference down later)
    sprintf(massLossFile, (folderLoc+"massLoss"+variant+".dat").c_str());
    //END

    //Export for diagnostics (Points Grain 0 and 1)
    
    sprintf(wallTFile,(folderLoc+"wallpositionsTop"+variant+".dat").c_str());
    sprintf(wallBFile,(folderLoc+"wallpositionsBottom"+variant+".dat").c_str());
    sprintf(wallLFile,(folderLoc+"wallpositionsLeft"+variant+".dat").c_str());
    sprintf(wallRFile,(folderLoc+"wallpositionsRight"+variant+".dat").c_str());

    
    
    FILE * positions           = fopen(posRotFile,"w");
    FILE * velocities          = fopen(velocityFile,"w");
    FILE * forces              = fopen(forceFile,"w");
    FILE * moments             = fopen(momentFile,"w");


    FILE * temper          = fopen(temperFile,"w");
    FILE * thick          = fopen(thickFile,"w");

    FILE * sampleLS       = fopen(sampleLSFile,"w");
    FILE * pointsout       = fopen(pointsoutFile,"w");
    FILE * pointsoutLS       = fopen(pointsoutLSFile,"w");
    FILE * normArea       = fopen(normAreaFile,"w");
    FILE * lossArea       = fopen(lossAreaFile,"w");
    FILE * normConc     = fopen(normConcFile,"w");
    FILE * fluxArea       = fopen(fluxAreaFile,"w");
    FILE * numbergout       = fopen(numbergoutFile,"w");
    FILE * Diametersout    = fopen(DiametersoutFile,"w");
    FILE * nperg           = fopen(npointspergrainFile,"w");
    FILE * oceanTGrid       = fopen(oceanTGridFile,"w");
    FILE * oceanTGridDim       = fopen(oceanTGridDimFile,"w");
    FILE * oceanTGridCoord       = fopen(oceanTGridCoordFile,"w");
    FILE * testParams            = fopen(testParamsFile,"w");
    FILE * finepositions          = fopen(posFineFile,"w");
    FILE * numberfines           = fopen(numberFineFile,"w");
    
    //Aug-18-2022 Change
    //BEGIN (reference down later)
    FILE * massLoss      =   fopen(massLossFile,"w");
    //END

    FILE * wallPosT             = fopen(wallTFile,"w");
    FILE * wallPosB             = fopen(wallBFile,"w");
    FILE * wallPosL             = fopen(wallLFile,"w");
    FILE * wallPosR             = fopen(wallRFile,"w");
    
    //For contact output
    string outDir2 = folderLoc2;
    string outDir3 = folderLoc3;
    
    //For thickness output
    string outDir4 = folderLoc5;
    
    // Global parameters
    double gDamping              = 2;       // global damping (Originally 0.2)
    double dt                    = 0.000007;         // time step (Originally 5e-6) //TODO: ADJUST FOR CORRECT TIMEFRAME, right now it is for contact mechanics, heat effects are stepped at 1 second

    //Should I replace it with   dt = 0.2*(2*sqrt(.5*200/kn));
//    double dt = 0.2*(2*sqrt(.5*200/kn));   0.004
//    double dt                    = 0.00005;
    
    unsigned int nT;
    if (year_case == 1)
    {
        nT = 2160000; //nT = 2073600; // 24 days are 2073600 s //300000 case 4  //1000000;//4000000; // 4000000  //    // total time steps compression (Originally 40000) // //TODO: ADJUST FOR CORRECT TIMEFRAME
        //!!!CHECK NEW WARNING!!!
    }
    else
    {
        nT = 4233600; //49 d  !!!CHECK NEW WARNING!!!
        //nT = 4147200; //48 d
        //nT = 414720*2.0;
        //nT = 414720*3.0;
        //nT = 4147200; //For 2018, 48 days
    }
    //If seconds use 1d  = 86400 s

    Vector2d offset;
    double offset_val; //= 4000; //**For larger domain
    if (year_case== 1)
    {
        offset_val = 380; //2020
    }
    else
    {
        offset_val = 400; //2018

    }
    
    offset << offset_val, offset_val;        //Watch out for Output extent, increase perhaps to 5000,2000 and sorround all over with walls ????????
                                 //If walls do not move in this case, should I export wallPos anyway and use world.getWalls ??????????

    //Create fluid grid and assign speeds

    //Fluid Field Grid
    size_t cell_sizex = offset_val/offset_val; //1 //Square cells //For performance size cant be too small.
    size_t cell_sizey = offset_val/offset_val; //1

    // size_t cell_sizex = offset_val/400; //1 //Square cells //For performance size cant be too small.
    // size_t cell_sizey = offset_val/400; //1
    
    size_t x_cells = size_t(offset(0))/(cell_sizex) + 1;
    size_t y_cells = size_t(offset(1))/(cell_sizey) + 1;


    //Initialize Fluid Grid
    vector<Vector2d> fluid_coord(x_cells * y_cells);
    for (size_t i = 0; i < y_cells; i++) {
        for (size_t j = 0; j < x_cells; j++) {
            fluid_coord[j+i*x_cells] <<  double(j)*double(cell_sizex),  double(i)*double(cell_sizey)  ;  //x , y coords
        }   
    }   

    //Initialize Water Velocity and Heat Grid
    vector<Vector2d> Uwg(x_cells * y_cells); //X and Y water velocities, U and V
    vector<Vector2d> Uwg0(x_cells * y_cells); //X and Y water velocities, U and V
    vector<double> oceanTemp(x_cells * y_cells); //Global Ocean Temp Grid
    vector<double> oceanTemp0(x_cells * y_cells);
    vector<double> alpha_ice_inter(x_cells * y_cells); //Proportion of Ice and Water (1 is only ice, 0 is only ocean, fraction in between)
    
    //Initialize flux vector for reference and other useful variables
    vector<double> netFlux(x_cells * y_cells);
    vector<double> netFlux2(x_cells * y_cells);
    vector<double> latFlux(x_cells * y_cells);
    vector<double> vertFlux(x_cells * y_cells);
    vector<double> lvRatio(x_cells * y_cells);
    double concf_free = initial_fine / (1-initial_coarse);
    double Tf = -1.8; //Melting temperature of sea ice
    double Qatm0 = QATM; //310 W/m2
    double Aatm0 = 70.0; // 70   W/m2 //55.6 W/m2 ??
    double Batm0 = 10.0; // 10   W/(m2*K)
    double a_o0 = 0.3; //Albedo of open ocean
    double a_i0 = 0.7; //Albedo of sea ice
    //Flux vector addition for calculation
    
    bool ice_contact; //Indicator for ice contact with ocean
    for (size_t i = 0; i < y_cells; i++) {
        for (size_t j = 0; j < x_cells; j++) {
            
            
            //Modify Ocean Temp Grid (Simple)
            // if ( fluid_coord[j+i*x_cells](1) <= 200.00)
            // {
            //     oceanTemp[j+i*x_cells] = 18.4;
            // }
            // else
            // {
            //     oceanTemp[j+i*x_cells] = -1.8;
            // }

            //Modify Ocean Temp using grain intersection (bbox and isPenetration)
            ice_contact = false;
            alpha_ice_inter[j+i*x_cells] = 0.0;
            for (size_t gi = 0; gi < ngrains; gi++) {
                //Start with a simple BBox radius criterion for simplicity

                ice_contact = under_ice(grains[gi], fluid_coord[j+i*x_cells]);

                if (ice_contact)
                {
                    alpha_ice_inter[j+i*x_cells] = 1.0;
                    break;
                }
            }

            //Initialize different temperature depending on ice cover or not
            if (ice_contact)
            {
                oceanTemp[j+i*x_cells] = -1.78;
            }
            else
            {
                if (year_case== 1)
                {
                    oceanTemp[j+i*x_cells] = -1.37; //Diff for 2020
                }
                else
                {
                    oceanTemp[j+i*x_cells] = -1.5;
                }
                
            }
                                                                        //Net flux coarse covered                                                                    //Net flux uncovered by coarse, mixing fines and ocean influence           
            netFlux[j+i*x_cells] = alpha_ice_inter[j+i*x_cells] * ( 1.0*( Qatm0*(1-a_i0) - (Aatm0 + Batm0*Tf) ) + 0.0 )  +  (1 - alpha_ice_inter[j+i*x_cells]) * ( concf_free*( Qatm0*(1-a_i0) - (Aatm0 + Batm0*Tf) ) + (1-concf_free) * ( Qatm0*(1-a_o0) - (Aatm0 + Batm0*oceanTemp[j+i*x_cells]) ) );
            netFlux2[j+i*x_cells] = alpha_ice_inter[j+i*x_cells] * ( 0.0 )  +  (1 - alpha_ice_inter[j+i*x_cells]) * ( concf_free*( Qatm0*(1-a_i0) - (Aatm0 + Batm0*Tf) ) + (1-concf_free) * ( Qatm0*(1-a_o0) - (Aatm0 + Batm0*oceanTemp[j+i*x_cells]) ) );

            //No edge cells for lateral flux of course (simplify doing PBC)
            if ( (i > 0 && i < y_cells - 1) && (j > 0 && j < x_cells - 1) ){
                //Track lateral flux
                latFlux[j+i*x_cells] = 1030.0*3991.0*mixed_layerh_in*1e-6*(kappa_in)* (oceanTemp[j+(i+1)*x_cells]+oceanTemp[(j+1)+i*x_cells]+oceanTemp[j+(i-1)*x_cells]+oceanTemp[(j-1)+i*x_cells]-4*oceanTemp[j+i*x_cells]);
                vertFlux[j+i*x_cells] = ( VERTQ*(oceanTemp[j+i*x_cells] - Tf) );
                if (vertFlux[j+i*x_cells] == 0){
                    lvRatio[j+i*x_cells] = 1111; //To show no vertical rate (very unlikely, avoid division by zero)
                }
                else{
                    lvRatio[j+i*x_cells] = latFlux[j+i*x_cells] /  vertFlux[j+i*x_cells];
                }
            }

            double vortex_vel = 1000.00 * 0.7; //1000 //FACTOR OF INITIAL March 1, 2023  //0.9-1.0 ok for stable thickness
            //VORTEX-like field in 4 symmetric parts around offset
            if (fluid_coord[j+i*x_cells](1)<offset(1)*0.5 && fluid_coord[j+i*x_cells](0)<offset(0)*0.5)
            {
               Uwg[j+i*x_cells] <<  0.0 ,  vortex_vel ;  //x , y velocities
            }
            else
            {    
                if (fluid_coord[j+i*x_cells](1)>offset(1)*0.5 && fluid_coord[j+i*x_cells](0)<offset(0)*0.5)
                {
                    Uwg[j+i*x_cells] << vortex_vel , 0.0 ;  //x , y velocities
                }
                else
                {  
                    if (fluid_coord[j+i*x_cells](1)>offset(1)*0.5 && fluid_coord[j+i*x_cells](0)>offset(0)*0.5)
                    {
                      Uwg[j+i*x_cells] << 0.0 , -vortex_vel ;  //x , y velocities
                    }
                    else
                    {
                      Uwg[j+i*x_cells] << -vortex_vel , 0.0;  //x , y velocities
                    }
                }
            }               
        }   
    }   
    //Save initial
    Uwg0 = Uwg;

    //Check grid
    for (size_t i = 0; i < y_cells; i++) {
        for (size_t j = 0; j < x_cells; j++) {
            //cout << "X: " << fluid_coord[j+i*x_cells](0) << " Y: " << fluid_coord[j+i*x_cells](1)  <<  " oceanTemp: " << oceanTemp[j+i*x_cells] << endl;
        }
    }
    //exit(1);


    // Create world

    size_t stepup = 0;  //-->

    double slopedir = 0.01; //For sinusoidal variation over time if needed. Just for initialization only
    double flowforceRef0 = 4.0e9;  //SCALE ????????? //Initial: 4.0e11
    //double flowforceRef0 = 2.7e12;   //FORMER: 4.0e11; //NEW TRY 27.0E12 from 86400 s * 5 day * contact_dt * F = 8 km/day * 112,000 (domain area) km2 * ConcC * 0.10 (frac of largest floe) * Ave thick in km * rho ice in kg/km3 so F = 27.0e12 approx. Try and see. 
    double flowforceRef = flowforceRef0 * floeVel[0](1)/ave_floeVel ; //Ref.
    double flowforcemin = 1.0e11; //Ref. for velocity of MODIS
    // double flowangle = 0.25*3.1415926535897932384; //pi/4
    double flowangle = floeVel[0](2); //From floe velocity 

    //START_TEMP = 2;
    //World2d world(grains, walls, offset, dt, gDamping, stepup, slopedir, flowangle, flowforceRef, FracProps, grainsWall, fluid_coord, Uwg, x_cells, y_cells, oceanTemp, START_TEMP, THERM_COEFF, MELT_MASS, BREAK_PROB);   //ADD    FracProps           //Use walls instead of offset ?????????
    World2d world(grains, walls, offset, dt, gDamping, stepup, slopedir, flowangle, flowforceRef, FracProps, grainsWall, fluid_coord, Uwg, x_cells, y_cells, oceanTemp, START_TEMP, THERM_COEFF, MELT_MASS, BREAK_PROB, cell_sizex, cell_sizey, fluid_mode);   //Enhanced for drag

    // Clear up memory from grains since world now has a copy
    grains.clear();
    grainsWall.clear();



    // Start clock
    double duration;
    double start = omp_get_wtime();

    // Create a new vector of grain objects for output purposes
    vector<Grain2d> grainList;
    GrainState2d grainList2;  //Created to display forces and moments for each grain

    for (size_t i = 0; i < grainList.size();i++)
    {
        //cout << "Properties Grain: " << i << " Radius: " << grainList[i].getRadius() << endl;
    }


    double Ke = 0.;

    size_t tiso = 2; //Control when can fracture and melting start  //Adjust in Grain2d.h temp adjustment
    
    //==============================================================================================
    //==============================================================================================

    // Start time integration (compression)
    int nsnaps;
    if (year_case == 1){
        nsnaps = 25; //24  // 100; //Switch to 30?? ##NEW
    }
    else{
        nsnaps = 49; //48 //2018 ##NEW
    }
        
    size_t nTout= nT/nsnaps;

    //Temperature Applications
    //size_t nBreak = BREAK_STEP; //or 2000 (3) //20000  //Last run used 4000, reduced to 4000 for frequency
    size_t nBreak = BREAK_PROB; //PROB. CHANGE instead of step STAY REF
    size_t nTemp = START_TEMP; //250 and 0.005 therm coeff decrease 0.9 conc in about 50,000 steps, you can use this step to update shape in interval, which is more expensive than updating temperatures
    
    //TEMPORAL CHANGES
    
    //nTemp = 1;
    //nBreak = 300; //300; //150 //450 //400 seemed to be fine
    
    size_t dstep = nTemp; //Start of melt and step are kept same
    //const size_t dstep = 500; //2000 //500 //Smaller slower, larger less stable
    //nTemp = dstep;
    
    // limitMaxDiam = 8.0;  //10 //15
    // limitMinDiam = 0.005;
    size_t nStress = 1000;
    //size_t nTout= nT/100;  //Total time steps are the resolution and 100 is the amount of snapshots or states captured by the output

    double AccArea0; //Initialization of original area of grains before melt and all
    double AccArea;
    double fluxAccArea = 0.0;


    //Initialize Damage
    //world.InitialDamageMesher();
    world.DamageInit();
    double load_fac0 = 0.0; //To get slope of sinusoidal loading
    
    //Define domain for studying concentration
    double Domain;
    if (year_case == 1)
    {
        Domain = 259.0 * 339.66; //2020
    }
    else
    {
        Domain = 0.076011266724044 * 1474000; //For 2018 Example, see imageprocess2_label_v2018.py (1340 x, 11000 y in pixels)
    }
    //double Domain = offset_val * offset_val; //80000, 80000 //80 , 80 //Larger areas
    
    //Print Parameters
    cout << "Initial Parameters" << endl;
    if (seed_on == false){
        cout << "Melt Step: " << START_TEMP << " LimitMaxD: " << limitMaxDiam << " LimitMinD: " << limitMinDiam << " Break Step: " << BREAK_STEP << " qVert: " << VERTQ << " Satm: " << QATM <<  " Hlayer: " << mixed_layerh_in <<  " Kappa: " << kappa_in << endl;  //No Seed
    }
    else{
        cout << "Melt Step: " << START_TEMP << " LimitMaxD: " << limitMaxDiam << " LimitMinD: " << limitMinDiam << " Break Step: " << BREAK_STEP << " qVert: " << VERTQ << " Satm: " << QATM <<  " Hlayer: " << mixed_layerh_in <<  " Kappa: " << kappa_in << " Seed: " << seedNo << endl;  //Seed
    }

    double GlobalConc = initial_fine + initial_coarse;
    double GlobalFine = initial_fine;
    
    //Define vector of Vector2d that contains: Conc(Proportional to Area) and Thickness
    size_t fine_bins = 1; //Real 100 TEMPo
    vector<Vector4d> Fine_Grains(fine_bins);
    vector<Vector2d> Fine_Velocities(fine_bins);
    //Use Initial fine to distribute each 1% of GlobalFine and assign a Gaussian Thickness around 1m
    std::default_random_engine generator;
    std::default_random_engine generator_pos;
    //std::normal_distribution<double> distribution(1.0,0.3); //1m mean, 0.3 std. dev
    double min_thick = 0.01; //Avoid using negative thickness given Gaussian distrib.
    
    double mean_thick;
    if (year_case== 1)
    {
        mean_thick = 0.3;
    }
    else
    {
        mean_thick = 0.5;
    }

    std::normal_distribution<double> distribution(mean_thick,0.1); //0.5m mean, 0.1 std. dev

    //For random fine position
    double min_pos = 0.01; //Avoid using negative thickness given Gaussian distrib.
    std::normal_distribution<double> distribution_posx(offset(0)*0.5, offset(0)*0.5*0.8); //Half domain mean, 0.1 std. dev of that, put around center just random
    std::normal_distribution<double> distribution_posy(offset(1)*0.5, offset(1)*0.5*0.8); //Half domain mean, 0.1 std. dev of that, put around center just random
    
    for (size_t indf = 0; indf < fine_bins; indf++)
    {
        Fine_Grains[indf](0) = 0.0 * (GlobalFine*Domain)/double(fine_bins); //TEMPo Remove 0.0 *
        Fine_Grains[indf](1) = 0.0 * max(distribution(generator) , min_thick); //Random thickness //TEMPo Remove thickness to remove thin floes 0.0 *
        Fine_Grains[indf](2) = max(distribution_posx(generator_pos) , min_pos); //Random xpos
        Fine_Grains[indf](3) = max(distribution_posy(generator_pos) , min_pos); //Random ypos
        
        //Fix for Periodic PBCs if needed
        if (Fine_Grains[indf](2) < 0.0 ) {
           Fine_Grains[indf](2) += offset(0);
        }
        else if (Fine_Grains[indf](2) > offset(0)) {
           Fine_Grains[indf](2) -= offset(0);
        }
        if (Fine_Grains[indf](3) < 0.0 ) {
           Fine_Grains[indf](3) += offset(1);
        }
        else if (Fine_Grains[indf](3) > offset(1)) {
           Fine_Grains[indf](3) -= offset(1);
        }

        //cout <<"indf: " << indf << " " << Fine_Grains[indf](0) << " , " << Fine_Grains[indf](1) << endl;
        Fine_Velocities[indf](0) = 0.0; //Start at rest
        Fine_Velocities[indf](1) = 0.0;
    }
    
    //Define global heat quantities
    //double qvert = 20.0; //W/(m2*K) //200  (0.1 W/m2*K) (Thickness in meters)
    
    //December 20, 2022 MODIF (Temporal)
    double qvert = VERTQ; //OFF for using qv
    double qvertOc = qvert; //OFF 
    //double qvert = 25.0; //ON for using Kappa
    //double qvertOc = qvert; //ON    //qvert*0.001; //From m to km //Why??

    //Done up already, duplicate.
    // if (caseNo == 1077){
    //     qvert = 5.0;
    //     qvertOc = 5.0;
    // }    
    
    //December 20, 2022 MODIF (Temporal)
    //double Khor = 20.0; //OFF for qv   // og 400 //Horizontal heat coff. in water //100 //Range: 100-1000 m2/s  (1 m2/s actually) (apply 1e-6 km2/s later)
    double Khor = kappa_in;
    //double Khor = Vert_q; //ON for kappa
    
    double k_adjust = cell_sizex * cell_sizey * 1e-6; //given that cell size is in km and but all units in heat equation are in meters, for kHor also in metric units
    double Tice = -1.8; //Melting temperature of sea ice
    double Qatm = QATM; //450 //490 // 200 //210 before fine heat removal try  W/m2
    double Aatm = 36.6; // 70   W/m2          Afit=36.6 W/m2 and slope Bfit=4.8 W/m2/C versus 38 and 0
    double Batm = 4.8; // 10   W/(m2*K) 

    //Other constants
    double rho_ocean = 1030.00; //kg/m3
    double cLf = 3991.00; //Prior 3000 //J/(kg * K) //CHECK???
    //double H_layer = 25.00; //Layer of mixing 20-30 meters
    double H_layer = mixed_layerh_in;

    //Albedo terms
    double alpha_ice =  0.8; //Related to albedo (0.6-0.8)
    double alpha_ocean = 0.3; //Related to albedo 0.3
    double a_o = 0.3; //Albedo of open ocean
    double a_i = 0.7; //Albedo of sea ice

    //Fine Average Terms
    double a_x; //Albedo of ocean plus fines
    double T_x; //Ocean and fines averaged temperature for a specific cell given presence of fines averaged on all fine cells
    double T_fave; //Average Temperature for all fine cells (where alpha_intersect = 0 only)
    double T_xave; //Average of T_x for all fine cells found
    double c_f; //Average Fine concentration of ice on all fine cells
    double A_cell = cell_sizex * cell_sizey; //Area of cell

    //Initialize to get averages
    T_fave = 0.0; 
    T_xave = 0.0; 
    c_f = 0.0;
    size_t n_fine_cells = 0;
    for (size_t i = 0; i < y_cells; i++) {
        for (size_t j = 0; j < x_cells; j++) {
            if ( alpha_ice_inter[j+i*x_cells] == 0 )
            {
                T_fave += oceanTemp0[j+i*x_cells];
                n_fine_cells++;
            }
        }
    }

    for (size_t i = 0; i < Fine_Grains.size(); i++)
    {
       c_f += Fine_Grains[i](0); //Add all Areas (original and new)
    }

    //Avoid problems if the systems is completely devoid of fine cells
    if (n_fine_cells < 1)
    {
        n_fine_cells = 1;
    }

    //Average Global Temperature of Fine Cells
    T_fave /= n_fine_cells;
    //Average Fines in Fine Cells
    c_f /=  (n_fine_cells * A_cell);
    //Average T_x in Fine Cells
    for (size_t i = 0; i < y_cells; i++) {
        for (size_t j = 0; j < x_cells; j++) {
            if ( alpha_ice_inter[j+i*x_cells] == 0 )
            {
                T_xave += c_f * Tice + oceanTemp0[j+i*x_cells] * (1 - c_f);
            }
        }
    }
    T_xave /= n_fine_cells; 
    T_x = T_xave;

    //Update a_x as well
    a_x = c_f * a_i + a_o * (1 - c_f);

      //Bin Size Info    
//    double melt_slope;
//    if (year_case == 1)
//    {
//        melt_slope = -3.0081e-7*1.1; //*0.8 // -1.26e-15; //Should be 3.0081 e -7     //Melt slope for fines, model implies -1.263277 % / day  //Tune with melt step and time units //2020
//    }
//    else
//    {
//        melt_slope = 1.2 * 1e-2  *  (-1.03/86400.0); //*(1e-2)*1.2; //2018 for -1.03 % / day
//    }
    
    double GlobalMaxD; //For fixed bin size
    double TaveGlobal;
    double initialAveThick;
    double MAXDD = 55.000; //42.84503
    double MINDD = 1.70271;
    
    size_t bin_size = 20; //Default is 10....  //40
    vector<double> Diameters(bin_size + 1); //Will be updated below and kept constant
    vector<double> MeltBin(bin_size + 1); //Control mass flux due to melting (Area Lost due to process, if positive it will be Area gained)
    vector<double> BreakBin(bin_size + 1); //Control mass flux due to breakage (Area Lost due to process, if positive it will be Area gained)
    vector<double> GMeltBin(bin_size + 1); //Control mass flux due to melting (Global)
    vector<double> GBreakBin(bin_size + 1); //Control mass flux due to breakage (Global)
    
    //FSD specific counts (April 24, 2023)
    vector<size_t> NLMelt(bin_size + 1); //Number of floes lost per bin due to melt across snapshot
    vector<size_t> NLBreak(bin_size + 1); //Number of floes lost per bin due to break across snapshot
    vector<size_t> NGMelt(bin_size + 1); //Number of floes gained per bin due to melt across snapshot
    vector<size_t> NGBreak(bin_size + 1); //Number of floes gained per bin due to break across snapshot
    
    //Lateral rate by bin
    vector<double> LatRate(bin_size + 1); 
    vector<double> LatRate_step(bin_size + 1); 
    
    //Initialize global values
    for (size_t i = 0; i < MeltBin.size(); i++)
    {
        GMeltBin[i] = 0.0;
        GBreakBin[i] = 0.0;
        
        //FSD specific counts (April 24, 2023)
        NLMelt[i] = 0;
        NLBreak[i] = 0;
        NGMelt[i] = 0;
        NGBreak[i] = 0;
        LatRate[i] = 0.0;
        LatRate_step[i] = 0.0;
    }
    
    //For fine variation
    size_t step_b = 0;
    size_t step_e = 0;
        
    //Limit Diameters for Processes
//    double limitMinDiam = 6.00; //6.00 DF
//    double limitMaxDiam = 30.00;  //20.00 DF  maxmin 27.00 and 13.00 as examples
    
    //To avoid Nan alterations
    double prevConc, prevNormArea, prevConcB, prevConcF, prevSurfArea, prevThickF, prevThickC;
    
    //For statistics
    double Mean_Diam, Max_Diam, Min_Diam, prevMean_Diam, prevMax_Diam, prevMin_Diam;
    
    //For external variables
    size_t out_Tcount = 0;  //Counts days
    size_t waveH_count = 0; //Counts hours
    size_t vel_index = 0;

    //For random changes
    bool change_switch = false;
    
    //Aug 18,2022 Change
    //BEGIN
    //Initialize mass loss variables WILL CALCULATE JUST BEFORE PRINT
    double tmass, tmasscoarse, tmassA, tmasscoarseA, tmassfines, loss_mcl, loss_mcv, loss_mcv_solar, loss_mcv_ocean, gain_fines, loss_fines, MLAT, MBASAL; // No need for init. (most)
    loss_mcl = 0.0; loss_mcv = 0.0; gain_fines = 0.0; loss_fines = 0.0; //Only these do
    loss_mcv_solar = 0.0; loss_mcv_ocean = 0.0; MLAT = 0.0; MBASAL = 0.0;
    
    //Initialize area loss variables
     double area_loss_bkg, area_loss_basal, area_loss_lat, area_loss_lat2, ave_delta_r, ave_delta_r_rate, ave_delta_r_rate_step;
     area_loss_bkg = 0.0; area_loss_basal = 0.0; area_loss_lat = 0.0; area_loss_lat2 = 0.0; ave_delta_r = 0.0; ave_delta_r_rate = 0.0; ave_delta_r_rate_step = 0.0;
     double prev_coarse_area = initial_coarse*Domain;
    
    //Thickness snapshots choice
    size_t snp1 = nT * 0.15;
    size_t snp2 = nT * 0.25;
    size_t snp3 = nT * 0.3;
    //END
    
    //*********/*********/*********/*********/*********/*********/*********/*********/*********/*********/*********//
    //MAIN LOOP
    for (size_t step = 0; step < nT; ++step)
    {
        //cout << "Start Timestep: " << step << endl; //For debugging
        //cout << "Adjust Grains" << endl;
        //world.readjustGrains();
        //cout << "Compute World" << endl;
        
        //Initialize in case no melt or break occurs
        for (size_t i = 0; i < MeltBin.size(); i++)
        {
            MeltBin[i] = 0.0;
            BreakBin[i] = 0.0;
            //LatRate[i] = 0.0;
        }
        
        //Adjust periodic movement
        change_switch = false;
        srand (time(NULL));
        double slopedir0 = slopedir;
        
        //size_t t_freq = nT/20; //How many times will field change (More freq.)
        size_t t_freq = nT/5; //How many times will field change (Less cycles)
        
        double load_fac = -0.5 * ( sin( (3.141592653589793*float(step)/t_freq) + 0.5*3.141592653589793 ) ) + 0.5;  //Ranges from 0 to 1
        double slope_dir = (load_fac - load_fac0);  //WARNING MAKE SURE BEGINNING IS A BIT ABOVE FREQ
        load_fac0 = load_fac;
        world.changeSlopeDir(slope_dir);

        if (step % t_freq == 0 && step > 1)
        {
            cout << "Change Currents!*!" << endl;
            change_switch = true;
        }

        if(change_switch == true){
            //Change angle
            //Start random number
            srand (time(NULL));
            //Define random in radians (UNIF DIST)
            double ang_line = rand() % 1257; //From -2pi to 2pi              
            ang_line -= 628;
            ang_line *= 0.01;
            cout << "Change angle line to: " << ang_line << endl;
            //world.changeFlowAngle(ang_line);

            //Change force (GAUSSIAN DIST)
            std::default_random_engine generator_force;
            std::normal_distribution<double> distribution_force(flowforceRef0, flowforceRef0*0.5); //0.5m mean, 0.5 std. dev
            double flowforce = max(distribution_force(generator_force) , flowforcemin); 
            cout << "Change flow force to: " << flowforce << endl;
            //world.changeFlowForce(flowforce);
        }
        
        //Constant velocity based on field data (all floes move the same)
        world.changeFlowAngle(floeVel[vel_index](2));  //In radians
        double flowforce = flowforceRef0 *  floeVel[vel_index](1) / ave_floeVel;  //Multiply by zero to make floes suspended in a STATIC ocean  //1.00 to move
      	world.changeFlowForce( flowforce ); //Magnitude
      	
      	//Update fluid grid (Jan 9, 2023) //Scale initial vortex depending on field data
       for (size_t i = 0; i < y_cells; i++) {
            for (size_t j = 0; j < x_cells; j++) {
                
                //Constant
                Uwg[j+i*x_cells](0) = Uwg[j+i*x_cells](0) * 1.0; //1.0  //CAREFUL THIS ACCUMULATES EACH STEP!!!!!!!! Use Uwg0 instead!!!! March 1, 2023
                Uwg[j+i*x_cells](1) = Uwg[j+i*x_cells](1) * 1.0;
                // //Variable
                // Uwg[j+i*x_cells](0) = Uwg[j+i*x_cells](0) * (floeVel[vel_index](1) / ave_floeVel);
                // Uwg[j+i*x_cells](1) = Uwg[j+i*x_cells](1) * (floeVel[vel_index](1) / ave_floeVel);
            }    
        }   
        world.changeUw(Uwg);
        //Update fluid grid (Jan 9, 2023)

        //****************************** COMPUTE WORLD, MELTING AND FRACTURE WORLD FUNCTIONS ******************************
        //Update contact interactions and ocean currents
        //cout << "Start WorldState" << endl;
        world.computeWorldState();
        //cout << "End WorldState" << endl;
        
        //Start processes after initializing a bit
        if (step>=tiso){

            //Add State particular steps
            if (step % nStress == 0 ){
                //world.StressStateForceFinder(int(step/nTout));
            }   

            //WARNING: Initializing thickness and melt requires an intial step inside Temperature Modif in Grains2d
            //If heat PDE is okay you DO NOT need to time step
            //if (step > 1){
            
            //Here dstep is to advance the process faster, could be equal to 1 if computer was faster
            if (step % nTemp == 0 && step > 1){
                double outwaterTemp = -0.4; //Generic value for constant
                if (caseNo != 1077){ //For now avoid melt for fracture experiments            
                    cout << "Melt step: " << step << endl;
                    world.TempState(MeltBin, Diameters, Out_temperature[out_Tcount](1), limitMaxDiam, limitMinDiam, dstep, qvert, Khor, alpha_ice, Qatm, Aatm, Batm, Fine_Grains, a_i, a_x, T_fave, T_xave, loss_mcv, loss_fines, NLMelt, NLBreak, NGMelt, NGBreak, loss_mcv_solar, loss_mcv_ocean, area_loss_basal, area_loss_lat, ave_delta_r, LatRate, MLAT, MBASAL);  
                    //world.TempState(MeltBin, Diameters, Out_temperature[out_Tcount](1), limitMaxDiam, limitMinDiam, dstep, qvert, Khor, alpha_ice, Qatm, Aatm, Batm, Fine_Grains, a_i, a_x, T_fave, T_xave);   //Change Aug 22, 2022
                    //world.TempState(MeltBin, Diameters, outwaterTemp, limitMaxDiam, limitMinDiam, dstep, qvert, Khor, alpha_ice, Qatm, Aatm, Batm, Fine_Grains);
                    
                }
                else
                {
                    //qvert = 0.000;
                    cout << "Melt step: " << step << endl;
                    world.TempState(MeltBin, Diameters, Out_temperature[out_Tcount](1), limitMaxDiam, limitMinDiam, dstep, qvert, Khor, alpha_ice, Qatm, Aatm, Batm, Fine_Grains, a_i, a_x, T_fave, T_xave, loss_mcv, loss_fines, NLMelt, NLBreak, NGMelt, NGBreak, loss_mcv_solar, loss_mcv_ocean, area_loss_basal, area_loss_lat, ave_delta_r, LatRate, MLAT, MBASAL);  
                }
            }   
            
            //Get lateral rate
            //ave_delta_r_rate = ave_delta_r/nTemp; //Get rate average radial reduction for this process for all floes divided by nTemp
            ave_delta_r_rate = ave_delta_r;
            ave_delta_r_rate_step += ave_delta_r_rate;
        
            if (step % nTemp == 0 && step > 1){   
                cout << "Melt Fully step: " << step << endl;
                world.Melty(MeltBin, Diameters, limitMinDiam, Fine_Grains, Fine_Velocities, loss_mcl, loss_mcv, gain_fines, NLMelt, NLBreak, NGMelt, NGBreak, loss_mcv_solar, loss_mcv_ocean, area_loss_bkg); //Change Aug 22, 2022
                //world.Melty(MeltBin, Diameters, limitMinDiam, Fine_Grains, Fine_Velocities);
                world.FineMelty(Fine_Grains, Fine_Velocities, loss_fines, NLMelt, NLBreak, NGMelt, NGBreak);
                //world.FineMelty(Fine_Grains, Fine_Velocities);  //Change Aug 22, 2022
                
                double Temp_fines = 0.0; //Area sum
                for (size_t i = 0; i < Fine_Grains.size(); i++)
                {
                   Temp_fines += Fine_Grains[i](0); //Add all Areas (original and new)
                }
                Temp_fines /= Domain; //From area to conc
                //Update all the time global fines
                //step_e = step; // Measure time step size
                //GlobalFine += melt_slope * double(step_e-step_b) + (MeltBin.back()/Domain);  //Substract melted, add fine fraction if applies
                //Has to be more than zero
                GlobalFine = max(Temp_fines,0.0);
                //step_b = step_e; //Update where we measure
            }    
            
            //world.StressState();
            //world.DamageState();
            
            //cout << "Try Fracture Routine" << endl;
            //if (step % (nBreak+3) == 0 && step <= 2000000){

            //Since breakage is still very parametrized we can control with timestep
            //if ( step > 1 ){  //PROB CHANGE ON
            //if ( (step) % (1) == 0 && step > 1 ){ //Modif Sep 30, 2022  //PROB CHANGE (Check for breakage for for all steps since Tfail varies, just skip bkg if conditions not fulfilled) CHANGE OCT-18-2022
            //if ( (step) % (nBreak) == 0 && step > 1 ){ //Modif Sep 30, 2022  //PROB CHANGE
            if ( (step-2) % (nBreak) == 0 && step > 1 && fluid_mode == 1 ){ //PROB CHANGE OFF  //Parametrized breakage
                cout << "Bkg" << step << endl;
                double WaveHOut = 0.5; //Generic values
                double WaveLOut = 125.0; //Generic Values
                //world.FracRoutine(BreakBin, Diameters, limitMaxDiam, Wave_HD[out_Tcount](1), Wave_LD[out_Tcount](1)); //Day
                //world.FracRoutine(BreakBin, Diameters, limitMaxDiam, WaveH[waveH_count](2), WaveL[waveH_count](2));  //Hourly
                world.FracRoutine(BreakBin, Diameters, limitMaxDiam, WaveHOut, WaveLOut, BREAK_PROB, NLMelt, NLBreak, NGMelt, NGBreak); //Constant over time
                //Check if fracture is randomly induced or contact-based
            }              
            
            //Just for test
            if ( caseNo == 1077 && step % nTout == 0 )
            {
                double WaveHOut = 0.5; //Generic values
                double WaveLOut = 125.0; //Generic Values
                world.FracRoutine_Nova(BreakBin, Diameters, limitMaxDiam, WaveHOut, WaveLOut, BREAK_PROB, outDir4);
            }
            
            if (rank == 0)
            {    
                //world.FractureState();
            }
            
            bool day_change = false;
            //Update external variables based on multiplicity of how you often get them from sensors
            if (step % (nTout/24) == 0) //Hourly
            {
                waveH_count++;
            }
            if (step % nTout == 0) //Daily
            {
                out_Tcount++;
                day_change = true;
            }

            if ( day_change && vel_index < floeVel.size()-1 ) //Only if change day and still vector space
            {
               if (int(step / nTout) >= floeVel[vel_index+1](0))  //Since we have gaps in data
               {
                 vel_index++;
               } 
            }
        }    

        
        //****************************** UPDATE GLOBAL OCEAN TEMP GRID ******************************
        oceanTemp0 = oceanTemp;

        //Get mean for convenience
        //Average Global Temperature in Domain (use prior time step to get next time step info) (FOR ALL)
        double TaveGlobalg = 0.0;
        for (size_t i = 0; i < y_cells; i++) {
            for (size_t j = 0; j < x_cells; j++) {
                TaveGlobalg += oceanTemp[j+i*x_cells];
            }
        }
        TaveGlobalg /= (x_cells*y_cells);

        //Initialize to get averages
        T_fave = 0.0; 
        T_xave = 0.0; 
        c_f = 0.0;
        n_fine_cells = 0;
        for (size_t i = 0; i < y_cells; i++) {
            for (size_t j = 0; j < x_cells; j++) {
                if ( alpha_ice_inter[j+i*x_cells] == 0 )
                {
                    T_fave += oceanTemp0[j+i*x_cells];
                    n_fine_cells++;
                }
            }
        }

        for (size_t i = 0; i < Fine_Grains.size(); i++)
        {
           c_f += Fine_Grains[i](0); //Add all Areas (original and new)
        }

        //Avoid problems if the systems is completely devoid of fine cells
        if (n_fine_cells == 0)
        {
            n_fine_cells = 1;
        }

        //Average Global Temperature of Fine Cells
        T_fave /= n_fine_cells;
        //Average Fines in Fine Cells
        c_f /=  (n_fine_cells * A_cell);
        //Average T_x in Fine Cells
        for (size_t i = 0; i < y_cells; i++) {
            for (size_t j = 0; j < x_cells; j++) {
                if ( alpha_ice_inter[j+i*x_cells] == 0 )
                {
                    T_xave += c_f * Tice + oceanTemp0[j+i*x_cells] * (1 - c_f);
                }
            }
        }
        T_xave /= n_fine_cells; 
        T_x = T_xave;

        //Update a_x as well
        a_x = c_f * a_i + a_o * (1 - c_f);
            

        //Update Ocean Temp Grid
        double flux_control;
        double external_temp = 1.4; //Can be coupled with field data instead of using a constant value
        bool PBC_grid = true; //True uses periodic boundary conditions, false uses dirichlet with a constant external_temp
        
        if (step % nTemp == 0 && step > 1){
           
            cout << "Updating Ocean Temp Grid" << endl;
            for (size_t i = 0; i < y_cells; i++) {
                for (size_t j = 0; j < x_cells; j++) {
                    
                    //Obtain alpha coefficient for this step
                    ice_contact = false;
                    alpha_ice_inter[j+i*x_cells] = 0.0; //0.0 Fine Cell, 1.0 Coarse Cell
                    for (size_t gi = 0; gi < grainList.size(); gi++) {

                        //Function to check ice contact
                        ice_contact = under_ice(grainList[gi], fluid_coord[j+i*x_cells]);

                        if (ice_contact)
                        {
                            alpha_ice_inter[j+i*x_cells] = 1.0;
                            break;
                        }
                    }

                    //Fine Cell (update averages)    //Coarse Cell (no need to update averages)
                    if ( alpha_ice_inter[j+i*x_cells] == 0 )
                    {
                        T_x = c_f * Tice + oceanTemp0[j+i*x_cells] * (1 - c_f);
                    }

                    flux_control = cLf * rho_ocean * H_layer ; //Accounts for density and mixed layer thickness right now ==> 1 kg/m2 units or 1e6 kg/km2, just multiply by right H and rho when ready
                    
                    if (i == 0) //DOWN
                    {
                        //cout << "Ice cover: " << alpha_ice_inter[j+i*x_cells] << endl;
                        if (PBC_grid)
                        {
                            if (j == 0) //Lower left corner
                            {
                                oceanTemp[j+i*x_cells] = oceanTemp0[j+i*x_cells] + dstep*k_adjust*(Khor) * (oceanTemp0[j+(i+1)*x_cells]+oceanTemp0[(j+1)+i*x_cells]+oceanTemp0[j+(y_cells-1)*x_cells]+oceanTemp0[(x_cells-1)+i*x_cells]-4*oceanTemp0[j+i*x_cells]) + (dstep/flux_control) * alpha_ice_inter[j+i*x_cells] * ( - qvert*(oceanTemp0[j+i*x_cells] - Tice) ) + (dstep/flux_control) * (1-alpha_ice_inter[j+i*x_cells]) * ( Qatm*(1-a_x) - (Aatm + Batm*T_x) - qvert*( T_fave - Tice) );          //Combined expression using alpha for Coarse and 1 - alpha for fine components
                            }
                            else if (j == x_cells-1) //Lower right corner
                            {
                                oceanTemp[j+i*x_cells] = oceanTemp0[j+i*x_cells] + dstep*k_adjust*(Khor) * (oceanTemp0[j+(i+1)*x_cells]+oceanTemp0[0+i*x_cells]+oceanTemp0[j+(y_cells-1)*x_cells]+oceanTemp0[(j-1)+i*x_cells]-4*oceanTemp0[j+i*x_cells]) + (dstep/flux_control) * alpha_ice_inter[j+i*x_cells] * ( - qvert*(oceanTemp0[j+i*x_cells] - Tice) ) + (dstep/flux_control) * (1-alpha_ice_inter[j+i*x_cells]) * ( Qatm*(1-a_x) - (Aatm + Batm*T_x) - qvert*( T_fave - Tice) );          //Combined expression using alpha for Coarse and 1 - alpha for fine components
                            }
                            else //No corner
                            {
                                oceanTemp[j+i*x_cells] = oceanTemp0[j+i*x_cells] + dstep*k_adjust*(Khor) * (oceanTemp0[j+(i+1)*x_cells]+oceanTemp0[(j+1)+i*x_cells]+oceanTemp0[j+(y_cells-1)*x_cells]+oceanTemp0[(j-1)+i*x_cells]-4*oceanTemp0[j+i*x_cells]) + (dstep/flux_control) * alpha_ice_inter[j+i*x_cells] * ( - qvert*(oceanTemp0[j+i*x_cells] - Tice) ) + (dstep/flux_control) * (1-alpha_ice_inter[j+i*x_cells]) * ( Qatm*(1-a_x) - (Aatm + Batm*T_x) - qvert*( T_fave - Tice) );          //Combined expression using alpha for Coarse and 1 - alpha for fine components
                            }
                        }
                        else
                        {
                            oceanTemp[j+i*x_cells] = external_temp; //Keep constant external heat (D.B.C.)
                        }
                    }
                    else if (j == 0) //LEFT
                    {
                        if (PBC_grid)
                        {
                            if (i == 0) //Lower left corner
                            {    
                                oceanTemp[j+i*x_cells] = oceanTemp0[j+i*x_cells] + dstep*k_adjust*(Khor) * (oceanTemp0[j+(i+1)*x_cells]+oceanTemp0[(j+1)+i*x_cells]+oceanTemp0[j+(y_cells-1)*x_cells]+oceanTemp0[(x_cells-1)+i*x_cells]-4*oceanTemp0[j+i*x_cells]) + (dstep/flux_control) * alpha_ice_inter[j+i*x_cells] * ( - qvert*(oceanTemp0[j+i*x_cells] - Tice) ) + (dstep/flux_control) * (1-alpha_ice_inter[j+i*x_cells]) * ( Qatm*(1-a_x) - (Aatm + Batm*T_x) - qvert*( T_fave - Tice) );          //Combined expression using alpha for Coarse and 1 - alpha for fine components
                            }
                            if (i == y_cells-1) //Upper left corner
                            {    
                                oceanTemp[j+i*x_cells] = oceanTemp0[j+i*x_cells] + dstep*k_adjust*(Khor) * (oceanTemp0[j+0*x_cells]+oceanTemp0[(j+1)+i*x_cells]+oceanTemp0[j+(i-1)*x_cells]+oceanTemp0[(x_cells-1)+i*x_cells]-4*oceanTemp0[j+i*x_cells]) + (dstep/flux_control) * alpha_ice_inter[j+i*x_cells] * ( - qvert*(oceanTemp0[j+i*x_cells] - Tice) ) + (dstep/flux_control) * (1-alpha_ice_inter[j+i*x_cells]) * ( Qatm*(1-a_x) - (Aatm + Batm*T_x) - qvert*( T_fave - Tice) );          //Combined expression using alpha for Coarse and 1 - alpha for fine components
                            }
                            else //No corner
                            {    
                                oceanTemp[j+i*x_cells] = oceanTemp0[j+i*x_cells] + dstep*k_adjust*(Khor) * (oceanTemp0[j+(i+1)*x_cells]+oceanTemp0[(j+1)+i*x_cells]+oceanTemp0[j+(i-1)*x_cells]+oceanTemp0[(x_cells-1)+i*x_cells]-4*oceanTemp0[j+i*x_cells]) + (dstep/flux_control) * alpha_ice_inter[j+i*x_cells] * ( - qvert*(oceanTemp0[j+i*x_cells] - Tice) ) + (dstep/flux_control) * (1-alpha_ice_inter[j+i*x_cells]) * ( Qatm*(1-a_x) - (Aatm + Batm*T_x) - qvert*( T_fave - Tice) );          //Combined expression using alpha for Coarse and 1 - alpha for fine components
                            }
                        }
                        else
                        {
                            oceanTemp[j+i*x_cells] = external_temp; //Keep constant external heat (D.B.C.)
                        }
                    }
                    else if (i == y_cells-1) //UP
                    {
                        if (PBC_grid)
                        {
                            if (j == 0) //Upper left corner
                            {
                                oceanTemp[j+i*x_cells] = oceanTemp0[j+i*x_cells] + dstep*k_adjust*(Khor) * (oceanTemp0[j+0*x_cells]+oceanTemp0[(j+1)+i*x_cells]+oceanTemp0[j+(i-1)*x_cells]+oceanTemp0[(x_cells-1)+i*x_cells]-4*oceanTemp0[j+i*x_cells]) + (dstep/flux_control) * alpha_ice_inter[j+i*x_cells] * ( - qvert*(oceanTemp0[j+i*x_cells] - Tice) ) + (dstep/flux_control) * (1-alpha_ice_inter[j+i*x_cells]) * ( Qatm*(1-a_x) - (Aatm + Batm*T_x) - qvert*( T_fave - Tice) );          //Combined expression using alpha for Coarse and 1 - alpha for fine components
                            }
                            else if (j == x_cells-1) //Upper right corner
                            {
                                oceanTemp[j+i*x_cells] = oceanTemp0[j+i*x_cells] + dstep*k_adjust*(Khor) * (oceanTemp0[j+0*x_cells]+oceanTemp0[0+i*x_cells]+oceanTemp0[j+(i-1)*x_cells]+oceanTemp0[(j-1)+i*x_cells]-4*oceanTemp0[j+i*x_cells]) + (dstep/flux_control) * alpha_ice_inter[j+i*x_cells] * ( - qvert*(oceanTemp0[j+i*x_cells] - Tice) ) + (dstep/flux_control) * (1-alpha_ice_inter[j+i*x_cells]) * ( Qatm*(1-a_x) - (Aatm + Batm*T_x) - qvert*( T_fave - Tice) );          //Combined expression using alpha for Coarse and 1 - alpha for fine components
                            }
                            else //No corner
                            {
                                oceanTemp[j+i*x_cells] = oceanTemp0[j+i*x_cells] + dstep*k_adjust*(Khor) * (oceanTemp0[j+0*x_cells]+oceanTemp0[(j+1)+i*x_cells]+oceanTemp0[j+(i-1)*x_cells]+oceanTemp0[(j-1)+i*x_cells]-4*oceanTemp0[j+i*x_cells]) + (dstep/flux_control) * alpha_ice_inter[j+i*x_cells] * ( - qvert*(oceanTemp0[j+i*x_cells] - Tice) ) + (dstep/flux_control) * (1-alpha_ice_inter[j+i*x_cells]) * ( Qatm*(1-a_x) - (Aatm + Batm*T_x) - qvert*( T_fave - Tice) );          //Combined expression using alpha for Coarse and 1 - alpha for fine components
                            }                    
                        }
                        else
                        {
                            oceanTemp[j+i*x_cells] = external_temp; //Keep constant external heat (D.B.C.)
                        }
                    }
                    else if (j == x_cells-1) //RIGHT
                    {
                        if (PBC_grid)
                        {
                            if (i == 0) //Lower right corner
                            {    
                                oceanTemp[j+i*x_cells] = oceanTemp0[j+i*x_cells] + dstep*k_adjust*(Khor) * (oceanTemp0[j+(i+1)*x_cells]+oceanTemp0[0+i*x_cells]+oceanTemp0[j+(y_cells-1)*x_cells]+oceanTemp0[(j-1)+i*x_cells]-4*oceanTemp0[j+i*x_cells]) + (dstep/flux_control) * alpha_ice_inter[j+i*x_cells] * ( - qvert*(oceanTemp0[j+i*x_cells] - Tice) ) + (dstep/flux_control) * (1-alpha_ice_inter[j+i*x_cells]) * ( Qatm*(1-a_x) - (Aatm + Batm*T_x) - qvert*( T_fave - Tice) );          //Combined expression using alpha for Coarse and 1 - alpha for fine components
                            }
                            if (i == y_cells-1) //Upper right corner
                            {    
                                oceanTemp[j+i*x_cells] = oceanTemp0[j+i*x_cells] + dstep*k_adjust*(Khor) * (oceanTemp0[j+(i+1)*x_cells]+oceanTemp0[0+i*x_cells]+oceanTemp0[j+(i-1)*x_cells]+oceanTemp0[(j-1)+i*x_cells]-4*oceanTemp0[j+i*x_cells]) + (dstep/flux_control) * alpha_ice_inter[j+i*x_cells] * ( - qvert*(oceanTemp0[j+i*x_cells] - Tice) ) + (dstep/flux_control) * (1-alpha_ice_inter[j+i*x_cells]) * ( Qatm*(1-a_x) - (Aatm + Batm*T_x) - qvert*( T_fave - Tice) );          //Combined expression using alpha for Coarse and 1 - alpha for fine components
                            }
                            else //No corner
                            {    
                                oceanTemp[j+i*x_cells] = oceanTemp0[j+i*x_cells] + dstep*k_adjust*(Khor) * (oceanTemp0[j+(i+1)*x_cells]+oceanTemp0[0+i*x_cells]+oceanTemp0[j+(i-1)*x_cells]+oceanTemp0[(j-1)+i*x_cells]-4*oceanTemp0[j+i*x_cells]) + (dstep/flux_control) * alpha_ice_inter[j+i*x_cells] * ( - qvert*(oceanTemp0[j+i*x_cells] - Tice) ) + (dstep/flux_control) * (1-alpha_ice_inter[j+i*x_cells]) * ( Qatm*(1-a_x) - (Aatm + Batm*T_x) - qvert*( T_fave - Tice) );          //Combined expression using alpha for Coarse and 1 - alpha for fine components    
                            }                    
                        }
                        else
                        {
                            oceanTemp[j+i*x_cells] = external_temp; //Keep constant external heat (D.B.C.)
                        }
                    }
                    else
                    {
                        //PRIOR EXPRESSION: //oceanTemp[j+i*x_cells] = oceanTemp0[j+i*x_cells]+dstep*k_adjust*(Khor)*(oceanTemp0[j+(i+1)*x_cells]+oceanTemp0[(j+1)+i*x_cells]+oceanTemp0[j+(i-1)*x_cells]+oceanTemp0[(j-1)+i*x_cells]-4*oceanTemp0[j+i*x_cells]) + (dstep/flux_control)*( qvertOc*alpha_ice_inter[j+i*x_cells]*(Tice-oceanTemp0[j+i*x_cells]) + (1 - alpha_ice_inter[j+i*x_cells])*(Qatm*(1-alpha_ocean)-(Aatm+Batm*oceanTemp0[j+i*x_cells])) );  //Add fine energy exchange (cooling as well) + dt*(Af/totalArea)*qvert*(Tfreeze-Tocean)  //Try later:  - meltV*(_grainTemp2D0.getGridValue2(ii,jj)-meltTemp);
                        
                        //MAIN EXPRESSION
                        //oceanTemp[j+i*x_cells] = oceanTemp0[j+i*x_cells] + dstep*k_adjust*(Khor) * (oceanTemp0[j+(i+1)*x_cells]+oceanTemp0[(j+1)+i*x_cells]+oceanTemp0[j+(i-1)*x_cells]+oceanTemp0[(j-1)+i*x_cells]-4*oceanTemp0[j+i*x_cells]) + (dstep/flux_control) * alpha_ice_inter[j+i*x_cells] * ( - qvert*(oceanTemp0[j+i*x_cells] - Tice) ) + (dstep/flux_control) * (1-alpha_ice_inter[j+i*x_cells]) * ( Qatm*(1-a_x) - (Aatm + Batm*T_x) - qvert*( T_fave - Tice) );          //Combined expression using alpha for Coarse and 1 - alpha for fine components
                                                    //Prior Step                //Denominator and Khor                              //Laplacian Term                                                                                                           //Denom                //Coarse Cell Switch      //Coarse Cell Expression                                                           //Denom                  //Fine Cell Switch               //Fine Cell Expression
                        //PROPOSED EXPRESSION
                        oceanTemp[j+i*x_cells] = oceanTemp0[j+i*x_cells] + dstep*k_adjust*(Khor) * (oceanTemp0[j+(i+1)*x_cells]+oceanTemp0[(j+1)+i*x_cells]+oceanTemp0[j+(i-1)*x_cells]+oceanTemp0[(j-1)+i*x_cells]-4*oceanTemp0[j+i*x_cells]) + (dstep/flux_control) * alpha_ice_inter[j+i*x_cells] * ( - qvert*(oceanTemp0[j+i*x_cells] - Tice) ) + (dstep/flux_control) * (1-alpha_ice_inter[j+i*x_cells]) * ( Qatm*(1-a_o)*(1-c_f) - (1-c_f)*(Aatm + Batm*oceanTemp0[j+i*x_cells]) - qvert*c_f*(T_fave - Tice) );          //Combined expression using alpha for Coarse and 1 - alpha for fine components
                        //If relevant, apply for boundary cells too.
                        
                        //Track lateral flux
                        latFlux[j+i*x_cells] = flux_control*k_adjust*(Khor)* (oceanTemp0[j+(i+1)*x_cells]+oceanTemp0[(j+1)+i*x_cells]+oceanTemp0[j+(i-1)*x_cells]+oceanTemp0[(j-1)+i*x_cells]-4*oceanTemp0[j+i*x_cells]);
                        vertFlux[j+i*x_cells] = ( qvert*(oceanTemp0[j+i*x_cells] - Tice) );
                        if (vertFlux[j+i*x_cells] == 0){
                            lvRatio[j+i*x_cells] = 1111; //To show no vertical rate (very unlikely, avoid division by zero)
                        }
                        else{
                            lvRatio[j+i*x_cells] = latFlux[j+i*x_cells] /  vertFlux[j+i*x_cells];
                        }
                        
                    }
                    
                    //Update net fluxes
                    netFlux[j+i*x_cells]  = alpha_ice_inter[j+i*x_cells] * ( 1.0*( Qatm*(1-a_i) - (Aatm + Batm*Tice) ) + 0.0  )  +  (1 - alpha_ice_inter[j+i*x_cells]) * ( c_f*( Qatm*(1-a_i) - (Aatm + Batm*Tice) ) + ( 1-c_f ) * ( Qatm*(1-a_o) - (Aatm + Batm*oceanTemp[j+i*x_cells]) ) );
                    netFlux2[j+i*x_cells] = alpha_ice_inter[j+i*x_cells] * (              0.0                                 )  +  (1 - alpha_ice_inter[j+i*x_cells]) * ( c_f*( Qatm*(1-a_i) - (Aatm + Batm*Tice) ) + ( 1-c_f ) * ( Qatm*(1-a_o) - (Aatm + Batm*oceanTemp[j+i*x_cells]) ) );
                    
                }
            }
            
            //world.changeoceanTemp(oceanTemp); //CHECK THIS UPDATE!!! APRIL 26, 2023!!!!  //Only needed for very precise thickness variations //Compare A and B, 1826/1827 vs 1828/1829
            //Temperature is unaffected
            //Coarse Concentration results are very similar and error is very close as well
            //Fine Concentration results are the same
            //Coarse thickness is reduced a bit faster. 
            //Hence, use this change to explore breakage in more rigor. Concentration, Mass Loss and FSD keep the same trends.
        }

        //******************************UPDATE OUTPUT AND GSD ******************************
        //Update Global changes from Discrete time step changes from Worlds, accumulate until they are printed at output steps
        for (size_t bi = 0; bi < MeltBin.size(); bi++)
        {
            GMeltBin[bi] += MeltBin[bi];
            GBreakBin[bi] += BreakBin[bi];
        }
        
        // for (size_t bi = 0; bi < MeltBin.size(); bi++)
        // {
        //     LatRate_step[bi] += LatRate[bi];
        // }

        //Calculate Total area and Exposed Area
        double TotIce = 0.0;
        double addArea;
        double AccArea = 0.0;
        
        if (step == 1)
        {   
            double SurfArea = 0; 
            for (size_t i = 0; i < grainList.size(); i++) {

                double SurfAreap = 0;
                for (size_t j = 0; j < grainList[i].getnpoints()-1; j++) {
                    SurfAreap += (-grainList[i].getPointList()[j]+grainList[i].getPointList()[j+1]).norm();            
                } 
                SurfAreap += (grainList[i].getPointList()[grainList[i].getnpoints()-1]-grainList[i].getPointList()[0]).norm();
                SurfArea += SurfAreap;     

            }
            //cout <<"Access Text File" << endl;
            fprintf(normArea,  "%4.8f\n", SurfArea); 
            //cout <<"Access Text File" << endl;
            
            double normalArea = 0;
            for (size_t i = 0; i < grainList.size(); i++) {
                //AccArea += grainList[i].getMass();
                vector<Vector2d> forArea0 = grainList[i].getPointList();
                //addArea = PointsAreaM(forArea0);
                
                double Area = 0;
                size_t n = forArea0.size();

                for (size_t pp = 0; pp < n-1; pp++)
                {
                    Area = Area + ( forArea0[pp](0) * forArea0[pp+1](1) -  forArea0[pp](1) * forArea0[pp+1](0) ); 
                }
                Area = Area + (forArea0[n-1](0) * forArea0[0](1) -  forArea0[n-1](1) * forArea0[0](0) ); 

                addArea = 0.5*abs(Area);

                AccArea = AccArea + addArea;
            }    
            AccArea0 = AccArea;
            fprintf(normArea,  "%4.8f\n", AccArea0);
            cout << "Init Surf Area= " << SurfArea << endl;
            cout << "AccArea0= " << AccArea0 << endl;
            normalArea = AccArea/AccArea0;
            
            double thick_coarse = 0.0;
            for (size_t i = 0; i < grainList.size(); i++) {
                vector<Vector2d> forArea = grainList[i].getPointList();
                //addArea = PointsAreaM(forArea);

                double Area = 0;
                size_t n = forArea.size();

                for (size_t pp = 0; pp < n-1; pp++)
                {
                    Area = Area + ( forArea[pp](0) * forArea[pp+1](1) -  forArea[pp](1) * forArea[pp+1](0) ); 
                }
                Area = Area + (forArea[n-1](0) * forArea[0](1) -  forArea[n-1](1) * forArea[0](0) ); 

                addArea = 0.5*abs(Area);

                TotIce = TotIce + addArea;

                //thick_coarse +=  grainList[i].getThickness();  //Get average thickness //Not weighted
                thick_coarse +=  addArea * grainList[i].getThickness();  //Get average thickness //Weighted

                //cout << "Grain Area from Points: " << PointsAreaM(forArea) << endl;
            }
            
            double thick_fines = 0.0; //Area sum
            for (size_t i = 0; i < Fine_Grains.size(); i++)
            {
               //thick_fines += Fine_Grains[i](1); //Add all Areas (original and new)  //Not weighted
               thick_fines += Fine_Grains[i](0) * Fine_Grains[i](1); //Add all Areas (original and new)  //Weighted
            }
            //double concentration = (TotIce / Domain) + initial_fine;
            double concentration = GlobalConc;
            
            //double thick_tot = (thick_coarse + thick_fines) / (double (grainList.size() +  Fine_Grains.size())); //By number
            double thick_tot = (thick_coarse + thick_fines) / (concentration*Domain); //By area Coarse + Fine
            double thick_coarse_ave = (thick_coarse) / (concentration*Domain); //By total area to observe separately
            double thick_fine_ave = (thick_fines) / (concentration*Domain); //By total area to observe separately
            initialAveThick = thick_tot;

            //Average Global Temperature in Domain
            double aveFlux0 = 0.0;
            double aveFlux20 = 0.0;
            double avelatFlux0 = 0.0;
            double avevertFlux0 = 0.0;
            double avelvRatio0 = 0.0;
            TaveGlobal = 0.0;
            for (size_t i = 0; i < y_cells; i++) {
                for (size_t j = 0; j < x_cells; j++) {
                    TaveGlobal += oceanTemp[j+i*x_cells];
                    aveFlux0  += netFlux[j+i*x_cells];
                    aveFlux20 += netFlux2[j+i*x_cells];
                    //Only exists at center cells, no borders
                    if ( (i > 0 && i < y_cells - 1) && (j > 0 && j < x_cells - 1) ){
                        avelatFlux0 += latFlux[j+i*x_cells];
                        avevertFlux0 += vertFlux[j+i*x_cells];
                        avelvRatio0 +=  lvRatio[j+i*x_cells];
                    }
                }
            }
            aveFlux0 /= (x_cells*y_cells);
            aveFlux20 /= (x_cells*y_cells);
            avelatFlux0 /= ((x_cells-2)*(y_cells-2));
            avevertFlux0 /= ((x_cells-2)*(y_cells-2));
            avelvRatio0 /= ((x_cells-2)*(y_cells-2));
            TaveGlobal /= (x_cells*y_cells);
            

            //fprintf(normConc,  "%4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f\n", 0.000, normalArea, concentration, max((TotIce / Domain), 0.000), (max(GlobalFine,0.000)), SurfArea, thick_tot, TaveGlobal, thick_coarse_ave, thick_fine_ave);
            fprintf(normConc,  "%4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f\n", 0.000, normalArea, initial_coarse + initial_fine, initial_coarse, initial_fine, SurfArea, thick_tot, TaveGlobal, thick_coarse_ave, thick_fine_ave, aveFlux0, aveFlux20, avelatFlux0, avevertFlux0, avelvRatio0); //Just for initial data
            fprintf(normArea,  "%4.8f %4.8f %4.8f\n", 0.000, normalArea, concentration);
            fprintf(fluxArea,  "%4.8f %4.8f\n", 0.000, fluxAccArea);
            fprintf(lossArea,  "%4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f\n", Domain, initial_coarse*Domain, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); 
            cout << "Normarea= " << normalArea << " time= " << 0 << endl;
            cout << "Fluxarea= " << fluxAccArea << " time= " << 0  << endl;
            cout << "Concentration= " << concentration << " time= " << 0  << endl;
            cout << "Average Thickness= " << thick_tot << " time= " << 0  << endl;
            prevConc = concentration;
            prevConcB = max( (TotIce / Domain) , 0.0000);
            prevConcF = max(GlobalFine,0.000);
            prevNormArea = normalArea;
            prevSurfArea = SurfArea;
            prevThickF = thick_fine_ave; 
            prevThickC = thick_coarse_ave;

        }

        //Initialize to Normalize Area
        if (step > tiso && step % nTout == 0 && step > 1){
            
            double SurfArea = 0;
            for (size_t i = 0; i < grainList.size(); i++) {

                double SurfAreap = 0;
                for (size_t j = 0; j < grainList[i].getnpoints()-1; j++) {
                    SurfAreap += (-grainList[i].getPointList()[j]+grainList[i].getPointList()[j+1]).norm();
                }
                SurfAreap += (grainList[i].getPointList()[grainList[i].getnpoints()-1]-grainList[i].getPointList()[0]).norm();
                SurfArea += SurfAreap;

            }
            
            double thick_coarse2 = 0.0;
            double normalArea = 0;
            fluxAccArea = 0.0;
            for (size_t i = 0; i < grainList.size(); i++) {
                //AccArea += grainList[i].getMass();
                vector<Vector2d> forArea00 = grainList[i].getPointList();
                //addArea = PointsAreaM(forArea00);
                
                double Area = 0;
                size_t n = forArea00.size();

                for (size_t pp = 0; pp < n-1; pp++)
                {
                    Area = Area + ( forArea00[pp](0) * forArea00[pp+1](1) -  forArea00[pp](1) * forArea00[pp+1](0) ); 
                }
                Area = Area + (forArea00[n-1](0) * forArea00[0](1) -  forArea00[n-1](1) * forArea00[0](0) ); 

                addArea = 0.5*abs(Area);

                AccArea = AccArea + addArea;
                //cout << "Grain Area from Mass: " << grainList[i].getMass() << endl;

                //thick_coarse2 +=  grainList[i].getThickness();  //Get average thickness //Not weighted
                thick_coarse2 +=  addArea * grainList[i].getThickness();  //Get average thickness //Weighted


                if (grainList[i].getPosition()[1] < -grainList[i].getPosition()[0] + 6500 ) // Less is down for Y //6500 depends on geometry of BCs
                {
                    fluxAccArea += grainList[i].getMass();
                }    
            }
            
            double thick_fines2 = 0.0; //Area sum
            for (size_t i = 0; i < Fine_Grains.size(); i++)
            {
               //thick_fines2 += Fine_Grains[i](1); //Add all Areas (original and new)  //Not weighted
               thick_fines2 += Fine_Grains[i](0) * Fine_Grains[i](1); //Add all Areas (original and new)  //Weighted
            }
                        
            //cout << "AccArea0= " << AccArea0 << endl;
            normalArea = AccArea/AccArea0;
            //Print time step and Accumulated Area
            double timeplot = step/nTout;
            
            for (size_t i = 0; i < grainList.size(); i++) {
                vector<Vector2d> forArea = grainList[i].getPointList();
                //addArea = PointsAreaM(forArea);
                
                double Area = 0;
                size_t n = forArea.size();

                for (size_t pp = 0; pp < n-1; pp++)
                {
                    Area = Area + ( forArea[pp](0) * forArea[pp+1](1) -  forArea[pp](1) * forArea[pp+1](0) ); 
                }
                Area = Area + (forArea[n-1](0) * forArea[0](1) -  forArea[n-1](1) * forArea[0](0) ); 

                addArea = 0.5*abs(Area);

                TotIce = TotIce + addArea;
                //cout << "For points: " << endl;
                for (size_t iidx = 0; iidx < forArea.size(); iidx++)
                {
                    //cout << forArea[iidx](0)  << " " << forArea[iidx](1)  << endl;
                }
                //cout << "Grain Area from Points: " << addArea << endl;
            }  
            double concentration = max((TotIce / Domain), 0.000) + (max(GlobalFine,0.000));
            
            //To avoid problems due to point interpolation
            double ConcB, ConcF;
            
            ConcB = max ((TotIce / Domain), 0.000);
            ConcF = max(GlobalFine,0.000);
            if (isnan(concentration) == 1)
            {
                concentration = prevConc;
            }
            if (isnan(ConcB) == 1)
            {
                cout << "Nan Value for Coarse Concentration" << endl;
                ConcB = prevConcB;
            }
            if (isnan(ConcF) == 1)
            {
                cout << "Nan Value for Fine Concentration" << endl;
                ConcF = prevConcF;
            }
            if (isnan(normalArea) == 1)
            {
                cout << "Nan Value for Area" << endl;
                normalArea = prevNormArea;
            }
            if (isnan(SurfArea) == 1)
            {
                SurfArea = prevSurfArea;
            }

            //Average Surface Flux
            double aveFlux, aveFlux2;
            aveFlux = 0.0;
            aveFlux2 = 0.0;
            double avelatFlux = 0.0;
            double avevertFlux = 0.0;
            double avelvRatio = 0.0;
            
            //Average Global Temperature in Domain
            TaveGlobal = 0.0;
            for (size_t i = 0; i < y_cells; i++) {
                for (size_t j = 0; j < x_cells; j++) {
                    TaveGlobal += oceanTemp[j+i*x_cells];
                    aveFlux += netFlux[j+i*x_cells];
                    aveFlux2 += netFlux2[j+i*x_cells];
                    //Only exists at center cells, no borders
                    if ( (i > 0 && i < y_cells - 1) && (j > 0 && j < x_cells - 1) ){
                        avelatFlux += latFlux[j+i*x_cells];
                        avevertFlux += vertFlux[j+i*x_cells];
                        avelvRatio +=  lvRatio[j+i*x_cells];
                    }
                }
            }
            aveFlux /= (x_cells*y_cells);
            aveFlux2 /= (x_cells*y_cells);
            avelatFlux /= ((x_cells-2)*(y_cells-2));
            avevertFlux /= ((x_cells-2)*(y_cells-2));
            avelvRatio /= ((x_cells-2)*(y_cells-2));
            TaveGlobal /= (x_cells*y_cells);

            //Thickness output
            //double thick_tot2 = (thick_coarse2 + thick_fines2) / (double (grainList.size() +  Fine_Grains.size())); //By number
            double thick_tot2 = (thick_coarse2 + thick_fines2) / (concentration*Domain); //By area
            double thick_coarse_ave2 = (thick_coarse2) / (concentration*Domain); //By total area to observe separately
            double thick_fine_ave2 = (thick_fines2) / (concentration*Domain); //By total area to observe separately

            if (isnan(thick_coarse_ave2) == 1)
            {
                cout << "Nan Value for Coarse Thick" << endl;
                thick_coarse_ave2 = prevThickC;
            }
            if (isnan(thick_fine_ave2) == 1)
            {
                cout << "Nan Value for Fine Thick" << endl;
                thick_fine_ave2 = prevThickF;
            }
            if (isnan(thick_tot2) == 1)
            { 
                thick_tot2 = thick_fine_ave2 + thick_coarse_ave2;
            }

            
            //fprintf(normConc,  "%4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f\n", timeplot, normalArea, concentration, ConcB, ConcF, SurfArea, thick_tot2, TaveGlobal, thick_coarse_ave2, thick_fine_ave2);
            fprintf(normConc,  "%4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f\n", timeplot, normalArea, concentration, min(ConcB, initial_coarse), ConcF, SurfArea, thick_tot2, TaveGlobal, thick_coarse_ave2, thick_fine_ave2, aveFlux, aveFlux2, avelatFlux, avevertFlux, avelvRatio); //Allow declining thread
            fprintf(normArea,  "%4.8f %4.8f %4.8f\n", timeplot, normalArea, concentration);
            fprintf(fluxArea,  "%4.8f %4.8f\n", timeplot, fluxAccArea);
            //For validation
            area_loss_lat2 = min(ConcB, initial_coarse)*Domain - area_loss_bkg - area_loss_basal;  //AreaDomain, AreaCoarse, AreaLossBkg, AreaLossBasal, AreaLossLat, AreaLossLat2, AreaLossCoarse
            //ave_delta_r_rate_step /=  (double(nTout) / double(nTemp)); //Average over all time steps for each capture considering it is done every nTemp
            ave_delta_r_rate_step /=  86400;
            fprintf(lossArea,  "%4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f\n", Domain, min(ConcB, initial_coarse)*Domain, area_loss_bkg, area_loss_basal, area_loss_lat, area_loss_lat2, min(ConcB, initial_coarse)*Domain - prev_coarse_area, ave_delta_r_rate_step); 
            ave_delta_r_rate_step = 0.0; //Restart to add again and get new rate
            prev_coarse_area = min(ConcB, initial_coarse)*Domain;
            cout << "Normarea= " << normalArea << " time= " << timeplot  << endl;
            cout << "Concentration= " << concentration << " time= " << timeplot  << endl;
            //cout << "Fluxarea= " << fluxAccArea << " time= " << timeplot  << endl;
            prevConc = concentration;
            prevNormArea = normalArea;
            prevSurfArea = SurfArea;
            prevConcB = ConcB;
            prevConcF = ConcF;
            prevThickC = thick_coarse_ave2; 
            prevThickF = thick_fine_ave2;
        } 

        //Obtain Grain Size Distribution
        //Switched to Main Caliper Diameter or Equivalent Diameter to Grains Area
        // Area = 3.141592653589793 * D^2 / 4          Mean D = 2 * sqrt ( Area / 3.141592653589793 )
        if (rank == 0)
        {    
            //if (step == 1){
            //cout << "Start GSD Output" << endl;
            if (step > 0 && (step == 1 || step % nTout == 0)){
                double TotMass = 0.0;
                double TotArea = 0.0;
                for (size_t i = 0; i < grainList.size(); i++) {
                    vector<Vector2d> forArea2 = grainList[i].getRefPointList();


                    //TotMass += grainList[i].getMass();
                    //addArea = PointsAreaM(forArea2);

                    double Area = 0;
                    size_t n = forArea2.size();

                    for (size_t pp = 0; pp < n-1; pp++)
                    {
                        Area = Area + ( forArea2[pp](0) * forArea2[pp+1](1) -  forArea2[pp](1) * forArea2[pp+1](0) ); 
                    }
                    Area = Area + (forArea2[n-1](0) * forArea2[0](1) -  forArea2[n-1](1) * forArea2[0](0) ); 

                    addArea = 0.5*abs(Area);

                    TotMass = TotMass + addArea;                    
                    TotArea = TotArea + addArea;
                } 
                //cout << "Total Mass: " << TotMass << endl;
                //cout << "Total Area: " << TotArea << endl;

                double MaxD = 0.0;
                double MinD = 40000000000000000.00000;
                size_t max_index;
                for (size_t i = 0; i < grainList.size(); i++) {
                    vector<Vector2d> forArea3 = grainList[i].getRefPointList();
                    //addArea = PointsAreaM(forArea3);


                    double Area = 0;
                    size_t n = forArea3.size();

                    for (size_t pp = 0; pp < n-1; pp++)
                    {
                        Area = Area + ( forArea3[pp](0) * forArea3[pp+1](1) -  forArea3[pp](1) * forArea3[pp+1](0) ); 
                    }
                    Area = Area + (forArea3[n-1](0) * forArea3[0](1) -  forArea3[n-1](1) * forArea3[0](0) ); 

                    addArea = 0.5*abs(Area);


                    //cout << "Size grain i: " << i << " is: " << 2 * sqrt( addArea /3.141592653589793) << endl;
                    //if (2.0*grainList[i].getRadius() > MaxD)
                    //if ( 2 * sqrt( grainList[i].getMass() /3.141592653589793)  > MaxD ) 
                    if ( 2 * sqrt( addArea / 3.141592653589793)  > MaxD ) 
                    {
                        //MaxD = 2.0*grainList[i].getRadius();
                        //MaxD = 2 * sqrt( grainList[i].getMass() /3.141592653589793);
                        MaxD = 2 * sqrt( addArea / 3.141592653589793);
                        max_index = i;
                    }
                    //if (2.0*grainList[i].getRadius() < MinD)
                    //if ( 2 * sqrt( grainList[i].getMass() /3.141592653589793)  < MinD ) 
                    if ( 2 * sqrt( addArea / 3.141592653589793)  < MinD ) 
                    {
                        //MinD = 2.0*grainList[i].getRadius();
                        //MinD = 2 * sqrt( grainList[i].getMass() /3.141592653589793);
                        MinD = 2 * sqrt( addArea / 3.141592653589793);
                    }
                }
                
                //Automatic
                //double DistRange = MaxD - MinD;
                
                //To keep standard GSD bins fix Max and MinD
                if (step == 1){  //Only changes here and is retained for the rest of the simulation
                    //GlobalMaxD = MaxD;
                    GlobalMaxD = MAXDD;
                }
                MaxD = GlobalMaxD;
                MinD = MINDD;  //10             //Choose very small final range since MinD will get smaller with melting
                double DistRange = MaxD - MinD; //Thus DistRange is fixed and so Diameters vector as well, only count will change
                
                //cout << "MinD/MaxD: " << MinD << " / " << MaxD << endl;
                //cout << "max_index: " << max_index << endl;
                //cout << "sqrt of 4: " << sqrt(4) << endl;

                //Bin size and Diameters defined at the beginning
                vector<double> Passes(bin_size+1);
                vector<double> Count(bin_size+1);
                vector<double> BinConc(bin_size+1);

                //Define Diameters (NOTE: Descending order)
                for (size_t i = 0; i < bin_size+1; i++)
                {
                    Diameters[i] = MaxD - (double(i)/double(bin_size))*DistRange;
                }

                //NEW MODIF
                double Area_fines = 0.0;
                for (size_t jf = 0; jf < Fine_Grains.size(); jf++) {
                    Area_fines += Fine_Grains[jf](0);
                }

                //Define Percent Passes                
                for (size_t i = 0; i < bin_size+1; i++)
                {
                    double MassSum = 0.0; //Temporal Mass for Adding at %Passes
                    for (size_t j = 0; j < grainList.size(); j++) {
                        //if (2.0*grainList[j].getRadius() <= Diameters[i]){
                        //if ( 2 * sqrt( grainList[j].getMass() /3.141592653589793) <= Diameters[i]){
                        const vector<Vector2d> forArea4 = grainList[j].getPointList();
                        //addArea = PointsAreaM(forArea4);
                        double Area = 0;
                        size_t n = forArea4.size();

                        for (size_t pp = 0; pp < n-1; pp++)
                        {
                            Area = Area + ( forArea4[pp](0) * forArea4[pp+1](1) -  forArea4[pp](1) * forArea4[pp+1](0) ); 
                        }
                        Area = Area + (forArea4[n-1](0) * forArea4[0](1) -  forArea4[n-1](1) * forArea4[0](0) ); 

                        addArea = 0.5*abs(Area);
                        
                        //Potential additional validation
                        //if (addArea < M_PI * pow(grainList[j].getRadius(),2) ){}

                        //cout << "Size grain i: " << j << " is: " << 2 * sqrt( addArea /3.141592653589793) << endl;
                        if ( 2 * sqrt( addArea /3.141592653589793) <= Diameters[i]){
                            //MassSum += grainList[j].getMass();
                            MassSum = MassSum + addArea ;
                        }
                        else if( (2 * sqrt( addArea /3.141592653589793) > Diameters[i]) && (i == 0) ) {  //If exceeding just handled as maxD threshold
                            MassSum = MassSum + addArea ;
                        } 
                    }
                    
                    //Older
                    //Passes[i] = MassSum*100.0/TotMass;
                    //NEW MODIF
                    if (i == 0)
                    {
                        MassSum = TotArea; //To avoid non 100 percent values
                    }
                    
                    if (TotArea > 0.0){
                        Passes[i] = ((MassSum)*100.0)/(TotArea); 
                    }
                    else{
                        Passes[i] = 0.000; 
                    }
                    //Passes[i] = ((MassSum + Area_fines)*100.0)/(TotArea + Area_fines);  //Mass Coarse below diam + Mass Fine (mass = area)  / Total Sea Ice Area
                    
                    cout << "FOR GSD DEBUGGING" << endl;
                    cout << "Mass Sum: " << MassSum << endl;
                    cout << "Tot Area Coarse: " << TotArea << endl;
                    cout << "Tot Area Fines: " << Area_fines << endl;
                }  

                //Define Count
                size_t count_size;
                for (size_t i = 0; i < bin_size; i++)
                {
                    count_size = 0;
                    double MassSumConc = 0.0; //Temporal Mass for Adding at %Concentration in a bin
                    for (size_t j = 0; j < grainList.size(); j++) {
                        //if ( 2.0*grainList[j].getRadius() <= Diameters[i] &&  2.0*grainList[j].getRadius() >= Diameters[i+1] ) {
                        //if ( 2 * sqrt( grainList[j].getMass() /3.141592653589793) <= Diameters[i] &&  2 * sqrt( grainList[j].getMass() /3.141592653589793) >= Diameters[i+1] ) {
                        vector<Vector2d> forArea5 = grainList[j].getPointList();
                        //addArea = PointsAreaM(forArea5);
                        
                        double Area = 0;
                        size_t n = forArea5.size();

                        for (size_t pp = 0; pp < n-1; pp++)
                        {
                            Area = Area + ( forArea5[pp](0) * forArea5[pp+1](1) -  forArea5[pp](1) * forArea5[pp+1](0) ); 
                        }
                        Area = Area + (forArea5[n-1](0) * forArea5[0](1) -  forArea5[n-1](1) * forArea5[0](0) ); 

                        addArea = 0.5*abs(Area);

                        if ( 2 * sqrt( addArea  /3.141592653589793) <= Diameters[i] &&  2 * sqrt( addArea  /3.141592653589793) >= Diameters[i+1] ) {
                            count_size++;
                            MassSumConc = MassSumConc + addArea ;
                        }
                    }
                    BinConc[i] = (MassSumConc * 100.0) / (Domain);
                    Count[i] = count_size;
                }  
                Count[bin_size] = 0;
                BinConc.back() = (GlobalFine)*100;

                //Print Output to File
                string fnamegsd;
                if (step == 1) {
                    fnamegsd = outDir3 + "GrainSizeDist_iter" + std::to_string(0) + "_0.dat";
                }
                else {
                    fnamegsd = outDir3 + "GrainSizeDist_iter" + std::to_string(step/nTout) + "_0.dat";
                }
                FILE * docgsd = fopen(fnamegsd.c_str(), "w");
                
                for (size_t bi = 0; bi < MeltBin.size(); bi++)
                {
                    LatRate_step[bi] = LatRate[bi];
                }
                
                for (size_t i = 0; i < bin_size+1; i++) {
                    //Purge any Nan value for export
                    if (isnan(GMeltBin[i]))
                    {
                        GMeltBin[i] = 0.0000;
                    }
                    if (isnan(GBreakBin[i]))
                    {
                        GBreakBin[i] = 0.0000;
                    }
                    
                    for (size_t i = 0; i < MeltBin.size(); i++)
                    {
                        LatRate_step[i] /= double(1); //Calculate sum of average over floes km lost per bin which is in itself km lost per day at these snapshots. Then you clear to find new length lost in the next day. Hence divide by 1.
                    }
                    
                    //FSD specific counts (April 24, 2023)
                    fprintf(docgsd, "%4.8f %4.8f %d %4.8f %4.8f %4.8f %d %d %d %d %8.16f\n", Diameters[i], Passes[i], (int) Count[i], BinConc[i], GMeltBin[i], GBreakBin[i], (int) NLMelt[i], (int) NLBreak[i], (int) NGMelt[i], (int) NGBreak[i], LatRate_step[i] ); // GSD Print
                } 
                fclose(docgsd);
                
                //Clear for next snapshot
                for (size_t i = 0; i < NLMelt.size(); i++)
                {
                    //FSD specific counts (April 24, 2023)
                    NLMelt[i] = 0;
                    NLBreak[i] = 0;
                    NGMelt[i] = 0;
                    NGBreak[i] = 0;
                    LatRate_step[i] = 0.0;
                    LatRate[i] = 0.0;
                }
            
                
            }
            //cout << "End GSD Output" << endl;
        }    

        if (rank == 0)
        {    

            // Contact output
            if (step % nTout == 0) {
                
                // Contact output (grain-grain) 
                CData cState = world.computeCstate();
                //string fname = outDir + "cinfo_iter" + std::to_string(t/tOut) + "_" + std::to_string(rank) + ".dat";
                string fname = outDir2 + "cinfo_iter" + std::to_string(step/nTout) + "_0.dat";
                FILE * cinfo = fopen(fname.c_str(), "w");
                for (size_t i = 0; i < cState._clocs.size(); i++) {
                    fprintf(cinfo, "%d %d ", cState._cpairs[i](0), cState._cpairs[i](1) ); // grains
                    fprintf(cinfo, "%.8f %.8f ",cState._forces[i](0), cState._forces[i](1));// f
                    fprintf(cinfo, "%.8f %.8f ",cState._normals[i](0), cState._normals[i](1));// n
                    fprintf(cinfo, "%.8f %.8f\n",cState._clocs[i](0), cState._clocs[i](1));// loc
                }
                fclose(cinfo);
            
                // // Contact output (with GrainWall)
                // CData cStateGWalls = world.computeCstateGWalls();
                // string fname2 = outDir2 + "cgwallinfo_iter" + std::to_string(step/nTout) + "_0.dat";
                // FILE * cinfo2 = fopen(fname2.c_str(), "w");
                // for (size_t i = 0; i < cStateGWalls._clocs.size(); i++) {
                //     fprintf(cinfo2, "%d %d ", cStateGWalls._cpairs[i](0), cStateGWalls._cpairs[i](1) ); // grains
                //     fprintf(cinfo2, "%.8f %.8f ",cStateGWalls._forces[i](0), cStateGWalls._forces[i](1));// f
                //     fprintf(cinfo2, "%.8f %.8f ",cStateGWalls._normals[i](0), cStateGWalls._normals[i](1));// n
                //     fprintf(cinfo2, "%.8f %.8f\n",cStateGWalls._clocs[i](0), cStateGWalls._clocs[i](1));// loc
                // }
                // fclose(cinfo2);
            }

 
            //Fine Data export for Debugging
            if (step % nTout == 0)
            {
                //cout << "*** FINE RESULTS FOR: " << step/nTout << " percent for Nfines: " << Fine_Grains.size() << " ***" << endl;
                for (size_t i = 0; i < Fine_Grains.size(); i++)
                {
                    //cout << "Fine Grain #: " << i << " Area: " << Fine_Grains[i](0) << " Thick: " << Fine_Grains[i](1) << " X: " << Fine_Grains[i](2) << " Y: " << Fine_Grains[i](3) << endl;
                }
            }
            

            //cout << "Start Output" << endl;
            // Get output for compression
            if (step % nTout == 0)
            {
                //Aug-18-2022 Change
                //BEGIN
                double time_plot = double(step)/double(nTout);
                //END
                
                // Compute kinetic energy
                grainList = world.getGrains();
                grainList2 = world.getGrainState(); //For Moment and Force Output
                
                //START CHANGE OCT-18-2022
                if (step == 0)
                {
                    cout << "Time Failure update" << endl;
                    double time_fail = double(BREAK_PROB);
                    world.changeGrainTfail(time_fail);
                    // cout << "Failure times!!! for grain number: " << grainList.size() << endl;
                    // //Initialize Random failure times (normal distribution, inherit failure time and origin time to new floes for simplicity)
                    // double mean_failure = double(BREAK_PROB);
                    // double temp_sample_fail;
                    // double sample_fail;
                    // const double zero_fail = 0.000;
                    // std::default_random_engine time_failure;
                    // std::normal_distribution<double> distribution_fail(mean_failure, 0.1 * mean_failure); // mean, 0.1 std. dev mean
                    
                    // for (size_t ik = 0; ik < grainList.size(); ik++)
                    // {
                    //     temp_sample_fail = max( distribution_fail(time_failure) , 0.000 );
                    //     sample_fail = temp_sample_fail;
                    //     std::cout << "Sample fail for grain: " << ik << " is: " << sample_fail << std::endl;
                    //     grainList[ik].changeTfail(sample_fail);
                    //     grainList[ik].changeTorigin(zero_fail);
                    //     cout << "New grain fail time: " << grainList[ik].getTfail() << endl;
                    // }
                }
                //END CHANGE OCT-18-2022

                Ke = 0;
                for (size_t i = 0; i < grainList.size(); i++) {
                    Ke += grainList[i].computeKineticEnergy();

                    //cout << "Grain: " << i << " with position X: " << grainList[i].getPosition()(0) << " and Y: " << grainList[i].getPosition()(1) << endl;
                    //cout << "Grain: " << i << " with Mass: " << grainList[i].getMass() << " and Inertia: " << grainList[i].getMoI() << endl;
                }
                // Print positions and velocities to file // and temperature and thickness
                //cout << "Save Data" << endl;
                for (size_t i = 0; i < grainList.size(); i++) {
                    fprintf(positions,  "%4.8f %4.8f %4.8f\n", grainList[i].getPosition()(0),
                            grainList[i].getPosition()(1), grainList[i].getTheta());

                    //cout << "output position over" << endl;
                    fprintf(velocities, "%4.8f %4.8f %4.8f\n", grainList[i].getVelocity()(0),
                            grainList[i].getVelocity()(1), grainList[i].getOmega());



                    fprintf(temper, "%4.8f \n", grainList[i].getTemperature());


                    fprintf(thick, "%4.12f \n", grainList[i].getThickness());


                }

                //Print fine Positions for Debugging
                size_t fines_count = 0;
                for (size_t i = 0; i < Fine_Grains.size(); i++) {
                    if (Fine_Grains[i](1) > 0.000)
                    {
                        fprintf(finepositions,  "%4.8f %4.8f \n", Fine_Grains[i](2), Fine_Grains[i](3)); 
                        fines_count++;
                    }
                }

                //Print Number of fines for reference
                double nfines = fines_count;
                //double nfines = Fine_Grains.size();
                fprintf(numberfines, "%4.8f \n", nfines);  
                
                //Aug-18-2022 Change
                //BEGIN
                //Calculate Total Mass Coarse (Using Mass) 
                tmasscoarse = 0.0;
                for (size_t i = 0; i < grainList.size(); i++) {
                    if (grainList[i].getMass() > 0.000)
                    {
                        tmasscoarse += grainList[i].getMass(); //In kg = kg/km3 * km3  //No need for  grainsList[i].getThickness0() * 0.001 because it's already calculated at the beginning and updated as thickness changes or area changes //_mass already uses density rho
                    }
                }
                
                
                //Calculate Total Mass Coarse (Using Polygon Area) 
                string fnameThickC;
                
                size_t max_size_index = 0;
                double max_area = 0.00;
                
                fnameThickC = outDir4 + "TSD_Coarse_step_" + std::to_string(step/nTout) + ".dat";    
                FILE * docTSDC = fopen(fnameThickC.c_str(), "w");
                tmasscoarseA = 0.0;
                for (size_t i = 0; i < grainList.size(); i++) {
                    if (grainList[i].getMass() > 0.000)
                    {
                        vector<Vector2d> ptArea = grainList[i].getPointList();

                        double tempArea = 0.0;
                        size_t ntemp = ptArea.size();

                        for (size_t pp = 0; pp < ntemp-1; pp++)
                        {
                            tempArea = tempArea + ( ptArea[pp](0) * ptArea[pp+1](1) -  ptArea[pp](1) * ptArea[pp+1](0) ); 
                        }
                        tempArea = tempArea + (ptArea[ntemp-1](0) * ptArea[0](1) -  ptArea[ntemp-1](1) * ptArea[0](0) ); 

                        tempArea = 0.5*abs(tempArea);

                        if (tempArea > max_area){
                            max_area = tempArea;
                            max_size_index = i;
                        }
                        
                        //For TSD
                        fprintf(docTSDC,  "%4.8f %4.8f\n", tempArea, grainList[i].getThickness() ); 
                        
                        tmasscoarseA += tempArea * grainList[i].getThickness() * 0.001 * rho; //In kg = kg/km3 * km3  //We need grainsList[i].getThickness0() * 0.001 because points only give area in km2 //Needs density rho
                    }
                }
                fclose(docTSDC);
                
                
                //Calculate Total Mass Fines (Usinge Fine_Grains)
                tmassfines = 0.0;
                for (size_t i = 0; i < Fine_Grains.size(); i++) {
                    if (Fine_Grains[i](1) > 0.000)
                    {
                        tmassfines += Fine_Grains[i](0) * Fine_Grains[i](1) * 0.001 * rho; //In kg = kg/km3 * km3 //Doesn't have already uses density rho
                    }
                }

                //Calculate Total Mass (Coarse + Fines) (using Mass Cooarse)
                tmass = tmasscoarse + tmassfines;
                
                //Calculate Total Mass (Coarse + Fines) (Using Polygon Area) 
                tmassA = tmasscoarseA + tmassfines;
                
                //TEMPORARY (DO NOT CALCULATE LOSSED UNTIL CHECKING BEST METHOD FOR MASS) (ASSUME 0 Before assigning to TempModif and Meltys)
                //loss_mcl = 0.0; loss_mcv = 0.0; gain_fines = 0.0; loss_fines; //loss_mfv = 0.0;   //Change Aug 22, 2022
                
                
                //double total_melt_loss = loss_mcl + loss_mcv + loss_mfv; // Total melt loss for coarse and fines  //Can deduce melt fines from tmassfines
                double total_melt_coarse = loss_mcl + loss_mcv; // Total melt loss for coarse
                
                //Print Mass Loss information
                //fprintf(massLoss, "%4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f\n", time_plot, tmass, tmasscoarse, tmassfines, loss_mcl, loss_mcv, gain_fines, loss_fines, total_melt_coarse, tmassA, tmasscoarseA);   //Change Aug 22, 2022
                //NEW top/down (solar/ ocean decomposition)
                fprintf(massLoss, "%4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %8.20f %4.8f\n", time_plot, tmass, tmasscoarse, tmassfines, loss_mcl, loss_mcv, gain_fines, loss_fines, total_melt_coarse, tmassA, tmasscoarseA, loss_mcv_solar, loss_mcv_ocean, MLAT, MBASAL);   //Change Aug 22, 2022
                cout << "S melt component total: " << loss_mcv_solar << endl;
                cout << "O melt component total: " << loss_mcv_ocean << endl;
                //END
                
                //How to input specific wall values [??]
                fprintf(wallPosT,"%4.8f %4.8f\n", world.getWalls()[3].getPosition()(0), world.getWalls()[1].getPosition()(1));

                //Print to Control Level Set Evolution Output (Print Left to Right (X), Up to Down (Y))
                size_t iii = 0;  //Careful, depends on grain number, if < will generate error
                //size_t iii = 35;  //All of these for particle 0 or original
                // size_t xxdim =  grainList[iii].getLset().getXdim();  //To get Geometric LSET -->
                // size_t yydim =  grainList[iii].getLset().getYdim();
                // double xgrid =  grainList[iii].getLset().getXdim();
                // double ygrid =  grainList[iii].getLset().getYdim();

                size_t xxdim =  grainList[iii].getgrainTemp2D().getXdim();  // To get Heat LSET
                size_t yydim =  grainList[iii].getgrainTemp2D().getYdim();
                double xgrid =  grainList[iii].getgrainTemp2D().getXdim();
                double ygrid =  grainList[iii].getgrainTemp2D().getYdim();

                // size_t xxdim =  grainList[iii].getLset().getXdim();  //To get Stress LSET -->
                // size_t yydim =  grainList[iii].getLset().getYdim();
                // double xgrid =  grainList[iii].getLset().getXdim();
                // double ygrid =  grainList[iii].getLset().getYdim();

                // size_t xxdim =  grainList[iii].getLset().getXdim();  //To get Damage LSET -->
                // size_t yydim =  grainList[iii].getLset().getYdim();
                // double xgrid =  grainList[iii].getLset().getXdim();
                // double ygrid =  grainList[iii].getLset().getYdim();
     
                //Print level set for study
                if (step==0)
                {
                   fprintf(sampleLS, "%4.8f \n", xgrid);
                   fprintf(sampleLS, "%4.8f \n", ygrid);
                }    
                
                //Study Level Set
                for (size_t jj = 0; jj < yydim; jj++){
                    for (size_t ii = 0; ii < xxdim; ii++){ 
                        fprintf(sampleLS, "%4.8f \n", grainList[iii].getMthick().getLevelset()[jj*xxdim + ii]); //Thickness LS
                        //fprintf(sampleLS, "%4.8f \n", grainList[iii].getgrainTemp2D().getLevelset()[jj*xxdim + ii]);    //Temp LS -->
                        //fprintf(sampleLS, "%4.8f \n", grainList[iii].getLset().getLevelset()[jj*xxdim + ii]);  //Geo LS -->
                    }     
                }
                //Point Set for Output Generation for Level Set
                for (size_t j = 0; j<grainList[iii].getnpoints(); j++) {

                        fprintf(pointsoutLS,  "%4.8f %4.8f \n", grainList[iii].getPointList()[j](0), grainList[iii].getPointList()[j](1) );
                }
                
                
            //    //TURN ON FOR MORE INFO
            //    //August 26, 2022
            //     //Print Thickness Level Set Snapshots for Thickness Analysis
            //     //Print Output to File
            //     string fnameT; 
            //     size_t xxdimT, yydimT;
            //     double xgridT, ygridT;
            //     //if (step / nTout == 0 || step / nTout == 4 || step / nTout == 15 || step / nTout == 22 || step / nTout == 32){ //Just for specific steps
            //     //if (step % nTout == 0 || step % nTout == 4 || step % nTout == 15 || step % nTout == 22 || step % nTout == 32){ //All inclusive (all steps)
            //     if (0 > 1){  //All Exclusive to not print (not print)
            //         for (size_t fidx = 0; fidx < grainList.size(); fidx++){
            //             fnameT = outDir4 + "Thickness_" + std::to_string(step/nTout)  + "_" + std::to_string(fidx) + "_0.dat";    
            //             FILE * docT = fopen(fnameT.c_str(), "w");
                        
            //             xxdimT =  grainList[fidx].getMthick().getXdim();  // To get THICKNESS LSET
            //             yydimT =  grainList[fidx].getMthick().getYdim();
            //             xgridT =  grainList[fidx].getMthick().getXdim();
            //             ygridT =  grainList[fidx].getMthick().getYdim();
             
            //             //BBRadius and Dimensions
            //             fprintf(docT, "%4.8f\n", grainList[fidx].getRadius());
            //             fprintf(docT, "%4.8f\n", xgridT);
            //             fprintf(docT, "%4.8f\n", ygridT);
            //             for (size_t jj = 0; jj < yydimT; jj++){
            //                 for (size_t ii = 0; ii < xxdimT; ii++){ 
            //                     fprintf(docT, "%4.8f\n", grainList[fidx].getMthick().getLevelset()[jj*xxdimT + ii] ); //Thickness LS
            //                 }     
            //             }
                        
            //             fclose(docT);
            //         }
            //     }    
                
            //     //Print Temperature, Geometric and Thickness Level Set for a Sample Grain
                
            //     //Output grain properties for study
            //     size_t grain_index;
            //     if (grainList.size() == 493){
            //         grain_index = 342;
            //     }
            //     else{
            //         //grain_index = 0;
            //         grain_index = max_size_index;
            //     }
            //     string fnameThick, fnameGeo, fnameTemp, fnamePos, fnameCM, fnamePoints; 
            //     fnameThick = outDir4 + "ThicknessLS_step_" + std::to_string(step/nTout)  + "_g_" + std::to_string(grain_index) + ".dat";    
            //     FILE * docThick = fopen(fnameThick.c_str(), "w");
            //     fnameGeo = outDir4 + "GeoLS_step_" + std::to_string(step/nTout)  + "_g_" + std::to_string(grain_index) + ".dat";    
            //     FILE * docGeo = fopen(fnameGeo.c_str(), "w");
            //     fnameTemp = outDir4 + "TempLS_step_" + std::to_string(step/nTout)  + "_g_" + std::to_string(grain_index) + ".dat";    
            //     FILE * docTemp = fopen(fnameTemp.c_str(), "w"); 
            //     fnamePos = outDir4 + "Pos_step_" + std::to_string(step/nTout)  + "_g_" + std::to_string(grain_index) + ".dat";    
            //     FILE * docPos = fopen(fnamePos.c_str(), "w"); 
            //     fnameCM = outDir4 + "CenterMass_step_" + std::to_string(step/nTout)  + "_g_" + std::to_string(grain_index) + ".dat";    
            //     FILE * docCM = fopen(fnameCM.c_str(), "w"); 
            //     fnamePoints = outDir4 + "GPoints_step_" + std::to_string(step/nTout)  + "_g_" + std::to_string(grain_index) + ".dat";    
            //     FILE * docPoints = fopen(fnamePoints.c_str(), "w"); 
                
            //     size_t xxdimLS =  grainList[grain_index].getgrainTemp2D().getXdim();  // Useful dims
            //     size_t yydimLS =  grainList[grain_index].getgrainTemp2D().getYdim();
            //     double xgridLS =  grainList[grain_index].getgrainTemp2D().getXdim();
            //     double ygridLS =  grainList[grain_index].getgrainTemp2D().getYdim();
                
            //     //Extra variables
            //     double lowLx, lowLy, lslocx, lslocy;
            //     double szx = xgridLS;
            //     double szy = ygridLS;
            //     Vector2d pointOcxy;
            //     //Shift to global coordinates and move to left bottom corner of level set
	           // lowLx = grainList[grain_index].getPosition()(0) - 0.5*szx +1.0; //IFF LS has unit of 1 for grid cells
	           // lowLy = grainList[grain_index].getPosition()(1) - 0.5*szy +1.0; //IFF LS has unit of 1 for grid cells
                
                
            //     //Properties
            //     fprintf(docThick, "%4.8f\n", grainList[grain_index].getRadius());
            //     fprintf(docThick, "%4.8f\n", xgridLS);
            //     fprintf(docThick, "%4.8f\n", ygridLS);
            //     fprintf(docGeo, "%4.8f\n", grainList[grain_index].getRadius());
            //     fprintf(docGeo, "%4.8f\n", xgridLS);
            //     fprintf(docGeo, "%4.8f\n", ygridLS);
            //     fprintf(docTemp, "%4.8f\n", grainList[grain_index].getRadius());
            //     fprintf(docTemp, "%4.8f\n", xgridLS);
            //     fprintf(docTemp, "%4.8f\n", ygridLS);
            //     fprintf(docPos, "%4.8f %4.8f\n", grainList[grain_index].getPosition()(0), grainList[grain_index].getPosition()(1) );
            //     fprintf(docCM, "%4.8f\n", grainList[grain_index].getRadius());
            //     fprintf(docCM, "%4.8f\n", grainList[grain_index].getCmLset()(0));
            //     fprintf(docCM, "%4.8f\n", grainList[grain_index].getCmLset()(1));
                
            //     //Print Each sample level set
            //     for (size_t jj = 0; jj < yydimLS; jj++){
            //         for (size_t ii = 0; ii < xxdimLS; ii++){ 
            //             fprintf(docThick, "%4.8f\n", grainList[grain_index].getMthick().getLevelset()[jj*xxdimLS + ii] ); //Thickness LS
            //             fprintf(docGeo, "%4.8f\n", grainList[grain_index].getLset().getLevelset()[jj*xxdimLS + ii] ); //Geo LS
                        
    	       //     	lslocx = lowLx + 0.5 + ii; //IFF LS has unit of 1 for grid cells
	           // 	    lslocy = lowLy + 0.5 + jj; //IFF LS has unit of 1 for grid cells
	           // 	    pointOcxy << lslocx , lslocy;
            //             fprintf(docTemp, "%4.8f\n", round_Ocean(pointOcxy, oceanTemp, x_cells, y_cells, offset) ); //Temp LS
            //         }
            //     }  
                
            //     //Print points for meshing if needed
            //     for (size_t pti = 0; pti < grainList[grain_index].getPointList().size(); pti++){
            //         fprintf(docPoints, "%4.8f %4.8f\n", grainList[grain_index].getPointList()[pti](0), grainList[grain_index].getPointList()[pti](1) );
            //     }
                
            //     fclose(docThick);
            //     fclose(docGeo);
            //     fclose(docTemp);
            //     fclose(docPos);
            //     fclose(docCM);
            //     fclose(docPoints);
                
            //     //Print Fine Floe List for Fine Thickness Size Distribution
            //     string fnameThickF;
            //     fnameThickF = outDir4 + "TSD_Fine_step_" + std::to_string(step/nTout) + ".dat";    
            //     FILE * docTSDF = fopen(fnameThickF.c_str(), "w");                
                
            //     for (size_t fidx = 0; fidx < Fine_Grains.size(); fidx++) {
            //         fprintf(docTSDF,  "%4.8f %4.8f %4.8f %4.8f\n", Fine_Grains[fidx](0), Fine_Grains[fidx](1), Fine_Grains[fidx](2), Fine_Grains[fidx](3)); 
            //     }
            //     fclose(docTSDF);
            // //TURN ON FOR MORE INFO
                
//                //Print Damage File
//                char DamageFile[300];
//                for (size_t fidx = 0; fidx < grainList.size(); fidx++)
//                {
//                    std::string nngg = std::__cxx11::to_string(fidx);
//                    std::string ttt = std::__cxx11::to_string(step/nTout);
//                    sprintf(DamageFile,(folderLoc4+"Damage_g_"+nngg+"_step_"+ttt+".dat").c_str()); //1 file per grain per step
//                    FILE * damg_out   = fopen(DamageFile,"w");
//                    vector<Vector3d> damg_temp = grainList[fidx].getgrainDamage();
//                    for (size_t jj = 0; jj < damg_temp.size(); jj++)
//                    {
//                        fprintf( damg_out, "%4.8f %4.8f %4.8f\n", damg_temp[jj](0) ,damg_temp[jj](1) , damg_temp[jj](2) );
//                    }
//                    fclose(damg_out);
//                }    

                // //Print Evolution of Points for Output that shows changes over time since Points are NOT constant due to melting or fracture
                for (size_t i = 0; i < grainList.size(); i++) {

                    for (size_t j = 0; j<grainList[i].getnpoints(); j++) {

                            fprintf(pointsout,  "%4.8f %4.8f \n", grainList[i].getPointList()[j](0), grainList[i].getPointList()[j](1) );
                    }

                }
                
               //Print Intergranular forces and Moments
               for (size_t i = 0; i < grainList.size(); i++) {
                   fprintf(forces,  "%4.8f %4.8f \n", grainList2._grainForces[i](0),
                            grainList2._grainForces[i](1));

                    //cout << "output position over" << endl;
                   fprintf(moments, "%4.8f \n", grainList2._grainMoments[i]);
      
                }

                
                Max_Diam = 0.0;
                Min_Diam = 40000000000000000.00000;
                double Sum_Diam = 0.0;
                double addAreaD;
                for (size_t i = 0; i < grainList.size(); i++)
                {
                    //vector<Vector2d> forAreaD = grainList[i].getRefPointList();
                    vector<Vector2d> forAreaD = grainList[i].getPointList();
                    //addArea = PointsAreaM(forArea3);


                    double AreaD = 0.0;
                    size_t nD = forAreaD.size();

                    for (size_t pp = 0; pp < nD-1; pp++)
                    {
                        AreaD = AreaD + ( forAreaD[pp](0) * forAreaD[pp+1](1) -  forAreaD[pp](1) * forAreaD[pp+1](0) );
                    }
                    AreaD = AreaD + (forAreaD[nD-1](0) * forAreaD[0](1) -  forAreaD[nD-1](1) * forAreaD[0](0) );

                    addAreaD = 0.5*abs(AreaD);
                    Sum_Diam += 2 * sqrt( addAreaD / 3.141592653589793);
                    
                    if ( 2 * sqrt( addAreaD / 3.141592653589793)  > Max_Diam )
                    {
                        Max_Diam = 2 * sqrt( addAreaD / 3.141592653589793);
                    }
                    if ( 2 * sqrt( addAreaD / 3.141592653589793)  < Min_Diam )
                    {
                        Min_Diam = 2 * sqrt( addAreaD / 3.141592653589793);
                    }
                }

                if (grainList.size() > 0)
                {
                    Mean_Diam = Sum_Diam / double(grainList.size());
                }
                else
                {
                    Mean_Diam = 0.0;
                }
                
                //Prevent Nans
                if (isnan(Max_Diam) == 1){
                    Max_Diam = prevMax_Diam;
                }
                if (isnan(Mean_Diam) == 1){
                    Mean_Diam = prevMean_Diam;
                }
                if (isnan(Min_Diam) == 1){
                    Min_Diam = prevMin_Diam;
                }
                
                fprintf(Diametersout, "%4.8f %4.8f %4.8f\n", Max_Diam, Mean_Diam, Min_Diam);
                prevMax_Diam  =  Max_Diam;
                prevMean_Diam =  Mean_Diam;
                prevMin_Diam  =  Min_Diam;
                    
                duration = omp_get_wtime() - start;
                printf( "Timestep %d of %d (%4.2f%% complete, %4.2f minutes)\n",
                        int(step), int(nT), 100*double(step+1)/double(nT), duration/60.);

                printf("Tot. Kin. E = %.16f\n", Ke);
                printf("Fine_Grain_Size = %d\n", Fine_Grains.size() );
                printf("Fine_Grain_Size 0 thick = %.8f\n", Fine_Grains[0](1) );
                //printf("Stress = %.16f\n", world.outputStress())
                //cout << "StressMacro = " << grainList2._stressVoigt << endl;  //?????????
                //cout << grainList[0].getPosition().transpose() <<endl; //Get position of a single grain for reference  -->
                cout << "Forces = " << grainList2._grainForces[0](1) << endl;
                cout << grainList[0].getPosition().transpose() << endl; //Get position of a single grain for reference  -->
                cout << "Number of grains: " << grainList.size() << endl;            // cout << "GM = " << grainList2.getGrainMoments() << endl;
                

                double nnnp = grainList.size();
                
                //PRINT NUMBER OF GRAINS TO CONTROL LATER ON PLOTTING OF POINTS
                //numbergout << grainList.size() << std::endl;
                fprintf(numbergout, "%4.8f \n", nnnp);

                //Print Number of Points to Define each grain at each timestep (in case breakage creates more grains)
                size_t numberPG;
                for (size_t i = 0; i < grainList.size(); i++) {                
                    
                    numberPG = grainList[i].getPointList().size();
                    fprintf(nperg, "%d\n", int(numberPG));
                }

                //Print All Relevant Test Parameters File
                if (step==0)
                {
                    fprintf(testParams, "%4.8f\n", double(nT));            //Total Time Steps
                    fprintf(testParams, "%4.8f\n", double(nTemp));         //nTemp or dstep  //Variation of ocean temperature stepping and melting
                    fprintf(testParams, "%4.8f\n", double(nBreak));        //Breakage Frequency
                    fprintf(testParams, "%4.8f\n", double(qvert));         //Vertical melting parameter for ice
                    fprintf(testParams, "%4.8f\n", double(Khor));          //Horizontal ocean temperature diffusivity
                    fprintf(testParams, "%4.8f\n", double(Qatm));          //Q atmospheric/solar forcing component
                    fprintf(testParams, "%4.8f\n", double(Aatm));          //A atmospheric/solar forcing component
                    fprintf(testParams, "%4.8f\n", double(Batm));          //B atmospheric/solar forcing component
                    fprintf(testParams, "%4.8f\n", double(alpha_ice));     //Ice albedo related term
                    fprintf(testParams, "%4.8f\n", double(alpha_ocean));   //Ocean albedo related term
                    fprintf(testParams, "%4.8f\n", double(limitMaxDiam));        //Bkg Melt Diameter Separator
                    fprintf(testParams, "%4.8f\n", double(limitMinDiam));        //Min Diameter Resolution
                }               


                //Output Ocean Temp Grid (Read Cols and then Rows, print one long list per time step)
                //GET GRID SIZES

                //RECENT CHANGE
                if (step==0)
                {
                   fprintf(oceanTGridDim, "%4.8f \n", double(x_cells));
                   fprintf(oceanTGridDim, "%4.8f \n", double(y_cells));
                }  
                //Get grid coordinates (fixed grid)
                if (step==0)
                {
                    for (size_t i = 0; i < y_cells; i++) {
                        for (size_t j = 0; j < x_cells; j++) {
                            fprintf(oceanTGridCoord, "%4.8f %4.8f ", fluid_coord[j+i*x_cells](0), fluid_coord[j+i*x_cells](1) ); 
                        }
                    }
                    fprintf(oceanTGridCoord, "\n"); 
                }    
                //Get ocean grid temperatures over time
                for (size_t i = 0; i < y_cells; i++) {
                    for (size_t j = 0; j < x_cells; j++) {
                        fprintf(oceanTGrid, "%4.8f ", oceanTemp[j+i*x_cells]); 
                    }
                }
                fprintf(oceanTGrid, "\n"); 


               // cout << "stress = " << grainList2.getStress() << endl;

                //cout << "stress = " << world.getStress() << endl;   //Sort of works

            //    cout << "stress = " << world.globalGrainState().stressVoigt.transpose()/volume << endl;   //??
            } //End rank 0 condition    
        } // end output step

        // Take a time step
        //cout << "Run Timestep" << endl;
        
        //world.grainTimestep(); //Update with force as usual LS-DEM
        
        double temp_speed = 0.0 * floeVel[vel_index](1)/86400.0;   //Remove 0.00 * TEMPo
        world.changeFlowSpeed(temp_speed);
        world.changeFlowAngle(floeVel[vel_index](2));
        world.grainTimestepVelocity(); //Update with reanalysis velocities, in km/s
        
        //world.fineTimestep(Fine_Grains, Fine_Velocities);
        //cout << "Finish Timestep:" << step << endl; //For debugging
        stepup++;
        world.changestepup(stepup);
    }

    fclose(positions);
    fclose(velocities);

    fclose(forces);
    fclose(moments);

    fclose(temper);
    fclose(thick);
    fclose(sampleLS);

    fclose(pointsout);
    fclose(pointsoutLS);
    
    fclose(numbergout);
    fclose(Diametersout);
    fclose(nperg);
    fclose(oceanTGrid);
    fclose(oceanTGridDim);
    fclose(oceanTGridCoord);
    fclose(testParams);
    fclose(finepositions);
    fclose(numberfines);
    
    //Aug 18, 2022 Change
    //BEGIN
    fclose(massLoss);
    //END
    fclose(lossArea);
    
    fclose(wallPosT);
    fclose(wallPosB);
    fclose(wallPosL);
    fclose(wallPosR);
    
    
    // Print simulation time
    if (rank ==0) {
        duration = omp_get_wtime() - start;
        printf("Time taken: %dh, %dm, %2.2fs\n",
                int(floor(duration/3600.)),
                -int(floor(duration/3600.))
                +int(floor(fmod(duration/60., 60.))),
                -floor(fmod(duration/60., 60.))*60.
                +fmod(duration, 3600.) );
    }

    // Close MPI
    MPI_Finalize();
    
    //Execute command to do FEM inside code (Just TEST)
//    std::string filename = "/Users/rigobertomoncadalopez/Dropbox/Caltech_PhD/Research/Code/LSDEM2D-SeaIceDense/FEMtest.edp -nw";
//    std::string command = "FreeFem++ ";
//    command += filename;
//    system(command.c_str());

    cout << "Program terminated" << endl;

    return 0;
}
