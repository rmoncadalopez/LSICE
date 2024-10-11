/*
 * World2d.h
 *
 *  Created on: October 25, 2016
 * Author: Reid Kawamoto - Konstantinos Karapiperis - Liuchi Li
 */

#ifndef WORLD2D_H_
#define WORLD2D_H_

#include "definitions.h"
#include "WorldStates.h"
#include "Grain2d.h"
#include "Wall2d.h"
//#include "Wall2dPeriodicx.h"
#include "Fracture.h"
#include "readInputFile.h"

class World2d {

public:
  World2d(const Vector2d & offset) :_offset(offset) {_gDamping = 0.; _dt = 0.; _ngrains = 0; _maxId = 0; _only_dem = false;}

  // Cone + wallGrains + grains + periodic bcs in cell of particular offset         //Walls instead of offset??
  World2d(const vector<Grain2d> & grains, vector<Wall2d> & walls,
      const Vector2d & offset, const double & dt, const double & gDamping, size_t & stepup, double & slopedir, double & flowangle, double & flowforce, const FracProps2D FracProps, const vector<Grain2d> & grainsWall,
            const vector<Vector2d> & fluid_coord, vector<Vector2d> & Uwg, const size_t & x_cells, const size_t & y_cells, vector<double> & oceanTemp,
            const size_t & START_TEMP, const double & THERM_COEFF, const double & MELT_MASS, const size_t & BREAK_PROB, size_t & cell_sizex, size_t & cell_sizey, size_t & fluid_mode):
    _grains(grains), _walls(walls),
    _offset(offset), _dt(dt), _gDamping(gDamping), _stepup(stepup), _slopedir(slopedir), _flowangle(flowangle), _flowforce(flowforce), _fracProps(FracProps), _grainsWall(grainsWall), 
        _fluid_coord(fluid_coord), _Uwg(Uwg), _x_cells(x_cells), _y_cells(y_cells), _oceanTemp(oceanTemp), 
        _START_TEMP(START_TEMP), _THERM_COEFF(THERM_COEFF), _MELT_MASS(MELT_MASS), _BREAK_PROB(BREAK_PROB), _cell_sizex(cell_sizex), _cell_sizey(cell_sizey), _fluid_mode(fluid_mode){  //ADD    _fracProps(FracProps)

    _ngrains = grains.size();

        _ngrainsWall = grainsWall.size();
            
        _nwalls = walls.size();

    _globalGrainState.resize(_ngrains,_nwalls);   _globalGrainState.reset();   //Don't add _ngrainsWall since they do not receive forces, only induce on grains TODO ????

  }

  void readjustGrains() {
    //Preprocess Particular Level Set FOR TIME STEP 0 AND EVERYTHING ELSE. THIS INVOLVES MODIFYING ALL FROM SCRATCH (CRITICAL STEP)  IT ONLY HAPPENS ONCE
    if (_stepup == 0)
    {
    //   cout << "Number of grains: " << _ngrains << endl;
    //   for (size_t i = 0; i < _ngrains; i++){
    //     cout << "GRAIN: " << i << endl;
    //     size_t xxdim =  _grains[i].getgrainTemp2D().getXdim();  // To get Heat LSET
    //     size_t yydim =  _grains[i].getgrainTemp2D().getYdim();
    //     vector<double> geovec = _grains[i].getLset().getLevelset();  //Level Set Vector
    //     vector<double> Tvec = _grains[i].getgrainTemp2D().getLevelset();    //Heat Vector
    //     double OuterTemp = -10.;

    //     //Initialize Damage Matrix
    //     MatrixXd newDamage = MatrixXd::Constant(yydim, xxdim, -1);  // is out of level set, D {0,1}
    //     MatrixXd newThick = _grains[i].getgrainThickness();
    //     //MatrixXd newTerrain = _grains[0].getglobalTerrain();   //Might be useful.

    //     //Delimit Ice extents 
    //     //size_t x0 = 18; // Or more is Ice
    //     //size_t x1 = 132; // Or Less is Ice (varies)
    //     //size_t y0 = 16; // Or More is Ice (varies)
    //     //size_t y1 = 46; // Or Less is Ice 


    //     // for (size_t jj = 0; jj < yydim; jj++)
    //     // {
    //     //   for (size_t ii = 0; ii < xxdim; ii++)
    //     //   { 
    //     //     //Eliminate positives under or left
    //     //     if (  (geovec[jj*xxdim + ii] < 0  &&    jj<y0) ||   (geovec[jj*xxdim + ii] < 0  &&    ii<x0) )
    //     //     {  
    //     //        geovec[jj*xxdim + ii] = 1; 
    //     //     }
    //     //   }     
    //     // }

    //     // //SOFTEN CORNERS LEFT (VERY PARTICULAR VALUES)
    //     // for (size_t jj = 16; jj < 47; jj++)
    //     // {
    //     //   for (size_t ii = 18; ii < 24; ii++)
    //     //   { 
    //     //     //Eliminate positives under or left
    //     //     if (  geovec[jj*xxdim + ii] >= 0 )
    //     //     {  
    //     //        geovec[jj*xxdim + ii] = -0.5 ; 
    //     //     }
    //     //   }     
    //     // }
        
    //     //BINARIZE GEOMETRY AND TEMPERATURE TOO (TRY)
    //     for (size_t jj = 0; jj < yydim; jj++)
    //     {
    //       for (size_t ii = 0; ii < xxdim; ii++)
    //       { 
    //         if (geovec[jj*xxdim + ii] <= 0)
    //         {  
    //            geovec[jj*xxdim + ii] = -1;
    //            newDamage(jj,ii) = 0;   //Inside start with 0 damage. Outside is -1; 
    //         }
    //         else if (geovec[jj*xxdim + ii] > 0) 
    //         { 
    //            geovec[jj*xxdim + ii] = 1; 
    //            Tvec[jj*xxdim + ii] = OuterTemp;
    //            newThick(jj,ii) = 0;
    //         }
    //       }     
    //     }

    //     //Go to Matrix
    //     MatrixXd Cls_Tag0 = MatrixDefVec(yydim, xxdim, geovec);

    //     //Soften Level Set
    //     size_t maxttries = 50;
    //     MatrixXd Geo_Matrix = LS_Soft(Cls_Tag0, maxttries); 

    //     //Go Back to Vector
    //     geovec = MatrixtoVec(Geo_Matrix);

    //     //BINARIZE SOFTENED GEOMETRY AND TEMPERATURE TOO (TRY)
    //     for (size_t jj = 0; jj < yydim; jj++)
    //     {
    //       for (size_t ii = 0; ii < xxdim; ii++)
    //       { 
    //         if (geovec[jj*xxdim + ii] <= 0)
    //         {  
    //            newDamage(jj,ii) = 0;   //Inside start with 0 damage. Outside is -1; 
    //         }
    //         else if (geovec[jj*xxdim + ii] > 0) 
    //         { 
    //            Tvec[jj*xxdim + ii] = OuterTemp;
    //            newDamage(jj,ii) = -1;
    //            newThick(jj,ii) = 0;
    //         }
    //       }     
    //     }

    //     //Update Modified Level Set and Temperature (Outer and Inner)
    //     //_grains[i].change_Full_lset(geovec, xxdim, yydim);
    //     //_grains[i].change_Full_lset0(geovec, xxdim, yydim);  //Save for future flow control of terrain
    //     //_grains[i].change_Full_Temp_lset(Tvec, xxdim, yydim);
    //     _grains[i].changegrainThickness(newThick);

    //     //BUT THESE NEED TO RE-UPDATE ALL INFO

    //     double massnew;
    //     Vector2d gnewcm;
    //     double Inew;
    //     double gdense = _grains[i].getDensity();

    //     //int Point_Densifier = 2; //Run program twice under these changes if needed. Watch out for Outputs! The denser the slower.


    //     _fracProps.findGrainProps(_grains[i].getLset(), massnew, gnewcm, Inew, gdense);


    //     //POINTS!!!!!
    //     //vector <Vector2d> gnewi(_grains[0].getnpoints()*Point_Densifier);   //Duplicate points from now on for higher resolution
    //     vector <Vector2d> gnewi(_grains[i].getnpoints());

    //     //LSET New Function to Create New Point Sets from LSET and Sorted
    //     cout<<"Point Generation Adjust for grain: "<< i <<endl;
    //     //_fracProps.PointsGen(_grains[0].getLset() , gnewi, _grains[0].getnpoints()*Point_Densifier, gnewcm );
    //     _fracProps.PointsGen(_grains[i].getLset() , gnewi, _grains[i].getnpoints(), gnewcm );
    //     cout<<"End Point Generation Adjust for grain: "<< i <<endl;

    //     cout<<"Position Adjust for grain: "<< i <<endl;
    //     const double cos1 = cos(_grains[i].getTheta());
    //     const double sin1 = sin(_grains[i].getTheta());
    //     double cosine1 = cos1;
    //     double sine1 = sin1;
         
    //     //ADJUST POSITION 
    //     Vector2d positionnew;
    //     Vector2d g1cm = _grains[i].getCmLset();
    //     positionnew = _grains[i].getPosition();

    //     //LAND WALL POSITIONING
    //     if (i == 0)
    //     {
    //         positionnew(0) += 0.;  //0
    //         positionnew(1) += -1.;  //-3
    //     }
    //     else if (i == 1)
    //     {
    //         positionnew(0) += -1;  //-3
    //         positionnew(1) += -48.; //-49
    //     } 
    //     else if (i == 2)
    //     {
    //         positionnew(0) += -88.;  //-85
    //         positionnew(1) += -44.;  //-45        
    //     } 
    //     else if (i == 3)
    //     {
    //         positionnew(0) += -40.; //-40
    //         positionnew(1) += -87.; //-87      
    //     } 
    //     else
    //     {
    //       //NOTHING
    //     }
    //     //positionnew = _grains[i].getPosition()-Vector2d((g1cm(0)-gnewcm(0))*cosine1 - (g1cm(1)-gnewcm(1))*sine1,
    //     //                                          (g1cm(0)-gnewcm(0))*sine1 + (g1cm(1)-gnewcm(1))*cosine1);
    //     cout << "Grain: " << i << " X: " << positionnew(0) << " Y: " << positionnew(1) << endl;

    //     //Shift a bit
    //     // positionnew(0) += 60;
    //     // positionnew(1) += 20;

    //     cout<<"Point Adjust for grain: "<< i <<endl;
    //     for (size_t kk = 0; kk < gnewi.size(); ++kk)
    //     {
    //       gnewi[kk] = gnewi[kk] + positionnew; // - gnewcm;
    //     }

    //     //NEW RADIUS
    //     cout<<"Radius Adjust for grain: "<< i <<endl;
    //     double maxRnew = _grains[i].getRadius(); //JUST POINTS
    //     // //for (size_t i = 0 ; i < _grains[0].getnpoints()*Point_Densifier ; i++)
    //     // for (size_t idx = 0 ; i < _grains[i].getnpoints(); idx++)
    //     // { 
    //     //   gnewi[idx] += -gnewcm;//-(g3cm-g1cm)-_grains[g].getPosition();
    //     //   //Find bbox radius
    //     //   if (gnewi[idx].norm() > maxRnew)
    //     //   {
    //     //     maxRnew = gnewi[idx].norm();
    //     //   }
    //     // }

    //     //int new_no_pts = (int)gnewi.size();   //OFF POINT DENSIFIER AND NEW NO PTS.


    //     //UPDATE ALL (Mass, CM, Inertia, Points, Position, Radius)
    //     cout<<"Prop Adjust for Grain: "<< i <<endl;
    //     //_grains[i].changeMass(massnew);
    //     //_grains[i].changeCM(gnewcm);
    //     //_grains[i].changeInertia(Inew); 
    //     _grains[i].changePoints(gnewi);
    //     //_grains[i].changePos(gnewcm);   //Move to reference center CURRENT USE 222
    //     //_grains[i].changePos(positionnew);  //Shifts in an odd shape
    //     //_grains[i].changeRadius(maxRnew);  //To other location congurent with later functions
    //     //_grains[0].changeNPoints(new_no_pts);  //Turn off if points Densifier is equal to 1.
    //     _grains[i].changeDamage(newDamage);
    //     cout<<"Prop Adjust End for Grain: "<< i <<endl;
    //   }  
    }
  
  }
    
    // update nearest neighbor list
  void updateNNLists(){
    int numprocessors, rank; 
    MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    //#pragma omp parallel default(none) shared(numprocessors,rank,cout) firstprivate(count, nnSize, nnList) // reduction(+:threadWallState,threadGrainState) //num_threads(1)
    #pragma omp parallel for default(none) shared(rank,numprocessors, cout) schedule(static,20) //num_threads(16) 
      for (size_t i = rank; i < _ngrains; i++){
        vector <size_t> nnList;
        for(size_t j=i+1; j < _ngrains; j+=numprocessors){
            if ( _grains[i].vRadiusCheck2(_grains[j]) || _grains[i].bcircleGrainIntersectionXOffset(_grains[j], _offset(0))  ||  _grains[i].bcircleGrainIntersectionYOffset(_grains[j], _offset(1))  )  {  //Check None or 2 for Vradius check??
              nnList.push_back(j);
            }

        }
        #pragma omp critical //WARNING INSPECT THIS CHANGE!!! TODO
        {
          _grains[i].changeNNList(nnList);
        }
      }
  }
    
  // computes the snapshot of the world (all grain forces/moments, all wall forces, total stress)
  // by updating _globalGrainState and _globalWallState
  void computeWorldState() {

    // reset the global states
    _globalGrainState.reset();


    // define temp variables
    Vector2d force;         // force on a grain
    force.fill(0);          // initialization
    double momenti = 0.;      // moment on grain i
    double momentj = 0.;      // momnet on grain j
    Vector2d cmvec;         // branch vector
    Vector3d stress;        // stress in the assembly

        // define parallel processing stuff (MPI is initialized at the mainfile)
    int numprocessors, rank;
      MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  #pragma omp parallel default(none) shared(numprocessors, rank, cout) firstprivate(force, momenti, momentj, cmvec, stress) //num_threads(4)
  {
        
        //Start Fluid Modification (Initialization of Function Parameters)
        //--------------------------------------------------------------------------------------------------------------------------------
        
        //Initialize constants for fluid interaction with ice
        double Cha = 1.7e-3;                         //Air-ice skin drag coefficient (unitless)
        double Cva = 1.3e-3;                         //Air-ice body drag coefficient (unitless)
        double Chw = 1.3e-6; //in km units  //1.3e-3;  //Water-ice skin drag coefficient (m/s)   Check units!!!
        double Cvw = 1.3e-6; //in km units  //1.3e-3;  //Water-ice body drag coefficient (m/s)   Check units!!!
        double rhoice =  0.91e9; // 0.910e-6;                    //Density of ice
        double rhoair =  1.23e6; //1.23e-9;                     //Density of air If rhoice is 0.91e-6;
        double rhowater = 1.025e9;  //1.025e-6;                  //Density of water
        double hice = 2/12;                          //Assume Unitary floe thickness, is thickness proportional to area?????
        
        
        //Constant over time but modify to vary over space ??function of _getPosition??
        Vector2d Ua = Vector2d (0.,0.);                    //Constant Wind Velocity over space and time
        Vector2d Uw = Vector2d (0.,0.);                    //Constant Water Current Velocity over space and time
        
    
        double hair = 0.;                           //Proportion of floe height above water surface
        double hwater = 0.;                         //Proportion of floe height below water surface
        Vector2d fluidForceha = Vector2d (0.,0.);
        Vector2d fluidForcehw = Vector2d (0.,0.);
        Vector2d fluidForceh = Vector2d (0.,0.);
        Vector2d fluidForceva = Vector2d (0.,0.);
        Vector2d fluidForcevw = Vector2d (0.,0.);
        Vector2d fluidForcev = Vector2d (0.,0.);
        Vector2d ppv = Vector2d (0.,0.);
        Vector2d ppvn = Vector2d (0.,0.);
        Vector2d midp = Vector2d (0.,0.);
        
        Vector2d fluidForce = Vector2d (0.,0.);
        double fluidMoment = 0.0 ;
        
        //--------------------------------------------------------------------------------------------------------------------------------
        //End of Fluid Modification


        // //Temperature constants
        // double Tair = -5.; //Degrees Celsius
        // double Twater = 20.; //Degrees Celsius
        // double dh = 0.;   //Size of space step
        // double alphaice=1.; //Thermal Coefficient of Ice (3)
        // double meltTemp=0.; //Melting point of Ice (confirm for salinity)
        // double meltid=0.;
        // double KIc = 50; //Ice Fracture Toughness Varies from 50-100 m^1/2 kPa (Fresh Water ice is Tougher than Sea Ice)
        // double afactor = 0.2; //Proportion of notch on slab for fracture initiation;
        // double fdim = 820; //Dimensional Constant for Adjustment
        // //vector<double> Tvec;
        // //vector<double> Tvec0;
        // double meltV = (-0.000001/20)*(Twater-meltTemp)*alphaice; //Positive grows, negative melts //abs(0.00001) is fine 0.000001   0.000005also
        //double mfactor=0.;
        //For the sake of stability of the explicit scheme dhsteps=80000
        //double Utemper[80000] = {};
        //double Utemper0[80000] = {};

    
      
    // Initialize local states in each thread to later add all together
    GrainState2d threadGrainState;

    threadGrainState.resize(_ngrains,_nwalls);
    threadGrainState.reset();

    size_t ncontacts = 0;

    // Each thread grubs 40 chunks of iterations dynamically
    #pragma omp for schedule(dynamic, 40)
    for (size_t i = rank; i < _ngrains; i+=numprocessors){
      //FOR LS-DEM
      if (_only_dem == false)
      {
          // one-way grain-grain contact detection (each grain against grains of higher index)
          for (size_t j = i+1; j < _ngrains; j++){
            // Normal contacts
            if (_grains[i].bcircleGrainIntersection(_grains[j])){
               if(_grains[i].findInterparticleForceMoment(_grains[j], _dt, force, momenti, momentj, ncontacts, Vector2d(0,0), _grainsWall.size() )) {
                            
                cmvec = _grains[i].getPosition() - _grains[j].getPosition();
                threadGrainState._grainForces[i] += force;
                            if (force.norm() > 1000 ){
                //if (i == 40 || i == 41 || i == 42 || i == 43 || i == 44 || i == 45){
                                //cout << "Contact on grain i: " << i << " and j: " << j << endl;    
                                //cout << "Force on grain i" << i << " x: " << force(0) << " and y: " << force(1) << " by: " << j <<endl;
                            }
                            threadGrainState._grainForces[j] -= force;
                threadGrainState._grainMoments[i] += momenti;
                threadGrainState._grainMoments[j] += momentj;
                threadGrainState._stressVoigt(0) += force(0)*cmvec(0);
                threadGrainState._stressVoigt(1) += force(1)*cmvec(1);
                threadGrainState._stressVoigt(2) += 0.5*(force(1)*cmvec(0) + force(0)*cmvec(1));
    
    
                            Vector3d contactStress = Vector3d(force(0)*cmvec(0),force(1)*cmvec(1),0.5*(force(1)*cmvec(0) + force(0)*cmvec(1)));
                            threadGrainState._grainStress[i] += contactStress;
                            threadGrainState._grainStress[j] += contactStress;
                            threadGrainState._grainCoord[i].push_back(_grains[j].getId());//lists contact ids
                            threadGrainState._grainCoord[j].push_back(_grains[i].getId());
                }
            }
                                   
                    
    //        // Contacts due to periodic bcs (Offset)
            else if (_grains[i].bcircleGrainIntersectionXOffset(_grains[j], _offset(0))){
               if(_grains[i].findInterparticleForceMoment(_grains[j], _dt, force, momenti, momentj, ncontacts, Vector2d(_offset(0),0), _grainsWall.size() )) {
                cmvec = _grains[i].getPosition() - _grains[j].getPosition();
    
                 if (_grains[j].getPosition()(0) > _grains[i].getPosition()(0))
                  cmvec = _grains[i].getPosition() + Vector2d(_offset(0),0) - _grains[j].getPosition();
                 else
                  cmvec = _grains[i].getPosition() - Vector2d(_offset(0),0) - _grains[j].getPosition();
     
                threadGrainState._grainForces[i] += force;
                threadGrainState._grainForces[j] -= force;
                threadGrainState._grainMoments[i] += momenti;
                threadGrainState._grainMoments[j] += momentj;
                threadGrainState._stressVoigt(0) += force(0)*cmvec(0);
                threadGrainState._stressVoigt(1) += force(1)*cmvec(1);
                threadGrainState._stressVoigt(2) += 0.5*(force(1)*cmvec(0) + force(0)*cmvec(1));
    
                            Vector3d contactStress = Vector3d(force(0)*cmvec(0),force(1)*cmvec(1),0.5*(force(1)*cmvec(0) + force(0)*cmvec(1)));
                            threadGrainState._grainStress[i] += contactStress;
                            threadGrainState._grainStress[j] += contactStress;
                            threadGrainState._grainCoord[i].push_back(_grains[j].getId());//lists contact ids
                            threadGrainState._grainCoord[j].push_back(_grains[i].getId());
                }
            }
            else if (_grains[i].bcircleGrainIntersectionYOffset(_grains[j], _offset(1))){
               if(_grains[i].findInterparticleForceMoment(_grains[j], _dt, force, momenti, momentj, ncontacts, Vector2d(0,_offset(1)), _grainsWall.size() )) {
                cmvec = _grains[i].getPosition() - _grains[j].getPosition();
    
                 if (_grains[j].getPosition()(1) > _grains[i].getPosition()(1))
                  cmvec = _grains[i].getPosition() + Vector2d(0,_offset(1)) - _grains[j].getPosition();
                 else
                  cmvec = _grains[i].getPosition() - Vector2d(0,_offset(1)) - _grains[j].getPosition();
    
                threadGrainState._grainForces[i] += force;
                threadGrainState._grainForces[j] -= force;
                threadGrainState._grainMoments[i] += momenti;
                threadGrainState._grainMoments[j] += momentj;
                threadGrainState._stressVoigt(0) += force(0)*cmvec(0);
                threadGrainState._stressVoigt(1) += force(1)*cmvec(1);
                threadGrainState._stressVoigt(2) += 0.5*(force(1)*cmvec(0) + force(0)*cmvec(1));
    
                            Vector3d contactStress = Vector3d(force(0)*cmvec(0),force(1)*cmvec(1),0.5*(force(1)*cmvec(0) + force(0)*cmvec(1)));
                            threadGrainState._grainStress[i] += contactStress;
                            threadGrainState._grainStress[j] += contactStress;
                            threadGrainState._grainCoord[i].push_back(_grains[j].getId());//lists contact ids
                            threadGrainState._grainCoord[j].push_back(_grains[i].getId());
                }
            }
                    
    
    
          } // end loop over grain j full LS-DEM
      }
      //For regular dem
      else
      {
          vector <size_t> nnList = _grains[i].getNNList();
        size_t nnSize = nnList.size();
        //cout << "nnList Size: " << nnSize << endl;
        // one-way grain-grain contact detection (each grain against grains of higher index) using NN
        for (size_t count = 0; count < nnSize; count++){
          size_t j = nnList[count];
            // Normal contacts
            if (_grains[i].vRadiusCheck2(_grains[j])){  //vRadius replaces bCircleGrain
               if(_grains[i].findInterparticleForceMomentDEM(_grains[j], _dt, force, momenti, momentj, ncontacts, Vector2d(0,0), _grainsWall.size(), _ntension, _nshear)) {
                            
                cmvec = _grains[i].getPosition() - _grains[j].getPosition();
                threadGrainState._grainForces[i] += force;
                            if (force.norm() > 1000 ){
                //if (i == 40 || i == 41 || i == 42 || i == 43 || i == 44 || i == 45){
                                //cout << "Contact on grain i: " << i << " and j: " << j << endl;    
                                //cout << "Force on grain i" << i << " x: " << force(0) << " and y: " << force(1) << " by: " << j <<endl;
                            }
                            threadGrainState._grainForces[j] -= force;
                threadGrainState._grainMoments[i] += momenti;
                threadGrainState._grainMoments[j] += momentj;
                threadGrainState._stressVoigt(0) += force(0)*cmvec(0);
                threadGrainState._stressVoigt(1) += force(1)*cmvec(1);
                threadGrainState._stressVoigt(2) += 0.5*(force(1)*cmvec(0) + force(0)*cmvec(1));
    
    
                            Vector3d contactStress = Vector3d(force(0)*cmvec(0),force(1)*cmvec(1),0.5*(force(1)*cmvec(0) + force(0)*cmvec(1)));
                            threadGrainState._grainStress[i] += contactStress;
                            threadGrainState._grainStress[j] += contactStress;
                            threadGrainState._grainCoord[i].push_back(_grains[j].getId());//lists contact ids
                            threadGrainState._grainCoord[j].push_back(_grains[i].getId());
                }
            }
    //        // Contacts due to periodic bcs (Offset)
            else if (_grains[i].bcircleGrainIntersectionXOffset(_grains[j], _offset(0))){ //Use radius so they work well, but check bonding PBC
               if(_grains[i].findInterparticleForceMomentDEM(_grains[j], _dt, force, momenti, momentj, ncontacts, Vector2d(_offset(0),0), _grainsWall.size(), _ntension, _nshear )) {
                cmvec = _grains[i].getPosition() - _grains[j].getPosition();
    
                 if (_grains[j].getPosition()(0) > _grains[i].getPosition()(0))
                  cmvec = _grains[i].getPosition() + Vector2d(_offset(0),0) - _grains[j].getPosition();
                 else
                  cmvec = _grains[i].getPosition() - Vector2d(_offset(0),0) - _grains[j].getPosition();
     
                threadGrainState._grainForces[i] += force;
                threadGrainState._grainForces[j] -= force;
                threadGrainState._grainMoments[i] += momenti;
                threadGrainState._grainMoments[j] += momentj;
                threadGrainState._stressVoigt(0) += force(0)*cmvec(0);
                threadGrainState._stressVoigt(1) += force(1)*cmvec(1);
                threadGrainState._stressVoigt(2) += 0.5*(force(1)*cmvec(0) + force(0)*cmvec(1));
    
                            Vector3d contactStress = Vector3d(force(0)*cmvec(0),force(1)*cmvec(1),0.5*(force(1)*cmvec(0) + force(0)*cmvec(1)));
                            threadGrainState._grainStress[i] += contactStress;
                            threadGrainState._grainStress[j] += contactStress;
                            threadGrainState._grainCoord[i].push_back(_grains[j].getId());//lists contact ids
                            threadGrainState._grainCoord[j].push_back(_grains[i].getId());
                }
            }
            else if (_grains[i].bcircleGrainIntersectionYOffset(_grains[j], _offset(1))){
               if(_grains[i].findInterparticleForceMomentDEM(_grains[j], _dt, force, momenti, momentj, ncontacts, Vector2d(0,_offset(1)), _grainsWall.size(), _ntension, _nshear )) {
                cmvec = _grains[i].getPosition() - _grains[j].getPosition();
    
                 if (_grains[j].getPosition()(1) > _grains[i].getPosition()(1))
                  cmvec = _grains[i].getPosition() + Vector2d(0,_offset(1)) - _grains[j].getPosition();
                 else
                  cmvec = _grains[i].getPosition() - Vector2d(0,_offset(1)) - _grains[j].getPosition();
    
                threadGrainState._grainForces[i] += force;
                threadGrainState._grainForces[j] -= force;
                threadGrainState._grainMoments[i] += momenti;
                threadGrainState._grainMoments[j] += momentj;
                threadGrainState._stressVoigt(0) += force(0)*cmvec(0);
                threadGrainState._stressVoigt(1) += force(1)*cmvec(1);
                threadGrainState._stressVoigt(2) += 0.5*(force(1)*cmvec(0) + force(0)*cmvec(1));
    
                            Vector3d contactStress = Vector3d(force(0)*cmvec(0),force(1)*cmvec(1),0.5*(force(1)*cmvec(0) + force(0)*cmvec(1)));
                            threadGrainState._grainStress[i] += contactStress;
                            threadGrainState._grainStress[j] += contactStress;
                            threadGrainState._grainCoord[i].push_back(_grains[j].getId());//lists contact ids
                            threadGrainState._grainCoord[j].push_back(_grains[i].getId());
                }
            }
    
          } // end loop over grain j regular DEM
      }     
      
            
            //Start Wall Modification
            //--------------------------------------------------------------------------------------------------------------------------------
      
            // walls
            for (size_t j = 0; j < _walls.size(); j++) {
              if (_walls[j].bcircleWallIntersection(_grains[i])) {
                  if (_walls[j].findWallForceMoment(_grains[i], force, momenti, ncontacts, stress, _dt)) {    //Watch out for stress  or Vector2d(0,0) ?????
                      threadGrainState._wallContacts[j] += ncontacts;
                      // compute things related to forces
                      threadGrainState._grainForces[i] += force;
                      threadGrainState._wallForces[j] -= force;
                      //cout << "Wall Force on grain i" << i << "x: " << force(0) << " and y: " << force(1) << " by: " << j <<endl;
                      // compute things related to moments
                      threadGrainState._grainMoments[i] += momenti;
                      threadGrainState._stressVoigt += stress;

                      threadGrainState._grainStress[i] += stress;
                      threadGrainState._grainCoord[i].push_back(_walls[j].getId());
                  }
              }
            }
      
            //End Wall Modification
            //--------------------------------------------------------------------------------------------------------------------------------

            //Start Land Wall Modification
            if (_grainsWall.size()>0)
            {    
                for (size_t j = 0; j < _ngrainsWall; j++) {
                    if (_grains[i].bcircleGrainIntersection(_grainsWall[j])){
                        //cout << "Contact of grain i: " << i << " and grain Wall j: " << j << endl;
                        if(_grains[i].findInterparticleForceMoment(_grainsWall[j], _dt, force, momenti, momentj, ncontacts, Vector2d(0,0), _grainsWall.size() )) {
                        //if(_grainsWall[j].findInterparticleForceMoment(_grains[i], _dt, force, momenti, momentj, ncontacts, Vector2d(0,0))) {
                                cmvec = _grains[i].getPosition() - _grainsWall[j].getPosition();
                                threadGrainState._grainForces[i] += force;
                                threadGrainState._grainMoments[i] += momenti;
                                //cout << "Force on grain i: " << i << " F: " << force << endl;
                                //cout << "Wall Grain Force on grain i" << i << " x: " << force(0) << " and y: " << force(1) << " by: " << j <<endl;
                                threadGrainState._stressVoigt(0) += force(0)*cmvec(0);
                                threadGrainState._stressVoigt(1) += force(1)*cmvec(1);
                                threadGrainState._stressVoigt(2) += 0.5*(force(1)*cmvec(0) + force(0)*cmvec(1));

                            Vector3d contactStress = Vector3d(force(0)*cmvec(0),force(1)*cmvec(1),0.5*(force(1)*cmvec(0) + force(0)*cmvec(1)));
                            threadGrainState._grainStress[i] += contactStress;
                            threadGrainState._grainCoord[i].push_back(_grainsWall[j].getId());//lists contact ids
                        }
                    }  
                }
            }    
               
            //End Land Wall Modification



       //Start Fluid Modification Function (add the additional fluid force to the grains)
       //--------------------------------------------------------------------------------------------------------------------------------
           
            //get force & moments of fluids
            //For Fluid Interaction Full
            //Cvw = 1.3e-27; // FORM Cvw = 1.3e-3; /Starts to have some effect
            //Chw = 1.3e-27; // SKIN Cvw = 1.3e-3; /Starts to have some effect
            
            //size_t nFluid = 10000; //Adjust fluid Time Step for Efficiency
            
            if (_stepup > 0 && i >= 0 )  //Fluid forces should be applied each time step even if they dont change
            //if (_stepup > 0 && i >= 0 && (_stepup % nFluid == 0))
            {
              //cout <<"Fluid Interaction for grain: " << i << " at step: " << _stepup <<  endl;
              //cout <<"Position X: " << _grains[i].getPosition()(0) << " Position Y: " << _grains[i].getPosition()(1) << endl;
              //_grains[i].fluidInteraction(Cha, Cva, Chw, Cvw, rhoice, rhoair, rhowater, hice, hair,  hwater, Ua, Uw, fluidForceha, fluidForcehw, fluidForceh,
              //                            fluidForceva, fluidForcevw, fluidForcev, ppv, ppvn, midp, force, fluidMoment, _stepup, _offset, _fluid_coord, _Uwg, _x_cells, _y_cells);

              if (_fluid_mode == 1){
                _grains[i].fluidInteractionSimple(Cha, Cva, Chw, Cvw, rhoice, rhoair, rhowater, hice, hair,  hwater, Ua, Uw, fluidForceha, fluidForcehw, fluidForceh, fluidForceva, fluidForcevw, fluidForcev, ppv, ppvn, midp, force, fluidMoment, _stepup, _slopedir, _flowangle, _flowforce, _offset);
              }
              else if (_fluid_mode == 2){
                if (i == 0)
                {
                    //cout << "Before Force Drag: " << endl;
                    //cout << "X: " << force(0) << " Y: " << force(1) << endl;
                    //cout << "Moment Drag: " << fluidMoment << endl;
                }
                _grains[i].fluidInteraction_Nova(Cha, Cva, Chw, Cvw, rhoice, rhoair, rhowater, hice, hair,  hwater, Ua, Uw, fluidForceha, fluidForcehw, fluidForceh, fluidForceva, fluidForcevw, fluidForcev, ppv, ppvn, midp, force, fluidMoment, _stepup, _slopedir, _flowangle, _flowforce, _offset, _cell_sizex, _cell_sizey, _Uwg, _x_cells, _y_cells, _fluid_coord, i); //Change fluid Jan 9, 2023
                if (i == 0)
                {
                    //cout << "After Force Drag: " << endl;
                    //cout << "X: " << force(0) << " Y: " << force(1) << endl;
                    //cout << "Moment Drag: " << fluidMoment << endl;
                }
                  
              }
              else if (_fluid_mode == 3){
                  if (_only_dem){
                      if (_ice_damp == false){
                          //_grains[i].fluidInteraction_NovaDEM(Cha, Cva, Chw, Cvw, rhoice, rhoair, rhowater, hice, hair,  hwater, Ua, Uw, fluidForceha, fluidForcehw, fluidForceh, fluidForceva, fluidForcevw, fluidForcev, ppv, ppvn, midp, force, fluidMoment, _stepup, _slopedir, _flowangle, _flowforce, _offset, _cell_sizex, _cell_sizey, _Uwg, _x_cells, _y_cells, _fluid_coord, i, _dt);
                          _grains[i].fluidInteraction_NovaDEM_Damp_big(Cha, Cva, Chw, Cvw, rhoice, rhoair, rhowater, hice, hair,  hwater, Ua, Uw, fluidForceha, fluidForcehw, fluidForceh, fluidForceva, fluidForcevw, fluidForcev, ppv, ppvn, midp, force, fluidMoment, _stepup, _slopedir, _flowangle, _flowforce, _offset, _cell_sizex, _cell_sizey, _Uwg, _x_cells, _y_cells, _fluid_coord, i, _dt, _dampMat, _curr_factor);
                      }
                      else{
                          //Big for floes whose Area > cell area, small for floes whose area < cell area
                          _grains[i].fluidInteraction_NovaDEM_Damp_big(Cha, Cva, Chw, Cvw, rhoice, rhoair, rhowater, hice, hair,  hwater, Ua, Uw, fluidForceha, fluidForcehw, fluidForceh, fluidForceva, fluidForcevw, fluidForcev, ppv, ppvn, midp, force, fluidMoment, _stepup, _slopedir, _flowangle, _flowforce, _offset, _cell_sizex, _cell_sizey, _Uwg, _x_cells, _y_cells, _fluid_coord, i, _dt, _dampMat, _curr_factor);
                          //_grains[i].fluidInteraction_NovaDEM_Damp_small(Cha, Cva, Chw, Cvw, rhoice, rhoair, rhowater, hice, hair,  hwater, Ua, Uw, fluidForceha, fluidForcehw, fluidForceh, fluidForceva, fluidForcevw, fluidForcev, ppv, ppvn, midp, force, fluidMoment, _stepup, _slopedir, _flowangle, _flowforce, _offset, _cell_sizex, _cell_sizey, _Uwg, _x_cells, _y_cells, _fluid_coord, i, _dt, _dampMat);
                      }
                  }
                  else{
                      _grains[i].fluidInteraction_Nova(Cha, Cva, Chw, Cvw, rhoice, rhoair, rhowater, hice, hair,  hwater, Ua, Uw, fluidForceha, fluidForcehw, fluidForceh, fluidForceva, fluidForcevw, fluidForcev, ppv, ppvn, midp, force, fluidMoment, _stepup, _slopedir, _flowangle, _flowforce, _offset, _cell_sizex, _cell_sizey, _Uwg, _x_cells, _y_cells, _fluid_coord, i);
                  }
              }
              else if (_fluid_mode == 4){
                  if ( _grains[i].getLoadGrain() ){ //Function applied only to loaded grains
                      if (_only_dem){
                          //Skip, instead do this Main Loop just for loadable grains.
                          //_grains[i].fluidInteractionLoad_NovaDEM(Cha, Cva, Chw, Cvw, rhoice, rhoair, rhowater, hice, hair,  hwater, Ua, Uw, fluidForceha, fluidForcehw, fluidForceh, fluidForceva, fluidForcevw, fluidForcev, ppv, ppvn, midp, force, fluidMoment, _stepup, _slopedir, _flowangle, _flowforce, _offset, _cell_sizex, _cell_sizey, _Uwg, _x_cells, _y_cells, _fluid_coord, i);
                      }
                      else{
                          _grains[i].fluidInteraction_Nova(Cha, Cva, Chw, Cvw, rhoice, rhoair, rhowater, hice, hair,  hwater, Ua, Uw, fluidForceha, fluidForcehw, fluidForceh, fluidForceva, fluidForcevw, fluidForcev, ppv, ppvn, midp, force, fluidMoment, _stepup, _slopedir, _flowangle, _flowforce, _offset, _cell_sizex, _cell_sizey, _Uwg, _x_cells, _y_cells, _fluid_coord, i);
                      }
                  }
              }  
              else
              {
                 cout << "WARNING: FLUID MODE is INCORRECT, no forces will be applied!!! EXIT PROCESS!!!" << endl; 
                 cout << "NOPE, NOPE" << endl;
                 exit(1);
              }
              

              threadGrainState._grainForces[i] += force;    //Contribute fluid interaction force
              threadGrainState._grainMoments[i] += fluidMoment; //Contribute fluid interaction moment
              //cout << "Force, Moment" << force(0) << " , " << force(1) << " , " <<fluidMoment << endl;
            }
            
//          momenti=0; //A constant field will not exert net moment

            //ANNULL ALL FORCES WHEN NEEDED ARBITRARILY
            if (i < 0)
            {  
              threadGrainState._grainForces[i] = force*0;
              threadGrainState._grainMoments[i] = momenti*0;
            }
        //-------------------------------------------------------------------------------------------------------------------------------
        //End Fluid Modification

        
           //Temperature Modification   
           //_grains[i].TemperatureModif(Tair, Twater, _dt, dh, alphaice, meltTemp, meltid, meltV, _stepup, KIc, afactor, fdim);    //Utemper, Utemper0
           
           //_grains[i].changeTemper(_grains[i].getMass()); 


          //Use to build new grain if needed otherwise they stay empty
          //if (i == 5){
          //    cout << "New mass/temper for grain: " << i << " is :" <<  _grains[i].getTemperature()  <<endl;
          //}

      } // end loop over grains
        // Assemble the global state
    #pragma omp critical
    {
      _globalGrainState += threadGrainState;

    }
  } // close opemp parallel section

  // MPI calls  sendbuf       recvbuff                                      count           type               op       comm
  MPI_Allreduce(MPI_IN_PLACE, _globalGrainState._grainForces[0].data(),     _ngrains*2,     MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, _globalGrainState._grainMoments.data(),       _ngrains,       MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, _globalGrainState._stressVoigt.data(),        3,              MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);
        
    //What about walls here???????
        

  } // end computeWorldState

  // Gets world state at the given snapshot for output purposes
  CData computeCstate() const {

    CData cDataRank;
    int numprocessors, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    #pragma omp parallel default(none) shared(numprocessors,rank,cout,cDataRank) // num_threads(1)
    {
      CData cDataThread;
      #pragma omp for schedule(dynamic, 5)
      for (size_t i = rank; i < _ngrains; i+=numprocessors) {
        if (_only_dem == false){
            //Grain-grain contact subloop
                    for (size_t j = i+1; j < _ngrains; j++) {
              // Normal contacts
              if (_grains[i].bcircleGrainIntersection(_grains[j])){
                CData cDataContact = _grains[i].findContactData(_grains[j], _dt, 0);
                if (cDataContact._clocs.size() > 0) {
                  cDataThread += cDataContact;
                }
              }
              // Contacts due to periodic bcs
              else if (_grains[i].bcircleGrainIntersectionXOffset(_grains[j], _offset(0))){
                CData cDataContact = _grains[i].findContactData(_grains[j], _dt, _offset(0));
                if (cDataContact._clocs.size() > 0) {
                  cDataThread += cDataContact;
                }
              }
            } // close grain subloop FULL LS-DEM
        }
        else{
            //Grain-grain contact subloop regular DEM
                    for (size_t j = i+1; j < _ngrains; j++) {
              // Normal contacts
              if (_grains[i].bcircleGrainIntersection(_grains[j])){
                CData cDataContact = _grains[i].findContactDataDEM(_grains[j], _dt, 0);
                if (cDataContact._clocs.size() > 0) {
                  cDataThread += cDataContact;
                }
              }
              // Contacts due to periodic bcs
              else if (_grains[i].bcircleGrainIntersectionXOffset(_grains[j], _offset(0))){
                CData cDataContact = _grains[i].findContactDataDEM(_grains[j], _dt, _offset(0));
                if (cDataContact._clocs.size() > 0) {
                  cDataThread += cDataContact;
                }
              }
            } // close grain subloop regular DEM
            
          }
                
                // //Grainwall contact subloop
                // for (size_t jj = 0; jj < _ngrainsWall; jj++) {
                //     // Normal contacts
                //     if (_grains[i].bcircleGrainIntersection(_grainsWall[jj])){
                //         CData cDataContact = _grains[i].findContactData(_grainsWall[jj], _dt, 0);
                //         if (cDataContact._clocs.size() > 0) {
                //             cDataThread += cDataContact;
                //         }
                //     }
                // } // close grainWall subloop
                
        // // TODO: Add here particle-wall contact data CORRECT????
    //             for (size_t j = 0; j < _walls.size(); j++) {
    //                     // Normal contacts
    //                     if (_walls[j].bcircleWallIntersection(_grains[i])){
    //                         CData cDataContact = _walls[j].findContactData(_grains[i], _dt, 0);
    //                         if (cDataContact._clocs.size() > 0) {
    //                             cDataThread += cDataContact;
    //                         }
    //                     }
    //             } // close wall subloop
           
      } // end loop over grains
      #pragma omp critical
      {
        cDataRank += cDataThread;
      }
    } // closes openmp parallel section
    return cDataRank;
  } // end computeCstate method


    //Contact State for Grain Walls
    CData computeCstateGWalls() const {

        CData cDataRank;
        int numprocessors, rank;
        MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        #pragma omp parallel default(none) shared(numprocessors,rank,cout,cDataRank) // num_threads(1)
        {
            CData cDataThread;
            #pragma omp for schedule(dynamic, 5)
            for (size_t i = rank; i < _ngrains; i+=numprocessors) {
                //Grain-grain contact subloop
                for (size_t j = 0; j < _ngrainsWall; j++) {
                    // Normal contacts
                    if (_grains[i].bcircleGrainIntersection(_grainsWall[j])){
                        CData cDataContact = _grains[i].findContactData(_grainsWall[j], _dt, 0);
                        if (cDataContact._clocs.size() > 0) {
                            cDataThread += cDataContact;
                        }
                    }
                } // close grain subloop
                
                // // TODO: Add here particle-wall contact data CORRECT????
    //             for (size_t j = 0; j < _walls.size(); j++) {
    //                     // Normal contacts
    //                     if (_walls[j].bcircleWallIntersection(_grains[i])){
    //                         CData cDataContact = _walls[j].findContactData(_grains[i], _dt, 0);
    //                         if (cDataContact._clocs.size() > 0) {
    //                             cDataThread += cDataContact;
    //                         }
    //                     }
    //             } // close wall subloop
           
            } // end loop over grains
            #pragma omp critical
            {
                cDataRank += cDataThread;
            }
        } // closes openmp parallel section
        return cDataRank;
    } // end computeCstate method

  // void TwoHighestForcesLocs(Vector2d & pF1,Vector2d & pF2,const size_t & g,const vector<PointInfo> & g1i){
  //   //takes grain index, g, and contact info, g1i, and finds locations of the 2 highest contact forces, pF1, pF2


  //   //ensure each contact is checked only once
  //   sort(_globalGrainState.grainCoord[g].begin(),_globalGrainState.grainCoord[g].end());
  //   unique(_globalGrainState.grainCoord[g].begin(),_globalGrainState.grainCoord[g].end());

  //   vector<CForce> cforces(_globalGrainState.grainCoord[g].size());
  //   //find points with 2 highest contact forces
  //   for (size_t c=0; c<_globalGrainState.grainCoord[g].size(); c++) { //for loop over contact grains

  //     size_t contact = _globalGrainState.grainCoord[g][c];
  //     double force=0;
  //     size_t contactpoint=0;

  //     //grain contacts
  //     if (contact < INT_MAX-5){// || contact == 2147483645){
  //       Grain2d other = _grains[FindGrainFromId(contact)];
  //       contactpoint = _grains[g].findInterparticleForceforFrac(other, force);
  //       //cout << "grain" << endl;
  //     }
  //     //wall contacts
  //     else if ( contact > INT_MAX-5){//wall contacts or jaw
  //       size_t w = FindWallFromId(contact);
  //       contactpoint = _walls[w].FindWallForceForFrac(_grains[g],force);
  //       //cout << "wall" << endl;
  //     }
  //     else {
  //       contactpoint = 0;
  //     }
  //     cforces[c].force = force;
  //     cforces[c].cpoint = g1i[contactpoint].point;
  //     cforces[c].cid = contact;

  //   }//end find points
  //   sort(cforces.begin(),cforces.end());
  //   pF1 = cforces[0].cpoint;
  //   pF2 = cforces[1].cpoint;


  // }

  //world.TempState();
  //world.StressState());
  //world.Damage();
  //world.Fracture();


  //If Global Grid is equal or less than local grid size we can use this round function
  double round_Ocean(Vector2d & pointxy, vector<double> & oceanTemp, const size_t & x_cells, const size_t & y_cells, const Vector2d & offset)
  {
    size_t idx1, idy1;  // j are x cols, i are y rows
    int cell_sizex = int(offset(0))/int(x_cells-1);
    int cell_sizey = int(offset(1))/int(y_cells-1);
    //cout << "Fines cell size x: " << cell_sizex << " cell_sizey: " << cell_sizey << endl;

    if ( isnan(pointxy(0)) ||  isnan(pointxy(1)) )
    {
      cout <<"WARNING NaN VALUES PROVIDED FOR INTERPOLATION, OUTPUT = 0" << endl;
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
    //cout << "Fines idx: " << idx1 << " idy: " << idy1 << endl;

    //Use this exact grid point rounded from the input point at local scale
    //cout << "Ocean Temp Fine" << endl;
    //cout << oceanTemp[idx1 + idy1*x_cells] << endl;
    return oceanTemp[idx1 + idy1*x_cells];
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

    //Verify if it's an exact grid point
    if ( ( fmod(input_p(0), int(cell_sizex)) == 0) && ( fmod(input_p(1), int(cell_sizey))  == 0  )  )
    {
       return oceanTemp[idx1 + idy1*x_cells];
    }

    //cout << "inputpx: " << input_p(0) << " inputpy: " << input_p(1) << endl;  
    //Find point values value based on location indices 
    size_t NNx_cells = size_t(offset(0))/(cell_sizex) + 1;
    size_t NNy_cells = size_t(offset(1))/(cell_sizey) + 1;
    //cout << "x_cells: " << NNx_cells << endl;
    //cout << "y_cells: " << NNy_cells << endl;

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

    if (abs(x) > _offset(0) ||  abs(y) > _offset(1) )
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

  //void TempState(vector<double> & MeltBin, vector<double> & Diameters, double & outer_TempW, double & limitMaxDiamV, double & limitMinDiamV, const size_t & dstep, double & qvert, double & Khor,  double & alpha_ice, double & Qatm, double & Aatm, double & Batm, vector<Vector4d> & Fine_Grains, double & a_i, double & a_x, double & T_fave, double & T_xave)
  //Change Aug 22, 2022
  void TempState(vector<double> & MeltBin, vector<double> & Diameters, double & outer_TempW, double & limitMaxDiamV, double & limitMinDiamV, const size_t & dstep, double & qvert, double & Khor,  double & alpha_ice, double & Qatm, double & Aatm, double & Batm, vector<Vector4d> & Fine_Grains, double & a_i, double & a_x, double & T_fave, double & T_xave,
                 double & loss_mcv, double & loss_fines, vector<size_t> & NLMelt, vector<size_t> & NLBreak, vector<size_t> & NGMelt, vector<size_t> & NGBreak, double & loss_mcv_solar, double & loss_mcv_ocean, double & area_loss_basal, double & area_loss_lat, double & ave_delta_r, vector<double> & LatRate, double & MLAT, double & MBASAL)
  {
    
    //Heat units will be in m, s, kg and K, re-scale if need w.r.t to global km units!!!

    //Temperature constants
    double dh = 0.;   //Size of space step
    double Twater = outer_TempW+1.8;
    double Tair = -5+1.8; //Degrees Celsius

    //Need to reset Melt base line (CAUTION)
    T_fave += 1.8;
    T_xave += 1.8;

    //double Twater = 20.; //Degrees Celsius //20 melt  //-50 cool
    
    double alphaice=1.; //Thermal Coefficient of Ice (3)
    double meltTemp= -1.8+1.8; //0 //Melting point of Ice (confirm for salinity) 0 melt freeze - 40  //You need to always use a zero temp for LS convenience, simply shift the external temperature!!!
    double meltid=0.;
    double KIc = 50;  //ADJUST //Ice Fracture Toughness Varies from 50-100 m^1/2 kPa (Fresh Water ice is Tougher than Sea Ice)
    double afactor = 0.2; //Proportion of notch on slab for fracture initiation;
    double fdim = 820; //Dimensional Constant for Adjustment

    double therm_dir = -1;  //-1 is melt, +1 is freeze
    double therm_coeff = _THERM_COEFF; //See below for ref values, how fast ice melts, adjustment factor

    therm_coeff = 2.5e-3; //2.5e-5*0.5; //in W/km*K   //2.5 W/m*K

    double rho_icem = 910; //kg/m3
    double Lf = 330000; //J/kg

    //vector<double> Tvec;
    //vector<double> Tvec0;
    //double meltV = therm_dir*(therm_coeff)*(Twater-meltTemp)*alphaice; //Same as qvert //-0.001 for freeze //Positive freezes, negative melts //abs(0.00001) is fine 0.000001   0.000005also
    
    double meltV = (qvert/(rho_icem*Lf)); //For coarse
    double fine_adjust_fac = 1.00; //fadjust //17.5
    double meltVf = (fine_adjust_fac*qvert/(rho_icem*Lf)); //For fines
    //Older//double meltVSun = 1.0 * ( (Qatm*(1-alpha_ice) - (Aatm + Batm * (meltTemp-1.8))) *   1.0  )   / (rho_icem*Lf) ; //- because it should remove thickness //Test
    //Newer
    double meltVSun = ( - Qatm*(1-a_i) + (Aatm + Batm *(meltTemp-1.8)) ) / (rho_icem*Lf) ;

    Vector2d pointOcxy; //Shifting point for each cell starting from lower left corner of LS
    double Tempval; //Value to use in Temp Local Grid interpolated from Global Grid
    
    //Update fine grains (only thickness)
    cout << "Update Fine for nfines: " << Fine_Grains.size() << endl;
    for (size_t i = 0; i < Fine_Grains.size(); i++){ 
      double Orig_thick = Fine_Grains[i](1);

      //pointOcxy << min(abs(Fine_Grains[i](2)), abs(_offset(0))) , min(abs(Fine_Grains[i](3)), abs(_offset(1))) ; //RN don't move for convenience sake
      //pointOcxy << min(abs(1.0), abs(_offset(0))) , min(abs(1.0), abs(_offset(1))) ;
      //pointOcxy << min(abs(_grains[i].getPosition()(0)), abs(_offset(0))) , min(abs(_grains[i].getPosition()(1)), abs(_offset(1))) ;
      //Bilinear interpolation of this point on Ocean Grid and get temp value (how fast will this be?)
      //Tempval = bilinear_Interp_Ocean(pointOcxy, _fluid_coord, _oceanTemp, _x_cells, _y_cells, _offset) + 1.8;
      
      // //For Debugging Purposes, remove grain later on and abort 
      // if (isnan(pointOcxy(0)) || isnan(pointOcxy(1)))
      // {
      //   cout << "Ill-positioned Grain, remove" << endl;
      //   //_grains[i].changeMass(0.0);
      //   Fine_Grains[i](1) = 0.00; 
      //   continue;
      // }

      // //Average Global Temperature in Domain
      // double TaveGlobal = 0.0;
      // for (size_t idx = 0; idx < _y_cells; idx++) {
      //     for (size_t jdx = 0; jdx < _x_cells; jdx++) {
      //         TaveGlobal += _oceanTemp[jdx+idx*_x_cells];
      //     }
      // }
      // TaveGlobal /= (_x_cells*_y_cells);
      // TaveGlobal += 1.8;

      //cout << "Fine Point Interpol" << endl;
      //Tempval = round_Ocean(pointOcxy, _oceanTemp, _x_cells, _y_cells, _offset) + 1.8;
      //Tempval = TaveGlobal;
      //cout << "Temp. compare fine #: " << i <<  " meltTemp: " << meltTemp << " Tempval: " << Tempval <<  " MeltV: " << meltVf << " MeltSun: " << meltVSun << endl;
      //cout << "Thick. compare fine #: " << i <<  " Before: " << Orig_thick << " After: " << Orig_thick + dstep*(meltVf)*(meltTemp - Tempval) + dstep*(meltVSun) << endl;
      //Older
      //Fine_Grains[i](1) = max(Orig_thick + dstep*(meltVf)*(meltTemp - Tempval) + dstep*(meltVSun) , 0.000); //Return zero if less than zero   //TODO ALPHA: Update using Global Grid Temp using Position and Bilinear Interpolation, don't forget 1.8 adjust
      //Newer fine average

      //MAIN EXPRESSION
      //Fine_Grains[i](1) = max( Orig_thick + dstep*(meltVf)*(meltTemp - T_fave) + (dstep/(rho_icem*Lf))*( -Qatm * (1-a_x) + (Aatm + Batm * (T_xave-1.8)) )   ,  0.000);
      //New thick             //Old thick     //qvert comp.  * (Tmelt-T_fave)        //Denom                      //Q, A, BT part
      
      //PROPOSED EXPRESSION
      Fine_Grains[i](1) = max( Orig_thick + dstep*(meltVf)*( (meltTemp-1.8) - (T_fave-1.8) ) + (dstep/(rho_icem*Lf))*( -Qatm * (1-a_i) + (Aatm + Batm * (meltTemp-1.8)) )   ,  0.000);
      //Check result
      
      //Calculate mass loss
      //Change Aug 22, 2022
      //BEGIN
      double loss_fines_temp = 0.0;
      loss_fines_temp = Fine_Grains[i](0) * max((Orig_thick -  Fine_Grains[i](1)) , 0.0) * 0.001 * rho_icem * 1000*1000*1000 ; //Adjusted to kg from km3
      loss_fines += loss_fines_temp;
      //END
    }
    
    //For lateral rate estimation on coarse grains (sum for averate)
    double delta_r_temp = 0.0;
    ave_delta_r = 0.0;
    double delta_r_vec = 0.0;
    //Lat. Rate for bins
    vector<double> deltaRV(LatRate.size());
    vector<size_t> nFloes(LatRate.size());
    for (size_t i = 0; i < LatRate.size(); i++)
    {
        //LatRate[i] = 0.0;
        deltaRV[i] = 0.0;
        nFloes[i] = 0;
    }
    
    double h_0, h_f;
    double hlat;

    //Update coarse grains (geometry and thickness)
    cout << "Update Coarse for ngrains: " << _ngrains << endl;
    for (size_t i = 0; i < _ngrains; i++){ 
      //if ( _grains[i].getMass() < _MELT_MASS && _grains[i].getMass() > 0.009*_MELT_MASS) //0.009 //All grains smaller than a threshold mass melt but not too small to avoid lset problems
      //if (i % 2 == 0)
      //if (i == 40 || i == 41 || i == 42 || i == 43 || i == 44 || i == 45)
      //{
          //cout << "For grain: " << i << endl;
          //Initial Area and Diameter before melting
          //double IAreaMelt = PointsArea(_grains[i].getPointList());
          //double IDiamMelt = MeanCaliperD(IAreaMelt);

          //cout << "Radius" << _grains[i].getRadius() << endl;
          //cout << "Point0" << _grains[i].getPointList()[0] << endl;
          

          // //Point Control for a big grain
          // if (i == 37)
          // {
          //    cout << "Points for Audit" << endl;
          //    for (size_t iik = 0; iik < _grains[i].getPointList().size(); iik++){ 
          //       cout << "Point: " << _grains[i].getPointList()[iik](0) << " " << _grains[i].getPointList()[iik](1) << endl;
          //    }
          // }

          //Avoid degenerate melting before change (check all)
          bool nanp_before = false;
          cout << "Check Nan Points before!!!!" << endl;
          for (size_t ip = 0; ip < _grains[i].getPointList().size(); ip++)
          {
            if (isnan(_grains[i].getPointList()[ip](0)) || isnan(_grains[i].getPointList()[ip](1)))
            {
              nanp_before = true;
              break;
            }
          }  

          if (nanp_before)
          {
            cout << "WARNING: Nan Points before!!!!" << endl;
            continue;  //Skip to next healthy grain
          }

          cout << "Get Area for floe Melt Classif." << endl;
          vector<Vector2d> VecOrigin = _grains[i].getPointList();
          double Area = 0.0;
          size_t n = VecOrigin.size();

          for (size_t ii = 0; ii < n-1; ii++)
          {
            Area += ( VecOrigin[ii](0) * VecOrigin[ii+1](1) -  VecOrigin[ii](1) * VecOrigin[ii+1](0) ); 
          }
          Area += (VecOrigin[n-1](0) * VecOrigin[0](1) -  VecOrigin[n-1](1) * VecOrigin[0](0) ); 

          //Avoid degenerate melting
          if (isnan(Area))
          {
              cout << "WARNING: Nan AREA!!!!" << endl;
              continue; //Skip to next healthy grain
          }

          double IAreaMelt = 0.5*abs(Area);
          double IDiamMelt = (2 * sqrt(IAreaMelt /3.141592653589793));

          h_0 = _grains[i].getThickness();
          //Basal Melting
          //if (i == 37000)
          //{_grains[i].TemperatureModifBasal(Tair, Twater, _dt, dh, alphaice, meltTemp, meltid, _stepup, _START_TEMP);}  //Not needed for new melt mode

          double limitMaxDiam = limitMaxDiamV;
          double limitMinDiam = limitMinDiamV;
          //Lateral Melting
          //if (i == 37 || i == 36) {

          size_t melt_flag = 1; //Big (1) or small (2) floe flag for Temperature function  //2 // > 5.00 works for limitMinDiam
          if (IDiamMelt <= limitMinDiam) 
          {
            melt_flag = 2;
          } 
          //if ( _grains[i].getMass() < _MELT_MASS && _grains[i].getMass() > 0.009*_MELT_MASS){ // && IDiamMelt < 15.00){
          //cout << "LS Thermal Modification" << endl;
          
          //cout << "IN Sample Mass for i: " << i << " is: " << _grains[i].getMass() << endl;

          if (i >= 0)
          //if (i == 242) //23
          {
             //cout << "TemperatureModif for grain " << i << endl;
             //Change Aug 22, 2022
             //BEGIN
             double loss_mcv_temp = 0.0;
             double loss_mcv_solar_temp = 0.0;
             double loss_mcv_ocean_temp = 0.0;
             hlat = 0.0;
             _grains[i].TemperatureModif(Tair, Twater, _dt, dh, alphaice, meltTemp, meltid, meltV, _stepup, KIc, afactor, fdim, _START_TEMP, Khor, meltVSun, melt_flag, dstep, _fluid_coord, _x_cells, _y_cells, _oceanTemp, _offset, loss_mcv_temp, loss_mcv_solar_temp, loss_mcv_ocean_temp, hlat); 
             loss_mcv += loss_mcv_temp;
             loss_mcv_solar += loss_mcv_solar_temp;
             loss_mcv_ocean += loss_mcv_ocean_temp;
             //END
             
             //_grains[i].TemperatureModif(Tair, Twater, _dt, dh, alphaice, meltTemp, meltid, meltV, _stepup, KIc, afactor, fdim, _START_TEMP, Khor, meltVSun, melt_flag, dstep, _fluid_coord, _x_cells, _y_cells, _oceanTemp, _offset); 
             //cout << "End TemperatureModif for grain " << i << endl;
          }   //Utemper, Utemper0 
          //_grains[i].changeTemper(_grains[i].getMass()); 

          //cout << "OUT Sample Mass for i: " << i << " is: " << _grains[i].getMass() << endl;
          
          //Avoid degenerate melting after change (check all)
          bool nanp_after = false;
          cout << "Check Nan Points after!!!!" << endl;
          for (size_t ip = 0; ip < _grains[i].getPointList().size(); ip++)
          {
            if (isnan(_grains[i].getPointList()[ip](0)) || isnan(_grains[i].getPointList()[ip](1)))
            {
              nanp_after = true;
              break;
            }
          }  

          if (nanp_after)
          {
            cout << "WARNING: Nan Points after!!!!" << endl;
            continue;  //Skip to next healthy grain
          }

          //else{
              //double meltVs = 0.05*meltV;
              //_grains[i].TemperatureModif(Tair, Twater, _dt, dh, alphaice, meltTemp, meltid, meltVs, _stepup, KIc, afactor, fdim, _START_TEMP);    //Utemper, Utemper0 
              //_grains[i].changeTemper(_grains[i].getMass()); 
          //}


          //Final Area and Diameter after melting
          //double FAreaMelt = PointsArea(_grains[i].getPointList());
          //double FDiamMelt = MeanCaliperD(FAreaMelt);

          cout << "Check Nan Area after!!!!" << endl;
          VecOrigin = _grains[i].getPointList();
          Area = 0.0;
          n = VecOrigin.size();
          for (size_t ii = 0; ii < n-1; ii++)
          {
            Area += ( VecOrigin[ii](0) * VecOrigin[ii+1](1) -  VecOrigin[ii](1) * VecOrigin[ii+1](0) ); 
          }
          Area += (VecOrigin[n-1](0) * VecOrigin[0](1) -  VecOrigin[n-1](1) * VecOrigin[0](0) ); 

          //Avoid degenerate melting after change
          if (isnan(Area))
          {
              cout << "WARNING: Nan AREA after!!!!" << endl;
              continue; //Skip to next healthy grain
          }

          cout << "Work with sizes after!!!!" << endl;
          double FAreaMelt = 0.5*abs(Area);
          double FDiamMelt = (2 * sqrt(FAreaMelt /3.141592653589793));


            
          double AreaMelt = abs(IAreaMelt - FAreaMelt);  //No freezing happening only melting I > F
          double AreaMeltcheck = IAreaMelt - FAreaMelt;
          double Area_MB = 0.5*(IAreaMelt + FAreaMelt);
          double deltaA = AreaMelt;
          
          //Get lateral radial calculate  //delta_r = sqrt(A_old/pi) - sqrt(A_new/pi) = r_old - r_new = 0.5*(d_old - d_new)
          if (AreaMeltcheck > 0){
            delta_r_temp += 0.5 * (IDiamMelt - FDiamMelt);
            delta_r_vec =  0.5 * (IDiamMelt - FDiamMelt);
          }
          else{
            delta_r_temp += 0.0; 
            delta_r_vec = 0.0;
          }
          
          if (AreaMeltcheck < 0){
              AreaMeltcheck = 0;
              Area_MB = 0.0;
              deltaA = 0.0;
          } 

          cout << "Area Melt for reference: " << AreaMelt << endl;

          //Lateral melt loss area (else thickness melt)
          if (_grains[i].getThickness() <= 0.0){
            area_loss_basal += abs(IAreaMelt);
          }
          else{
            area_loss_lat += abs(AreaMeltcheck);
          }
          
          //Basal melt cumulative mass
          h_f = _grains[i].getThickness();
          double deltahB = 0.00;
          if (h_0 > h_f){
              deltahB = h_0 - h_f;
          }
          MBASAL += Area_MB * deltahB * rho_icem * 1000000.0;     //(km2 * m * kg/m3) -->  (km2 * m * kg/m3) * 1000^2 m2/1km2 --> 1000000.0 * (kg)
          
          //Lateral melt cumulative mass
          
          if (hlat <= 0.0){ //Check this
              hlat = deltahB;
          }
          MLAT += deltaA * hlat * rho_icem * 1000000.0; 
          
          //Find bin that lost this area (assume minor melt, staying on same bin)
          cout << "Work with bins lower!!!!" << endl;
          size_t loc_index = 0;
          size_t loc_index2 = 0;
          //April 24, 2023
          size_t locl = 0;
          size_t locg = 0;
          for (size_t bi = 0; bi < Diameters.size()-1; bi++){  
              if (IDiamMelt > Diameters[0]){
                  locl = 1000;  //Out of bin range (not even fines) //Max threshold has no info (avoid errors)
                  break;
              }
              else if ( IDiamMelt < Diameters[Diameters.size()-1] )  //Min size threshold control
              {
                 loc_index = Diameters.size()-2;
                 locl = Diameters.size()-2;
                 break;
              }
              else{
                  if ( IDiamMelt < Diameters[bi] &&  IDiamMelt >= Diameters[bi+1] ) 
                  {
                     loc_index = bi;
                     locl = bi;
                     break;
                  }
              }
              // else if ( IDiamMelt > Diameters[0] )  //Max size threshold control
              // {
              //    loc_index = 0;
              //    break;
              // }
          }
          cout << "Work with bins upper!!!!" << endl;
          for (size_t bi = 0; bi < Diameters.size()-1; bi++){  
              if (FDiamMelt > Diameters[0]){
                  locg = 1000;  //Out of bin range (not even fines) //Max threshold has no info (avoid errors)
                  break;
              }
              else if ( FDiamMelt < Diameters[Diameters.size()-1] )  //Min size threshold control
              {
                 loc_index2 = Diameters.size()-2;
                 locg = Diameters.size()-2;
                 break;
              }
              else{
                  if ( FDiamMelt < Diameters[bi] &&  FDiamMelt >= Diameters[bi+1] ) 
                  {
                     loc_index2 = bi;
                     locg = bi;
                     break;
                  }
              }
              // else if ( FDiamMelt > Diameters[0] )  //Max size threshold control
              // {
              //    loc_index2 = 0;
              //    break;
              // }
          }
          
          //April 24, 2023
          if (locl == 1000){
                if (locg != 1000){
                    NGMelt[locg] = NGMelt[locg] + 1;
                }
          }
          else{
              if (locl < locg && locg != 1000) //Only if there is a size transition, otherwise leave to Melty or FineMelty or Break
              {
                     
                    NLMelt[locl] = NLMelt[locl] + 1;
                    NGMelt[locg] = NGMelt[locg] + 1;
              }
          }
          
          //Check if grain changed bin as well, if it did then substract all mass from initial bin and move it to final bin
          cout << "Assign area to BINS!!!!" << endl;
          cout << "MeltBin.size(): " << MeltBin.size() << " loc_index: " << loc_index << " loc_index2: " << loc_index2 << endl;  
          if (loc_index != loc_index2){
            MeltBin[loc_index] += -IAreaMelt;
            MeltBin[loc_index2] += FAreaMelt;
          } 
          else{
            MeltBin[loc_index] += -AreaMelt;    //Otherwise if grain stayed in bin remove melt area
          }

          //cout << "MELT: AreaOriginal: " << IAreaMelt << " IDiam: " << IDiamMelt  << " Area new: " << FAreaMelt << " FDiam: " << FDiamMelt  << " loc: " << loc_index << " loc2: " << loc_index2 << endl;
          
          //Save bin lateral rate
          if (locl != 1000){
              deltaRV[locl] = deltaRV[locl] + delta_r_vec;
              nFloes[locl] =  nFloes[locl] + 1;   
          }

      //}
      //else  //Consider melting larger floes moooooore slowly
      //{
          //_grains[i].TemperatureModif(Tair, Twater, _dt, dh, alphaice, meltTemp, meltid, 0.05*meltV, _stepup, KIc, afactor, fdim, _START_TEMP);    //Utemper, Utemper0 
          //_grains[i].changeTemper(_grains[i].getMass()); 
      //}
    } 
    
    if (_ngrains > 0){
        ave_delta_r = delta_r_temp / _ngrains;
    }
    
    for (size_t i = 0; i < LatRate.size(); i++)
    {
        if (nFloes[i] > 0){
            LatRate[i] = LatRate[i] + (double(deltaRV[i])/double(nFloes[i]));
        }
        else{
            LatRate[i] = LatRate[i] + 0.0;
        }
    }

    cout<<"End of TempState" << endl;
    return; 
  }

  //Create a unique mesh of nodes and elements for each grain, init damage and use with SVL when needed
  void InitialDamageMesher()
  {    
    for (size_t i = 0; i < _ngrains; i++)
    { 
        //Output points of Grains to a .dat file to convert into mesh
        string outDir = "./Output/SeaIce/SVL/";  //DIR
        string fname = outDir + "g_points_init.dat";
        FILE * out_points = fopen(fname.c_str(), "w");

        //Use position to reference wrt to local coord, not global

        for (size_t ip = 0; ip < _grains[i].getnpoints(); ip++) 
        {
            fprintf(out_points,  "%4.8f %4.8f\n", _grains[i].getPointList()[ip](0)-_grains[i].getPosition()(0), _grains[i].getPointList()[ip](1) - _grains[i].getPosition()(1) );  //For pydist2dmesh
        }
        fclose(out_points);

        //Generate Grain Number Tag and TimeStep
        string fname2 = outDir + "g_tag_init.dat"; 
        FILE * outr = fopen(fname2.c_str(), "w");
        fprintf(outr, "%d\n", int(i));// grain Number Tag
        fprintf(outr, "%4.8f\n", _grains[i].getRadius());  // Grain Radius to Scale Mesh Adequately
        fclose(outr);   

        //Generate point node file from LSDEM wrt to Local coord syste
        std::string filename = "points2mesh.py"; 
        std::string command = "python3 ";
        command += filename;
        system(command.c_str());

    //     //Import text
    //     char tempfname[300];
    //     char tempfname2[300];
    //     std::string nngg = std::__cxx11::to_string(i); //int(i);
    //     //tempfname = "./Output/SeaIce/Mesh/PointMesh_" + nngg + ".dat";
    //     //tempfname2 = "./Output/SeaIce/Mesh/npoints_" + nngg + ".dat";

    //     sprintf(tempfname, ("./Output/SeaIce/Mesh/PointMesh_" + nngg + ".dat").c_str());  //DIR
    //     sprintf(tempfname2, ("./Output/SeaIce/Mesh/npoints_" + nngg + ".dat").c_str());   //DIR
        
    //     std::string file_point = tempfname;
    //     std::string file_num = tempfname2;

    //     //Generate MatrixXd Object for flexibility and to export to other functions
    //     vector<Vector2d> point_Mesh = readPointFile(file_point, file_num);
    //     vector<Vector3d> initDamage(point_Mesh.size());

    //     for (size_t ip = 0; ip < point_Mesh.size(); ip++) 
    //     {
    //         initDamage[ip](0) = point_Mesh[ip](0);
    //         initDamage[ip](1) = point_Mesh[ip](1);
    //         initDamage[ip](2) = 0.0;
    //     }

    //     _grains[i].changeDamage(initDamage);
    }

  }

  //Zero Damage Init
  void DamageInit()
  {
    for (size_t i = 0; i < _ngrains; i++)
    { 

        //Import text
        char tempfname[300];
        char tempfname2[300];
        std::string nngg = std::to_string(i); //int(i);
        //tempfname = "./Output/SeaIce/Mesh/PointMesh_" + nngg + ".dat";
        //tempfname2 = "./Output/SeaIce/Mesh/npoints_" + nngg + ".dat";

        sprintf(tempfname, ("./Output/SeaIce/Mesh/PointMesh_" + nngg + ".dat").c_str());  //DIR
        sprintf(tempfname2, ("./Output/SeaIce/Mesh/npoints_" + nngg + ".dat").c_str());   //DIR
        
        std::string file_point = tempfname;
        std::string file_num = tempfname2;

        //Generate MatrixXd Object for flexibility and to export to other functions
        vector<Vector2d> point_Mesh = readPointFile(file_point, file_num);
        vector<Vector3d> initDamage(point_Mesh.size());

        for (size_t ip = 0; ip < point_Mesh.size(); ip++) 
        {
            initDamage[ip](0) = point_Mesh[ip](0);
            initDamage[ip](1) = point_Mesh[ip](1);
            initDamage[ip](2) = 0.0;
        }
        _grains[i].changeDamage(initDamage);
    }
  }

  void StressStateForceFinder(const int & step_time)
  {

    //Important Initializations
    double rho_ice = 910;                    //Density of ice
    double rho_air = 1.23e-3;                //Density of air If rhoice is 0.91e-6;
    double rho_water = 1025;                  //Density of water

    double Tair = -5; //Degrees Celsius
    double Twater = 20.; //Degrees Celsius //50
    double Ticeedge = -20.; //Degrees Celcius

    //For Internal Stresses 
    double A_rheo = 0.001;

    // //For flow velocity adjust
    size_t flow_step;  //TODO: control at main.cpp  //60  //3 //6

    if (_stepup > 0 ) //Only 1????
    {
      flow_step = 1; //3
    } 

    //Generate a CData vector for each ice grain using other grains and shape walls
    //If grain ith CData component has loc.size() > 0, then we have forces applied and we use all these force to get internal stress, else stress is aprox. zero
    //We can also filter the two largest forces
    //Use stress to damage and break
    //Combine with thickness later on
 
    cout <<"Begin Stress Analysis" << endl;   
    for (size_t i = 0; i < _ngrains; i++){ 
        CData grainContactF;  //Temporary CData for all forces on grain
        vector<size_t> ptCIdx; //Vector that specifies exact point of force application (for Convenience)
        vector<size_t> ptCIdx2; //Vector that specifies exact point of force application (for Convenience)
        vector<size_t> ptCIdx3; //Vector that specifies exact point of force application (for Convenience)
        vector<size_t> ptCIdxF; //Vector that specifies exact point of force application (for Convenience)
        for (size_t j = i+1; j < _ngrains; j++) {
            // Normal contacts
            if (_grains[i].bcircleGrainIntersection(_grains[j])){
                CData cDataContact = _grains[i].findContactDataIDX(_grains[j], _dt, 0, ptCIdx);
                if (cDataContact._clocs.size() > 0) {
                    grainContactF += cDataContact;
                }
            }
            // Contacts due to periodic bcs
            else if (_grains[i].bcircleGrainIntersectionXOffset(_grains[j], _offset(0))){
                CData cDataContact = _grains[i].findContactDataIDX(_grains[j], _dt, _offset(0), ptCIdx2);
                if (cDataContact._clocs.size() > 0) {
                     grainContactF += cDataContact;
                }
            }
        }
        
        if (_grainsWall.size()>0)
        {    
            for (size_t j = 0; j < _ngrainsWall; j++) {
                // Normal contacts
                if (_grains[i].bcircleGrainIntersection(_grainsWall[j])){
                    CData cDataContact = _grains[i].findContactDataIDX(_grainsWall[j], _dt, 0, ptCIdx3);
                    if (cDataContact._clocs.size() > 0) {
                         grainContactF += cDataContact;
                    }
                }
            } // close grain subloop
        }
    
        //Join all three index vectors
        for (size_t k = 0; k < ptCIdx.size(); k++){
            ptCIdxF.push_back(ptCIdx[k]);
        }
        for (size_t k = 0; k < ptCIdx2.size(); k++){
            ptCIdxF.push_back(ptCIdx2[k]);
        }
        for (size_t k = 0; k < ptCIdx3.size(); k++){
            ptCIdxF.push_back(ptCIdx3[k]);
        }


        cout << "Contact Location #: " << grainContactF._clocs.size() << endl;
        if (_grains[i].getVelocity().norm() < 10000){  //Laxer Vel. Req.  //Modify based on BCs
        //if (grainContactF._clocs.size()>2 && _grains[i].getVelocity().norm() < 1000){  //Laxer Vel. Req.  //
        //if (grainContactF._clocs.size()>2 && _grains[i].getVelocity().norm() < 1000 && i >5 && i <11){  //Laxer Vel. Req.  //
        //if (grainContactF._clocs.size()>2 && _grains[i].getVelocity().norm() < 1000 ){  //Laxer Vel. Req.  //
        //if (grainContactF._clocs.size()>2 && _grains[i].getVelocity().norm() < 0.01 ){

            string outDir2 = "./Output/SeaIce/ContactInfo/Detail/";  //DIR
            string fname3 = outDir2 + "cinfoALL_iter_" + std::to_string(_stepup/70000) + "_grain_" + std::to_string(i) + "_0.dat";  //TODO: CHANGE 70000
            FILE * cinfo3 = fopen(fname3.c_str(), "w");
            for (size_t ii = 0; ii < grainContactF._clocs.size(); ii++) {
                //cout << "GrainForces on grain i: " << i << " FX: " <<  grainContactF._forces[ii](0) << " FY: " << grainContactF._forces[ii](0) << " exerted by grain j: " << grainContactF._cpairs[ii](1) << " at Coord X : " << grainContactF._clocs[ii](0) << " and Y: " << grainContactF._clocs[ii](1) << " or Point idx: " << ptCIdxF[ii] <<endl;
            
                fprintf(cinfo3, "%d %d ", grainContactF._cpairs[ii](0), grainContactF._cpairs[ii](1) ); // grains in contact
                fprintf(cinfo3, "%.8f %.8f ", grainContactF._forces[ii](0), grainContactF._forces[ii](1));// force x and y
                //fprintf(cinfo3, "%.8f %.8f ",cStateAll._normals[i](0), cStateAll._normals[i](1));// n
                fprintf(cinfo3, "%.8f %.8f ", grainContactF._clocs[ii](0), grainContactF._clocs[ii](1));// loc x and y in grain
                fprintf(cinfo3, "%d\n" , int(ptCIdxF[ii]) ); // Pointlist location
            }
            for (size_t jj = 0; jj < _grains[i].getnpoints(); jj++) 
            {
                fprintf(cinfo3,  "%d %4.8f %4.8f \n", int(jj), _grains[i].getPointList()[jj](0), _grains[i].getPointList()[jj](1) );
            }
            fclose(cinfo3);    

            //Generate Forces
            string outDir3 = "./Output/SeaIce/SVL/";  //DIR
            string fname4 = outDir3 + "g_forceP.dat"; 
            FILE * cinfo4 = fopen(fname4.c_str(), "w");
            for (size_t ii = 0; ii < grainContactF._clocs.size(); ii++) {            
                fprintf(cinfo4, "%.8f %.8f ", grainContactF._forces[ii](0), grainContactF._forces[ii](1));// force x and y
                fprintf(cinfo4, "%.8f %.8f\n", grainContactF._clocs[ii](0), grainContactF._clocs[ii](1));// loc x and y in grain
                //fprintf(cinfo3, "%d\n" , int(ptCIdxF[ii]) ); // Pointlist location
            }
            fclose(cinfo4);    

            //Generate Grain Number Tag and TimeStep
            string fname5 = outDir3 + "g_tag.dat"; 
            FILE * cinfo5 = fopen(fname5.c_str(), "w");
          
            fprintf(cinfo5, "%d\n", int(i));// grain Number Tag
            fprintf(cinfo5, "%d\n", int(step_time));// Time Step
            fprintf(cinfo5, "%4.8f\n", _grains[i].getRadius());// Grain Radius to Scale Mesh Adequately
            fprintf(cinfo5, "%4.8f\n", _grains[i].getPosition()(0) ); //Position X
            fprintf(cinfo5, "%4.8f\n", _grains[i].getPosition()(1) ); //Position Y 
            fprintf(cinfo5, "%4.8f\n", _grains[i].getTheta() ); //Rotation Angle too           
            fclose(cinfo5);                

            //Generate Point List and Generate SVL Code
            
            //Output points of Grains to a .dat file to convert into mesh
            string fname6 = outDir3 + "g_points.dat";
            FILE * out_points = fopen(fname6.c_str(), "w");

            for (size_t ip = 0; ip < _grains[i].getnpoints(); ip++) 
            {
                fprintf(out_points,  "%4.8f %4.8f\n", _grains[i].getPointList()[ip](0), _grains[i].getPointList()[ip](1) );  //For pydist2dmesh
            }
            fclose(out_points);

            //Use pygmsh script dat2gmsh.py to convert .dat file to geo and gmsh file
            //std::string filename = "/Users/rigobertomoncadalopez/Dropbox/Caltech_PhD/Research/Code/LSDEM2D-SeaIceDenseFrac/dat2gmsh.py";
            //std::string filename2 = "/Users/rigobertomoncadalopez/Dropbox/Caltech_PhD/Research/Code/LSDEM2D-SeaIceDenseFrac/Output/SeaIce/g_points.dat";
            // std::string filename = "dat2gmsh.py ";
            // std::string filename2 = "g_points.dat";
            // std::string command = "python3 ";
            // command += filename;
            // command += filename2;
            // system(command.c_str());

            //Generate py SVL file from LSDEM
            std::string filename = "points2svl_tag.py"; 
            std::string command = "python3 ";
            command += filename;
            system(command.c_str());

            //Set up Python Export Path
            std::string commandex = "export PYTHONPATH=/Users/rigobertomoncadalopez/Dropbox/Caltech_PhD/Research/Code/LSDEM2D-SeaIceDenseFrac/SeismoVLAB_V1.15_SERIAL/01-Pre_Process";
            system(commandex.c_str());

            //Run Py file to Prepare .json and then use py file itself to execute in SVLAB and get Solution file with Stresses (also vtk)    
            std::string nngg = std::to_string(i);//int(i);
            std::string ttss = std::to_string(step_time);//int(step_time)
            //cout << "Grain Numba: " << nngg << " and time_step: " << ttss << endl;
            std::string filename2 = outDir3 + "Grain_" + nngg + "_iter_"+ ttss +".py";
            std::string command2 = "python3 ";
            command2 += filename2;
            system(command2.c_str());

            //Cut solution file into a pure X, Y, Sigma XX, Sigma YY, Sigma XY version
            //std::string filenamecut = "cut_stress.py"; 
            std::string filenamecut = "cut_stress2.py";  //Faster
            std::string commandcut = "python3 ";
            commandcut += filenamecut;
            system(commandcut.c_str());

            //Convert FEM output file into Stress Object that can be fed to damage info and into LS-DEM in general 
            //Import text
            char tempfname[300];
            char tempfname2[300];
            sprintf(tempfname, "./Output/SeaIce/SVL/Solution/IceTria3/stress_results.dat");  //DIR
            sprintf(tempfname2, "./Output/SeaIce/SVL/Solution/IceTria3/stress_nelem.dat");   //DIR
            std::string file_stress = tempfname;
            std::string file_num = tempfname2;

            //Generate MatrixXd Object for flexibility and to export to other functions
            //MatrixXd stress_M = readStressFile(file_stress, file_num);
            MatrixXd stress_M = readStressFile2(file_stress, file_num);

            //Print for Debug
            cout << "Stress matrix for Trial DEBUG" << endl;
            cout << stress_M << endl;


            cout <<"Damage Evolution" << endl;
            //Damage Evolution
            double Dc = 0.9; //Critical Damage
            //origin B = 0.003
            //Fracture toughness 100-140 kPa m^1/2 (Timco & Frederking, 1982)
            double B_damage = 0.00030;    //MULTI 0.006;    //v1 0.002 // 0.006 //0.00006   //0.0002 4 pieces  //0.00002     //Effective damage rate 65 MPa^-r a^-1     0.001 is slow enough for 4900 taucrit, 40 steps, 7000 show gradual fail 1  ///  0.006 is fast for tauc 2700 and w_factor 600 //Orig 0.003 for 100 time steps  //0.001
            //delta t will be one for the mean time !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
            //Let's use tau_crit as Damage Threshold     0.17 MPa 

            double r_damage = 0.43;  //Damage Law Exponent 0.43
            double k_damage = 0.5;  //Other Damage Law Exponent, look for reference value, Mercenier 2018 does not specify it.
            double numer_damage; //Useful
            double denom_damage;  //Useful
            double taucrit = 10; //Critical value for Damage Onset

            vector<Vector3d> Old_Damage = _grains[i].getgrainDamage();
            vector<Vector3d> New_Damage = Old_Damage;

            size_t SizeM = New_Damage.size();
            if (New_Damage.size() == stress_M.rows())
            {
                cout <<"GOOD: Stress and Damage Size MATCH!!!" << endl;
            }
            else if (New_Damage.size() > stress_M.rows())
            {
                cout << "DS: " <<  New_Damage.size() <<  " SS: " << stress_M.rows()  << endl;
                cout <<"WARNING: Stress and Damage Size do not match, Not all stress will be represented" << endl;
                SizeM = stress_M.rows();
            }
            else
            {
                cout << "DS: " <<  New_Damage.size() <<  " SS: " << stress_M.rows()  << endl;
                cout <<"WARNING: Stress and Damage Size do not match, Segmentation Fault Inminent unless..." << endl;
                SizeM = New_Damage.size();
            }    

            for (size_t jj = 0; jj < SizeM; ++jj)   
            {
                if (Old_Damage[jj](2) >= 0.0)
                {  
                    //We will NOT assume healing. Damage can only stay the same or increase
                    //if ( abs(stress_M(jj,4)) > taucrit  &&  Old_Damage[jj](2) < Dc )  //EVALUATE XY Stress
                    if ( (stress_M(jj,2)) > taucrit  &&  Old_Damage[jj](2) < Dc )  //EVALUATE XY Stress //abs?
                    {
                         //numer_damage = abs(stress_M(jj,4)) - taucrit;
                         numer_damage = (stress_M(jj,2)) - taucrit; //abs?
                         denom_damage = 1 -  Old_Damage[jj](2);

                         New_Damage[jj](2) =  Old_Damage[jj](2) + (  B_damage * pow (numer_damage , r_damage)  / pow (denom_damage, (k_damage + r_damage)) );
                        
                         //Bound damage to Critical Damage (afterwards it does not really matter because, for visualization purposes )
                         if (New_Damage[jj](2) > Dc)
                         {
                            New_Damage[jj](2) = Dc;
                         }
                    }
                    else
                    {
                         //Nothing happens (damage stays the same as before unless the cell melts)
                    }
                }  
            }      
            _grains[i].changeDamage(New_Damage);

        }

        //Calculate Stress using all these Forces for grain i using Data Above!!!!
        //vector<double> newstress = _grains[i].InnerStress(grainContactF); //Output is a level set vector
        //_grains[i].changeInternalStress(newstress, _grains[i].getgrainStress2D().getXdim(), _grains[i].getgrainStress2D().getYdim())
    }  
    cout <<"End Stress Analysis" << endl;  


    // for (size_t i = 0; i < _ngrains; i++){ 
    //   if (i == 0)
    //   {
    //      //Always soften level set before Stress Calculation to assure correct BCs. Weird cells can arise from breakage or melting
    //     size_t xxdim =  _grains[i].getgrainTemp2D().getXdim();  // To get Heat LSET
    //     size_t yydim =  _grains[i].getgrainTemp2D().getYdim();
    //     vector<double> geovec = _grains[0].getLset().getLevelset();  //Level Set Vector
    //     //Go to Matrix
    //     MatrixXd Cls_Tag0 = MatrixDefVec(yydim, xxdim, geovec);
    //     //Soften Level Set
    //     size_t maxttries = 50;
    //     MatrixXd Geo_Matrix = LS_Soft(Cls_Tag0, maxttries); 
    //     //Go Back to Vector
    //     geovec = MatrixtoVec(Geo_Matrix);
    //     _grains[i].change_Full_lset(geovec, xxdim, yydim);

    //     //CALCULATE INNER STRESS AND FLOW/CHANGE IN GEOMETRY
    //     if ( (_stepup + 1) % flow_step == 0 && _stepup > 10) //28   //34
    //     {
    //       //cout << "Finding Flow Stress for grain: " << 0 << " at step: " << _stepup <<endl; 
    //       //_grains[i].StressModifV(rho_ice, rho_water, Tair, Twater, Ticeedge, A_rheo, _yw, 1, _stepup);
    //     } 

    //     //RECALCULATE ONLY STRESS BUT NOT FLOW IN ORDER TO GO TO FRACTURE/DAMAGE
    //     cout << "Finding Fixed Stress for grain: " << i << " at step: " << _stepup <<endl; 
    //     _grains[i].StressModifElastic(rho_ice, rho_water, Tair, Twater, Ticeedge, A_rheo, 0, _stepup);
    //   }
    //   else
    //   {

    //   }
    // } 
  }

  void DamageState()
  {
    for (size_t i = 0; i < _ngrains; i++){ 
      if (i == 0){
        
      }
      else
      {

      }
    } 
  }

  void FractureState()
  {
    
    size_t try_index = 92;

    cout << "Original Position for Testing for grain: " << try_index << endl;
    cout << "X: "  << _grains[try_index].getPosition()[0] << " Y: " << _grains[try_index].getPosition()[1] << endl;

    Vector2d velOrig(_grains[try_index].getVelocity()[0] , _grains[try_index].getVelocity()[1]);
    double tempOrig  = _grains[try_index].getTemperature();
    double thickOrig = _grains[try_index].getThickness();
    double thetaOrig = _grains[try_index].getTheta();
    double omegaOrig = _grains[try_index].getOmega();

    //Print temporal point file
    string outDir = "./Output/SeaIce/Temporal/";
    string fname = outDir + "Temp_Points.dat";
    string fname2 = outDir + "Temp_CM.dat";
    string fname3 = outDir + "Temp_POS.dat";
    FILE * points_export = fopen(fname.c_str(), "w");
    
    for (size_t i = 0; i < _grains[try_index].getPointList().size(); i++){ 
        fprintf( points_export,  "%4.8f %4.8f\n", _grains[try_index].getPointList()[i](0), _grains[try_index].getPointList()[i](1) );
    }   
    fclose(points_export);

    FILE * CM_export = fopen(fname2.c_str(), "w");
    fprintf(CM_export,  "%4.8f %4.8f\n", _grains[try_index].getCmLset()(0), _grains[try_index].getCmLset()(1) );
    fclose(CM_export);

    FILE * POS_export = fopen(fname3.c_str(), "w");
    fprintf(POS_export,  "%4.8f %4.8f\n", _grains[try_index].getPosition()(0), _grains[try_index].getPosition()(1) );
    fclose(POS_export);

    //Convert point file to new temporal morphology //This prints a temporal grainproperty0 only for point processing.
    std::string filename = "/Users/rigobertomoncadalopez/Dropbox/Caltech_PhD/Research/Code/LSDEM2D-SeaIceDense/Points_to_LS.py";
    std::string command = "python3 ";
    command += filename;
    system(command.c_str());


    //Grains  
    char tempfname[300];  
    sprintf(tempfname, "./Output/SeaIce/Temporal/Temp_morphologies.dat");
    string file_morph = tempfname;
    sprintf(tempfname, "./Output/SeaIce/Temporal/Temp_positions.dat"); //Should come from point interpolation centroid and respective position file 
    string file_pos = tempfname;
    sprintf(tempfname, "../Output/SeaIce/Temporal/Temp_velocities.dat");  //RE-UPDATE
    string file_vel = tempfname;
    sprintf(tempfname, "./Output/SeaIce/Temporal/");
    string morph_dir = tempfname;
    //Get temperature-related quantitities
    sprintf(tempfname, "./Output/SeaIce/Temporal/Temp_temper.dat"); //RE-UPDATE
    string file_temper = tempfname;
    sprintf(tempfname, "./Output/SeaIce/Temporal/Temp_thick.dat"); //RE-UPDATE
    string file_thick = tempfname;


    //Convert temporal morphology to fractured new grain
    vector<Grain2d> grainsTemp = generateGrainsFromFiles(file_morph, morph_dir, file_pos, file_vel, file_temper, file_thick);
    _grains[try_index] = grainsTemp[0]; //Update grain except for:

    //Recover Velocity, Thickness and Temperature, Theta Rotation and Omega Angular Velocity
    _grains[try_index].changeVelocity(velOrig);
    _grains[try_index].changeTemper(tempOrig);
    _grains[try_index].changeThickness(thickOrig);
    _grains[try_index].changeTheta(thetaOrig);
    _grains[try_index].changeOmega(omegaOrig);

    cout << "Properties New Grain: " << try_index << endl;
    cout << "Position X: " << _grains[try_index].getPosition()[0] << " Y: " << _grains[try_index].getPosition()[1] << endl;
    cout << "Velocity X: " << _grains[try_index].getVelocity()[0] << " Y: " << _grains[try_index].getVelocity()[1] << endl;
    cout << "CM X: " << _grains[try_index].getCmLset()[0] << " Y: " << _grains[try_index].getCmLset()[1] << endl;
    cout << "Radius: " << _grains[try_index].getRadius() << endl;


    for (size_t i = 0; i < _ngrains; i++){ 
      if (i == 0){
        
      }
      else
      {

      }
    } 
  }  


  //Function that eliminates grains if they melt
  //Change Aug 22, 2022
  void Melty(vector<double> & MeltBin, vector<double> & Diameters, double & limitMinDiamV, vector<Vector4d> & Fine_Grains, vector<Vector2d> & Fine_Velocities, double & loss_mcl, double & loss_mcv, double & gain_fines, vector<size_t> & NLMelt, vector<size_t> & NLBreak, vector<size_t> & NGMelt, vector<size_t> & NGBreak, double & loss_mcv_solar, double & loss_mcv_ocean, double & area_loss_bkg)
  //void Melty(vector<double> & MeltBin, vector<double> & Diameters, double & limitMinDiamV, vector<Vector4d> & Fine_Grains, vector<Vector2d> & Fine_Velocities)
  {
    bool coarse_enough = true; //Assume ok to continue
    if (_ngrains < 2)
    {
      coarse_enough = false; //Make melty unusable to avoid problems
    }

    double min_mass = _MELT_MASS*0.01;   //10e-2 of melting order (try 1) or 1/30 of melt_mass (2)    
    
    if (coarse_enough)
    {
      for (size_t i = 0; i < _ngrains; i++){  

        vector<Vector2d> VecOrigin = _grains[i].getPointList();
        double Area = 0.0;
        size_t n = VecOrigin.size();

        for (size_t ii = 0; ii < n-1; ii++)
        {
          Area += ( VecOrigin[ii](0) * VecOrigin[ii+1](1) -  VecOrigin[ii](1) * VecOrigin[ii+1](0) ); 
        }
        Area += (VecOrigin[n-1](0) * VecOrigin[0](1) -  VecOrigin[n-1](1) * VecOrigin[0](0) ); 

        double IAreaMelt = 0.5*abs(Area);
        double IDiamMelt = (2 * sqrt(IAreaMelt /3.141592653589793));

        //if (_grains[i].getMass() < min_mass && _ngrains>=2){ //Safeguard to prevent segmentation failure
        //if (IDiamMelt <= limitMinDiamV || _grains[i].getMass() == 0){

        //Melty will process grains in 2 different ways: Deleting very small grains and sending them to fines OR deleting completely a grain that reaches zero thickness regardless of size
        
        //First we check thickness (in this case nothing is sent to fines)
        if (_grains[i].getThickness() <= 0.0 || _grains[i].getMass() == 0.0){  //Use thickness to decide melt

          cout << "Grain L0ST due to thickness" << endl;
          //Find area lost for each melted grain if it melts and then in which bin it belongs
          //double AreaMelt = PointsArea(_grains[i].getPointList());
          //double DiamMelt = MeanCaliperD(AreaMelt);

          vector<Vector2d> VecOrigin = _grains[i].getPointList();
          double Area = 0.0;
          size_t n = VecOrigin.size();

          for (size_t ii = 0; ii < n-1; ii++)
          {
            Area += ( VecOrigin[ii](0) * VecOrigin[ii+1](1) -  VecOrigin[ii](1) * VecOrigin[ii+1](0) ); 
          }
          Area += (VecOrigin[n-1](0) * VecOrigin[0](1) -  VecOrigin[n-1](1) * VecOrigin[0](0) ); 

          if (isnan(Area) == 1)
          {
            cout << "WARNING: Nan in Melty algorithm for grain: " << i << endl;
            Area = 1.0; //Just to debug
          }

          double AreaMelt = 0.5*abs(Area);
          double DiamMelt = (2 * sqrt(AreaMelt /3.141592653589793));

          size_t loc_index;
          size_t locC_index; //April 24, 2023 FSD
          //Find bin that lost this area
          for (size_t bi = 0; bi < Diameters.size()-1; bi++){  
               if ( DiamMelt > Diameters[0] )  //Max size threshold control
               {
                  loc_index = 0;
                  locC_index = 1000;
                  break;
               }
              else if ( DiamMelt < Diameters[Diameters.size()-1] )  //Min size threshold control
              {
                 loc_index = Diameters.size()-2;
                 locC_index = Diameters.size()-2;
                 break;
              }
              else{
                  if ( DiamMelt < Diameters[bi] &&  DiamMelt >= Diameters[bi+1] ) 
                  {
                     loc_index = bi;
                     locC_index = bi;
                     break;
                  }
              }
          }
          
          if (locC_index != 1000){
            NLMelt[locC_index] = NLMelt[locC_index] + 1;
          }
          
          //Substract from lost bin
          MeltBin[loc_index] += -AreaMelt;
          
          //Accumulate basal area loss
          //area_loss_basal += AreaMelt; //No need already accounted for in TempState after basal melt made use declare _thickness = 0, if so in Melty we only remove the grain. Different because we are calculating area, not mass.
          
          //Change Aug 22, 2022
          //BEGIN
          loss_mcv += AreaMelt * max(_grains[i].getThickness(), 0.0) * 0.001 * _grains[i].getDensity();
          loss_mcv_solar += 0.5* AreaMelt * max(_grains[i].getThickness(), 0.0) * 0.001 * _grains[i].getDensity(); //Small mass, will redistribute
          loss_mcv_ocean += 0.5* AreaMelt * max(_grains[i].getThickness(), 0.0) * 0.001 * _grains[i].getDensity(); //Small mass, will redistribute
          //cout << " melt component melty: " << loss_mcv << endl;
          //cout << "S melt component melty: " << loss_mcv_solar << endl;
          //cout << "O melt component melty: " << loss_mcv_ocean << endl;
          //END
          
          //Add to fines (last bin)  //Not added because thickness is zero
          //MeltBin.back() += AreaMelt;

          cout << "MELTY: AreaLOST: " << AreaMelt << " DiamLost: " << DiamMelt << " Thickness: " << _grains[i].getThickness() << " Mass: " << _grains[i].getMass() <<   " loc: " << loc_index << endl;

          //cout << "Surely Delete Grain" << i <<endl; 
          _grains.erase(_grains.begin()+i);
          _ngrains--;
          //cout << "Delete Grain end" << endl; 

          _globalGrainState.reset();
          _globalGrainState.resize(_ngrains, _nwalls);
          //cout << "begin compute world again to be sure" << endl;
          computeWorldState();
          //cout << "end compute world again really sure " << endl;
          //return;
        }
   
        //Second: We check if the grains is small or if NaN Values (under LS-DEM Level Set Stable resolution capacity) and will be sent to fines
        //else if (_grains[i].getRemove()) //We use the convenient remove flag
        else if (_grains[i].getRemove() || isnan(_grains[i].getPosition()(0)) || isnan(_grains[i].getPosition()(1)) || isnan(_grains[i].getPointList()[0](0)) ) //We use the convenient remove flag
        {
          cout << "Grain L0ST due to Small Size" << endl;
          vector<Vector2d> VecOrigin = _grains[i].getPointList();
          double Area = 0.0;
          size_t n = VecOrigin.size();

          for (size_t ii = 0; ii < n-1; ii++)
          {
            Area += ( VecOrigin[ii](0) * VecOrigin[ii+1](1) -  VecOrigin[ii](1) * VecOrigin[ii+1](0) ); 
          }
          Area += (VecOrigin[n-1](0) * VecOrigin[0](1) -  VecOrigin[n-1](1) * VecOrigin[0](0) ); 

          if (isnan(Area) == 1)
          {
            cout << "WARNING: Nan in Melty algorithm for grain: " << i << endl;
            Area = 1.0; //Just to debug
          }

          double AreaMelt = 0.5*abs(Area); 
          double DiamMelt = (2 * sqrt(AreaMelt /3.141592653589793));

          size_t loc_index;
          size_t locC_index; //April 24, 2023
          //Find bin that lost this area
          for (size_t bi = 0; bi < Diameters.size()-1; bi++){  
              if ( DiamMelt >= Diameters[0] )  //Max size threshold control
              {
                  loc_index = 0;
                  locC_index = 1000;
                  break;
              }
              else if ( DiamMelt < Diameters[Diameters.size()-1] )  //Min size threshold control
              {
                 loc_index = Diameters.size()-2;
                 locC_index = Diameters.size()-2;
                 break;
              }
              else{
                  if ( DiamMelt < Diameters[bi] &&  DiamMelt >= Diameters[bi+1] ) 
                  {
                     loc_index = bi;
                     locC_index = bi;
                     break;
                  }
              }
          }
          
          if (locC_index != 1000){
            NLMelt[locC_index] =  NLMelt[locC_index] + 1;
          }
          NGMelt[Diameters.size()-1] = NGMelt[Diameters.size()-1] + 1;
          
          //Substract from lost bin
          MeltBin[loc_index] += -AreaMelt;
          
          //Accumulate breakage area loss (will assume melty coarse to fine, is for sufficiently small broken coarse grains. Lateral melt will exclusively come from areal differences)
          area_loss_bkg += AreaMelt;
          
          //Change Aug 22, 2022
          //BEGIN
          loss_mcl += AreaMelt * max(_grains[i].getThickness(), 0.0) * 0.001 * _grains[i].getDensity();
          gain_fines += AreaMelt * max(_grains[i].getThickness(), 0.0) * 0.001 * _grains[i].getDensity();
          //END
          
          //Add to fines (last bin)
          //MeltBin.back() += AreaMelt; //Not anymore since MeltBin will be assembled with Fine_Grains

          //So add to fine grains
          //double loss_factor = 0.01; //How much mass/area remains when breaking //CHECK!!!
          Vector4d New_Fine;
          Vector2d New_Vel;
          cout << "Inherited Thickness from Coarse: " << _grains[i].getThickness() << endl;
          New_Fine << AreaMelt, _grains[i].getThickness(), _grains[i].getPosition()(0), _grains[i].getPosition()(1);
          New_Vel << _grains[i].getVelocity()(0) , _grains[i].getVelocity()(1); //Transfer same velocity from original grain
          Fine_Grains.push_back(New_Fine); //Add to end of the list our new small grain
          Fine_Velocities.push_back(New_Vel); //Add to end of the list our new velocity for this small grain


          cout << "MELTY: AreaLOST: " << AreaMelt << " DiamLost: " << DiamMelt << " Thickness: " << _grains[i].getThickness() << " Mass: " << _grains[i].getMass() <<   " loc: " << loc_index << endl;

          //cout << "Surely Delete Grain" << i <<endl; 
          _grains.erase(_grains.begin()+i);
          _ngrains--;
          //cout << "Delete Grain end" << endl; 

          _globalGrainState.reset();
          _globalGrainState.resize(_ngrains, _nwalls);
          //cout << "begin compute world again to be sure" << endl;
          computeWorldState();
          //cout << "end compute world again really sure " << endl;
          //return;
        }   
        
        //Otherwise Melty does nothing and the grain stays intact, for now.
      }  
    }  
    // //Now let's repeat the process for Fine_Grains which will only melt by Thickness 
    // vector<size_t> delete_flag;
    // for (size_t i = 0; i < Fine_Grains.size(); i++)
    // {  
    //   if (Fine_Grains[i](1) <= 0.01)  //If thickness is zero  //0.000001
    //   {
    //      // Fine_Grains[i](0) = 0.0; //Make its mass go down to zero (CHECK PERFORMANCE OF THIS)
    //      // Fine_Grains[i](1) = 0.0;
    //      delete_flag.push_back(i); //Delete flag vector
    //   }
    // }
    // size_t del_counter = 0; //To shift as we erase positions, only applies if Fines reach critical low thickness
    // for (size_t i = 0; i < delete_flag.size(); i++) //Now loop through flag vector to modify your original Fines (faster)
    // { 
    //   Fine_Grains.erase(Fine_Grains.begin()+i-del_counter); //Delete objects in ascending order of index, so shift as more objects are deleted
    //   del_counter++; //Shift as our Fine_Grains vector shrinks to delete the right elements
    // }

    // //Alternative Method if Deleting becomes a problem (slower)
    // vector<Vector4d> New_Fine_Grains; //Save remaining fines here
    // Vector4d temporalV; //Temporal vector for saving fines
    // for (size_t i = 0; i < Fine_Grains.size(); i++)
    // {
    //   if (Fine_Grains[i](1) > 0.01)  //If thickness is zero  //0.000001 //Then the fine is not saved anymore
    //   {
    //      temporalV = Fine_Grains[i]; //We will preserve this fine
    //      New_Fine_Grains.push_back(temporalV); //Add to new fine list
    //   }
    // } 
    // if (New_Fine_Grains.size() < 1) // Add an artifical fine to avoid problems if list is empty
    // {
    //    temporalV[0] = 0.5;  temporalV[0] = 0.015;  temporalV[0] = 0.01;  temporalV[0] = 0.01;
    //    New_Fine_Grains.push_back(temporalV); //Add to new fine list 
    // }
    // Fine_Grains = New_Fine_Grains; //Re-update fine list

    return;  //Only exit after looping ALL the grains
  }

  //Change Aug 22, 2022
  void FineMelty(vector<Vector4d> & Fine_Grains, vector<Vector2d> & Fine_Velocities, double & loss_fines, vector<size_t> & NLMelt, vector<size_t> & NLBreak, vector<size_t> & NGMelt, vector<size_t> & NGBreak)
  //void FineMelty( vector<Vector4d> & Fine_Grains, vector<Vector2d> & Fine_Velocities  )
  {
    //Change Aug 22, 2022
    double rho_icem = 910; //kg/m3
    
    //Just make area and thickness equal to zero to stop considering and refilter
    double SumArea = 0.0;
    double SumThick = 0.0;
    for (size_t i = 0; i < Fine_Grains.size(); i++)
    {  
      //April 24, 2023
      if ( (Fine_Grains[i](1) <= 0.01 && Fine_Grains[i](1) != 0.0) || (Fine_Grains[i](1) < 0.0) )  //If thickness is zero  //0.000001
      {
         //Change Aug 22, 2022
         loss_fines += Fine_Grains[i](0) * max(Fine_Grains[i](1), 0.0) * 0.001 * rho_icem * 1000*1000*1000 ;  //Mass loss in kg using km3 units
          
         Fine_Grains[i](0) = 0.0; //Make its mass go down to zero (CHECK PERFORMANCE OF THIS)
         Fine_Grains[i](1) = 0.0;
         Fine_Velocities[i](0) = 0.0; //Velocity also goes down to avoid problems
         Fine_Velocities[i](1) = 0.0;
         NLMelt[NLMelt.size()-1] = NLMelt[NLMelt.size()-1] + 1;
      }
      SumArea += Fine_Grains[i](0);
      SumThick += Fine_Grains[i](1);
    }

    //April 24, 2023
    if (SumArea < 0.1 || SumThick < 0.01) // Add an artifical thickness and area to avoid problems if list is only zeros
    {
        Fine_Grains[0](0) = 1; //Only for the first fine
        Fine_Grains[0](1) = 0.001; //Only for the first fine
    }

    // //Alternative Method if Deleting becomes a problem (slower)
    // vector<Vector4d> New_Fine_Grains; //Save remaining fines here
    // Vector4d temporalV; //Temporal vector for saving fines
    // for (size_t i = 0; i < Fine_Grains.size(); i++)
    // {
    //   if (Fine_Grains[i](1) > 0.01)  //If thickness is zero  //0.000001 //Then the fine is not saved anymore
    //   {
    //      temporalV = Fine_Grains[i]; //We will preserve this fine
    //      New_Fine_Grains.push_back(temporalV); //Add to new fine list
    //   }
    // } 
    // if (New_Fine_Grains.size() < 1) // Add an artifical fine to avoid problems if list is empty
    // {
    //    temporalV[0] = 0.5;  temporalV[0] = 0.015;  temporalV[0] = 0.01;  temporalV[0] = 0.01;
    //    New_Fine_Grains.push_back(temporalV); //Add to new fine list 
    // }
    // Fine_Grains = New_Fine_Grains; //Re-update fine list
    return;
  }


  void TwoHighestForcesLocs(Vector2d & pF1,Vector2d & pF2,const size_t & g,const vector<PointInfo> & g1i){
        //takes grain index, g, and contact info, g1i, and finds locations of the 2 highest contact forces, pF1, pF2


        //ensure each contact is checked only once
        sort(_globalGrainState._grainCoord[g].begin(),_globalGrainState._grainCoord[g].end());
        unique(_globalGrainState._grainCoord[g].begin(),_globalGrainState._grainCoord[g].end());

        vector<CForce> cforces(_globalGrainState._grainCoord[g].size());
        //find points with 2 highest contact forces
        for (size_t c=0; c<_globalGrainState._grainCoord[g].size(); c++) { //for loop over contact grains

            size_t contact = _globalGrainState._grainCoord[g][c];
            double force=0;
            size_t contactpoint=0;

            size_t ID_Diff = 100000; //For ID differentiation grains from grainWalls, depends on size of grains !!!!!

            //grain contacts
            if (contact < INT_MAX-5 && contact < ID_Diff){// || contact == 2147483645){
                Grain2d other = _grains[FindGrainFromId(contact)];
                contactpoint = _grains[g].findInterparticleForceforFrac(other, force);
                //cout << "grain" << endl;
            }

            //Grain Wall Contacts
            else if (contact > ID_Diff && contact < INT_MAX-5){
                if (_grainsWall.size()>0)
                { 
                    Grain2d other = _grainsWall[FindGrainWallFromId(contact)];
                    contactpoint = _grains[g].findInterparticleForceforFrac(other, force); 
                }
            }


            //wall contacts
            else if ( contact > INT_MAX-5){//wall contacts or jaw
                size_t w = FindWallFromId(contact);
                contactpoint = _walls[w].FindWallForceForFrac(_grains[g],force);
                //cout << "wall" << endl;
            }
            else {
                contactpoint = 0;
            }
            cforces[c].force = force;
            cforces[c].cpoint = g1i[contactpoint].point;
            cforces[c].cid = contact;

        }//end find points
        sort(cforces.begin(),cforces.end());
        pF1 = cforces[0].cpoint;
        pF2 = cforces[1].cpoint;


    }

    void ClosestPoint(Vector2d & pF1,const vector<Vector2d> & gp0) {
        size_t np1 = gp0.size();
        //find closest point to center
        Vector2d newp = Vector2d(DBL_MAX,DBL_MAX);
        for (size_t p = 0; p<np1; p++ ){
            if (gp0[p].norm() < newp.norm()){
                newp = gp0[p];
            }
            pF1 = newp;
        }
    }

    //Find closest point to specific point
    void ClosestPoint_line(Vector2d & pF1, const vector<Vector2d> & gp0, Vector2d & LineP) {
        size_t np1 = gp0.size();
        //find closest point to line point
        double dist = DBL_MAX;
        for (size_t p = 0; p<np1; p++ ){
            if (  (gp0[p]-LineP).norm() < dist) {
                dist = (gp0[p]-LineP).norm();
                pF1 = gp0[p];
            }
        }
    }
    
    void changeGrainTfail(double & BREAK_PROB)
    {
        //cout << "Failure times!!! for grain number: " << _grains.size() << endl;
        //Initialize Random failure times (normal distribution, inherit failure time and origin time to new floes for simplicity)
        double mean_failure = double(BREAK_PROB);
        double temp_sample_fail;
        double sample_fail;
        double zero_fail = 0.000;
        std::default_random_engine time_failure;
        std::normal_distribution<double> distribution_fail(mean_failure, 0.1 * mean_failure); // mean, 0.1 std. dev mean
        
        for (size_t ik = 0; ik < _grains.size(); ik++)
        {
            temp_sample_fail = max( distribution_fail(time_failure) , 0.000 );
            sample_fail = temp_sample_fail;
            //std::cout << "Sample fail for grain: " << ik << " is: " << sample_fail << std::endl;
            _grains[ik].changeTfail(sample_fail);
            _grains[ik].changeTorigin(zero_fail);
            //cout << "New grain fail time: " << _grains[ik].getTfail() << endl;
        }
    }

    // Checks every grain if it is flagged for fracture, then cuts grain with a line connecting the 2 highest contact force points
    void FracRoutine(vector<double> & BreakBin, vector<double> & Diameters, double & limitMaxDiamV, double & WaveH, double & WaveL, int & PROB_B, vector<size_t> & NLMelt, vector<size_t> & NLBreak, vector<size_t> & NGMelt, vector<size_t> & NGBreak) {
        //cout << "GO BREAK!" << endl;
        vector<size_t> iter(_ngrains);
        for (size_t g = 0; g < _ngrains; g++){
            iter[g] = g;
        }
        cout <<"Random Shuffle" << endl;
        random_shuffle(iter.begin(),iter.end()); ////CHECK!!!
        cout << "ngrains control before loop: " << _ngrains << endl;
      for (size_t g1 = 0; g1 < _ngrains; g1++){
             size_t g = iter[g1]; //CRITICAL For 1 random breakage (no problem since 1 grain chosen per step) MODIF_2 STAB   //Sep 30, 2022 Change for Prob. Break Try off 
             //cout << "ngrains control in  loop: " << _ngrains << " grain in process: " << g1 <<  endl;
             //size_t g = g1; //For thickness control (allows looping all grains) MODIF_2 STAB                                   //Sep 30, 2022 Change for Prob. Break Try on

            //cout <<"Grain to Try out for break i = " << g << endl;  //BUGGGGG
            //check if grain should fracture
            //cout << "Get Volume" << endl;
            double Volume = _grains[g].getMass()/_grains[g].getDensity();
          //cout << "Get Stress" << endl; 
            Vector3d grainStress = _globalGrainState._grainStress[g]/Volume;  //EDIT
            //Vector3d grainStress;
            //grainStress << 0.0, 0.0, 0.0;

            if ( (fabs(grainStress(0)) > 0 || fabs(grainStress(1)) > 0 || fabs(grainStress(2)) > 0) && g >= 0 ){
                //cout << "Volume: " << Volume << " grainStress: " << grainStress(0) << " " << grainStress(1) << " " << grainStress(2) << endl;
            }    

            //check grain yield stress
            double Sig1 = (grainStress(0)+grainStress(1))/2. + sqrt(pow((grainStress(0)-grainStress(1))/2.,2)+pow(grainStress(2),2));
            double Sig2 = (grainStress(0)+grainStress(1))/2. - sqrt(pow((grainStress(0)-grainStress(1))/2.,2)+pow(grainStress(2),2));

            // Von Mises criterion
//          double VonMises = sqrt(grainStress(0)*grainStress(0)+grainStress(1)*grainStress(1)
//                  -grainStress(0)*grainStress(1)+3.*grainStress(2)*grainStress(2));
//          _grains[g].changeCritStress(VonMises);

            // Tresca criterion
//          double Tresca  = fabs(Sig1-Sig2);
//          _grains[g].changeCritStress(Tresca);

            // maximum principal stress criterion
            _grains[g].changeCritStress(Sig1);  //EDIT


            //Diam Control
            //cout << "Get Diameter " << endl;
            vector<Vector2d> VecOriginZ = _grains[g].getPointList();
            double AreaZ = 0.0;
            size_t nZ = VecOriginZ.size();

            for (size_t i = 0; i < nZ-1; i++)
            {
              AreaZ += ( VecOriginZ[i](0) * VecOriginZ[i+1](1) -  VecOriginZ[i](1) * VecOriginZ[i+1](0) ); 
            }
            AreaZ += (VecOriginZ[nZ-1](0) * VecOriginZ[0](1) -  VecOriginZ[nZ-1](1) * VecOriginZ[0](0) ); 

            double IAreaZ = 0.5*abs(AreaZ);
            double IDiamZ = (2 * sqrt(IAreaZ /3.141592653589793));

            
            if (g >= 0){
            //cout << "CritSIG: " << Sig1 << " grainYield: " << _grains[g].getYield() << " for grain #: " << g << endl;
            }

            if (g >= 0){

                //if ( (Sig1 > _grains[g].getYield() || fabs(Sig2) > 100.*_grains[g].getYield()) && (_grains[g].getMass() >= _MELT_MASS) && (_grains[g].getMass() > 0.0)  ){ 
                if ((Sig1 > _grains[g].getYield() || fabs(Sig2) > 100.*_grains[g].getYield()) && (IDiamZ >= 20.00 && IDiamZ > 0.0)) {
                    //cout << "We HAVE Breakage for Grain: " << g << " with position X: " << _grains[g].getPosition()(0) << " and Y:" << _grains[g].getPosition()(1) << endl;
                    //cout << "Testing Mass: " << _grains[g].getMass() << " Melt Mass: " << _MELT_MASS;
                    //cout << "Contact points #: " << _globalGrainState._grainCoord[g].size() << endl;
                    _grains[g].changeFracFlag(true);
                }
                else
                {
                    //cout << "NO Breakage for Grain: " << g << " with position X: " << _grains[g].getPosition()(0) << " and Y:" << _grains[g].getPosition()(1) << endl;
                }

            }

            //flag for open arbitrary break
            bool open_flag = false;

            srand (time(NULL));

            //size_t og = g; //For thickness criterion

            int B_rand_index = rand() % PROB_B + 1;  //Rand prob intended inverted + 1  //Sep 30, 2022 Change for Prob. Break Try Stay
            //int B_rand_index = rand() % 100;
            //int prob = _BREAK_PROB;  //30  //Sep 30, 2022 Change for Prob. Break Try 
            Vector2d bk_point;
            int prob = PROB_B;

            double limitMaxDiam = limitMaxDiamV;
            double limitMaxDiamT = limitMaxDiam*0.5;


            double adjustIbr = 0.1;
            double Ibr = WaveH * _grains[g].getThickness() * 5.5e9 * 0.5 * pow(0.27e6,-1.0) * pow(WaveL, -2.0);
            Ibr *= adjustIbr;
            //cout << "Break Criterion from Dumont using meters and seconds as units (> 0.014 to break): " << Ibr << endl;
            
            //Obtain times for failure
            const double Tfail = _grains[g].getTfail();
            const double originTfail = _grains[g].getTorigin();
            //cout << "CONFIRM Tfail: " << _grains[g].getTfail() << endl;
            //cout << "CONFIRM Tfailorigin: " << _grains[g].getTorigin() << endl;
            const double curr_step = _stepup;
            //cout << "CONFIRM stepup: " << curr_step << endl;
            
            
            //FLOE BREAKAGE CRITERIA
            //Ibr = 0.00001; //ON BREAK ONLY USING THICKNESS REDUCTION DONT USE DUMONT OFF:Use Dumont as well
            //Use of random break, random grain, random angle, WAVE RELATED CRITERION
            //if ( B_rand_index < prob && (_grains[g].getMass() >= _MELT_MASS) && (_grains[g].getMass() > 0.0) ) //Then break
            if ( (IDiamZ >= limitMaxDiam) && (_grains[g].getMass() > 0.0) ) //Then break  //Using Ibr criterion //METHOD 2a Skip Ibr //MODIF_2 PREC   Sep 30, 2022 Change Off PROB CHANGE OFF
            //if ( Ibr > 0.014 && (IDiamZ >= limitMaxDiam) && (_grains[g].getMass() > 0.0) ) //Then break  //Using Ibr criterion //METHOD 2 MODIF_2 PREC
            //if ( B_rand_index >= prob - 2 && (IDiamZ >= limitMaxDiam) && (_grains[g].getMass() > 0.0) ) //Then break  //Using random occurrence  //METHOD 1  Sep 30, 2022 Change On   //Probability related criteria
            //if ( curr_step - originTfail >= Tfail )   //Oct 18, 2022 Time Failure Criterion Oct 18, 2022 PROB CHANGE ON
            {
              if (g == 0)
              {
                cout << "Tfail: " << Tfail << " at step: " << _stepup << endl;
                cout << "originTfail: " << originTfail << endl;
              }
              //cout << "Ibr BREAK" << endl;
              //size_t ng = rand() % (_ngrains);
              //cout << "To WAVE break: " << ng << " out of ngrains total: " << _ngrains << endl;
              //g = ng; ///CHECK!!! //ALREADY SHUFFLE NO NEED

             // //SIZE CONDITION
             // if (IDiamZ >= 20.00)
             // {
              
              //NEW OFF
              _grains[g].changeFracFlag(true);
              open_flag = true;
              //NEW OFF

             // }

             
              //To avoid breaking at center chose any point from 0 to Perc_Dist% of BBRox Radius randomly
              int perc_dist = 85;  //100
              //// 50% (rand() % 50) + 50
              bk_point(0) = (float(rand() % perc_dist))*(1.0/100.0)*_grains[g].getRadius();
              bk_point(1) = (float(rand() % perc_dist))*(1.0/100.0)*_grains[g].getRadius();

              //Fix Dist, also angle  (TEMPORAL)
              //bk_point(0) = 25*(1.0/100.0)*_grains[og].getRadius();
              //bk_point(1) = 25*(1.0/100.0)*_grains[og].getRadius();


              // Ice Breakage Criterion and Waves
              // Ibr = H_s  * h * Y * 0.5 * sigma_c^-1 * lambda ^ -2 > 0.014   Hs = sign. wave Height m, h = thickness m, lambda = peak wavelength, sigma_c = 0.27 MPa, Y = 5.5 GPa, lambda = cT or c/f, use Period and Stokes Wave Velocity to get lambda
              // Source Li et al. 2021, using Dumont et al. 2011

            }  

            //Thickness condition (if wave does not do it)
            double crit_thick = 0.005;  //0.1 Prior cases 1-78  Sep 30, 2022 Change on
            //double crit_thick =  0.1; //Try finer 1 cm average for all //0.2 m or 2e-4 km    Sep 30, 2022 Change off
            
            if (_grains[g].getThickness() < crit_thick && IDiamZ >= limitMaxDiam && open_flag == false && _grains[g].getMass() > 0.0) //Be more lenient for breakage? WONT Happen for only B
            //if (_grains[og].getThickness() < crit_thick && IDiamZ >= 20.00 && ng != og )
            {
              cout << "THICK BREAK for thickness: " << _grains[g].getThickness() << endl;
              //cout << "To THICK break: " << og << " out of ngrains total: " << _ngrains << endl;
              _grains[g].changeFracFlag(true);
              open_flag = true;

              //To avoid breaking at center chose any point from 0 to Perc_Dist% of BBRox Radius randomly
              int perc_dist = 85;  //100
              //// 50% (rand() % 50) + 50
              bk_point(0) = (float(rand() % perc_dist))*(1.0/100.0)*_grains[g].getRadius();
              bk_point(1) = (float(rand() % perc_dist))*(1.0/100.0)*_grains[g].getRadius();

              //Fix Dist, also angle  (TEMPORAL)
              //bk_point(0) = 25*(1.0/100.0)*_grains[og].getRadius();
              //bk_point(1) = 25*(1.0/100.0)*_grains[og].getRadius();

              //g = og; ///CHECK!!!
            }

            double ang_line = rand() % 1257; //From -2pi to 2pi              
            ang_line -= 628;
            ang_line *= 0.01;
            
            //ang_line = 0.5*3.1415926535897932384; //Fix to pi/4 for the mean time also distance (TEMPORAL)

            // cout << "Breakage Real Start for grain: " << g << endl;
            //if ( (_grains[g].getFracFlag() && _globalGrainState._grainCoord[g].size() >= 2) ||  (_grains[g].getFracFlag() && open_flag)      ) {
            if (    (_grains[g].getFracFlag() && open_flag)     ) { //Only thickness prob., no stress for now              
            //cout << "Breakage Real Start for grain IN: " << g << endl;
//              cout << "coord # " << _globalGrainState._grainCoord[g].size() << endl;
                /*
                 * Get Grain 1 Info
                 */
                //cout << "Breakage Real Start for grain define: " << g << endl;
                Vector2d pF1;
                Vector2d pF2;
                Vector2d g1cm = _grains[g].getCmLset();
                vector<Vector2d> g1p0 = _grains[g].getRefPointList();
                size_t np1 = g1p0.size();
                size_t gid = _grains[g].getId();
                Levelset2d g1lset = _grains[g].getLset();
                vector<double> g1Lset = g1lset.getLevelset();
                size_t Xdim = g1lset.getXdim();
                size_t Ydim = g1lset.getYdim();

                vector<PointInfo> g1i(np1);

                const double cos1 = cos(_grains[g].getTheta());
                const double sin1 = sin(_grains[g].getTheta());
                
                //cout << "Breakage Real Start for grain LOOP Point contact " << g << endl;
                for (size_t i = 0; i<np1;i++){
                    //must unrotate and bring to origin
//                  g1i[i].point << (g1p[i](0)-_grains[g].getPosition()(0))*cos1 + (g1p[i](1)-_grains[g].getPosition()(1))*sin1,
//                                            -(g1p[i](0)-_grains[g].getPosition()(0))*sin1 + (g1p[i](1)-_grains[g].getPosition()(1))*cos1;
                    g1i[i].point = g1p0[i];
                    //
                    g1i[i].shear = _grains[g].getNodeShears()[i];
                    g1i[i].contact = _grains[g].getNodeContact()[i];
                }

                //cout << "Testing Mass: " << _grains[g].getMass() << " Melt Mass: " << _MELT_MASS << " Diameter Size: " << IDiamZ;

                /*
                 * Find Splitter Shape
                 */
                if (open_flag == false)
                {
                    //cout << "TwoHighestForcesLocs" << endl;
                    TwoHighestForcesLocs(pF1,pF2,g,g1i);
                }

                
                //Splitter for two arbitrary points not based on forces
                if (open_flag == true)
                {
                   //cout << "Arbitrary Cut" << endl;
                   //replacing pF1, pF2
                   //Point_Slicer(pF1a, pF2a)
                   Vector2d pF1a;
                   Vector2d pF2a;
                   Vector2d LPA; 
                   Vector2d LPB;
                   LPA(0) = bk_point(0) +  _grains[g].getRadius()*cos(ang_line);
                   LPA(1) = bk_point(1) + _grains[g].getRadius()*sin(ang_line);
                   LPB(0) = bk_point(0) -_grains[g].getRadius()*cos(ang_line);
                   LPB(1) = bk_point(1) -_grains[g].getRadius()*sin(ang_line);
                   //cout << "ClosestPointLine" << endl;
                   ClosestPoint_line(pF1a, g1p0, LPA);
                   ClosestPoint_line(pF2a, g1p0, LPB);
                   pF1 = pF1a;
                   pF2 = pF2a;
                   cout <<  "bk_point: " << bk_point(0) << " " << bk_point(1) << endl;
                   cout <<  "pF1: " << pF1(0) << " " << pF1(1) << endl;
                   cout <<  "pF2: " << pF2(0) << " " <<  pF2(1) << endl;
                }  

                //cout << "Get pF2" << endl;
                double bigV = M_PI*_grains[g].getRadius()*_grains[g].getRadius();
                double g1V  = _grains[g].getMass()/_grains[g].getDensity();
                if (bigV/g1V > _fracProps.getRoMratio() || pF2 == pF1){
                    ClosestPoint(pF1,g1p0);
                    pF2 = Vector2d(0.,0.);
                }

                //cout << "pF1: " << pF1 << endl;
                //cout << "pF2: " << pF2 << endl;
                /*
                 * Create Splitter
                 */

                //cout << "Make Splitter" << endl;
                //get points that define split line in the level set reference frame
                Vector2d splitpt1 = pF1 + g1cm;
                Vector2d splitpt2 = pF2 + g1cm;

                //cout << "Making Splitter and Surface Lines" << endl;
                //make splitter surface points
                vector<Vector2d> splpts0 = _fracProps.MakeSplitter(_grains[g],splitpt1,splitpt2); //A. Make splitter line
                //cout << "Splts Pts: " << splpts0[0] <<  endl;
                //make splitter level set
                vector<double> splitset = g1lset.makeFracSurfaceLine(splitpt1,splitpt2); //B. Make LS split
                //cout << "Splts Set: " << splitset[0] << endl;
                //output splitset if desired
                // _fracProps.SaveSplitter(splpts0,splitset,Xdim,Ydim);

                /*
                 * Make New Grains (using Lset and Surface Points)
                 */


                //make new Lsets  //C. Filter new LS
                //cout << "Make New LSets" << endl;
                vector<double> g2Lset(splitset.size());
                vector<double> g3Lset(splitset.size());
                for(size_t i=0; i<splitset.size();i++){
                    g2Lset[i] = max(splitset[i],g1Lset[i]);
                    g3Lset[i] = max(splitset[i]*-1.,g1Lset[i]);
                }
                //cout << "Splts Set2: " << g2Lset[0] << endl;
                //cout << "Splts Set3: " << g3Lset[0] << endl;

                //reinit due to imperfect set ops
                //g2Lset = reinit(g2Lset); //need to build in Levelset2d.h
                //g3Lset = reinit(g3Lset);

                Levelset2d g2lset(g2Lset, Xdim, Ydim);
                Levelset2d g3lset(g3Lset, Xdim, Ydim);

                //2)pts
                vector<PointInfo> g2i, g3i;

                //cout << "Work with Splitter Points" << endl;
                //decides which splitter points stay
                //Input: splpts0, g2lset, g1cm;  Output: g2i, g3i
                _fracProps.KeepSplitterPoints(splpts0,g1lset,g1cm,g2i,g3i); //g2i and g3i start the same //D. Confirm splitter points penetrating LS, use PointInfo format
                //decides which previous grain points stay
                //Input: splitpt1, splitpt2, g1i, g1cm;  Output: g2i, g3i
                _fracProps.KeepGrainPoints(splitpt1,splitpt2,g1i,g1cm,g2i,g3i);  //g2i and g3i are differentiated  E. Add splitter line and points for each respective half of split

                //make sure potential grains has a reasonable amount of surface points
                if (g2i.size()<5 || g3i.size()<5 ){ 
                //Adjust for size
                //if (g2i.size()<3 || g3i.size()<3 ){   
                    //cout << "Break WARNING: TOO SMALL NaN" << endl;
                    //return; //to ensure no NaNs in the crumbles
                    continue;
                }


                //save number of points for grains
                size_t np2 = g2i.size();
                size_t np3 = g3i.size();

                /*
                 * Find Physical Properties
                 */
                //cout << "Find New Properties" << endl;
                double mass2;
                Vector2d g2cm;
                double I2;
                _fracProps.findGrainProps0(g2lset, mass2, g2cm, I2);

                double mass3;
                Vector2d g3cm;
                double I3;
                _fracProps.findGrainProps0(g3lset, mass3, g3cm, I3);
                /*
                 * Position grains correctly, removing remnants of the grain 1 reference frame
                 *
                 */

                Vector2d position2 = _grains[g].getPosition()-Vector2d((g1cm(0)-g2cm(0))*cos1 - (g1cm(1)-g2cm(1))*sin1,
                                                              (g1cm(0)-g2cm(0))*sin1 + (g1cm(1)-g2cm(1))*cos1);
                double maxR2 = 0;
                vector<Vector2d> g2p0, g3p0;

                //shift points to correct position
                for (size_t i = 0; i<np2;i++){

                    g2i[i].point += g1cm-g2cm;
                    if (g2i[i].point.norm() > maxR2){
                        maxR2 = g2i[i].point.norm();
                    }
                }

                // reorder points to counterclockwise
                sort(g2i.begin(), g2i.end(),pointsort);

                for (size_t i=0;i<np2;i++){
                    g2p0.push_back(g2i[i].point);
                }


                Vector2d position3 = _grains[g].getPosition()-Vector2d((g1cm(0)-g3cm(0))*cos1 - (g1cm(1)-g3cm(1))*sin1,
                                                                  (g1cm(0)-g3cm(0))*sin1 + (g1cm(1)-g3cm(1))*cos1);
                double maxR3 = 0;
                //shift points to correct position
                for (size_t i = 0; i<np3;i++){
                    g3i[i].point += g1cm-g3cm;

                    //Find bbox radius
                    if (g3i[i].point.norm() > maxR3){
                        maxR3 = g3i[i].point.norm();
                    }
                }
                // sort indexes based on comparing values in g3p0 b
                sort(g3i.begin(), g3i.end(),pointsort);

                for (size_t i=0;i<np3;i++){
                    g3p0.push_back(g3i[i].point);
                }


                //check thinness of potential grains in case they could slide into other grains  !!!!TODO: ADD TO MELT
                Vector2d g2min,g3min; //closest surface point to centroid
                double g2minR,g3minR; //minimum distance from centoid to surface
                ClosestPoint(g2min,g2p0);
                ClosestPoint(g3min,g3p0);
                g2minR = g2min.norm();
                g3minR = g3min.norm();
                if (g2minR < .5|| g3minR < .5) {
                    //cout << "Break Radius NaN" << endl;
                    //return;
                    continue;
                }


//              if (M_PI*maxR2*maxR2/(mass2) > _fracProps.getMaxRoM() || M_PI*maxR3*maxR3/(mass3) > _fracProps.getMaxRoM() ){
//                  return; //to ensure no NaNs in the crumbles
//              }

                /*
                 * Create grains
                 */
                //Data Which Will NOT depend on fracture
                vector<double> g2tempV =  _grains[g].getgrainTemp2D().getLevelset();
                vector<double> g3tempV =  _grains[g].getgrainTemp2D().getLevelset(); //g2lset g3lset

                // vector<double> g2tempV =  g2lset.getLevelset();
                // vector<double> g3tempV =  g3lset.getLevelset();

                vector<double> g2ThickV =  _grains[g].getMthick().getLevelset();
                vector<double> g3ThickV =  _grains[g].getMthick().getLevelset(); //g2lset g3lset

                //Pre-correct if needed for thickness
                vector<double> g2Geo =  g2lset.getLevelset();
                vector<double> g3Geo =  g3lset.getLevelset(); //g2lset g3lset

                for (size_t ii = 0; ii<g2lset.getXdim(); ii++) 
                {            
                  for (size_t jj = 0; jj<g2lset.getYdim(); jj++)
                  {
                    if ( g2Geo[(jj*g2lset.getXdim())+ii]>0.0 && g2ThickV[(jj*g2lset.getXdim())+ii] > 0.0 ) //You need to redefine thickness in case LS is changed by breakage or other reasons, to avoid bugs
                    {
                          //cout << "Correcting Thickness after break Origin!!!" << endl;
                          g2ThickV[(jj*g2lset.getXdim())+ii] = -_grains[g].getThickness(); 
                    }
                  }
                } 
                for (size_t ii = 0; ii<g3lset.getXdim(); ii++) 
                {            
                  for (size_t jj = 0; jj<g3lset.getYdim(); jj++)
                  {
                    if ( g3Geo[(jj*g3lset.getXdim())+ii]>0.0 && g3ThickV[(jj*g3lset.getXdim())+ii] > 0.0 ) //You need to redefine thickness in case LS is changed by breakage or other reasons, to avoid bugs
                    {
                          //cout << "Correcting Thickness after break Origin!!!" << endl;
                          g3ThickV[(jj*g3lset.getXdim())+ii] = -_grains[g].getThickness(); 
                    }
                  }
                }       



                Levelset2d g2temp(g2tempV, g2lset.getXdim(), g2lset.getYdim());
                Levelset2d g3temp(g3tempV, g3lset.getXdim(), g3lset.getYdim());

                Levelset2d g2ThickLS(g2ThickV, g2lset.getXdim(), g2lset.getYdim());
                Levelset2d g3ThickLS(g3ThickV, g3lset.getXdim(), g3lset.getYdim());                

                //cout << "Creating grain2 with npoints: " << g2p0.size() << endl;
                //cout << "Grain 2 Points:" << endl;
                for (size_t i = 0; i < g2p0.size() ; i++ )
                {
                    //cout << g2p0[i](0) << " " << g2p0[i](1) << endl;
                }

                //cout << "Making New Grains Split 2" << endl;
                bool tempflag = false;
                int temploc = 0;
                const double current_step = _stepup;
                //Grain2d g2 = Grain2d(mass2, position2, _grains[g].getVelocity(), I2, _grains[g].getTheta(), _grains[g].getOmega(),
                //                        g2cm, g2p0, maxR2, g2lset, _grains[g].getKn(), _grains[g].getKs(), _grains[g].getMu(), _maxId+1);
                Grain2d g2 = Grain2d(mass2, position2, _grains[g].getVelocity(), 
                 mass2, _grains[g].getTemperature(), _grains[g].getThickness(), g2ThickLS, _grains[g].getTemperature(), _grains[g].getThickness0(), _grains[g].getUtemper(), _grains[g].getUtemper(), g2temp, g2temp, 
                 I2, _grains[g].getTheta(), _grains[g].getOmega(), g2cm, g2p0, g2p0.size(), maxR2, 
                 g2lset, g2lset, 0, _grains[g].getId()+1, _grains[g].getKn(), _grains[g].getKs(), _grains[g].getMu(), tempflag, temploc,
               _grains[g].getgrainStress2D(), _grains[g].getgrainDamage() , _grains[g].getgrainThickness(), g2p0, _grains[g].getTfail(), current_step); //_grains[g].getId()+1

                //cout << "For grain2 Tfail: " << g2.getTfail() << " at step: " << _stepup << endl;
                //cout << "For grain2 originTfail: " << g2.getTorigin() << endl;

                // double new_mass;
                // double new_I;
                // Vector2d new_cm;
                // Levelset2d ref_lset;

                // ref_lset = g2.getLset();
                // _fracProps.findGrainProps0(ref_lset, new_mass, new_cm, new_I);
                // g2.changeInertia(new_I);
                // g2.changeMass(new_mass*(10e6)); 
                // g2.changeCM(new_cm);

                g2.changeDensity(_grains[g].getDensity());
                for (size_t i=0;i<np2;i++){
                    g2.getNodeShearsNonConst()[i]  = g2i[i].shear;
                    g2.getNodeContactNonConst()[i] = g2i[i].contact;
                }

                //cout << "Save New Grains Split 2: " << endl;
                //_fracProps.SaveGrain(g2,g2p0);


                //Grain2d g3 = Grain2d(mass3, position3, _grains[g].getVelocity(), I3, _grains[g].getTheta(), _grains[g].getOmega(),
                //                        g3cm, g3p0, maxR3, g3lset, _grains[g].getKn(), _grains[g].getKs(), _grains[g].getMu(), _maxId+2);
                

                //cout << "Creating grain3 with npoints: " << g3p0.size() << endl;
                //cout << "Grain 3 Points:" << endl;
                for (size_t i = 0; i < g3p0.size() ; i++ )
                {
                    //cout << g3p0[i](0) << " " << g3p0[i](1) << endl;
                }

                //cout << "Making New Grains Split 3" << endl;
                Grain2d g3 = Grain2d(mass3, position3, _grains[g].getVelocity(), 
                             mass3, _grains[g].getTemperature(), _grains[g].getThickness(), g3ThickLS, _grains[g].getTemperature(), _grains[g].getThickness0(), _grains[g].getUtemper(), _grains[g].getUtemper(), g3temp, g3temp, 
                             I3, _grains[g].getTheta(), _grains[g].getOmega(), g3cm, g3p0, g3p0.size(), maxR3, 
                             g3lset, g3lset, 1, _grains[g].getId()+2, _grains[g].getKn(), _grains[g].getKs(), _grains[g].getMu(), tempflag, temploc,
                           _grains[g].getgrainStress2D(), _grains[g].getgrainDamage() , _grains[g].getgrainThickness(), g3p0, _grains[g].getTfail(), current_step);  //_grains[g].getId()+2
                g3.changeDensity(_grains[g].getDensity());
                
                //cout << "For grain3 Tfail: " << g3.getTfail() << " at step: " << _stepup << endl;
                //cout << "For grain3 originTfail: " << g3.getTorigin() << endl;


                // ref_lset = g3.getLset();
                // _fracProps.findGrainProps0(ref_lset, new_mass, new_cm, new_I);
                // g3.changeInertia(new_I);
                // g3.changeMass(new_mass*(10e6));
                // g3.changeCM(new_cm);

//              for (size_t i=0;i<np3;i++){
//                  g3.getNodeShearsNonConst()[i]  = g3i[i].shear;
//                  g3.getNodeContactNonConst()[i] = g3i[i].contact;
//              }

                //send shears to "master" grains
                //cout << "send shears to master grains" << endl;
//              double dist = 50;
////                cout << "g3isize " << g3i.size() << " " << np3 << endl;
//              vector<Vector2d> g3_moved_pts = g3.getPointList();
                for (size_t i=0;i<np3;i++){
//                  if (g3i[i].contact > _maxId){ //contact with wall; shear stays
                        g3.getNodeShearsNonConst()[i] = g3i[i].shear;
                        g3.getNodeContactNonConst()[i] = g3i[i].contact;

////                        cout << "\nwall " << g3i[i].contact << "\n"<<endl;
//                  }
//                  else if (g3i[i].contact>0){ //contact with grain; shear goes to other grain
//                      //loop over each possible point and search for closest one
//                      size_t pt_id;
//
//                      vector<Vector2d> target_grain_pts = _grains[FindGrainFromId(g3i[i].contact)].getPointList();
//                      for (size_t j = 0;j<target_grain_pts.size();j++){
//                          if ((g3_moved_pts[i]-target_grain_pts[j]).norm() < dist){
//                              dist = (g3_moved_pts[i]-target_grain_pts[j]).norm();
//                              pt_id = j;
//                          }
//                      }
////                        cout << "\ndistance " << dist << "\n"<<endl;
//                      _grains[FindGrainFromId(g3i[i].contact)].getNodeContactNonConst()[pt_id] = g3.getId();
//                      _grains[FindGrainFromId(g3i[i].contact)].getNodeShearsNonConst()[pt_id] = g3i[i].shear;
//                  }
                }

                //cout << "Save New Grains Split 3: " << endl;
                //_fracProps.SaveGrain(g3,g3p0);

                double yield = _grains[g].getYield();
                g3.changeYield(yield*pow(_grains[g].getDiameter()/g3.getDiameter(),3./3.));
                g2.changeYield(yield*pow(_grains[g].getDiameter()/g2.getDiameter(),3./3.));


                //find shears where g was slave and move to master
                for (size_t i = 0; i<g;i++){
                    //check if grain was a master grain on g
                    vector<size_t> master_contacts = _grains[i].getNodeContact();
                    bool master = std::binary_search(master_contacts.begin(),master_contacts.end(),gid); //true if master grain was contacting fractured grain
                    if (master) {
//                      cout << "Position " << g2.getPosition().transpose() << " " << g3.getPosition().transpose() << " " << _grains[i].getPosition().transpose() << endl;
                        //move shears to g2,g3
                        //find node
                        vector<size_t> nodes;
                        size_t n=0;
                        for (size_t j = 0; j<master_contacts.size();j++){
                            if (master_contacts[j] == gid){
                                nodes.push_back(j);
                                n++;
                            }
                        }

                        //get point location
                        vector<Vector2d> master_pts = _grains[i].getPointList();

                        vector<Vector2d> g2_pts = g2.getPointList();
                        vector<Vector2d> g3_pts = g3.getPointList();

                        for(size_t j=0;j<n;j++){
                            //find closest pt on either g2 or g3
                            Vector2d master_pt = master_pts[nodes[j]];
                            size_t kmin = 0;
                            double dist = 50.;
                            for (size_t k = 0; k<np2;k++){
                                double dist_k = (master_pt-g2_pts[k]).norm();
//                              cout << "dist_k " << dist_k << endl;
                                if (dist_k < dist){
                                    dist = dist_k;
                                    kmin = k;
                                }
                            }
                            for (size_t k = np2; k<np2+np3;k++){
                                double dist_k = (master_pt-g3_pts[k-np2]).norm();
                                if (dist_k < dist){
                                    dist = dist_k;
                                    kmin = k;
                                }
                            }
//                          cout << "dist " << dist << endl;
                            if (dist < 1) {
                                //send shears to fragments
                                if (kmin>=np2) {
                                    kmin -= np2; //on grain 3
                                    g3.getNodeShearsNonConst()[kmin]  = _grains[i].getNodeShears()[nodes[j]];
                                    g3.getNodeContactNonConst()[kmin] = _grains[i].getNodeContact()[nodes[j]];
                                } else { //on grain 2
                                    g2.getNodeShearsNonConst()[kmin]  = _grains[i].getNodeShears()[nodes[j]];
                                    g2.getNodeContactNonConst()[kmin] = _grains[i].getNodeContact()[nodes[j]];
                                }
                                _grains[i].getNodeShearsNonConst()[nodes[j]] = 0;
                                _grains[i].getNodeContactNonConst()[nodes[j]] = 0;
                            }
                        }
                    }
                }


                //Bin Mass changes due to Breakage
                //Initial Area and Diameter before Breaking
                //double IArea = PointsArea(_grains[g].getPointList());
                //double IDiam = MeanCaliperD(IArea);

                //cout << "Find Area 2" << endl;
                vector<Vector2d> VecOrigin = _grains[g].getPointList();
                double Area = 0.0;
                size_t n = VecOrigin.size();

                for (size_t i = 0; i < n-1; i++)
                {
                  Area += ( VecOrigin[i](0) * VecOrigin[i+1](1) -  VecOrigin[i](1) * VecOrigin[i+1](0) ); 
                }
                Area += (VecOrigin[n-1](0) * VecOrigin[0](1) -  VecOrigin[n-1](1) * VecOrigin[0](0) ); 

                double IArea = 0.5*abs(Area);
                double IDiam = (2 * sqrt(IArea /3.141592653589793));
                double F1Area, F2Area;

                //Avoid degenerate breakage
                if (isnan(IArea))
                {
                    //cout << "Breakage Area NaN" << endl;
                    //return;
                    continue;
                }

                //Final Area and Diameter after breaking in 2
                //double F1Area = PointsArea(g2.getPointList());
                //double F1Diam = MeanCaliperD(F1Area);
                //double F2Area = PointsArea(g3.getPointList());
                //double F2Diam = MeanCaliperD(F2Area);

                VecOrigin = g2.getPointList();
                Area = 0.0;
                n = VecOrigin.size();

                for (size_t i = 0; i < n-1; i++)
                {
                  Area += ( VecOrigin[i](0) * VecOrigin[i+1](1) -  VecOrigin[i](1) * VecOrigin[i+1](0) ); 
                }
                Area += (VecOrigin[n-1](0) * VecOrigin[0](1) -  VecOrigin[n-1](1) * VecOrigin[0](0) ); 

                if (isnan(Area))
                {
                    F1Area = IArea * (g2.getMass()/_grains[g].getMass());
                }
                else
                {
                    F1Area = 0.5*abs(Area);
                }

                
                //cout << "Find Area 3" << endl;
                VecOrigin = g3.getPointList();
                Area = 0.0;
                n = VecOrigin.size();

                for (size_t i = 0; i < n-1; i++)
                {
                  Area += ( VecOrigin[i](0) * VecOrigin[i+1](1) -  VecOrigin[i](1) * VecOrigin[i+1](0) ); 
                }
                Area += (VecOrigin[n-1](0) * VecOrigin[0](1) -  VecOrigin[n-1](1) * VecOrigin[0](0) ); 

                if (isnan(Area))
                {
                    F2Area = IArea * (g3.getMass()/_grains[g].getMass());
                }
                else
                {
                    F2Area = 0.5*abs(Area);
                }
                
                
                if ((F2Area + F1Area) > (1.2*IArea) ||  (F2Area + F1Area) < (0.8*IArea))   //20% error respected otherwise force using mass proportions
                {
                    F1Area = IArea * (g2.getMass()/_grains[g].getMass());
                    F2Area = IArea * (g3.getMass()/_grains[g].getMass());
                }
                  

                double F1Diam = (2 * sqrt(F1Area /3.141592653589793));
                double F2Diam = (2 * sqrt(F2Area /3.141592653589793));

                //Find bin that lost this area (assume minor break, one staying on same bin other changes)
                size_t loc_index, loc_index2, loc_index3;
                size_t locC_index, locC_index2, locC_index3; //April 24, 2023 for FSD
                for (size_t bi = 0; bi < Diameters.size()-1; bi++){  
                    if ( IDiam < Diameters[bi] &&  IDiam >= Diameters[bi+1] ) 
                    {
                       loc_index = bi;
                       locC_index = bi;
                       break;
                    }
                    else if ( IDiam < Diameters[Diameters.size()-1] )  //Min size threshold control
                    {
                       loc_index = Diameters.size()-2;
                       locC_index = Diameters.size()-2;
                    }
                    else if ( IDiam >= Diameters[0] )  //Max size threshold control
                    {
                       loc_index = 0;
                       locC_index = 1000;
                       break;
                    }
                }
                for (size_t bi = 0; bi < Diameters.size()-1; bi++){  
                    if ( F1Diam < Diameters[bi] &&  F1Diam >= Diameters[bi+1] ) 
                    {
                       loc_index2 = bi;
                       locC_index2 = bi;
                       break;
                    }
                    else if ( F1Diam < Diameters[Diameters.size()-1] )  //Min size threshold control
                    {
                       loc_index2 = Diameters.size()-2;
                       locC_index2 = Diameters.size()-2;
                    }
                    else if ( F1Diam >= Diameters[0] )  //Max size threshold control
                    {
                       loc_index2 = 0;
                       locC_index2 = 1000;
                       break;
                    }
                }
                for (size_t bi = 0; bi < Diameters.size()-1; bi++){  
                    if ( F2Diam < Diameters[bi] &&  F2Diam >= Diameters[bi+1] ) 
                    {
                       loc_index3 = bi;
                       locC_index3 = bi;
                       break;
                    }
                    else if ( F2Diam < Diameters[Diameters.size()-1] )  //Min size threshold control
                    {
                       loc_index3 = Diameters.size()-2;
                       locC_index3 = Diameters.size()-2;
                    }
                    else if ( F2Diam >= Diameters[0] )  //Max size threshold control
                    {
                       loc_index3 = 0;
                       locC_index3 = 1000;
                       break;
                    }
                }
                
                //April 24, 2023 FSD
                if (locC_index2 != 1000){
                    NGBreak[locC_index2] = NGBreak[locC_index2] + 1;
                }
                if (locC_index3 != 1000){
                    NGBreak[locC_index3] = NGBreak[locC_index3] + 1;
                }
                if (locC_index != 1000){
                    NLBreak[locC_index] = NLBreak[locC_index] + 1;
                }

                //cout << "BKG: AreaOriginal: " << IArea << " Area2: " << F1Area << " Area3: " << F2Area << " loc: " << loc_index << " loc2: " << loc_index2 << " loc3: " << loc_index3 << endl;
                //cout << "BKG: DiamOriginal: " << IDiam << " Diam2: " << F1Diam << " Diam3: " << F2Diam << " loc: " << loc_index << " loc2: " << loc_index2 << " loc3: " << loc_index3 << endl;

                //Avoid degenerate breakage, skip if it happens
                if (isnan(IDiam))
                {
                    //cout << "WARNING: BreakBinsSkipped due to NaN" << endl;
                    //return;
                    continue;
                }


                //Check if grain changed bin as well, if it did then substract all mass from initial bin and move it to final bin
                if (loc_index == loc_index2 && loc_index == loc_index3){ //No change
                } 
                else if (loc_index == loc_index2 && loc_index != loc_index3)  //Some Loss
                {
                    BreakBin[loc_index] += -(F2Area);
                    BreakBin[loc_index3] += (F2Area);
                }
                else if (loc_index == loc_index3 && loc_index != loc_index2) //Some Loss
                {
                    BreakBin[loc_index] += -(F1Area);
                    BreakBin[loc_index2] += (F1Area);
                }
                else{ //All lost
                    BreakBin[loc_index] += -(IArea);  //or (F1Area + F2Area), but what if bkg is not mass conservative, thus lose IArea just in case
                    BreakBin[loc_index2] += (F1Area);
                    BreakBin[loc_index3] += (F2Area);
                }

                /*
                 * Add/Remove grains from world
                 */
                //remove fractured grain
                //cout << "Number of grains BEFORE: " << _ngrains << endl;
                _grains.erase(_grains.begin()+g);
                _ngrains--;
                //Add new grains if not problematic  _2 is for a smaller resolution of grains
                if (_fracProps.ApplyFilter_2(g2)){
                    //cout << "Insert 2" << endl;
                    _grains.insert(_grains.begin(),g2);
                    _ngrains++;
                }
                if (_fracProps.ApplyFilter_2(g3)) {
                    //cout << "Insert 3" << endl;
                    _grains.insert(_grains.begin(), g3);
                    _ngrains++;
                }
                //cout << "Number of grains NOW: " << _ngrains << endl;

//              _grains.insert(_grains.begin(),g2);
//              _grains.insert(_grains.begin(),g3);
//
//              _ngrains++;
//              _ngrains++;

                _maxId++;_maxId++;
                //cout << "Reset Global State" << endl;
                _globalGrainState.reset();
                _globalGrainState.resize(_ngrains, _nwalls);
                computeWorldState();
                
                //NEW CHANGE OFF
                //cout << "Breakage Real End for grain: " << g << endl;
                return; //Sep 30, 2022 Change for Prob. Break Try OFF //TURN OFF FOR ALL GRAINS TO BE ANALYZED //ON FOR ONLY ONE GRAIN PER STEP
                //NEW CHANGE OFF
                
                //TWO OPTIONS: THICKNESS ONLY AND ALL GRAINS BREAK BUT FIX BUG ORRRR ONLY ONE BREAK AND TWO CRITERIA ?? A or B??
                //cout << "Breakage Real End for grain: " << g << endl;
            }//if frac flag
        }//grains loop
        //return; //CHECK!!!  //Sep 30, 2022 Change for Prob. Break Try  ON
    }
    
    //##################################################################################################################################################################################################################
    //##################################################################################################################################################################################################################
    
    void exportGrainDataBreak(size_t & grain_index, string & outDirT){
    
        size_t nTout = 4233600/49; //Just for 2018!!!!!!!
        size_t step_time = _stepup/nTout;
        //After identifying grain to break, export data for external FEM
        //Print Temperature, Geometric and Thickness Level Set for a Sample Grain
        
        //Output grain properties for study
        //string outDir4 = "./Output/SeaIce2018_Fluid_1077/Vertical_Thickness/";   //WARNING MAKE MORE GENERAL IN THE FUTURE
        string outDir4 = outDirT;   //MORE GENERAL IN THE FUTURE
        
        string fnameThick, fnameGeo, fnameTemp, fnamePos, fnameCM, fnamePoints, fContacts, fTractions, fStress; 
        fnameThick = outDir4 + "ThicknessLS_Bkg_g_" + std::to_string(grain_index) + "_t_" + std::to_string(step_time) + ".dat";    
        FILE * docThick = fopen(fnameThick.c_str(), "w");
        fnameGeo = outDir4 + "GeoLS_Bkg_g_" + std::to_string(grain_index)  + "_t_" + std::to_string(step_time) + ".dat";      
        FILE * docGeo = fopen(fnameGeo.c_str(), "w");
        fnameTemp = outDir4 + "TempLS_Bkg_g_" + std::to_string(grain_index)  + "_t_" + std::to_string(step_time) + ".dat";   
        FILE * docTemp = fopen(fnameTemp.c_str(), "w"); 
        fnamePos = outDir4 + "Pos_Bkg_g_" + std::to_string(grain_index)  + "_t_" + std::to_string(step_time) + ".dat";    
        FILE * docPos = fopen(fnamePos.c_str(), "w"); 
        fnameCM = outDir4 + "CenterMass_Bkg_g_" + std::to_string(grain_index)  + "_t_" + std::to_string(step_time) + ".dat";      
        FILE * docCM = fopen(fnameCM.c_str(), "w"); 
        fnamePoints = outDir4 + "GPoints_Bkg_g_" + std::to_string(grain_index)  + "_t_" + std::to_string(step_time) + ".dat";    
        FILE * docPoints = fopen(fnamePoints.c_str(), "w"); 
        fContacts = outDir4 + "NodesF_Bkg_g_" + std::to_string(grain_index)  + "_t_" + std::to_string(step_time) + ".dat";       
        FILE * docContacts = fopen(fContacts.c_str(), "w"); 
        fTractions = outDir4 + "TractionF_Bkg_g_" + std::to_string(grain_index)  + "_t_" + std::to_string(step_time) + ".dat";     
        FILE * docTractions = fopen(fTractions.c_str(), "w"); 
        fStress = outDir4 + "Stress_Bkg_g_" + std::to_string(grain_index)  + "_t_" + std::to_string(step_time) + ".dat";     
        FILE * docStress = fopen(fStress.c_str(), "w"); 
        
        size_t xxdimLS =  _grains[grain_index].getgrainTemp2D().getXdim();  // Useful dims 
        size_t yydimLS =  _grains[grain_index].getgrainTemp2D().getYdim();
        double xgridLS =  _grains[grain_index].getgrainTemp2D().getXdim();
        double ygridLS =  _grains[grain_index].getgrainTemp2D().getYdim();
        
        //Extra variables
        double lowLx, lowLy, lslocx, lslocy;
        double szx = xgridLS;
        double szy = ygridLS;
        Vector2d pointOcxy;
        //Shift to global coordinates and move to left bottom corner of level set
        lowLx = _grains[grain_index].getPosition()(0) - 0.5*szx +1.0; //IFF LS has unit of 1 for grid cells
        lowLy = _grains[grain_index].getPosition()(1) - 0.5*szy +1.0; //IFF LS has unit of 1 for grid cells
        
        
        //Properties
        fprintf(docThick, "%4.8f\n", _grains[grain_index].getRadius());
        fprintf(docThick, "%4.8f\n", xgridLS);
        fprintf(docThick, "%4.8f\n", ygridLS);
        fprintf(docGeo, "%4.8f\n", _grains[grain_index].getRadius());
        fprintf(docGeo, "%4.8f\n", xgridLS);
        fprintf(docGeo, "%4.8f\n", ygridLS);
        fprintf(docTemp, "%4.8f\n", _grains[grain_index].getRadius());
        fprintf(docTemp, "%4.8f\n", xgridLS);
        fprintf(docTemp, "%4.8f\n", ygridLS);
        fprintf(docPos, "%4.8f %4.8f %4.8f\n", _grains[grain_index].getPosition()(0), _grains[grain_index].getPosition()(1), _grains[grain_index].getTheta() );
        fprintf(docCM, "%4.8f\n", _grains[grain_index].getRadius());
        fprintf(docCM, "%4.8f\n", _grains[grain_index].getCmLset()(0));
        fprintf(docCM, "%4.8f\n", _grains[grain_index].getCmLset()(1));
        
        //Print Each sample level set
        for (size_t jj = 0; jj < yydimLS; jj++){
            for (size_t ii = 0; ii < xxdimLS; ii++){ 
                fprintf(docThick, "%4.8f\n", _grains[grain_index].getMthick().getLevelset()[jj*xxdimLS + ii] ); //Thickness LS
                fprintf(docGeo, "%4.8f\n", _grains[grain_index].getLset().getLevelset()[jj*xxdimLS + ii] ); //Geo LS
                
              lslocx = lowLx + 0.5 + ii; //IFF LS has unit of 1 for grid cells
              lslocy = lowLy + 0.5 + jj; //IFF LS has unit of 1 for grid cells
              pointOcxy << lslocx , lslocy;
                fprintf(docTemp, "%4.8f\n", round_Ocean(pointOcxy, _oceanTemp, _x_cells, _y_cells, _offset) ); //Temp LS
            }
        }  
        
        //Print points for meshing if needed
        for (size_t pti = 0; pti < _grains[grain_index].getPointList().size(); pti++){
            fprintf(docPoints, "%4.8f %4.8f\n", _grains[grain_index].getPointList()[pti](0), _grains[grain_index].getPointList()[pti](1) );
        }
        
        fclose(docThick);
        fclose(docGeo);
        fclose(docTemp);
        fclose(docPos);
        fclose(docCM);
        fclose(docPoints);
        
        //Output forces acting on grain
        size_t total_c_nodes = 0;
        //Node Forces (Contacts)
        for (size_t i = 0; i < _ngrains; i++){ 
            CData grainContactF;  //Temporary CData for all forces on grain
            vector<size_t> ptCIdx; //Vector that specifies exact point of force application (for Convenience)
            vector<size_t> ptCIdx2; //Vector that specifies exact point of force application (for Convenience)
            vector<size_t> ptCIdx3; //Vector that specifies exact point of force application (for Convenience)
            vector<size_t> ptCIdxF; //Vector that specifies exact point of force application (for Convenience)
            for (size_t j = i+1; j < _ngrains; j++) {
                // Normal contacts
                if (_grains[i].bcircleGrainIntersection(_grains[j])){
                    CData cDataContact = _grains[i].findContactDataIDX(_grains[j], _dt, 0, ptCIdx);
                    if (cDataContact._clocs.size() > 0) {
                        grainContactF += cDataContact;
                    }
                }
                // Contacts due to periodic bcs
                else if (_grains[i].bcircleGrainIntersectionXOffset(_grains[j], _offset(0))){
                    CData cDataContact = _grains[i].findContactDataIDX(_grains[j], _dt, _offset(0), ptCIdx2);
                    if (cDataContact._clocs.size() > 0) {
                         grainContactF += cDataContact;
                    }
                }
            }
            
            if (_grainsWall.size()>0)
            {    
                for (size_t j = 0; j < _ngrainsWall; j++) {
                    // Normal contacts
                    if (_grains[i].bcircleGrainIntersection(_grainsWall[j])){
                        CData cDataContact = _grains[i].findContactDataIDX(_grainsWall[j], _dt, 0, ptCIdx3);
                        if (cDataContact._clocs.size() > 0) {
                             grainContactF += cDataContact;
                        }
                    }
                } // close grain subloop
            }
        
            //Join all three index vectors
            for (size_t k = 0; k < ptCIdx.size(); k++){
                ptCIdxF.push_back(ptCIdx[k]);
            }
            for (size_t k = 0; k < ptCIdx2.size(); k++){
                ptCIdxF.push_back(ptCIdx2[k]);
            }
            for (size_t k = 0; k < ptCIdx3.size(); k++){
                ptCIdxF.push_back(ptCIdx3[k]);
            }
            
            //total_c_nodes = grainContactF._clocs.size();

            for (size_t ii = 0; ii < grainContactF._clocs.size(); ii++) {       
                if (grainContactF._cpairs[ii](0) == grain_index || grainContactF._cpairs[ii](1) == grain_index){ //Grain in question is in contact
                    fprintf(docContacts, "%d %d ", grainContactF._cpairs[ii](0), grainContactF._cpairs[ii](1) ); //grain contact 0 grain contact 1    
                    fprintf(docContacts, "%.8f %.8f ", grainContactF._forces[ii](0), grainContactF._forces[ii](1) ); // force x and y
                    fprintf(docContacts, "%.8f %.8f\n", grainContactF._clocs[ii](0), grainContactF._clocs[ii](1) );  // loc x and y in grain
                    //fprintf(cinfo3, "%d\n" , int(ptCIdxF[ii]) ); // Pointlist location
                    total_c_nodes += 1;
                }
            }
        }
        
        
        //Traction Forces
        
        double Chw = 1.3e-6; //Same as above
        double rhowater = 1.025e9; //Same as above
        Vector2d Uwf;
        Vector2d fluid_drag_force;
        double adjust_fff =  1e-4 * (1/0.000007); //1e3 few, 1e7 too much, even 5e3 too much Same as in grain!!! //Initially 1e2 seems okay but still tooooo much , 5e-2 still too high converge, 5e-4 still a bitty too high but improving
        Vector2d floeLocation = _grains[grain_index].getPosition();
        for (size_t i = 0; i < _y_cells; i++) {
            for (size_t j = 0; j < _x_cells; j++) {
               
                Vector2d pt_compare = _fluid_coord[j+i*_x_cells];
                bool contact_p = _grains[grain_index].under_ice(pt_compare);
                
                if (contact_p)
                {
                    //Assign current velocity
                    Uwf = _Uwg[j+i*_x_cells];
                    //SKIN DRAG SIMPLE
                fluid_drag_force = rhowater*Chw*(Uwf-_grains[grain_index].getVelocity()).norm()*(Uwf-_grains[grain_index].getVelocity()) * (_cell_sizex * _cell_sizey); 
                    
                    
                    //Based on centroid location
                    double dist_x = pt_compare(0) - floeLocation(0);
                    double dist_y = pt_compare(1) - floeLocation(1);
                    
                    fprintf(docTractions, "%4.8f %4.8f %4.8f %4.8f\n", dist_x, dist_y, fluid_drag_force(0) * adjust_fff, fluid_drag_force(1) * adjust_fff ); //pos x, pos y, fx, fy (pos wrt. to position or grain center) 
                } 
            }    
        } 
        
        //Save stress state for validation in FEM
        //Save stress info due to interfloe contact nodal forces, traction drag stress must be added later in FEM code since it is easier to assign loads in complex geometries this way
        double Volume = _grains[grain_index].getMass()/_grains[grain_index].getDensity();
        Vector3d grainStress; 
        
        cout << "Find grain node stress with node #: " << total_c_nodes << endl;
        if ( total_c_nodes < 1 ){
            grainStress << 0.0, 0.0, 0.0;
        }
        else{
             grainStress = _globalGrainState._grainStress[grain_index]/Volume;
        }
        
        //check grain yield stress
        double Sig1 = (grainStress(0)+grainStress(1))/2. + sqrt(pow((grainStress(0)-grainStress(1))/2.,2)+pow(grainStress(2),2));
        double Sig2 = (grainStress(0)+grainStress(1))/2. - sqrt(pow((grainStress(0)-grainStress(1))/2.,2)+pow(grainStress(2),2));
        double yield_prop = _grains[grain_index].getYield();
        
        fprintf(docStress, "%4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f\n", grainStress(0), grainStress(1), grainStress(2), Sig1, Sig2, yield_prop, Volume);    //Key:  Stress(0), Stress(1), Stress(2), Sigma Main Component 1, Sigma Main Component 2, Material Property Yield Stress, Volume
        
        
        fclose(docContacts);
        fclose(docTractions);
        fclose(docStress);
         
        return;
    }

    //New function for breakage
    // Checks every grain if it is flagged for fracture, then cuts grain with a line connecting the 2 highest contact force points
    void FracRoutine_Nova(vector<double> & BreakBin, vector<double> & Diameters, double & limitMaxDiamV, double & WaveH, double & WaveL, int & PROB_B, string & outDirT) {
        
        //Temporary use for testing
        for (size_t g = 0; g < _ngrains; g++){
            if (g == 0)
            {
                //Save force and LS Stats and stress state
                exportGrainDataBreak(g, outDirT);
            }
        }
        
        size_t nTout = 4233600/49; //Just for 2018!!!!!!!
        size_t step_time = _stepup/nTout;
        
        return;  //De-activate for now   
        
        //Just for a sample try
        if (step_time != 37){
            return;
        }

        cout << "GO BREAK ARBITRARY!" << endl;
        
        //cout << "GO BREAK!" << endl;
        vector<size_t> iter(_ngrains);
        for (size_t g = 0; g < _ngrains; g++){
            iter[g] = g;
        }
        cout <<"Random Shuffle" << endl;
        random_shuffle(iter.begin(),iter.end()); ////CHECK!!!
        cout << "ngrains control before loop: " << _ngrains << endl;
      for (size_t g1 = 0; g1 < _ngrains; g1++){
             //size_t g = iter[g1]; //CRITICAL For 1 random breakage (no problem since 1 grain chosen per step) MODIF_2 STAB   //Sep 30, 2022 Change for Prob. Break Try off 
             cout << "ngrains control in  loop: " << _ngrains << " grain in process: " << g1 <<  endl;
             size_t g = g1; //For thickness control (allows looping all grains) MODIF_2 STAB                                   //Sep 30, 2022 Change for Prob. Break Try on

            //cout <<"Grain to Try out for break i = " << g << endl;  //BUGGGGG
            //check if grain should fracture
            //cout << "Get Volume" << endl;
            double Volume = _grains[g].getMass()/_grains[g].getDensity();
          //cout << "Get Stress" << endl; 
            Vector3d grainStress = _globalGrainState._grainStress[g]/Volume;  //EDIT
            //Vector3d grainStress;
            //grainStress << 0.0, 0.0, 0.0;

            if ( (fabs(grainStress(0)) > 0 || fabs(grainStress(1)) > 0 || fabs(grainStress(2)) > 0) && g >= 0 ){
                //cout << "Volume: " << Volume << " grainStress: " << grainStress(0) << " " << grainStress(1) << " " << grainStress(2) << endl;
            }    

            //check grain yield stress
            double Sig1 = (grainStress(0)+grainStress(1))/2. + sqrt(pow((grainStress(0)-grainStress(1))/2.,2)+pow(grainStress(2),2));
            double Sig2 = (grainStress(0)+grainStress(1))/2. - sqrt(pow((grainStress(0)-grainStress(1))/2.,2)+pow(grainStress(2),2));

            // Von Mises criterion
//          double VonMises = sqrt(grainStress(0)*grainStress(0)+grainStress(1)*grainStress(1)
//                  -grainStress(0)*grainStress(1)+3.*grainStress(2)*grainStress(2));
//          _grains[g].changeCritStress(VonMises);

            // Tresca criterion
//          double Tresca  = fabs(Sig1-Sig2);
//          _grains[g].changeCritStress(Tresca);

            // maximum principal stress criterion
            _grains[g].changeCritStress(Sig1);  //EDIT


            //Diam Control
            //cout << "Get Diameter " << endl;
            vector<Vector2d> VecOriginZ = _grains[g].getPointList();
            double AreaZ = 0.0;
            size_t nZ = VecOriginZ.size();

            for (size_t i = 0; i < nZ-1; i++)
            {
              AreaZ += ( VecOriginZ[i](0) * VecOriginZ[i+1](1) -  VecOriginZ[i](1) * VecOriginZ[i+1](0) ); 
            }
            AreaZ += (VecOriginZ[nZ-1](0) * VecOriginZ[0](1) -  VecOriginZ[nZ-1](1) * VecOriginZ[0](0) ); 

            double IAreaZ = 0.5*abs(AreaZ);
            double IDiamZ = (2 * sqrt(IAreaZ /3.141592653589793));

            
            if (g >= 0){
            //cout << "CritSIG: " << Sig1 << " grainYield: " << _grains[g].getYield() << " for grain #: " << g << endl;
            }

            if (g >= 0){

                //if ( (Sig1 > _grains[g].getYield() || fabs(Sig2) > 100.*_grains[g].getYield()) && (_grains[g].getMass() >= _MELT_MASS) && (_grains[g].getMass() > 0.0)  ){ 
                if ((Sig1 > _grains[g].getYield() || fabs(Sig2) > 100.*_grains[g].getYield()) && (IDiamZ >= 20.00 && IDiamZ > 0.0)) {
                    //cout << "We HAVE Breakage for Grain: " << g << " with position X: " << _grains[g].getPosition()(0) << " and Y:" << _grains[g].getPosition()(1) << endl;
                    //cout << "Testing Mass: " << _grains[g].getMass() << " Melt Mass: " << _MELT_MASS;
                    //cout << "Contact points #: " << _globalGrainState._grainCoord[g].size() << endl;
                    _grains[g].changeFracFlag(true);
                }
                else
                {
                    //cout << "NO Breakage for Grain: " << g << " with position X: " << _grains[g].getPosition()(0) << " and Y:" << _grains[g].getPosition()(1) << endl;
                }

            }
            _grains[g].changeFracFlag(true);
            
            //After identifying grain to break, export data for external FEM
            
            // if ( _grains[g].getFracFlag() ){
            //     exportGrainDataBreak(g);
            // }
            
            //NOTE: This is repeated right now!!!
            

            //flag for open arbitrary break
            bool open_flag = false;
            open_flag = true;

            srand (time(NULL));

            //size_t og = g; //For thickness criterion

            int B_rand_index = rand() % PROB_B + 1;  //Rand prob intended inverted + 1  //Sep 30, 2022 Change for Prob. Break Try Stay
            //int B_rand_index = rand() % 100;
            //int prob = _BREAK_PROB;  //30  //Sep 30, 2022 Change for Prob. Break Try 
            Vector2d bk_point;
            int prob = PROB_B;

            double limitMaxDiam = limitMaxDiamV;
            double limitMaxDiamT = limitMaxDiam*0.5;


            double adjustIbr = 0.1;
            double Ibr = WaveH * _grains[g].getThickness() * 5.5e9 * 0.5 * pow(0.27e6,-1.0) * pow(WaveL, -2.0);
            Ibr *= adjustIbr;
            //cout << "Break Criterion from Dumont using meters and seconds as units (> 0.014 to break): " << Ibr << endl;
            
            //Obtain times for failure
            const double Tfail = _grains[g].getTfail();
            const double originTfail = _grains[g].getTorigin();
            //cout << "CONFIRM Tfail: " << _grains[g].getTfail() << endl;
            //cout << "CONFIRM Tfailorigin: " << _grains[g].getTorigin() << endl;
            const double curr_step = _stepup;
            //cout << "CONFIRM stepup: " << curr_step << endl;
            
            
            //FLOE BREAKAGE CRITERIA
            //Ibr = 0.00001; //ON BREAK ONLY USING THICKNESS REDUCTION DONT USE DUMONT OFF:Use Dumont as well
            //Use of random break, random grain, random angle, WAVE RELATED CRITERION
            //if ( B_rand_index < prob && (_grains[g].getMass() >= _MELT_MASS) && (_grains[g].getMass() > 0.0) ) //Then break
            if ( (IDiamZ >= limitMaxDiam) && (_grains[g].getMass() > 0.0) ) //Then break  //Using Ibr criterion //METHOD 2a Skip Ibr //MODIF_2 PREC   Sep 30, 2022 Change Off PROB CHANGE OFF
            //if ( Ibr > 0.014 && (IDiamZ >= limitMaxDiam) && (_grains[g].getMass() > 0.0) ) //Then break  //Using Ibr criterion //METHOD 2 MODIF_2 PREC
            //if ( B_rand_index >= prob - 2 && (IDiamZ >= limitMaxDiam) && (_grains[g].getMass() > 0.0) ) //Then break  //Using random occurrence  //METHOD 1  Sep 30, 2022 Change On   //Probability related criteria
            //if ( curr_step - originTfail >= Tfail )   //Oct 18, 2022 Time Failure Criterion Oct 18, 2022 PROB CHANGE ON
            {
              if (g == 0)
              {
                cout << "Tfail: " << Tfail << " at step: " << _stepup << endl;
                cout << "originTfail: " << originTfail << endl;
              }
              //cout << "Ibr BREAK" << endl;
              //size_t ng = rand() % (_ngrains);
              //cout << "To WAVE break: " << ng << " out of ngrains total: " << _ngrains << endl;
              //g = ng; ///CHECK!!! //ALREADY SHUFFLE NO NEED

             // //SIZE CONDITION
             // if (IDiamZ >= 20.00)
             // {
              
              //NEW OFF
              _grains[g].changeFracFlag(true);
              open_flag = true;
              //NEW OFF

             // }

             
              //To avoid breaking at center chose any point from 0 to Perc_Dist% of BBRox Radius randomly
              int perc_dist = 85;  //100
              //// 50% (rand() % 50) + 50
              bk_point(0) = (float(rand() % perc_dist))*(1.0/100.0)*_grains[g].getRadius();
              bk_point(1) = (float(rand() % perc_dist))*(1.0/100.0)*_grains[g].getRadius();

              //Fix Dist, also angle  (TEMPORAL)
              //bk_point(0) = 25*(1.0/100.0)*_grains[og].getRadius();
              //bk_point(1) = 25*(1.0/100.0)*_grains[og].getRadius();


              // Ice Breakage Criterion and Waves
              // Ibr = H_s  * h * Y * 0.5 * sigma_c^-1 * lambda ^ -2 > 0.014   Hs = sign. wave Height m, h = thickness m, lambda = peak wavelength, sigma_c = 0.27 MPa, Y = 5.5 GPa, lambda = cT or c/f, use Period and Stokes Wave Velocity to get lambda
              // Source Li et al. 2021, using Dumont et al. 2011

            }  

            //Thickness condition (if wave does not do it)
            double crit_thick = 0.005;  //0.1 Prior cases 1-78  Sep 30, 2022 Change on
            //double crit_thick =  0.1; //Try finer 1 cm average for all //0.2 m or 2e-4 km    Sep 30, 2022 Change off
            
            if (_grains[g].getThickness() < crit_thick && IDiamZ >= limitMaxDiam && open_flag == false && _grains[g].getMass() > 0.0) //Be more lenient for breakage? WONT Happen for only B
            //if (_grains[og].getThickness() < crit_thick && IDiamZ >= 20.00 && ng != og )
            {
              cout << "THICK BREAK for thickness: " << _grains[g].getThickness() << endl;
              //cout << "To THICK break: " << og << " out of ngrains total: " << _ngrains << endl;
              _grains[g].changeFracFlag(true);
              open_flag = true;

              //To avoid breaking at center chose any point from 0 to Perc_Dist% of BBRox Radius randomly
              int perc_dist = 85;  //100
              //// 50% (rand() % 50) + 50
              bk_point(0) = (float(rand() % perc_dist))*(1.0/100.0)*_grains[g].getRadius();
              bk_point(1) = (float(rand() % perc_dist))*(1.0/100.0)*_grains[g].getRadius();

              //Fix Dist, also angle  (TEMPORAL)
              //bk_point(0) = 25*(1.0/100.0)*_grains[og].getRadius();
              //bk_point(1) = 25*(1.0/100.0)*_grains[og].getRadius();

              //g = og; ///CHECK!!!
            }

            double ang_line = rand() % 1257; //From -2pi to 2pi              
            ang_line -= 628;
            ang_line *= 0.01;
            
            //ang_line = 0.5*3.1415926535897932384; //Fix to pi/4 for the mean time also distance (TEMPORAL)

            // cout << "Breakage Real Start for grain: " << g << endl;
            //if ( (_grains[g].getFracFlag() && _globalGrainState._grainCoord[g].size() >= 2) ||  (_grains[g].getFracFlag() && open_flag)      ) {
            if (    (_grains[g].getFracFlag() && open_flag)     ) { //Only thickness prob., no stress for now              
            //cout << "Breakage Real Start for grain IN: " << g << endl;
//              cout << "coord # " << _globalGrainState._grainCoord[g].size() << endl;
                /*
                 * Get Grain 1 Info
                 */
                //cout << "Breakage Real Start for grain define: " << g << endl;
                Vector2d pF1;
                Vector2d pF2;
                Vector2d g1cm = _grains[g].getCmLset();
                vector<Vector2d> g1p0 = _grains[g].getRefPointList();
                size_t np1 = g1p0.size();
                size_t gid = _grains[g].getId();
                Levelset2d g1lset = _grains[g].getLset();
                vector<double> g1Lset = g1lset.getLevelset();
                size_t Xdim = g1lset.getXdim();
                size_t Ydim = g1lset.getYdim();

                vector<PointInfo> g1i(np1);

                const double cos1 = cos(_grains[g].getTheta());
                const double sin1 = sin(_grains[g].getTheta());
                
                //cout << "Breakage Real Start for grain LOOP Point contact " << g << endl;
                for (size_t i = 0; i<np1;i++){
                    //must unrotate and bring to origin !!!KEEP AN EYE ON THIS!!!!
//                  g1i[i].point << (g1p[i](0)-_grains[g].getPosition()(0))*cos1 + (g1p[i](1)-_grains[g].getPosition()(1))*sin1,
//                                            -(g1p[i](0)-_grains[g].getPosition()(0))*sin1 + (g1p[i](1)-_grains[g].getPosition()(1))*cos1;
                    g1i[i].point = g1p0[i];
                    //
                    g1i[i].shear = _grains[g].getNodeShears()[i];
                    g1i[i].contact = _grains[g].getNodeContact()[i];
                }

                //cout << "Testing Mass: " << _grains[g].getMass() << " Melt Mass: " << _MELT_MASS << " Diameter Size: " << IDiamZ;

                //Arbitrary fracture pre-liminaries
                //Use a stress-based and thickness and yield criterion to decide if fracture or not.
                //If decided Use loads and thickness data to get internal stress, use points to get mesh, find damage in external FEM from this loads and stress
                //Use damage map to identify location of fracture, use fracture location to create an arbitrary splitter with python script
                
                //Run Python script here:
                
                // REMEMBER TO TURN ON
                //std::string filenameArb = "arbitrary_fracture_testv2.py;";  //Might require an HPC shell to run python3 arbitrary_fracture_testv2.py and a pause as a result, remember to adapt script for HPC (NEW)
                //std::string command = "python3 ";
                //command += filenameArb;
                //system(command.c_str());
                // //Optional (requires making .sh HPC file)
                // std::string commandshell = "sbatch -A geomechanics HPCArbPython.sh;";
                // system(commandshell.c_str());
                // std::string commandp = "pause 120s;"
                // system(commandp.c_str());
                
                
                //Arbitrary fracture step 1:
                //Import the points of the arbitrary splitter as a text file for Vector2d, validate if there are only 2 points (very unlikely), also import if orientation is X (1) or Y (0) from another small file created by the python script
                cout << "Import external splitter" << endl;
                char tempfnameS[300];
                //string outDir4 = "./Output/SeaIce2018_1077/Vertical_Thickness/";  //WARNING, python script must aim here and generalize in the future.
                string outDir4 = "./Output/";  //WARNING, python script must aim here and generalize in the future.
                sprintf(tempfnameS, (outDir4 + "arb_splitter_points.dat").c_str() );
                string file_arbSplitter = tempfnameS;
                size_t xdom = 1; // Define as 1 as default for xorientation
                cout << "Read external splitter file" << endl;
                vector<Vector2d> splpts0_arb = arbSplitterFile(file_arbSplitter, xdom); //Read arbitrary splitter file and convert to 2D Pointlist in Ref. Configuration and indicate if xdom = true or false (NEW)
                vector<Vector2d> splpts0_arb_LS = splpts0_arb; //For level set comparisons
                
                /*
                 * Find Splitter Shape
                 */
                if (open_flag == false)
                {
                    //cout << "TwoHighestForcesLocs" << endl;
                    TwoHighestForcesLocs(pF1,pF2,g,g1i);
                }

                
                //Splitter for two arbitrary points not based on forces
                if (open_flag == true)
                {
                   //cout << "Arbitrary Cut" << endl;
                   //replacing pF1, pF2
                   //Point_Slicer(pF1a, pF2a)
                   Vector2d pF1a;
                   Vector2d pF2a;
                   Vector2d LPA; 
                   Vector2d LPB;
                   LPA(0) = bk_point(0) +  _grains[g].getRadius()*cos(ang_line);
                   LPA(1) = bk_point(1) + _grains[g].getRadius()*sin(ang_line);
                   LPB(0) = bk_point(0) -_grains[g].getRadius()*cos(ang_line);
                   LPB(1) = bk_point(1) -_grains[g].getRadius()*sin(ang_line);
                   //cout << "ClosestPointLine" << endl;
                   ClosestPoint_line(pF1a, g1p0, LPA);
                   ClosestPoint_line(pF2a, g1p0, LPB);
                   pF1 = pF1a;
                   pF2 = pF2a;
                   cout <<  "bk_point: " << bk_point(0) << " " << bk_point(1) << endl;
                   cout <<  "pF1: " << pF1(0) << " " << pF1(1) << endl;
                   cout <<  "pF2: " << pF2(0) << " " <<  pF2(1) << endl;
                }  

                //cout << "Get pF2" << endl;
                double bigV = M_PI*_grains[g].getRadius()*_grains[g].getRadius();
                double g1V  = _grains[g].getMass()/_grains[g].getDensity();
                if (bigV/g1V > _fracProps.getRoMratio() || pF2 == pF1){
                    ClosestPoint(pF1,g1p0);
                    pF2 = Vector2d(0.,0.);
                }

                //cout << "pF1: " << pF1 << endl;
                //cout << "pF2: " << pF2 << endl;
                
                /*
                 * Create Splitter
                 */

                //cout << "Make Splitter" << endl;
                //get points that define split line in the level set reference frame
                Vector2d splitpt1 = pF1 + g1cm;
                Vector2d splitpt2 = pF2 + g1cm;

                //cout << "Making Splitter and Surface Lines" << endl;
                //make splitter surface points
                vector<Vector2d> splpts0 = _fracProps.MakeSplitter(_grains[g],splitpt1,splitpt2);
                cout << "Print Splitter Point Output as reference of normal method: " << endl;
                for(size_t i=0; i<splpts0.size();i++){
                    cout << splpts0[i](0) << " , " << splpts0[i](1) << endl;
                }
                //Fracture nova in multi mode will replace this function with a python script that generates splitter points  and then we read this point list file above (check for correct coordinate reference) and call it splpts0_arb
                
        
                //cout << "Splts Pts: " << splpts0[0] <<  endl;
                //make splitter level set
                //vector<double> splitset = g1lset.makeFracSurfaceLine(splitpt1,splitpt2);
                //cout << "Splts Set: " << splitset[0] << endl;
                //output splitset if desired
                // _fracProps.SaveSplitter(splpts0,splitset,Xdim,Ydim);
                
                //Arbitrary fracture step 2:
                cout << "Print Splitter Point Output ARBITRARY TO DEBUG: " << endl;
                for(size_t i=0; i<splpts0_arb.size();i++){
                    cout << splpts0_arb[i](0) << " , " << splpts0_arb[i](1) << endl;
                }
                //Reference using cm for LS_Set comparison (if line is centered in 0,0)  #WARNING CONFIRM IF YOU NEED ROTATION FOR THIS or for splpts0_arb!!!!
                for(size_t ia=0; ia<splpts0_arb_LS.size();ia++){
                    splpts0_arb_LS[ia] += g1cm;
                }
                
                //Replace splitset function with new function that use arbitrary shape pointList and x-y direction to create split level set (NEW)
                cout << "makeFracSurfaceLine_arb" << endl;
                vector<double> splitset_arb = g1lset.makeFracSurfaceLine_arb(splpts0_arb_LS, xdom);
                

                /*
                 * Make New Grains (using Lset and Surface Points)
                 */


                //make new Lsets
                //cout << "Make New LSets" << endl;
               // vector<double> g2Lset(splitset.size());
                //vector<double> g3Lset(splitset.size());
                // for(size_t i=0; i<splitset.size();i++){
                //     g2Lset[i] = max(splitset[i],g1Lset[i]);
                //     g3Lset[i] = max(splitset[i]*-1.,g1Lset[i]);
                // }
                
                vector<double> g2Lset(splitset_arb.size());
                vector<double> g3Lset(splitset_arb.size());
                
                //This function doesn't change for arbitrary fracture. If your splitSet is correct, everything is fine and the new level sets come from this processing.
                for(size_t i=0; i<splitset_arb.size();i++){
                    g2Lset[i] = max(splitset_arb[i],g1Lset[i]);
                    g3Lset[i] = max(splitset_arb[i]*-1.,g1Lset[i]);
                }
                
                
                //cout << "Splts Set2: " << g2Lset[0] << endl;
                //cout << "Splts Set3: " << g3Lset[0] << endl;

                //reinit due to imperfect set ops
                //g2Lset = reinit(g2Lset); //need to build in Levelset2d.h
                //g3Lset = reinit(g3Lset);

                Levelset2d g2lset(g2Lset, Xdim, Ydim);
                Levelset2d g3lset(g3Lset, Xdim, Ydim);

                //2)pts
                vector<PointInfo> g2i, g3i;

                //cout << "Work with Splitter Points" << endl;
                //decides which splitter points stay
                //Input: splpts0, g2lset, g1cm;  Output: g2i, g3i
                //_fracProps.KeepSplitterPoints(splpts0,g1lset,g1cm,g2i,g3i); //g2i and g3i start the same
                //This function aims to keep points that are only inside the grain, the python functions does this, but it can be used again for safety. It is not supposed to affect anything and can be kept ALMOST the same.
                
                //Arbitrary fracture Step 3: Keep correct splitter points and make the two separate grain point sets (NEW)
                cout << "KeepSplitterPoints_arb" << endl;
                _fracProps.KeepSplitterPoints_arb(splpts0_arb, g1lset, g1cm, g2i, g3i);
                
                
                //decides which previous grain points stay
                //Input: splitpt1, splitpt2, g1i, g1cm;  Output: g2i, g3i
                //_fracProps.KeepGrainPoints(splitpt1,splitpt2,g1i,g1cm,g2i,g3i);  //g2i and g3i are differentiated
                
                cout << "KeepGrainPoints_arb" << endl;
                //Include new function for 2 and 3 grain points using arbitrary pointlist and x-y orientation (NEW)
                _fracProps.KeepGrainPoints_arb(splpts0_arb_LS, g1i, g1cm, g2i, g3i, xdom);  //g2i and g3i are differentiated

                //FROM HERE FRACTURE STAYS THE SAME. You use LS2, LS3, g2i and g3i to make new grains as usual.
                
                
                //make sure potential grains has a reasonable amount of surface points
                if (g2i.size()<5 || g3i.size()<5 ){ 
                //Adjust for size
                //if (g2i.size()<3 || g3i.size()<3 ){   
                    //cout << "Break WARNING: TOO SMALL NaN" << endl;
                    //return; //to ensure no NaNs in the crumbles
                    continue;
                }


                //save number of points for grains
                size_t np2 = g2i.size();
                size_t np3 = g3i.size();

                /*
                 * Find Physical Properties
                 */
                //cout << "Find New Properties" << endl;
                double mass2;
                Vector2d g2cm;
                double I2;
                //_fracProps.findGrainProps0(g2lset, mass2, g2cm, I2);
                _fracProps.findGrainProps0D(g2lset, mass2, g2cm, I2, _grains[g].getDensity());

                double mass3;
                Vector2d g3cm;
                double I3;
                //_fracProps.findGrainProps0(g3lset, mass3, g3cm, I3);
                _fracProps.findGrainProps0D(g3lset, mass3, g3cm, I3, _grains[g].getDensity());
                /*
                 * Position grains correctly, removing remnants of the grain 1 reference frame
                 *
                 */

                Vector2d position2 = _grains[g].getPosition()-Vector2d((g1cm(0)-g2cm(0))*cos1 - (g1cm(1)-g2cm(1))*sin1,
                                                              (g1cm(0)-g2cm(0))*sin1 + (g1cm(1)-g2cm(1))*cos1);
                double maxR2 = 0;
                vector<Vector2d> g2p0, g3p0;

                //shift points to correct position
                for (size_t i = 0; i<np2;i++){

                    g2i[i].point += g1cm-g2cm;
                    if (g2i[i].point.norm() > maxR2){
                        maxR2 = g2i[i].point.norm();
                    }
                }

                // reorder points to counterclockwise
                sort(g2i.begin(), g2i.end(),pointsort);

                for (size_t i=0;i<np2;i++){
                    g2p0.push_back(g2i[i].point);
                }


                Vector2d position3 = _grains[g].getPosition()-Vector2d((g1cm(0)-g3cm(0))*cos1 - (g1cm(1)-g3cm(1))*sin1,
                                                                  (g1cm(0)-g3cm(0))*sin1 + (g1cm(1)-g3cm(1))*cos1);
                double maxR3 = 0;
                //shift points to correct position
                for (size_t i = 0; i<np3;i++){
                    g3i[i].point += g1cm-g3cm;

                    //Find bbox radius
                    if (g3i[i].point.norm() > maxR3){
                        maxR3 = g3i[i].point.norm();
                    }
                }
                // sort indexes based on comparing values in g3p0 b
                sort(g3i.begin(), g3i.end(),pointsort);

                for (size_t i=0;i<np3;i++){
                    g3p0.push_back(g3i[i].point);
                }


                //check thinness of potential grains in case they could slide into other grains  !!!!TODO: ADD TO MELT
                Vector2d g2min,g3min; //closest surface point to centroid
                double g2minR,g3minR; //minimum distance from centoid to surface
                ClosestPoint(g2min,g2p0);
                ClosestPoint(g3min,g3p0);
                g2minR = g2min.norm();
                g3minR = g3min.norm();
                if (g2minR < .5|| g3minR < .5) {
                    //cout << "Break Radius NaN" << endl;
                    //return;
                    continue;
                }


//              if (M_PI*maxR2*maxR2/(mass2) > _fracProps.getMaxRoM() || M_PI*maxR3*maxR3/(mass3) > _fracProps.getMaxRoM() ){
//                  return; //to ensure no NaNs in the crumbles
//              }

                /*
                 * Create grains
                 */
                //Data Which Will NOT depend on fracture
                vector<double> g2tempV =  _grains[g].getgrainTemp2D().getLevelset();
                vector<double> g3tempV =  _grains[g].getgrainTemp2D().getLevelset(); //g2lset g3lset

                // vector<double> g2tempV =  g2lset.getLevelset();
                // vector<double> g3tempV =  g3lset.getLevelset();

                vector<double> g2ThickV =  _grains[g].getMthick().getLevelset();
                vector<double> g3ThickV =  _grains[g].getMthick().getLevelset(); //g2lset g3lset

                //Pre-correct if needed for thickness
                vector<double> g2Geo =  g2lset.getLevelset();
                vector<double> g3Geo =  g3lset.getLevelset(); //g2lset g3lset

                for (size_t ii = 0; ii<g2lset.getXdim(); ii++) 
                {            
                  for (size_t jj = 0; jj<g2lset.getYdim(); jj++)
                  {
                    if ( g2Geo[(jj*g2lset.getXdim())+ii]>0.0 && g2ThickV[(jj*g2lset.getXdim())+ii] > 0.0 ) //You need to redefine thickness in case LS is changed by breakage or other reasons, to avoid bugs
                    {
                          //cout << "Correcting Thickness after break Origin!!!" << endl;
                          g2ThickV[(jj*g2lset.getXdim())+ii] = -_grains[g].getThickness(); 
                    }
                  }
                } 
                for (size_t ii = 0; ii<g3lset.getXdim(); ii++) 
                {            
                  for (size_t jj = 0; jj<g3lset.getYdim(); jj++)
                  {
                    if ( g3Geo[(jj*g3lset.getXdim())+ii]>0.0 && g3ThickV[(jj*g3lset.getXdim())+ii] > 0.0 ) //You need to redefine thickness in case LS is changed by breakage or other reasons, to avoid bugs
                    {
                          //cout << "Correcting Thickness after break Origin!!!" << endl;
                          g3ThickV[(jj*g3lset.getXdim())+ii] = -_grains[g].getThickness(); 
                    }
                  }
                }    
                
                //Get average thickness values
                double thick_2 = 0.000;
                double thick_3 = 0.000;
                size_t ct_2 = 0;
                size_t ct_3 = 0;
                
                //Count for average thickness    
                for (size_t ii = 0; ii<g2lset.getXdim(); ii++) 
                {            
                  for (size_t jj = 0; jj<g2lset.getYdim(); jj++)
                  {
                    if ( g2ThickV[(jj*g2lset.getXdim())+ii] > 0.0 ) //Count
                    {
                        thick_2 += g2ThickV[(jj*g2lset.getXdim())+ii];
                        ct_2 += 1;
                    }
                  }
                }
                for (size_t ii = 0; ii<g3lset.getXdim(); ii++) 
                {            
                  for (size_t jj = 0; jj<g3lset.getYdim(); jj++)
                  {
                    if ( g3ThickV[(jj*g3lset.getXdim())+ii] > 0.0 ) //Count
                    {
                        thick_3 += g3ThickV[(jj*g3lset.getXdim())+ii];
                        ct_3 += 1;
                    }
                  }
                }

                //Find value to average new thickness
                if (ct_2 > 0)
                {
                    thick_2 /= double(ct_2);
                }
                else{
                    thick_2 = 0.000;
                }
                if (ct_3 > 0)
                {
                    thick_3 /= double(ct_3);
                }
                else{
                    thick_3 = 0.000;
                }
                mass2 *= thick_2*0.001;
                I2 *= thick_2;
                mass3 *= thick_3*0.001;
                I3 *= thick_3;                

                Levelset2d g2temp(g2tempV, g2lset.getXdim(), g2lset.getYdim());
                Levelset2d g3temp(g3tempV, g3lset.getXdim(), g3lset.getYdim());

                Levelset2d g2ThickLS(g2ThickV, g2lset.getXdim(), g2lset.getYdim());
                Levelset2d g3ThickLS(g3ThickV, g3lset.getXdim(), g3lset.getYdim());                

                //cout << "Creating grain2 with npoints: " << g2p0.size() << endl;
                //cout << "Grain 2 Points:" << endl;
                for (size_t i = 0; i < g2p0.size() ; i++ )
                {
                    //cout << g2p0[i](0) << " " << g2p0[i](1) << endl;
                }

                //cout << "Making New Grains Split 2" << endl;
                bool tempflag = false;
                int temploc = 0;
                const double current_step = _stepup;
                //Grain2d g2 = Grain2d(mass2, position2, _grains[g].getVelocity(), I2, _grains[g].getTheta(), _grains[g].getOmega(),
                //                        g2cm, g2p0, maxR2, g2lset, _grains[g].getKn(), _grains[g].getKs(), _grains[g].getMu(), _maxId+1);
                Grain2d g2 = Grain2d(mass2, position2, _grains[g].getVelocity(), 
                 mass2, _grains[g].getTemperature(), _grains[g].getThickness(), g2ThickLS, _grains[g].getTemperature(), _grains[g].getThickness0(), _grains[g].getUtemper(), _grains[g].getUtemper(), g2temp, g2temp, 
                 I2, _grains[g].getTheta(), _grains[g].getOmega(), g2cm, g2p0, g2p0.size(), maxR2, 
                 g2lset, g2lset, 0, _grains[g].getId()+1, _grains[g].getKn(), _grains[g].getKs(), _grains[g].getMu(), tempflag, temploc,
               _grains[g].getgrainStress2D(), _grains[g].getgrainDamage() , _grains[g].getgrainThickness(), g2p0, _grains[g].getTfail(), current_step); //_grains[g].getId()+1
               g2.changeDensity(_grains[g].getDensity());
                //cout << "For grain2 Tfail: " << g2.getTfail() << " at step: " << _stepup << endl;
                //cout << "For grain2 originTfail: " << g2.getTorigin() << endl;

                // double new_mass;
                // double new_I;
                // Vector2d new_cm;
                // Levelset2d ref_lset;

                // ref_lset = g2.getLset();
                // _fracProps.findGrainProps0(ref_lset, new_mass, new_cm, new_I);
                // g2.changeInertia(new_I);
                // g2.changeMass(new_mass*(10e6)); 
                // g2.changeCM(new_cm);

                for (size_t i=0;i<np2;i++){
                    g2.getNodeShearsNonConst()[i]  = g2i[i].shear;
                    g2.getNodeContactNonConst()[i] = g2i[i].contact;
                }

                //cout << "Save New Grains Split 2: " << endl;
                //_fracProps.SaveGrain(g2,g2p0);


                //Grain2d g3 = Grain2d(mass3, position3, _grains[g].getVelocity(), I3, _grains[g].getTheta(), _grains[g].getOmega(),
                //                        g3cm, g3p0, maxR3, g3lset, _grains[g].getKn(), _grains[g].getKs(), _grains[g].getMu(), _maxId+2);
                

                //cout << "Creating grain3 with npoints: " << g3p0.size() << endl;
                //cout << "Grain 3 Points:" << endl;
                for (size_t i = 0; i < g3p0.size() ; i++ )
                {
                    //cout << g3p0[i](0) << " " << g3p0[i](1) << endl;
                }

                //cout << "Making New Grains Split 3" << endl;
                Grain2d g3 = Grain2d(mass3, position3, _grains[g].getVelocity(), 
                             mass3, _grains[g].getTemperature(), _grains[g].getThickness(), g3ThickLS, _grains[g].getTemperature(), _grains[g].getThickness0(), _grains[g].getUtemper(), _grains[g].getUtemper(), g3temp, g3temp, 
                             I3, _grains[g].getTheta(), _grains[g].getOmega(), g3cm, g3p0, g3p0.size(), maxR3, 
                             g3lset, g3lset, 1, _grains[g].getId()+2, _grains[g].getKn(), _grains[g].getKs(), _grains[g].getMu(), tempflag, temploc,
                           _grains[g].getgrainStress2D(), _grains[g].getgrainDamage() , _grains[g].getgrainThickness(), g3p0, _grains[g].getTfail(), current_step);  //_grains[g].getId()+2
                g3.changeDensity(_grains[g].getDensity());
                
                //cout << "For grain3 Tfail: " << g3.getTfail() << " at step: " << _stepup << endl;
                //cout << "For grain3 originTfail: " << g3.getTorigin() << endl;


                // ref_lset = g3.getLset();
                // _fracProps.findGrainProps0(ref_lset, new_mass, new_cm, new_I);
                // g3.changeInertia(new_I);
                // g3.changeMass(new_mass*(10e6));
                // g3.changeCM(new_cm);

//              for (size_t i=0;i<np3;i++){
//                  g3.getNodeShearsNonConst()[i]  = g3i[i].shear;
//                  g3.getNodeContactNonConst()[i] = g3i[i].contact;
//              }

                //send shears to "master" grains
                //cout << "send shears to master grains" << endl;
//              double dist = 50;
////                cout << "g3isize " << g3i.size() << " " << np3 << endl;
//              vector<Vector2d> g3_moved_pts = g3.getPointList();
                for (size_t i=0;i<np3;i++){
//                  if (g3i[i].contact > _maxId){ //contact with wall; shear stays
                        g3.getNodeShearsNonConst()[i] = g3i[i].shear;
                        g3.getNodeContactNonConst()[i] = g3i[i].contact;

////                        cout << "\nwall " << g3i[i].contact << "\n"<<endl;
//                  }
//                  else if (g3i[i].contact>0){ //contact with grain; shear goes to other grain
//                      //loop over each possible point and search for closest one
//                      size_t pt_id;
//
//                      vector<Vector2d> target_grain_pts = _grains[FindGrainFromId(g3i[i].contact)].getPointList();
//                      for (size_t j = 0;j<target_grain_pts.size();j++){
//                          if ((g3_moved_pts[i]-target_grain_pts[j]).norm() < dist){
//                              dist = (g3_moved_pts[i]-target_grain_pts[j]).norm();
//                              pt_id = j;
//                          }
//                      }
////                        cout << "\ndistance " << dist << "\n"<<endl;
//                      _grains[FindGrainFromId(g3i[i].contact)].getNodeContactNonConst()[pt_id] = g3.getId();
//                      _grains[FindGrainFromId(g3i[i].contact)].getNodeShearsNonConst()[pt_id] = g3i[i].shear;
//                  }
                }

                //cout << "Save New Grains Split 3: " << endl;
                //_fracProps.SaveGrain(g3,g3p0);

                double yield = _grains[g].getYield();
                g3.changeYield(yield*pow(_grains[g].getDiameter()/g3.getDiameter(),3./3.));
                g2.changeYield(yield*pow(_grains[g].getDiameter()/g2.getDiameter(),3./3.));


                //find shears where g was slave and move to master
                for (size_t i = 0; i<g;i++){
                    //check if grain was a master grain on g
                    vector<size_t> master_contacts = _grains[i].getNodeContact();
                    bool master = std::binary_search(master_contacts.begin(),master_contacts.end(),gid); //true if master grain was contacting fractured grain
                    if (master) {
//                      cout << "Position " << g2.getPosition().transpose() << " " << g3.getPosition().transpose() << " " << _grains[i].getPosition().transpose() << endl;
                        //move shears to g2,g3
                        //find node
                        vector<size_t> nodes;
                        size_t n=0;
                        for (size_t j = 0; j<master_contacts.size();j++){
                            if (master_contacts[j] == gid){
                                nodes.push_back(j);
                                n++;
                            }
                        }

                        //get point location
                        vector<Vector2d> master_pts = _grains[i].getPointList();

                        vector<Vector2d> g2_pts = g2.getPointList();
                        vector<Vector2d> g3_pts = g3.getPointList();

                        for(size_t j=0;j<n;j++){
                            //find closest pt on either g2 or g3
                            Vector2d master_pt = master_pts[nodes[j]];
                            size_t kmin = 0;
                            double dist = 50.;
                            for (size_t k = 0; k<np2;k++){
                                double dist_k = (master_pt-g2_pts[k]).norm();
//                              cout << "dist_k " << dist_k << endl;
                                if (dist_k < dist){
                                    dist = dist_k;
                                    kmin = k;
                                }
                            }
                            for (size_t k = np2; k<np2+np3;k++){
                                double dist_k = (master_pt-g3_pts[k-np2]).norm();
                                if (dist_k < dist){
                                    dist = dist_k;
                                    kmin = k;
                                }
                            }
//                          cout << "dist " << dist << endl;
                            if (dist < 1) {
                                //send shears to fragments
                                if (kmin>=np2) {
                                    kmin -= np2; //on grain 3
                                    g3.getNodeShearsNonConst()[kmin]  = _grains[i].getNodeShears()[nodes[j]];
                                    g3.getNodeContactNonConst()[kmin] = _grains[i].getNodeContact()[nodes[j]];
                                } else { //on grain 2
                                    g2.getNodeShearsNonConst()[kmin]  = _grains[i].getNodeShears()[nodes[j]];
                                    g2.getNodeContactNonConst()[kmin] = _grains[i].getNodeContact()[nodes[j]];
                                }
                                _grains[i].getNodeShearsNonConst()[nodes[j]] = 0;
                                _grains[i].getNodeContactNonConst()[nodes[j]] = 0;
                            }
                        }
                    }
                }


                //Bin Mass changes due to Breakage
                //Initial Area and Diameter before Breaking
                //double IArea = PointsArea(_grains[g].getPointList());
                //double IDiam = MeanCaliperD(IArea);

                //cout << "Find Area 2" << endl;
                vector<Vector2d> VecOrigin = _grains[g].getPointList();
                double Area = 0.0;
                size_t n = VecOrigin.size();

                for (size_t i = 0; i < n-1; i++)
                {
                  Area += ( VecOrigin[i](0) * VecOrigin[i+1](1) -  VecOrigin[i](1) * VecOrigin[i+1](0) ); 
                }
                Area += (VecOrigin[n-1](0) * VecOrigin[0](1) -  VecOrigin[n-1](1) * VecOrigin[0](0) ); 

                double IArea = 0.5*abs(Area);
                double IDiam = (2 * sqrt(IArea /3.141592653589793));
                double F1Area, F2Area;

                //Avoid degenerate breakage
                if (isnan(IArea))
                {
                    //cout << "Breakage Area NaN" << endl;
                    //return;
                    continue;
                }

                //Final Area and Diameter after breaking in 2
                //double F1Area = PointsArea(g2.getPointList());
                //double F1Diam = MeanCaliperD(F1Area);
                //double F2Area = PointsArea(g3.getPointList());
                //double F2Diam = MeanCaliperD(F2Area);

                VecOrigin = g2.getPointList();
                Area = 0.0;
                n = VecOrigin.size();

                for (size_t i = 0; i < n-1; i++)
                {
                  Area += ( VecOrigin[i](0) * VecOrigin[i+1](1) -  VecOrigin[i](1) * VecOrigin[i+1](0) ); 
                }
                Area += (VecOrigin[n-1](0) * VecOrigin[0](1) -  VecOrigin[n-1](1) * VecOrigin[0](0) ); 

                if (isnan(Area))
                {
                    F1Area = IArea * (g2.getMass()/_grains[g].getMass());
                }
                else
                {
                    F1Area = 0.5*abs(Area);
                }

                
                //cout << "Find Area 3" << endl;
                VecOrigin = g3.getPointList();
                Area = 0.0;
                n = VecOrigin.size();

                for (size_t i = 0; i < n-1; i++)
                {
                  Area += ( VecOrigin[i](0) * VecOrigin[i+1](1) -  VecOrigin[i](1) * VecOrigin[i+1](0) ); 
                }
                Area += (VecOrigin[n-1](0) * VecOrigin[0](1) -  VecOrigin[n-1](1) * VecOrigin[0](0) ); 

                if (isnan(Area))
                {
                    F2Area = IArea * (g3.getMass()/_grains[g].getMass());
                }
                else
                {
                    F2Area = 0.5*abs(Area);
                }
                
                
                if ((F2Area + F1Area) > (1.2*IArea) ||  (F2Area + F1Area) < (0.8*IArea))   //20% error respected otherwise force using mass proportions
                {
                    F1Area = IArea * (g2.getMass()/_grains[g].getMass());
                    F2Area = IArea * (g3.getMass()/_grains[g].getMass());
                }
                  

                double F1Diam = (2 * sqrt(F1Area /3.141592653589793));
                double F2Diam = (2 * sqrt(F2Area /3.141592653589793));

                //Find bin that lost this area (assume minor break, one staying on same bin other changes)
                size_t loc_index, loc_index2, loc_index3;
                for (size_t bi = 0; bi < Diameters.size()-1; bi++){  
                    if ( IDiam <= Diameters[bi] &&  IDiam >= Diameters[bi+1] ) 
                    {
                       loc_index = bi;
                    }
                    else if ( IDiam < Diameters[Diameters.size()-1] )  //Min size threshold control
                    {
                       loc_index = Diameters.size()-2;
                    }
                    else if ( IDiam > Diameters[0] )  //Max size threshold control
                    {
                       loc_index = 0;
                       break;
                    }
                }
                for (size_t bi = 0; bi < Diameters.size()-1; bi++){  
                    if ( F1Diam <= Diameters[bi] &&  F1Diam >= Diameters[bi+1] ) 
                    {
                       loc_index2 = bi;
                    }
                    else if ( F1Diam < Diameters[Diameters.size()-1] )  //Min size threshold control
                    {
                       loc_index2 = Diameters.size()-2;
                    }
                    else if ( F1Diam > Diameters[0] )  //Max size threshold control
                    {
                       loc_index2 = 0;
                       break;
                    }
                }
                for (size_t bi = 0; bi < Diameters.size()-1; bi++){  
                    if ( F2Diam <= Diameters[bi] &&  F2Diam >= Diameters[bi+1] ) 
                    {
                       loc_index3 = bi;
                    }
                    else if ( F2Diam < Diameters[Diameters.size()-1] )  //Min size threshold control
                    {
                       loc_index3 = Diameters.size()-2;
                    }
                    else if ( F2Diam > Diameters[0] )  //Max size threshold control
                    {
                       loc_index3 = 0;
                       break;
                    }
                }

                //cout << "BKG: AreaOriginal: " << IArea << " Area2: " << F1Area << " Area3: " << F2Area << " loc: " << loc_index << " loc2: " << loc_index2 << " loc3: " << loc_index3 << endl;
                //cout << "BKG: DiamOriginal: " << IDiam << " Diam2: " << F1Diam << " Diam3: " << F2Diam << " loc: " << loc_index << " loc2: " << loc_index2 << " loc3: " << loc_index3 << endl;

                //Avoid degenerate breakage, skip if it happens
                if (isnan(IDiam))
                {
                    //cout << "WARNING: BreakBinsSkipped due to NaN" << endl;
                    //return;
                    continue;
                }


                //Check if grain changed bin as well, if it did then substract all mass from initial bin and move it to final bin
                if (loc_index == loc_index2 && loc_index == loc_index3){ //No change
                } 
                else if (loc_index == loc_index2 && loc_index != loc_index3)  //Some Loss
                {
                    BreakBin[loc_index] += -(F2Area);
                    BreakBin[loc_index3] += (F2Area);
                }
                else if (loc_index == loc_index3 && loc_index != loc_index2) //Some Loss
                {
                    BreakBin[loc_index] += -(F1Area);
                    BreakBin[loc_index2] += (F1Area);
                }
                else{ //All lost
                    BreakBin[loc_index] += -(IArea);  //or (F1Area + F2Area), but what if bkg is not mass conservative, thus lose IArea just in case
                    BreakBin[loc_index2] += (F1Area);
                    BreakBin[loc_index3] += (F2Area);
                }



                /*
                 * Add/Remove grains from world
                 */
                //remove fractured grain
                //cout << "Number of grains BEFORE: " << _ngrains << endl;
                _grains.erase(_grains.begin()+g);
                _ngrains--;
                //Add new grains if not problematic  _2 is for a smaller resolution of grains
                if (_fracProps.ApplyFilter_2(g2)){
                    //cout << "Insert 2" << endl;
                    _grains.insert(_grains.begin(),g2);
                    _ngrains++;
                }
                if (_fracProps.ApplyFilter_2(g3)) {
                    //cout << "Insert 3" << endl;
                    _grains.insert(_grains.begin(), g3);
                    _ngrains++;
                }
                //cout << "Number of grains NOW: " << _ngrains << endl;

//              _grains.insert(_grains.begin(),g2);
//              _grains.insert(_grains.begin(),g3);
//
//              _ngrains++;
//              _ngrains++;

                _maxId++;_maxId++;
                //cout << "Reset Global State" << endl;
                _globalGrainState.reset();
                _globalGrainState.resize(_ngrains, _nwalls);
                computeWorldState();
                
                //NEW CHANGE OFF
                //cout << "Breakage Real End for grain: " << g << endl;
                return; //Sep 30, 2022 Change for Prob. Break Try OFF //TURN OFF FOR ALL GRAINS TO BE ANALYZED //ON FOR ONLY ONE GRAIN PER STEP
                //NEW CHANGE OFF
                
                //TWO OPTIONS: THICKNESS ONLY AND ALL GRAINS BREAK BUT FIX BUG ORRRR ONLY ONE BREAK AND TWO CRITERIA ?? A or B??
                //cout << "Breakage Real End for grain: " << g << endl;
            }//if frac flag
        }//grains loop
        return; //CHECK!!!  //Sep 30, 2022 Change for Prob. Break Try  ON
    }
    
    static bool pointsort(const PointInfo & a,const PointInfo & b){
        if (a.point(0) >= 0. && b.point(0) < 0.){
            return true;
        }
        if (a.point(0) < 0. && b.point(0) >=0.) {
            return false;
        }

        double det = a.point(0)*b.point(1)-b.point(0)*a.point(1);

        return det<0;
    }


    void addGravityForce(const double & grav) {
        for (size_t grainid = 0; grainid < _grains.size(); grainid++) {
            _globalGrainState._grainForces[grainid](1) -= grav*_grains[grainid].getMass();
        }
    }




  // Checks every grain if it is flagged for fracture, then cuts grain with a line connecting the 2 highest contact force points
  double randf(double low,double high){
      return (rand()/(double)(RAND_MAX))*abs(low-high)+low;
  }


  void FracRoutineModif() {

    //double rhoice = 0.910e-6; //or 0.910

    vector<size_t> iter(_ngrains);
    for (size_t g = 0; g < _ngrains; g++){
      iter[g] = g;
    }


    //To avoid conflic no shuffling for now
    //random_shuffle(iter.begin(),iter.end());
    
    if ( _ngrains <2 ) {
      //Choose a random time to break that random grain
      // if (_stepup>10000 && _stepup % 400 == 0){
      //   double trigger = randf(-1,1);
      //   cout << "trigger: " << trigger <<endl;
      //   if (trigger>0.6){
      //      //Choose Random Grain to Break
      //     size_t rg = randf(1,_ngrains-1);
      //     cout << "rg: " << rg <<endl;
      //     if (_grains[rg].getFracLoc()==0){
      //         _grains[rg].changeFracFlag(false); 
      //     }
      //     else{
      //       _grains[rg].changeFracFlag(true);
      //       _grains[rg].changeFracLoc(_grains[rg].getFracLoc());
      //     }   //To be modified by a Random Factor Later instead of simply a midpoint
      //   }
      // }
    }
    else if (_ngrains>=2){      
      if ( _stepup>4000 && _stepup % 500 == 0){ //400 vs 250
        double trigger = randf(-1,1);
        cout << "trigger: " << trigger <<endl;
        if (trigger>-0.8){ //0 versus 0.5  //-0.3
           //Choose Random Grain to Break
           size_t rg = randf(0,_ngrains-1);
           cout << "rg: " << rg <<endl;
          if (_grains[rg].getFracLoc()==0){  //Avoid change if grain is too small
            _grains[rg].changeFracFlag(false); 
          }
          else{
            _grains[rg].changeFracFlag(true);
            _grains[rg].changeFracLoc(_grains[rg].getFracLoc());
          }   //To be modified by a Random Factor Later instead of simply a midpoint
        }
      }
    }
    else{
    }  

    for (size_t g1 = 0; g1 < _ngrains; g1++){
      size_t g = iter[g1];
      //check if grain should fracture
      //double Volume = _grains[g].getMass()/rhoice;


      //HERE WE WILL NOT USE STRESS BUT A GRAIN BREAK POINT COORDINATE FOR VERTICAL PLANE FRACTURE
      //Vector3d grainStress = _globalGrainState.grainStress[g]/Volume;

      //check grain yield stress
      //double Sig1 = (grainStress(0)+grainStress(1))/2. + sqrt(pow((grainStress(0)-grainStress(1))/2.,2)+pow(grainStress(2),2));
      //double Sig2 = (grainStress(0)+grainStress(1))/2. - sqrt(pow((grainStress(0)-grainStress(1))/2.,2)+pow(grainStress(2),2));

      // Von Mises criterion
//      double VonMises = sqrt(grainStress(0)*grainStress(0)+grainStress(1)*grainStress(1)
//          -grainStress(0)*grainStress(1)+3.*grainStress(2)*grainStress(2));
//      _grains[g].changeCritStress(VonMises);

      // Tresca criterion
//      double Tresca  = fabs(Sig1-Sig2);
//      _grains[g].changeCritStress(Tresca);

      // maximum principal stress criterion
      //_grains[g].changeCritStress(Sig1);

      //THIS IS RIGHT NOW CHANGED ON THE TEMPERATURE MODIF SUBROUTINE
      //if (Sig1 > _grains[g].getYield() || fabs(Sig2) > 100.*_grains[g].getYield()){
      //  _grains[g].changeFracFlag(true);
      //}

      //if (_grains[g].getFracFlag()
      //    && _globalGrainState._grainCoord[g].size() >= 2) 
      //Temporaly and arbitrarily break grains
      //if ( (_stepup==10000 && g == _ngrains-1) || (_stepup==16000 && g == _ngrains-1) || (_stepup == 20000 && g == _ngrains-1) ) {   //} || (_stepup == 32000  && i== 3) ){
      

      // if (_stepup>10000 && _stepup % 400 == 0){
      //   double trigger = randf(-1,1);
      //   cout << "trigger: " << trigger <<endl;
      //   if (trigger>0.6){
      

      //if ( (_stepup==10001 && g == _ngrains-1) || (_stepup==16001 && g == _ngrains-1) || (_stepup == 20001 && g == _ngrains-1) || (_stepup == 24001  &&  g == _ngrains-1) || (_stepup == 28001  && g == _ngrains-1) || (_stepup == 30001  && g == _ngrains-1)  || (_stepup == 35001  && g == 1) || (_stepup == 35101  && g == 2) ) { //} || (_stepup == 32000  && i== 3) ){
      // if ( (_stepup==6001 && g == _ngrains-1) || (_stepup==8001 && g == _ngrains-1) || (_stepup == 9001 && g == _ngrains-1)  ) { //} || (_stepup == 32000  && i== 3) ){
      //   _grains[g].changeFracFlag(true);
      //   _grains[g].changeFracLoc(_grains[g].getFracLoc());   //To be modified by a Random Factor Later instead of simply a midpoint
      // }

        if (_grains[g].getFracFlag() == true ){
            cout << "We really have fracture for grain: " << g << endl;
  //      cout << "coord # " << _globalGrainState._grainCoord[g].size() << endl;
          /*
           * Get Grain 1 Info
           */
          
          //Vector2d pF1;
          //Vector2d pF2;
          
          //Get Source Grain Properties
          Vector2d g1cm = _grains[g].getCmLset();

          vector<Vector2d> g1p0 = _grains[g].getPointList();

          size_t np1 = g1p0.size();
          size_t gid = _grains[g].getId();
          Levelset2d g1lset = _grains[g].getLset();
          vector<double> g1Lset = g1lset.getLevelset();

          Levelset2d g1lsettemp = _grains[g].getgrainTemp2D();
          vector<double> g1Lsettemp = g1lsettemp.getLevelset();
          
          size_t Xdim = g1lset.getXdim();
          size_t Ydim = g1lset.getYdim();

          vector<Vector2d>  g1ipoint(np1);
          vector<double>  g1ishear(np1);
          vector<size_t>  g1icontact(np1);

          double gdense = _grains[g].getDensity();

          


          //vector<PointInfo> g1i(np1);

          const double cos1 = cos(_grains[g].getTheta());
          const double sin1 = sin(_grains[g].getTheta());
          for (size_t i = 0; i<np1;i++){
            //must unrotate and bring to origin
  //          g1i[i].point << (g1p[i](0)-_grains[g].getPosition()(0))*cos1 + (g1p[i](1)-_grains[g].getPosition()(1))*sin1,
  //                        -(g1p[i](0)-_grains[g].getPosition()(0))*sin1 + (g1p[i](1)-_grains[g].getPosition()(1))*cos1;
            //g1ipoint = g1p0[i];
            //
            //g1ishear = _grains[g].getNodeShears()[i];
            //g1icontact = _grains[g].getNodeContact()[i];
          }

          /*
           * Find Splitter Shape
           */
          //TwoHighestForcesLocs(pF1,pF2,g,g1i);

          //double bigV = M_PI*_grains[g].getRadius()*_grains[g].getRadius();
          //double g1V  = _grains[g].getMass()/rhoice;
          
          //if (bigV/g1V > _fracProps.getRoMratio() || pF2 == pF1){
          //  ClosestPoint(pF1,g1p0);
          //  pF2 = Vector2d(0.,0.);
          //}



          /*
           * Create Splitter
           */


          //get points that define split line in the level set reference frame
          //Vector2d splitpt1 = pF1 + g1cm;
          //Vector2d splitpt2 = pF2 + g1cm;
         
          //REDEFINE WITH NEW DATA
          //Vector2d splitpt1
          //Vector2d splitpt2

          //make splitter surface points
          //vector<Vector2d> splpts0 = _fracProps.MakeSplitter(_grains[g],splitpt1,splitpt2);
          //make splitter level set
          //vector<double> splitset = g1lset.makeFracSurfaceLine(splitpt1,splitpt2);
          //output splitset if desired
          // _fracProps.SaveSplitter(splpts0,splitset,Xdim,Ydim);



          size_t pointbreak = _grains[g].getFracLoc();





          /*
           * Make New Grains (using Lset and Surface Points)
           */


          //make new Lsets
          //vector<double> g2Lset(splitset.size());
          //vector<double> g3Lset(splitset.size());
          //for(size_t i=0; i<splitset.size();i++){
          //  g2Lset[i] = max(splitset[i],g1Lset[i]);
          //  g3Lset[i] = max(splitset[i]*-1.,g1Lset[i]);
          //}

          //reinit due to imperfect set ops
          //g2Lset = reinit(g2Lset); //need to build in Levelset2d.h
          //g3Lset = reinit(g3Lset);

          
          // cout<<"LSET TO SPLIT"<<endl; 
          // for (size_t j = 0; j < g1lset.getYdim(); j++) {
          //     for (size_t i = 0; i < g1lset.getXdim(); i++) {
          //          cout<<"Value LSETSPLIT at i: "<< i <<" j: " << j <<" is: " << g1Lset[(j*g1lset.getXdim())+i] << endl;
          //     }
          // }

          //Declare new level set vectors
          vector<double> g2Lset = g1Lset;
          vector<double> g3Lset = g1Lset;

          //Temperature BCs of this new level sets will be updated in ComputeWorld()/TemperaModif() in which GeoLS define B.C. for Temperature
          vector<double> g2tempV = g1Lsettemp;
          vector<double> g3tempV = g1Lsettemp;


          //LSET New Function for New Level Set
          // cout<<"LevelSetSplit"<<endl;
          // cout<<"pointbreak: "<< pointbreak<< endl;
          // cout<<"Length: "<<_grains[g].getLset().getXdim()<<endl;
          // _fracProps.VertSplitLSETS(_grains[g],pointbreak,g2Lset,g3Lset);
          // WITH TEMPERATURE
          _fracProps.VertSplitLSETS(_grains[g],pointbreak,g2tempV,g3tempV);
          _fracProps.VertSplitLSETSG(_grains[g],pointbreak,g2Lset,g3Lset);

          //New level sets
          Levelset2d g2lset(g2Lset, Xdim, Ydim);
          Levelset2d g3lset(g3Lset, Xdim, Ydim);

            Levelset2d g2temp(g2tempV, Xdim, Ydim);
            Levelset2d g3temp(g3tempV, Xdim, Ydim);



          /*
           * Find Physical Properties
           */
          double mass2;
          Vector2d g2cm;
          double I2;
          cout << "frac density" << gdense << endl;
          _fracProps.findGrainProps(g2temp, mass2, g2cm, I2, gdense);

          double mass3;
          Vector2d g3cm;
          double I3;
          _fracProps.findGrainProps(g3temp, mass3, g3cm, I3, gdense);
          /*
           * Position grains correctly, removing remnant
           -//s of the grain 1 reference frame
           *
           */
          cout << "New mass for grain 1 gen is : " <<  mass2 <<endl;
          cout << "New mass for grain 2 gen is : " <<  mass3 <<endl;
          //mass2 = mass2*0.91e-6;
          //mass3 = mass3*0.91e-6;


          //POINTS!!!!!
          //vector<PointInfo> g2i, g3i;
          vector <Vector2d> g2i(_grains[g].getnpoints()); 
          vector <Vector2d> g3i(_grains[g].getnpoints()); 

          //decides which splitter points stay
          //Input: splpts0, g2lset, g1cm;  Output: g2i, g3i
          //_fracProps.KeepSplitterPoints(splpts0,g1lset,g1cm,g2i,g3i);
          //decides which previous grain points stay
          //Input: splitpt1, splitpt2, g1i, g1cm;  Output: g2i, g3i
          //_fracProps.KeepGrainPoints(splitpt1,splitpt2,g1i,g1cm,g2i,g3i);

          //make sure potential grains has a reasonable amount of surface points
          //if (g2i.size()<5 || g3i.size()<5 ){
          //  return; //to ensure no NaNs in the crumbles
          //}


          //LSET New Function to Create New Point Sets from LSET and Sorted
          //cout<<"Point Generation"<<endl;
          bool fail_ind2 = false;
          bool fail_ind3 = false;
          cout << "New Point Generation Frac" << endl;
          _fracProps.PointsGen(g2temp,g2i, _grains[g].getnpoints(), g2cm, fail_ind2 );
          _fracProps.PointsGen(g3temp,g3i, _grains[g].getnpoints(), g3cm, fail_ind3 );

          if (fail_ind2 || fail_ind3)
          {
            cout <<"WARNING: Wrong Point INTERPOLATION contact problems may occur in fracture" << endl;
            return; //Leave points the same to avoid mess and don't break
          }

          
           //cout<<"Point Generation End"<<endl;

          //save number of points for grains
          size_t np2 = g2i.size();
          size_t np3 = g3i.size();



          double cosine1 = cos1;
          double sine1 = sin1;


          Vector2d position2;
          //position2 = _grains[g].getPosition()-Vector2d((g1cm(0)-g2cm(0))*cos1 - (g1cm(1)-g2cm(1))*sin1,
          //                                          (g1cm(0)-g2cm(0))*sin1 + (g1cm(1)-g2cm(1))*cos1);
          position2 = _grains[g].getPosition()-Vector2d((g1cm(0)-g2cm(0))*cosine1 - (g1cm(1)-g2cm(1))*sine1,
                                                    (g1cm(0)-g2cm(0))*sine1 + (g1cm(1)-g2cm(1))*cosine1);


          double maxR2 = 0;
          //vector<Vector2d> g2p0, g3p0;

          //shift points to correct position
          //cout<<"Radius Loop 2"<<endl;
          //cout<<"FOR DEBUG 2 ***"<<endl;
          //cout<<"g1cm x: "<< g1cm(0) << " , y: " << g1cm(1) <<  endl;
          //cout<<"g2cm x: "<< g2cm(0) << " , y: " << g2cm(1) <<  endl;
          //cout<<"Origin Grain Position X: "<< _grains[g].getPosition()(0) << " , y: " << _grains[g].getPosition()(1) << endl;
          //cout<<"Position 2 X: "<< position2(0) << " , y: " << position2(1) << endl;
          for (size_t i = 0; i<np2; i++){

            //g2i[i] += -g1cm+g2cm-_grains[g].getPosition();

            
            //g2i[i] += -(position2-g1cm)-_grains[g].getPosition();   //ADJUST A BIT
            g2i[i] += -g2cm;//-(g2cm-g1cm)-_grains[g].getPosition();

            if (g2i[i].norm() > maxR2){
              maxR2 = g2i[i].norm();
            }
          }


          //POINT REORDERING (need PointInfo)
          // reorder points to counterclockwise. ///!!!!!!!!!!!
          // sort(g2i.begin(), g2i.end(),pointsort);

          // for (size_t i=0;i<np2;i++){
          //   g2p0.push_back(g2i[i].point);
          // }




          //Vector2d position3 = _grains[g].getPosition()-Vector2d((g1cm(0)-g3cm(0))*cos1 - (g1cm(1)-g3cm(1))*sin1,
          //                                                        (g1cm(0)-g3cm(0))*sin1 + (g1cm(1)-g3cm(1))*cos1);
          Vector2d position3 = _grains[g].getPosition()-Vector2d((g1cm(0)-g3cm(0))*cosine1 - (g1cm(1)-g3cm(1))*sine1,
                                                                  (g1cm(0)-g3cm(0))*sine1 + (g1cm(1)-g3cm(1))*cosine1);
          double maxR3 = 0;
          //shift points to correct position
          //cout<<"Radius Loop 3"<<endl;

          //cout<<"FOR DEBUG 3 ***"<<endl;
          //cout<<"g1cm x: "<< g1cm(0) << " , y: " << g1cm(1) <<  endl;
          //cout<<"g3cm x: "<< g3cm(0) << " , y: " << g3cm(1) <<  endl;
          //cout<<"Origin Grain Position X: "<< _grains[g].getPosition()(0) << " , y: " << _grains[g].getPosition()(1) << endl;
          //cout<<"Position 3 X: "<< position3(0) << " , y: " << position3(1) << endl;
          for (size_t i = 0; i<np3; i++){
            //g3i[i] += g1cm-g3cm-_grains[g].getPosition();
            
            g3i[i] += -g3cm;//-(g3cm-g1cm)-_grains[g].getPosition();

            //Find bbox radius
            if (g3i[i].norm() > maxR3){
              maxR3 = g3i[i].norm();
            }
          }

          //POINT REORDERING (need PointInfo)
          // sort indexes based on comparing values in g3p0 b
          //sort(g3i.begin(), g3i.end(),pointsort);

          //for (size_t i=0;i<np3;i++){
          //  g3p0.push_back(g3i[i].point);
          //}




          //check thinness of potential grains in case they could slide into other grains
          Vector2d g2min,g3min; //closest surface point to centroid
          double g2minR,g3minR; //minimum distance from centoid to surface
          ClosestPoint(g2min,g2i);
          ClosestPoint(g3min,g3i);
          g2minR = g2min.norm();
          g3minR = g3min.norm();

          //Check distance tolerance
          if (g2minR < .005|| g3minR < .005) {
            return;
          }


  //        if (M_PI*maxR2*maxR2/(mass2) > _fracProps.getMaxRoM() || M_PI*maxR3*maxR3/(mass3) > _fracProps.getMaxRoM() ){
  //          return; //to ensure no NaNs in the crumbles
  //        }

          /*
           * Create grains
           */

          //Grain2d g2 = Grain2d(mass2, position2, _grains[g].getVelocity(), I2, _grains[g].getTheta(), _grains[g].getOmega(),
          //            g2cm, g2p0, maxR2, g2lset, _grains[g].getKn(), _grains[g].getKs(), _grains[g].getMu(), _maxId+1);

          //cout<<"New grains are created"<<endl;
          bool tempflag = false;
          int temploc = 0;
          double time_fail = 1.000;
          double orig_fail = 0.000;
          Grain2d g2 = Grain2d(mass2, position2, _grains[g].getVelocity(), 
                           mass2, _grains[g].getTemperature(), _grains[g].getThickness(), _grains[g].getMthick(), _grains[g].getTemperature(), _grains[g].getThickness(), _grains[g].getUtemper(), _grains[g].getUtemper(), g2temp, g2temp, 
                           I2, _grains[g].getTheta(), _grains[g].getOmega(), g2cm, g2i, _grains[g].getnpoints(), maxR2, 
                           g2lset, g2lset, 0, _grains[g].getId()+1, _grains[g].getKn(), _grains[g].getKs(), _grains[g].getMu(), tempflag, temploc,
                           _grains[g].getgrainStress2D(), _grains[g].getgrainDamage() , _grains[g].getgrainThickness(), g2i, time_fail, orig_fail); //_grains[g].getId()+1
      


          //g2.changeDensity(_grains[g].getDensity());
          for (size_t i=0;i<np2;i++){
            //g2.getNodeShearsNonConst()[i]  = g2i[i].shear;
            //g2.getNodeContactNonConst()[i] = g2i[i].contact;
          }

          //This is to serve the new grain as a file
          //_fracProps.SaveGrain(g2,g2p0);


          //Grain2d g3 = Grain2d(mass3, position3, _grains[g].getVelocity(), I3, _grains[g].getTheta(), _grains[g].getOmega(),
          //            g3cm, g3p0, maxR3, g3lset, _grains[g].getKn(), _grains[g].getKs(), _grains[g].getMu(), _maxId+2);

          Grain2d g3 = Grain2d(mass3, position3, _grains[g].getVelocity(), 
                           mass3, _grains[g].getTemperature(), _grains[g].getThickness(), _grains[g].getMthick(), _grains[g].getTemperature(), _grains[g].getThickness(), _grains[g].getUtemper(), _grains[g].getUtemper(), g3temp, g3temp, 
                           I3, _grains[g].getTheta(), _grains[g].getOmega(), g3cm, g3i, _grains[g].getnpoints(), maxR3, 
                           g3lset, g3lset, 1, _grains[g].getId()+2, _grains[g].getKn(), _grains[g].getKs(), _grains[g].getMu(), tempflag, temploc,
                           _grains[g].getgrainStress2D(), _grains[g].getgrainDamage() , _grains[g].getgrainThickness(), g3i, time_fail, orig_fail);  //_grains[g].getId()+2
          


          //g3.changeDensity(_grains[g].getDensity());
  //        for (size_t i=0;i<np3;i++){
  //          g3.getNodeShearsNonConst()[i]  = g3i[i].shear;
  //          g3.getNodeContactNonConst()[i] = g3i[i].contact;
  //        }

          //send shears to "master" grains
  //        double dist = 50;
  ////        cout << "g3isize " << g3i.size() << " " << np3 << endl;
  //        vector<Vector2d> g3_moved_pts = g3.getPointList();
          for (size_t i=0;i<np3;i++){
  //          if (g3i[i].contact > _maxId){ //contact with wall; shear stays

              //g3.getNodeShearsNonConst()[i] = g3i[i].shear;
              //g3.getNodeContactNonConst()[i] = g3i[i].contact;


  ////            cout << "\nwall " << g3i[i].contact << "\n"<<endl;
  //          }
  //          else if (g3i[i].contact>0){ //contact with grain; shear goes to other grain
  //            //loop over each possible point and search for closest one
  //            size_t pt_id;
  //
  //            vector<Vector2d> target_grain_pts = _grains[FindGrainFromId(g3i[i].contact)].getPointList();
  //            for (size_t j = 0;j<target_grain_pts.size();j++){
  //              if ((g3_moved_pts[i]-target_grain_pts[j]).norm() < dist){
  //                dist = (g3_moved_pts[i]-target_grain_pts[j]).norm();
  //                pt_id = j;
  //              }
  //            }
  ////            cout << "\ndistance " << dist << "\n"<<endl;
  //            _grains[FindGrainFromId(g3i[i].contact)].getNodeContactNonConst()[pt_id] = g3.getId();
  //            _grains[FindGrainFromId(g3i[i].contact)].getNodeShearsNonConst()[pt_id] = g3i[i].shear;
  //          }
          }

          //_fracProps.SaveGrain(g3,g3p0);

          //double yield = _grains[g].getYield();
          //g3.changeYield(yield*pow(_grains[g].getDiameter()/g3.getDiameter(),3./3.));
          //g2.changeYield(yield*pow(_grains[g].getDiameter()/g2.getDiameter(),3./3.));


          //find shears where g was slave and move to master
          //cout << "Master" << endl; 
          for (size_t i = 0; i<g;i++){
            //check if grain was a master grain on g
            vector<size_t> master_contacts = _grains[i].getNodeContact();
            bool master = std::binary_search(master_contacts.begin(),master_contacts.end(),gid); //true if master grain was contacting fractured grain
            if (master) {
  //            cout << "Position " << g2.getPosition().transpose() << " " << g3.getPosition().transpose() << " " << _grains[i].getPosition().transpose() << endl;
              //move shears to g2,g3
              //find node
              vector<size_t> nodes;
              size_t n=0;
              for (size_t j = 0; j<master_contacts.size();j++){
                if (master_contacts[j] == gid){
                  nodes.push_back(j);
                  n++;
                }
              }

              //get point location
              vector<Vector2d> master_pts = _grains[i].getPointList();

              vector<Vector2d> g2_pts = g2.getPointList();
              vector<Vector2d> g3_pts = g3.getPointList();

              for(size_t j=0;j<n;j++){
                //find closest pt on either g2 or g3
                Vector2d master_pt = master_pts[nodes[j]];
                size_t kmin = 0;
                double dist = 50.;
                for (size_t k = 0; k<np2;k++){
                  double dist_k = (master_pt-g2_pts[k]).norm();
  //                cout << "dist_k " << dist_k << endl;
                  if (dist_k < dist){
                    dist = dist_k;
                    kmin = k;
                  }
                }
                for (size_t k = np2; k<np2+np3;k++){
                  double dist_k = (master_pt-g3_pts[k-np2]).norm();
                  if (dist_k < dist){
                    dist = dist_k;
                    kmin = k;
                  }
                }
  //              cout << "dist " << dist << endl;
                if (dist < 1) {
                  //send shears to fragments
                  if (kmin>=np2) {
                    kmin -= np2; //on grain 3
                    g3.getNodeShearsNonConst()[kmin]  = _grains[i].getNodeShears()[nodes[j]];
                    g3.getNodeContactNonConst()[kmin] = _grains[i].getNodeContact()[nodes[j]];
                  } else { //on grain 2
                    g2.getNodeShearsNonConst()[kmin]  = _grains[i].getNodeShears()[nodes[j]];
                    g2.getNodeContactNonConst()[kmin] = _grains[i].getNodeContact()[nodes[j]];
                  }
                  _grains[i].getNodeShearsNonConst()[nodes[j]] = 0;
                  _grains[i].getNodeContactNonConst()[nodes[j]] = 0;
                }
              }
            }
          }

          /*
           * Add/Remove grains from world
           */
          //remove fractured grain
          //cout << "Delete original" << endl; 
          _grains.erase(_grains.begin()+g);
          _ngrains--;
          //cout << "Delete original end" << endl; 
          
          //Add new grains if not problematic
          //cout << "Add new" << endl; 
          //Our ordering goes left right

          for (size_t i = 0; i < _ngrains; i++){
            if (i<g){
                 _grains[i].changeId(_grains[i].getId()+2);
              }
              else{
                 _grains[i].changeId(_grains[i].getId()+1); 
              }
          }


          if (_fracProps.ApplyFilter(g3)) {
            _grains.insert(_grains.begin(),g3);
            //_grains.push_back(g2);    
            _ngrains++;
          }
          if (_fracProps.ApplyFilter(g2)){
            _grains.insert(_grains.begin(),g2);
            //_grains.push_back(g3);  
            _ngrains++;
          }
          //cout << "Add new end" << endl; 

  //        _grains.insert(_grains.begin(),g2);
  //        _grains.insert(_grains.begin(),g3);
  //
  //        _ngrains++;
  //        _ngrains++;

          //_maxId++;_maxId++;

          //_globalGrainState._grainForces[0] = Vector2d(0.,0.);
          //_globalGrainState._grainMoments[0] = 0.0;
          //_globalGrainState._wallForces[0] = Vector2d(0.,0.);
          //_globalGrainState._grainForces[_ngrains] = _globalGrainState._grainForces[0];
          //_globalGrainState._grainMoments[_ngrains] = _globalGrainState._grainMoments[0];
          //_globalGrainState._wallForces[_ngrains] = _globalGrainState._wallForces[0];

          _globalGrainState.reset();
          _globalGrainState.resize(_ngrains, _nwalls);

         cout << "Force Grain 2" << endl; 

          cout << "begin compute world again" << endl;
          computeWorldState();
          cout << "end compute world again" << endl;
          //return; 
        }//if frac flag
    }//grains loop
    return; 
  }
  
  
  //Information Output for Processing
  void outputGrainData (size_t & idx, string & folderLoc, size_t & now_time){
      
    size_t try_index = idx;

    cout << "Original Position for Testing for grain: " << try_index << endl;
    cout << "X: "  << _grains[try_index].getPosition()[0] << " Y: " << _grains[try_index].getPosition()[1] << endl;

    Vector2d velOrig(_grains[try_index].getVelocity()[0] , _grains[try_index].getVelocity()[1]);
    double tempOrig  = _grains[try_index].getTemperature();
    double thickOrig = _grains[try_index].getThickness();
    double thetaOrig = _grains[try_index].getTheta();
    double omegaOrig = _grains[try_index].getOmega();

    //Print temporal point file
    string outDir = folderLoc;
    string fname = outDir + "Grain_Test_Points_step_" + std::to_string(now_time) + "_g_" + std::to_string(try_index) + ".dat";           
    string fname2 = outDir + "Test_CM_step_" + std::to_string(now_time) + "_g_" + std::to_string(try_index) + ".dat";
    string fname3 = outDir + "Test_POS_step_" + std::to_string(now_time) + "_g_" + std::to_string(try_index) + ".dat";
    string fname4 = outDir + "Geo_LS_step_" + std::to_string(now_time) + "_g_" + std::to_string(try_index) + ".dat";
    string fname5 = outDir + "Thick_LS_step_" + std::to_string(now_time) + "_g_" + std::to_string(try_index) + ".dat";
    
    FILE * points_export = fopen(fname.c_str(), "w");
    
    for (size_t i = 0; i < _grains[try_index].getPointList().size(); i++){ 
        fprintf( points_export,  "%4.8f %4.8f %4.8f\n", _grains[try_index].getPointList()[i](0), _grains[try_index].getPointList()[i](1), 0.0000 );
    }   
    fclose(points_export);

    FILE * CM_export = fopen(fname2.c_str(), "w");
    fprintf(CM_export,  "%4.8f\n", _grains[try_index].getRadius() );
    fprintf(CM_export,  "%4.8f\n", _grains[try_index].getCmLset()(0) );
    fprintf(CM_export,  "%4.8f\n", _grains[try_index].getCmLset()(1) );
    fclose(CM_export);

    FILE * POS_export = fopen(fname3.c_str(), "w");
    fprintf(POS_export,  "%4.8f %4.8f\n", _grains[try_index].getPosition()(0), _grains[try_index].getPosition()(1) );
    fclose(POS_export);
    
    //Level Set Export
    size_t xxdim =  _grains[try_index].getLset().getXdim();  // To get Heat LSET
    size_t yydim =  _grains[try_index].getLset().getYdim();
    double xgrid =  _grains[try_index].getLset().getXdim();
    double ygrid =  _grains[try_index].getLset().getYdim();
    
    FILE * geo_LS = fopen(fname4.c_str(), "w");

    fprintf(geo_LS, "%4.8f \n", xgrid);
    fprintf(geo_LS, "%4.8f \n", ygrid);

    for (size_t jj = 0; jj < yydim; jj++){
        for (size_t ii = 0; ii < xxdim; ii++){ 
            fprintf(geo_LS, "%4.8f \n", _grains[try_index].getLset().getLevelset()[jj*xxdim + ii]); //Geo LS
        }     
    }
    fclose(geo_LS);
    
    FILE * thick_LS = fopen(fname5.c_str(), "w");

    fprintf(thick_LS, "%4.8f \n", xgrid);
    fprintf(thick_LS, "%4.8f \n", ygrid);

    for (size_t jj = 0; jj < yydim; jj++){
        for (size_t ii = 0; ii < xxdim; ii++){ 
            fprintf(thick_LS, "%4.8f \n", _grains[try_index].getMthick().getLevelset()[jj*xxdim + ii]); //Thick LS
        }     
    }
    fclose(thick_LS);
    

    //Execute Points to gmsh, ls etc.
    // //Convert point file to new temporal morphology //This prints a temporal grainproperty0 only for point processing.
    // std::string filename = "/Users/rigobertomoncadalopez/Dropbox/Caltech_PhD/Research/Code/LSDEM2D-SeaIceDense/Points_to_LS.py";
    // std::string command = "python3 ";
    // command += filename;
    // system(command.c_str());


    // //Grains  
    // char tempfname[300];  
    // sprintf(tempfname, "./Output/SeaIce/Temporal/Temp_morphologies.dat");
    // string file_morph = tempfname;
    // sprintf(tempfname, "./Output/SeaIce/Temporal/Temp_positions.dat"); //Should come from point interpolation centroid and respective position file 
    // string file_pos = tempfname;
    // sprintf(tempfname, "../Output/SeaIce/Temporal/Temp_velocities.dat");  //RE-UPDATE
    // string file_vel = tempfname;
    // sprintf(tempfname, "./Output/SeaIce/Temporal/");
    // string morph_dir = tempfname;
    // //Get temperature-related quantitities
    // sprintf(tempfname, "./Output/SeaIce/Temporal/Temp_temper.dat"); //RE-UPDATE
    // string file_temper = tempfname;
    // sprintf(tempfname, "./Output/SeaIce/Temporal/Temp_thick.dat"); //RE-UPDATE
    // string file_thick = tempfname;


    // //Convert temporal morphology to fractured new grain
    // vector<Grain2d> grainsTemp = generateGrainsFromFiles(file_morph, morph_dir, file_pos, file_vel, file_temper, file_thick);
    // _grains[try_index] = grainsTemp[0]; //Update grain except for:

    // //Recover Velocity, Thickness and Temperature, Theta Rotation and Omega Angular Velocity
    // _grains[try_index].changeVelocity(velOrig);
    // _grains[try_index].changeTemper(tempOrig);
    // _grains[try_index].changeThickness(thickOrig);
    // _grains[try_index].changeTheta(thetaOrig);
    // _grains[try_index].changeOmega(omegaOrig);

    // cout << "Properties New Grain: " << try_index << endl;
    // cout << "Position X: " << _grains[try_index].getPosition()[0] << " Y: " << _grains[try_index].getPosition()[1] << endl;
    // cout << "Velocity X: " << _grains[try_index].getVelocity()[0] << " Y: " << _grains[try_index].getVelocity()[1] << endl;
    // cout << "CM X: " << _grains[try_index].getCmLset()[0] << " Y: " << _grains[try_index].getCmLset()[1] << endl;
    // cout << "Radius: " << _grains[try_index].getRadius() << endl;


    // for (size_t i = 0; i < _ngrains; i++){ 
    //   if (i == 0){
        
    //   }
    //   else
    //   {

    //   }
    // } 
      
      
    return;
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
  MatrixXd LS_Soft(MatrixXd & Cls_Tag0, const size_t & maxttries) //const size_t & SZY, const size_t & SZX)
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
          //weird_ind = 0; 
       } 

    }

    return Cls_Tag;      

  }


  void AddGrain ( Grain2d & grain ){
    //grain.changeId(_maxId+1);
    //_maxId++;
    _grains.push_back(grain);
    _ngrains++;
    _globalGrainState.resize(_ngrains,_nwalls);
    _globalGrainState.reset();
  }

  void RemoveGrain(const size_t & grainid){
        size_t grain = FindGrainFromId(grainid);
        _grains[grain] = _grains.back();
        _grains.pop_back();
        _ngrains--;
        _globalGrainState.resize(_ngrains, _nwalls);
  }


  void RemoveGrain2(const size_t & grainid){
    size_t grain = FindGrainFromId(grainid);
    //_grains[grain] = _grains.back();
    //_grains.pop_back();
    //_ngrains--;
    //_globalGrainState.resize(_ngrains,_nwalls);

    _grains.erase(_grains.begin()+grain);
    _ngrains--;
    _globalGrainState.reset();
    _globalGrainState.resize(_ngrains, _nwalls);
  }

  size_t FindGrainFromId(const size_t & id) {//need to fix
    for (size_t g = 0; g<_ngrains; g++){
      size_t gid = _grains[g].getId();
      if (gid == id){
        return g;
      }
    }
    cout << "ERROR" << endl;
    return 0;
  }

   size_t FindWallFromId(const size_t & id) {//need to fix
        for (size_t w = 0; w<_walls.size(); w++){
            size_t wid = _walls[w].getId();
            if (wid == id){
                return w;
            }
        }
        cout << "ERROR" << endl;
        return 0;
    }


    size_t FindGrainWallFromId(const size_t & id) {//need to fix
        for (size_t w = 0; w<_grainsWall.size(); w++){
            size_t wid = _grainsWall[w].getId();
            if (wid == id){
                return w;
            }
        }
        cout << "ERROR" << endl;
        return 0;
    }

  void applyBodyForce(Vector2d bodyForce) {
    for (size_t i = 0; i < _ngrains; i++) {
      _globalGrainState._grainForces[i] += bodyForce;
    }
  }

  // Need to prescribe negative acceleration for gravity
  void applyAcceleration(Vector2d acceleration) {
      //double mfactor=0.;  //Temper

    for (size_t i = 0; i < _ngrains; i++) {
            //_grains[i].changeMassf(mfactor); //Temper
      _globalGrainState._grainForces[i] += _grains[i].getMass()*acceleration; //*mfactor;
    }
  }




  // take a timestep for each grain based off of the world's _globalGrainState
  void grainTimestep() {

    int numprocessors, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    #pragma omp parallel for default(none) schedule(static,1) //num_threads(12)
    for (size_t i = 0; i < _ngrains; i++) {
      _grains[i].takeTimestep(_globalGrainState._grainForces[i], _globalGrainState._grainMoments[i], _gDamping, _dt);
      if (_grains[i].getPosition()(0) < 0 ) {
        _grains[i].moveGrain(Vector2d(_offset(0), 0.) );
      }
      else if (_grains[i].getPosition()(0) > _offset(0)) {
        _grains[i].moveGrain(Vector2d(-_offset(0), 0.) );
      }
      if (_grains[i].getPosition()(1) < 0 ) {
        _grains[i].moveGrain(Vector2d(0,_offset(1)) );
      }
      else if (_grains[i].getPosition()(1) > _offset(1)) {
        _grains[i].moveGrain(Vector2d(0,-_offset(1)) );
      }

    }
  }
  
  // take a timestep for each grain based off of the world's _globalGrainState
  void grainTimestepVelocity() {

    int numprocessors, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    #pragma omp parallel for default(none) schedule(static,1) //num_threads(12)
    for (size_t i = 0; i < _ngrains; i++) {
      _grains[i].takeTimestepVelocity(_flowspeed, _flowangle, _globalGrainState._grainMoments[i], _gDamping, _dt);
      if (_grains[i].getPosition()(0) < 0 ) {
        _grains[i].moveGrain(Vector2d(_offset(0), 0.) );
      }
      else if (_grains[i].getPosition()(0) > _offset(0)) {
        _grains[i].moveGrain(Vector2d(-_offset(0), 0.) );
      }
      if (_grains[i].getPosition()(1) < 0 ) {
        _grains[i].moveGrain(Vector2d(0,_offset(1)) );
      }
      else if (_grains[i].getPosition()(1) > _offset(1)) {
        _grains[i].moveGrain(Vector2d(0,-_offset(1)) );
      }

    }
  }
  
  //Using fixed node list, make a binary list for timestepping
  //Only works if _ngrains is stable and has a constant size (TODO) 
  //If floes are removed we need to also shift movable_nodes
  void grainTimestepMov() {
    int numprocessors, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    #pragma omp parallel for default(none) schedule(static,1) //num_threads(12)
    for (size_t i = 0; i < _ngrains; i++) {
        if (_movable_floes[i] == 1){
          _grains[i].takeTimestep(_globalGrainState._grainForces[i], _globalGrainState._grainMoments[i], _gDamping, _dt);
          if (_grains[i].getPosition()(0) < 0 ) {
            _grains[i].moveGrain(Vector2d(_offset(0), 0.) );
          }
          else if (_grains[i].getPosition()(0) > _offset(0)) {
            _grains[i].moveGrain(Vector2d(-_offset(0), 0.) );
          }
          if (_grains[i].getPosition()(1) < 0 ) {
            _grains[i].moveGrain(Vector2d(0,_offset(1)) );
          }
          else if (_grains[i].getPosition()(1) > _offset(1)) {
            _grains[i].moveGrain(Vector2d(0,-_offset(1)) );
          }
        }
        //Special for loading grains
        else if(_grains[i].getLoadGrain()) {
            _grains[i].takeTimestepLoad(_globalGrainState._grainForces[i], _globalGrainState._grainMoments[i], _gDamping, _dt, _loadvel);
          if (_grains[i].getPosition()(0) < 0 ) {
            _grains[i].moveGrain(Vector2d(_offset(0), 0.) );
          }
          else if (_grains[i].getPosition()(0) > _offset(0)) {
            _grains[i].moveGrain(Vector2d(-_offset(0), 0.) );
          }
          if (_grains[i].getPosition()(1) < 0 ) {
            _grains[i].moveGrain(Vector2d(0,_offset(1)) );
          }
          else if (_grains[i].getPosition()(1) > _offset(1)) {
            _grains[i].moveGrain(Vector2d(0,-_offset(1)) );
          }
        }
    }
  }
  
  //Load on both extremes
  void grainTimestepMovUD() {
    int numprocessors, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    #pragma omp parallel for default(none) schedule(static,1) //num_threads(12)
    for (size_t i = 0; i < _ngrains; i++) {
        if (_movable_floes[i] == 1){
          _grains[i].takeTimestep(_globalGrainState._grainForces[i], _globalGrainState._grainMoments[i], _gDamping, _dt);
          if (_grains[i].getPosition()(0) < 0 ) {
            _grains[i].moveGrain(Vector2d(_offset(0), 0.) );
          }
          else if (_grains[i].getPosition()(0) > _offset(0)) {
            _grains[i].moveGrain(Vector2d(-_offset(0), 0.) );
          }
          if (_grains[i].getPosition()(1) < 0 ) {
            _grains[i].moveGrain(Vector2d(0,_offset(1)) );
          }
          else if (_grains[i].getPosition()(1) > _offset(1)) {
            _grains[i].moveGrain(Vector2d(0,-_offset(1)) );
          }
        }
        //Special for loading grains Up
        else if(_grains[i].getLoadGrain() && _movable_floes[i] == 0) {
            _grains[i].takeTimestepLoad(_globalGrainState._grainForces[i], _globalGrainState._grainMoments[i], _gDamping, _dt, _loadvel);
          if (_grains[i].getPosition()(0) < 0 ) {
            _grains[i].moveGrain(Vector2d(_offset(0), 0.) );
          }
          else if (_grains[i].getPosition()(0) > _offset(0)) {
            _grains[i].moveGrain(Vector2d(-_offset(0), 0.) );
          }
          if (_grains[i].getPosition()(1) < 0 ) {
            _grains[i].moveGrain(Vector2d(0,_offset(1)) );
          }
          else if (_grains[i].getPosition()(1) > _offset(1)) {
            _grains[i].moveGrain(Vector2d(0,-_offset(1)) );
          }
        }
        //Special for loading grains Down (for fixed under base case)
        else if(_grains[i].getLoadGrain() == false && _movable_floes[i] == 0) {
            _grains[i].takeTimestepLoad(_globalGrainState._grainForces[i], _globalGrainState._grainMoments[i], _gDamping, _dt, _loadvelD);
          if (_grains[i].getPosition()(0) < 0 ) {
            _grains[i].moveGrain(Vector2d(_offset(0), 0.) );
          }
          else if (_grains[i].getPosition()(0) > _offset(0)) {
            _grains[i].moveGrain(Vector2d(-_offset(0), 0.) );
          }
          if (_grains[i].getPosition()(1) < 0 ) {
            _grains[i].moveGrain(Vector2d(0,_offset(1)) );
          }
          else if (_grains[i].getPosition()(1) > _offset(1)) {
            _grains[i].moveGrain(Vector2d(0,-_offset(1)) );
          }
        }       
    }
  }
  
  void change_movable_floes(vector<size_t> & new_movable_floes) {
    _movable_floes = new_movable_floes;
  }
  
  void change_loadvel(Vector2d & newloadvel) {
    _loadvel = newloadvel;
  }

  void change_loadvelD(Vector2d & newloadvel) {
    _loadvelD = newloadvel;
  }

  // move grain
  void moveGrain(const size_t & grainid, const Vector2d & amount) {
    _grains[grainid].moveGrain(amount);
  }

  // rotate grain
  void rotateGrain(const size_t & grainid, const double & amount) {
    _grains[grainid].changeRot(amount);
  }


  //************ START FINE TIME STEPPING  *******************************//
  vector<Vector2d> finefluidinteraction(vector<Vector4d> & Fine_Grains, vector<Vector2d> & Fine_Velocities) //Handle as simple point mass with only skin drag and no moment
  {
    vector<Vector2d> forceFines(Fine_Grains.size());
    Vector2d Uw = Vector2d (0.,0.); //Water velocity

    //Useful constants
    // //Initialize constants for fluid interaction with ice
    // double Cha = 1.7e-3;                         //Air-ice skin drag coefficient (unitless)
    // double Cva = 1.3e-3;                         //Air-ice body drag coefficient (unitless)
    // double Chw = 1.3e-6; //in km units  //1.3e-3;  //Water-ice skin drag coefficient (m/s)   Check units!!!
    // double Cvw = 1.3e-6; //in km units  //1.3e-3;  //Water-ice body drag coefficient (m/s)   Check units!!!
    // double rhoice =  0.91e9; // 0.910e-6;                    //Density of ice
    // double rhoair =  1.23e6; //1.23e-9;                     //Density of air If rhoice is 0.91e-6;
    // double rhowater = 1.025e9;  //1.025e-6;                  //Density of water
    // //??//double hice = 2/12;                          //Assume Unitary floe thickness, is thickness proportional to area?????
    // Vector2d Ua = Vector2d (0.,0.);                    //Constant Wind Velocity over space and time
    // Vector2d Uw = Vector2d (0.,0.);                    //Constant Water Current Velocity over space and time
        
    //VORTEX-like field  
    // double field_vel = 4.0e11;  //8.0e10 * 0.0 //8.0e1; //2.0e12    //-- 8.0e11; //2.0e12
    // //double field_vel = 1.0 * 0.00006; //8 //300 //2500.0;  //8.0 //400 bk  //TODO: Adjust to scale
    // if (_slopedir > 0)
    // {
    //     Uw << field_vel , field_vel  ; 
    // }
    // else if (_slopedir < 0)
    // {
    //     Uw << -field_vel , -field_vel  ; 
    // }
    // else
    // {
    //    cout << "WARNING SLOPE ZERO, INSPECT!!!" << endl;
    // }

    //Random direction field  
    double field_force = _flowforce; //4.0e11  //8.0e10 * 0.0 //8.0e1; //2.0e12    //-- 8.0e11; //2.0e12
    //double field_vel = 1.0 * 0.00006; //8 //300 //2500.0;  //8.0 //400 bk  //TODO: Adjust to scale

    Uw << field_force * cos(_flowangle) , field_force * sin(_flowangle) ; 


   for (size_t i = 0; i < Fine_Grains.size(); i++) 
   {
     //No surf. area adjust (drag the same for now)
     forceFines[i] = Uw;
     //double adjust_factor = 1; //100 //1000000
     //double moment_adjust_factor = 0.00000; //1000000
     //SKIN DRAG SIMPLE (Import Constant from Above)
     //fluidForcehw = rhowater*Chw*((_mass*0.91/adjust_factor)/_density)*(Uw-Fine_Velocities[i]).norm()*(Uw-Fine_Velocities[i]);   
     //fluidForceha = rhoair*Cha*((floeArea*_thick)/_density)*Ua.norm()*Ua; //use mass0/density to get volumne or in this case to get 2D area
     //fluidForcehw = rhowater*Chw*((floeArea*_thick)/_density)*(Uw-_velocity).norm()*(Uw-_velocity);  //???????????
     //fluidForceha = rhoair*Cha*((_mass*0.91/adjust_factor)/_density)*Ua.norm()*Ua; //use mass0/density to get volumne or in this case to get 2D area
     //fluidForcehw = rhowater*Chw*((_mass*0.91/1000000)/_density)*(Uw-_velocity).norm()*(Uw-_velocity);  //???????????
     //fluidForceh = fluidForceha + fluidForcehw;
   }   

    return forceFines;
  } 


  //Define a function(s) to move fines per timestep
  void fineTimestep(vector<Vector4d> & Fine_Grains, vector<Vector2d> & Fine_Velocities)
  {
    //Obtain forces from Fluid Drag (ignore intergranular contact for now)
    vector<Vector2d> forceFines(Fine_Grains.size()); 
    forceFines = finefluidinteraction(Fine_Grains, Fine_Velocities);  //Contribute fluid interaction force
    double fine_mass;
    
    //Update Velocities with Forces and Positions with Velocities (Verlet)
    for (size_t i = 0; i < Fine_Grains.size(); i++) 
    {
      fine_mass = _grains[0].getDensity() * Fine_Grains[i](0) * Fine_Grains[i](1) * 0.001; //Assuming flat cylinder fine for simplicity
      
      if (fine_mass > 0.0001)  //No force or position update for deleted grains
      {
        Fine_Velocities[i] = 1/(1+_gDamping*_dt/2)*( (1-_gDamping*_dt/2)*Fine_Velocities[i] + _dt*forceFines[i]/fine_mass );
        Fine_Grains[i](2) = Fine_Grains[i](2) + _dt*Fine_Velocities[i](0);
        Fine_Grains[i](3) = Fine_Grains[i](3) + _dt*Fine_Velocities[i](1);
      }   
    }

    //Account for PBCs
    for (size_t i = 0; i < Fine_Grains.size(); i++) {
      if (Fine_Grains[i](2) < 0.0 ) {
         Fine_Grains[i](2) += _offset(0);
      }
      else if (Fine_Grains[i](2) > _offset(0)) {
        Fine_Grains[i](2) -= _offset(0);
      }
      if (Fine_Grains[i](3) < 0.0 ) {
         Fine_Grains[i](3) += _offset(1);
      }
      else if (Fine_Grains[i](3) > _offset(1)) {
        Fine_Grains[i](3) -= _offset(1);
      }
    }
  }


  //************END FINE TIME STEPPING  *******************************//            


 //    //Remove Grain due to Melting (Fit to 2D)    //Verify that Grain Identity is Correct over Time!!!!
  // void RemoveGrain(const size_t & grainid){
  // //_grains[grainid] = _grains.back();
  // //_grains.pop_back();
  // //_ngrains--;
  
  // _grains.erase(_grains.begin()+grainid);
  // _ngrains=_grains.size();
  // _globalGrainState.resize(_ngrains, _nwalls);
  // _globalGrainState.reset();
  // }


  // change friction for particles
  void changeMuGrains(const double & newmu) {
    for (size_t i = 0; i < _ngrains; ++i) {
      _grains[i].changeMu(newmu);
    }
  }

  void changeKnGrains(const double & newkn) {
    for (size_t i = 0; i < _ngrains; ++i) {
      _grains[i].changeKn(newkn);
    }
  }


  void changeOffset(const Vector2d & offset){
    _offset = offset;
  }
  void changeDt(const double & dt) {
    _dt = dt;
  }
  void changeGdamping(const double & gDamping) {
    _gDamping = gDamping;
  }

   void changestepup(size_t & stepup) {
    _stepup = stepup;
  }

 void changeSlopeDir(double & slopedir) {
    _slopedir = slopedir;
  }

 void changeFlowAngle(double & flowang) {
    _flowangle = flowang;
  }
  
   void changeFlowSpeed(double & flowspeed) {
    _flowspeed = flowspeed;
  }

   void changeFlowForce(double & flowfor) {
    _flowforce = flowfor;
  }
  
  ///Change for fluid grid (January 9, 2023)
  void changeUw (vector<Vector2d> & newUw)
  {
      _Uwg = newUw;
  }
  
  //For Realistic Data
  void changeoceanTemp (vector<double> & newoceanTemp)
  {
      _oceanTemp = newoceanTemp;
  }

  double PointsArea ( const vector<Vector2d> & VecOrigin ) 
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

  double MeanCaliperD (double & Area)
  {
     return (2 * sqrt(Area /3.141592653589793));
  }


    // const get methods
  const vector<Grain2d> & getGrains() const {
    return _grains;
        
  }
    
    const vector<Grain2d> & getGrainsWall() const {
        return _grainsWall;
        
    }
    
    
    const vector<Wall2d> & getWalls() const {
        return _walls;
    }

  const GrainState2d & getGrainState() const {
    return _globalGrainState;
  }

  //   non const get methods (members can be changed such as wall positions)
  vector<Grain2d> & getGrainsNonConst() {
    return _grains;
  }
    
    
    
    
    
    vector<Wall2d> & getWallsNC() {
        return _walls;
    }
    
    
  const vector<Vector2d> & getGrainForces() const {
    return _globalGrainState._grainForces;
  }
  const vector<double> & getGrainMoments() const {
    return _globalGrainState._grainMoments;
  }
  const Vector3d & getStress() const {
    return _globalGrainState._stressVoigt;
  }


  void outputStress() {
        //cout << "Fx = " << _globalGrainState._grainForces(0) << endl;   //??
    //cout << "Fy= " << _globalGrainState._grainForces(1) << endl;     //??
    //cout << "M = " << _globalGrainState._grainMoments << endl;      //??

    cout << "sigxx = " << _globalGrainState._stressVoigt(0) << endl;
    cout << "sigyy = " << _globalGrainState._stressVoigt(1) << endl;
    cout << "sigxy = " << _globalGrainState._stressVoigt(2) << endl;
  }
  
  //Bonded Particle Method Properties (May 7, 2023)
  void createBonds() {
        Vector3i bondInfoIJ;
        size_t grainGrainBonds = 0;
        
        
        //Shuffle bond order before creating bonds in regular DEM
        bool shuffle_grains = false; 
        vector<size_t> iter(_ngrains);
        if (shuffle_grains){
            cout << "Shuffle Grain Index" << endl;
            for (size_t g = 0; g < _ngrains; g++){
                iter[g] = _ngrains - 1 - g;  //Invert just for to check
                cout << "Shuffle: " << g << " into: " << iter[g] << endl;
            }
            //random_shuffle(iter.begin(),iter.end());
            //std::random_device rd;
            //std::mt19937 grrr(rd());
            //std::shuffle(iter.begin(), iter.end(), grrr);
        }
        
        // Create bonds between grains FULL LS-DEM
        if (_only_dem == false){
            for (size_t i = 0; i<_ngrains;i++){
                for (size_t j = i+1; j < _ngrains; j++){
                    bondInfoIJ = _grains[i].createInterparticleBond(_grains[j]);
                    grainGrainBonds++;
                }//Bonded Particle Method Properties (May 7, 2023)
            }
        }
        //Regular DEM
        else{
            if (shuffle_grains == false){
                for (size_t i = 0; i<_ngrains;i++){
                    for (size_t j = i+1; j < _ngrains; j++){
                        //Avoid wasting time comparing far grains
                        if ( (_grains[i].getPosition() - _grains[j].getPosition()).norm() < (_grains[i].getRadius() + _grains[j].getRadius() + _grains[i].getCohesiveDistance()) ){
                            bondInfoIJ = _grains[i].createInterparticleBondDEM3(_grains[j]);  //2 is multiple surface, 3 is only 1 for efficiency
                            grainGrainBonds++;
                        }
                    }//Bonded Particle Method Properties (May 7, 2023)
                }
            }
            else{
                for (size_t i = 0; i<_ngrains;i++){
                    for (size_t j = i+1; j < _ngrains; j++){
                        //Avoid wasting time comparing far grains
                        if ( (_grains[iter[i]].getPosition() - _grains[iter[j]].getPosition()).norm() < (_grains[iter[i]].getRadius() + _grains[iter[j]].getRadius() + _grains[iter[i]].getCohesiveDistance()) ){
                            bondInfoIJ = _grains[iter[i]].createInterparticleBondDEM3(_grains[iter[j]]);  //2 is multiple surface, 3 is only 1 for efficiency
                            grainGrainBonds++;
                        }
                    }//Bonded Particle Method Properties (May 7, 2023)
                }                
            }
        }
        
        printf("Bonds Created\n");
    }

    void outputBonds() {
        vector<Vector3i> bondInfoIJ;
        size_t grainGrainBonds = 0;

        for (size_t i = 0; i < _ngrains; i++){
            bondInfoIJ = _grains[i].getBondInfo();
            for (size_t j = 0; j < bondInfoIJ.size(); ++j) {
                if (bondInfoIJ[j](0) == 1) {
                    if (size_t(bondInfoIJ[j](2)) < (_ngrains+1)) {
                        grainGrainBonds++;
                    }
                }
            }
        }
        printf("Grain-grain bonds: %lu\n", grainGrainBonds);
        //Print detail of bonds
        // for (size_t i = 0; i < _ngrains; i++){
        //     bondInfoIJ = _grains[i].getBondInfo();
        //     for (size_t j = 0; j < bondInfoIJ.size(); ++j) {
        //         if (bondInfoIJ[j](0) == 1) {
        //             if (size_t(bondInfoIJ[j](2)) < (_ngrains+1)) {
        //                 cout << "Bond Data Grain i: " << i << " data: " <<  bondInfoIJ[j](0) << " " << bondInfoIJ[j](1) << " "  << bondInfoIJ[j](2) << endl;
        //             }
        //         }
        //     }
        // }
    }
    
    void outputFullBonds(vector<Vector3i> & RbondInfoIJexp, vector<size_t> & Rexp_bondpts, vector<Vector2d> & Rexp_bondForceNormal, vector<Vector2d> & Rexp_bondForceShear, 
                          vector<double> & Rexp_bondThisMomentNormal, vector<Vector2d> & Rexp_bondNormals, vector<double> & Rexp_fSig, vector<double> & Rexp_fTau, vector<double> & Rexp_sigrF, vector<double> & Rexp_taurF ) {
        vector<Vector3i> bondInfoIJ; //Temporary per grain to filter active

        //Data to get bond info for all grains
        vector<Vector3i>    bondInfoIJexp;
        vector<size_t>    exp_bondpts;
      vector<Vector2d>  exp_bondForceNormal;
      vector<Vector2d>  exp_bondForceShear;
      vector<double>    exp_bondThisMomentNormal;
      vector<Vector2d>  exp_bondNormals;
      vector<double>    exp_fSig;
      vector<double>    exp_fTau;
      vector<double>    exp_sigrF;
      vector<double>    exp_taurF;        

        for (size_t i = 0; i < _ngrains; i++){
            bondInfoIJ = _grains[i].getBondInfo();
            for (size_t j = 0; j < bondInfoIJ.size(); ++j) {
                if (bondInfoIJ[j](0) == 1) { //Only save active bonds
                    if (size_t(bondInfoIJ[j](2)) < (_ngrains+1)) { //Additional restriction, please check
                        bondInfoIJexp.push_back(bondInfoIJ[j]);
                        //Now use all bond get functions to add to bond data
                        exp_bondpts.push_back(_grains[i].getbondpts()[j]);
                        exp_bondForceNormal.push_back(_grains[i].getbondForceNormal()[j]);
                        exp_bondForceShear.push_back(_grains[i].getbondForceShear()[j]);
                        exp_bondThisMomentNormal.push_back(_grains[i].getbondThisMomentNormal()[j]);
                        exp_bondNormals.push_back(_grains[i].getbondNormals()[j]);
                        exp_fSig.push_back(_grains[i].getfSig()[j]);
                        exp_fTau.push_back(_grains[i].getfTau()[j]);
                        exp_sigrF.push_back(_grains[i].getsigrF()[j]);
                        exp_taurF.push_back(_grains[i].gettaurF()[j]);
                    }    
                }
            }
        }
        
        //Save output to send to main
        RbondInfoIJexp        =  bondInfoIJexp;
        Rexp_bondpts         =  exp_bondpts;
      Rexp_bondForceNormal =  exp_bondForceNormal;
      Rexp_bondForceShear  =  exp_bondForceShear;
      Rexp_bondThisMomentNormal  =    exp_bondThisMomentNormal;
      Rexp_bondNormals =  exp_bondNormals;
      Rexp_fSig      =  exp_fSig;
      Rexp_fTau      =  exp_fTau;
      Rexp_sigrF       =  exp_sigrF;
      Rexp_taurF       =  exp_taurF;
    }
    
    void changeOnlyDEM(bool & new_only_dem) {
    _only_dem = new_only_dem;
  }
  
    //Function to create defects or notches in BPM
  void weakenBonds(double & xWeak, double & yWeak){
     size_t ct_break = 0;
     vector<Vector3i> bondInfoIJ;
       for (size_t i = 0; i < _ngrains; i++){
            bondInfoIJ = _grains[i].getBondInfo();
            for (size_t j = 0; j < bondInfoIJ.size(); ++j) {
                if (bondInfoIJ[j](0) == 1) {
                    if (size_t(bondInfoIJ[j](2)) < (_ngrains+1)) {
                        //Inspection
                        //cout << "grain: " << i << " bond ptx: " << _grains[i].getPointList()[_grains[i].getbondpts()[j]](0) << " bond pty: " << _grains[i].getPointList()[_grains[i].getbondpts()[j]](1) << " bond Normals x: " << _grains[i].getbondNormals()[j](0) << " bonds Normals y: " << _grains[i].getbondNormals()[j](1) << endl; 
                        if (   (_grains[i].getPointList()[_grains[i].getbondpts()[j]](1) > yWeak &&  _grains[i].getPointList()[_grains[i].getbondpts()[j]](1) < yWeak + _grains[i].getRadius())  &&  (abs(_grains[i].getbondNormals()[j](1)) > 0.0) &&  (_grains[i].getPointList()[_grains[i].getbondpts()[j]](0) < xWeak)   ) {  //Horizontal Left NOTCH   //Use vertically normal bonds, horizontals are: -1.0, 0.0 and within yWeak location and grain size as reference to bound defect to one line
              //if (   (_grains[i].getPointList()[_grains[i].getbondpts()[j]](0) > xWeak &&  _grains[i].getPointList()[_grains[i].getbondpts()[j]](0) < xWeak + _grains[i].getRadius())  &&  (abs(_grains[i].getbondNormals()[j](0)) > 0.0) &&  (_grains[i].getPointList()[_grains[i].getbondpts()[j]](1) < yWeak)   ) {  //Vertical Down NOTCH 
              //if (   (_grains[i].getPointList()[_grains[i].getbondpts()[j]](0) > xWeak &&  _grains[i].getPointList()[_grains[i].getbondpts()[j]](0) < xWeak + _grains[i].getRadius())  &&  (abs(_grains[i].getbondNormals()[j](0)) > 0.0) &&  (_grains[i].getPointList()[_grains[i].getbondpts()[j]](1) > yWeak)   ) {  //Vertical Up NOTCH 
                   cout << "Bond broken due to generation of defect or weakness!!!" << endl;
                   cout << "xpos: " << _grains[i].getPointList()[_grains[i].getbondpts()[j]](0) <<  " ypos: " << _grains[i].getPointList()[_grains[i].getbondpts()[j]](1) << endl; 
                   _grains[i].deleteInterparticleBondDEM(j);
                   ct_break++;
                        }
                    }
                }
            }
        }
        cout << "Bonds removed due to defect: " << ct_break << endl;
        return;
  }
  
  double bond_logistic(double & tinput){
      double th_temp = 0.1;
      double coeff_temp = 3.6;
      //f(T) = ( (1 - th_temp) / (1+exp(3.6 * T)) ) + th_temp
      return ( (1 - th_temp) / (1 + exp(coeff_temp * tinput)) ) + th_temp;
  }

  double bond_qvODE(double & origrF, double & tinput, const size_t & dstep, double & qvert, double & Qatm, double & Aatm, double & Batm){
        double newrF;
        double th_temp = 0.00001; //Minimum weakness value 
        
        //Material values
        double meltTemp= -1.8; //0 for LS, -1.8 degree Celsius for normal implementations //Melting point of Ice (confirm for salinity)
        double rho_icem = 910; //kg/m3
        double Lf = 330000; //J/kg
        double a_i = 0.7; //Albedo of sea ice
    
        //Thickness reduction terms from ODE
        double meltV = (qvert/(rho_icem*Lf)); //For coarse ocean melt
        double meltVSun = 0.00 * ( - Qatm*(1-a_i) + (Aatm + Batm *(meltTemp)) ) / (rho_icem*Lf) ;  //For coarse sun, atmosphere melt, multiply by zero to only consider ocean effect

        //Thickness ODE used to weaken bonds
        newrF =      max(  (origrF + dstep*(meltV)*(meltTemp - tinput) + dstep * meltVSun) ,  th_temp);
        //New thick     //Old thick    //qvert comp.  * (Tmelt-T_fave)        //Denom                      //Q, A, BT part
      newrF =      min(origrF, newrF); //Avoid healing    
     
      return newrF;
  }
  
  //Include logistic bond weakening based on temperature
  void updateBondTemp(const size_t & dstep, double & qvert, double & Qatm, double & Aatm, double & Batm, double & og_mass){ 
        bool irr_damage = true; //Apply irreversible damage
      vector<Vector3i> bondInfoIJ;
      double temp_bond, bond_redux, prevsigrF, prevtaurF, temp_G;
      bool grain_thick = true; //Turn on to also reduce grain thickness and change their mass
      double new_thick, old_thick, mult_loss, old_mass;
      Vector2d posG;
      
        for (size_t i = 0; i < _ngrains; i++){
            bondInfoIJ = _grains[i].getBondInfo();
            
            if (grain_thick){
                old_thick = _grains[i].getThickness();
                posG(0) = _grains[i].getPosition()(0);
                posG(1) = _grains[i].getPosition()(1);
                temp_G = round_Ocean(posG, _oceanTemp, _x_cells, _y_cells, _offset);
                if (temp_G > -1.8){
                    new_thick = bond_qvODE(old_thick, temp_G, dstep, qvert, Qatm, Aatm, Batm);
                    _grains[i].changeThickness(new_thick);
                    mult_loss = new_thick/_grains[i].getThickness0();
                    _grains[i].changeMass(og_mass*mult_loss);
                }
                
            }
            
            for (size_t j = 0; j < bondInfoIJ.size(); ++j) {
                if (bondInfoIJ[j](0) == 1) {
                    if (size_t(bondInfoIJ[j](2)) < (_ngrains+1)) {
                        //Find position of an admissible bond
                        Vector2d posBond;
                        posBond(0) = _grains[i].getPointList()[_grains[i].getbondpts()[j]](0);
                        posBond(1) = _grains[i].getPointList()[_grains[i].getbondpts()[j]](1);
                        
                        //Find temperature based on position
                        temp_bond = round_Ocean(posBond, _oceanTemp, _x_cells, _y_cells, _offset);
                        
                        if (temp_bond > -1.8){   //Only change if above fusion point of sea ice
                            prevsigrF = _grains[i].getsigrF()[j];
                            prevtaurF = _grains[i].gettaurF()[j];
                            
                            // //Find reduction factor using logistic function
                            //bond_redux = bond_logistic(temp_bond);
                            // Use thickness reduction ODE for reference
                            bond_redux = bond_qvODE(prevsigrF, temp_bond, dstep, qvert, Qatm, Aatm, Batm); //Modify if tau is reduced differently!!!
        
                            //Modify bond information
                            double new_sigrF = bond_redux;
                            double new_taurF = bond_redux;
                            
                            if (irr_damage){
                                 if (prevsigrF > new_sigrF){
                                      _grains[i].changesigrF(j, new_sigrF);
                                 }
                                 if (prevtaurF > new_taurF){
                                      _grains[i].changetaurF(j, new_taurF);
                                 }
                            }
                            else{
                                //Update bond info in grains (recovery for cold waters, reversible damage)
                                _grains[i].changesigrF(j, new_sigrF);
                                _grains[i].changetaurF(j, new_taurF);
                            }
                        }
                  }
                }
            }
        }
  }
  
  
  
  
  
  
  //BEGIN DAMPING FUNCTIONS
  //Clustering functions
  
  //Exist for a particular cluster
  bool exists_list(vector<size_t> & listV, size_t & elementV){
        for (size_t i = 0; i < listV.size(); ++i) {
            if (elementV == listV[i]){
                return true;
            }
        }
        return false;
  }
  
  //Exist in all clusters or segments
  bool exists_seg(vector<vector<size_t>> & segment, size_t & elementV){
        for (size_t i = 0; i < segment.size(); ++i) {
            if ( exists_list(segment[i],elementV) ){
                return true;
            }
        }
        return false;
  }

    //No need
    // size_t exists_seg_index(vector<vector<size_t>> & segment, size_t & elementV){
    //     for (size_t i = 0; i < segment.size(); ++i) {
    //         if ( exists_list(segment[i],elementV) ){
    //             return i;
    //         }
    //     }
    //     return 10000000000; //If badly used we'll get an error
    // }
    
    //Fixed floe detection functions from floe list
  vector<size_t> fixed_detect(vector<size_t> & floeList, vector<size_t> & fixed_floes){
      size_t idx;
      size_t lenList = floeList.size();
      bool fixg; //Assume not fixed
      vector<size_t> fixedList(lenList);
      for (size_t i = 0; i < lenList; ++i) {
          fixg = false; //Assume not fixed everytime and then try to prove wrong.
          idx = floeList[i];
          for (size_t j = 0; j < fixed_floes.size(); ++j) {
              if (idx == fixed_floes[j]){
                  fixg = true;
                  break;
              }
          }
         //Define fix cluster list per floe in list
          if(fixg){
              fixedList[i] = 1;
          }
          else{
              fixedList[i] = 0;
          }
      }
      return fixedList;
  }
  
  //Distance and fixed floe detection functions (using cluster fixed list)
  bool fixed_find(vector<size_t> & floeList, vector<size_t> & fixedList, size_t & edge_index){
      size_t idx;
      for (size_t i = 0; i < floeList.size(); ++i) {
          idx = floeList[i];
          if(idx == edge_index && fixedList[i] == 1){
              return true;
          }
      }
      return false;
  }
  
  bool cluster_intersects(vector<size_t> & floeList, Vector2d & PointFind){
      size_t idx;
      for (size_t i = 0; i < floeList.size(); ++i) {
          idx = floeList[i];
            if ( (_grains[idx].getPosition() - PointFind).norm() <= _grains[idx].getRadius() ){
                return true;
            }
      }
      return false;
  }
  
  size_t min_distFx(vector<size_t> & floeList, Vector2d & PointFind){
      double minD = 100000000000.000000;
      size_t idx, edge_index;
      for (size_t i = 0; i < floeList.size(); ++i) {
          idx = floeList[i];
            if ( (_grains[idx].getPosition() - PointFind).norm() < minD ){
                minD = (_grains[idx].getPosition() - PointFind).norm();
                edge_index = idx;
            }
      }
      return edge_index;
  }
  
  
  double max_distFx(vector<size_t> & floeList){
      double maxD = 0.000000;
      size_t idx, jdx;
      for (size_t i = 0; i < floeList.size(); ++i) {
          idx = floeList[i];
          for (size_t j = i+1; j < floeList.size(); ++j) {
              jdx = floeList[j];
              if ( (_grains[idx].getPosition() - _grains[jdx].getPosition()).norm() > maxD ){
                  maxD = (_grains[idx].getPosition() - _grains[jdx].getPosition()).norm();  
              }
          }
      }
      return maxD;
  }

  //Sum of included bonds to keep looping list  
  size_t sum_included(vector<size_t> & bond_included){
      size_t sumV = 0;
      for (size_t i = 0; i < bond_included.size(); i++){
          sumV += bond_included[i];
      }
      return sumV;
  }
  
  //Average of a vector (self)
  double average_vector(vector<double> & vecV){
      size_t sizeV = vecV.size();
      double sumV = 0.0;
      if (sizeV == 0){
          return 0.0;
      }
      for (size_t i = 0; i < sizeV; i++){
          sumV += vecV[i];
      }
      return sumV/double(sizeV);
  }
  
  vector<size_t> exhaust_options(vector<Vector3i> & bondV, size_t & rootg, vector<size_t> & bond_included){
    vector<double> ready_list; //List to track linked grains covered
    ready_list.push_back(0.0);
    vector<size_t> temp_seg; //Temporary segment of full floe
    temp_seg.push_back(rootg);
    size_t idg  = 0; //Placeholder of grain to inspect
    size_t floeL; //For follower
    size_t floeF; //For leader
    while (average_vector(ready_list) < 1){
        for (size_t i = 0; i < bondV.size(); i++){
            if (bond_included[i] == 0){ 
                floeL = size_t(bondV[i](1));
                floeF = size_t(bondV[i](2));
                if (  (temp_seg[idg] == floeL) || (temp_seg[idg] == floeF)  ){   //Avoid repeats
                    bond_included[i] = 1;
                    if (temp_seg[idg] == floeL){
                        if (exists_list(temp_seg, floeF) == false){
                            temp_seg.push_back(floeF);
                            ready_list.push_back(0.0);
                        }
                    }
                    else{
                        if (exists_list(temp_seg, floeL) == false){
                            temp_seg.push_back(floeL);
                            ready_list.push_back(0.0);
                        }
                    }
                }
            }
        }
        ready_list[idg] = 1;
        idg += 1;
    }
    return temp_seg;
  }

  //Clustering of Floes depending on Bond Connectivity (Critical Function)
  bool bond_clustering(vector<vector<size_t>> & clusterFloes, vector<vector<size_t>> & clusterFixed, vector<size_t> & fixed_floes){
      vector<vector<size_t>> segments;
      vector<vector<size_t>> fixed_segments;
      vector<size_t> temp_segment;
      vector<size_t> temp_fixed;
      
      //Generate Bond Connectivity list using grain Bond Data
      //Data to get bond info for all grains
        vector<Vector3i>    bondInfoIJ;  //Place holder to go all over
        vector<Vector3i>    bondV;       //List of all 0: on/off 1: leader index 2: follower index Bond Connectivity Data
        size_t floeL; //For follower
        size_t floeF; //For leader
         
        for (size_t i = 0; i < _ngrains; i++){
            bondInfoIJ = _grains[i].getBondInfo();
            for (size_t j = 0; j < bondInfoIJ.size(); ++j) {
                if (bondInfoIJ[j](0) == 1) { //Only save active bonds
                    if (size_t(bondInfoIJ[j](2)) < (_ngrains+1)) { //Additional restriction, please check
                        bondV.push_back(bondInfoIJ[j]);
                    }    
                }
            }
        }
        //Bond included to speed up process
        vector<size_t> bond_included(bondV.size()); //
        for (size_t i = 0; i < bondV.size(); i++){
            bond_included[i] = 0; //Initialized all zeros or not included. 1 means already included in cluster and must be skipped.
        }
        
        size_t counter_seek = 0; //In case you need iterations for large ensembles.
        size_t max_iters = 22; //To achieve full connectivity this methods needs a lot of iterations for larger floes. Works for a 8000 floe, 25K bond scale. Less values work for smaller floes. Check damping map.  //17-->
        size_t min_size = 4; //Minimum cluster size to even consider for damping
        
        if (bondV.size() > 0){  //Else do nothing
            while( sum_included(bond_included) < bondV.size() ){ //Only leave if all bonds are added to clusters
                for (size_t i = 0; i < bondV.size(); i++){   //Loop over all bonds to cluster
                    if (bond_included[i] == 0){  //Only analyze if the bonds is not added to the cluster already
                        floeL = size_t(bondV[i](1));
                        floeF = size_t(bondV[i](2));
                        if (temp_segment.size() == 0){ //First consider empty temp_segment or new segment
                            if ( exists_seg(segments, floeL) == false && exists_seg(segments, floeF) == false ){  //None present in other parts of all the clusters, hence add to list
                                bond_included[i] = 1;
                                temp_segment.push_back(floeL);
                                temp_segment.push_back(floeF);
                            }
                            else{  //Repeated or disconnect bonds, just mark as done. Since a b already belong to the prior cluster. why would a or b be repeated in another?
                                bond_included[i] = 1;
                            }
                        }
                        else{  //Now for a non-empty new temp_segment things need to be checked as well
                            //if ( ( exists_list(temp_segment, floeL) || exists_list(temp_segment, floeF) ) && ( ( exists_list(temp_segment, floeL) && exists_list(temp_segment, floeF) ) == false) ){ //At least one has to exist on the temp_Segment list for connection and not both at the same time to avoid repeats
                            if ( exists_list(temp_segment, floeL) || exists_list(temp_segment, floeF)  ){
                                bond_included[i] = 1;
                                if ( exists_list(temp_segment, floeL) == false and exists_seg(segments, floeL) == false){ //Only add if new to segment cluster and to current temp_segment
                                    temp_segment.push_back(floeL);
                                } 
                                if ( exists_list(temp_segment, floeF) == false and exists_seg(segments, floeF) == false){ //Only add if new to segment cluster and to current temp_segment
                                    temp_segment.push_back(floeF);
                                }
                            }
                        }
                    }
                }
                
                counter_seek += 1;
                if (counter_seek > max_iters){
                    if (temp_segment.size() >= min_size){
                        segments.push_back(temp_segment);
                        temp_fixed = fixed_detect(temp_segment, fixed_floes); //Make the same for fixed floe cluster list for convenience
                        fixed_segments.push_back(temp_fixed);
                    }
                    temp_segment.clear(); //Clear to start again
                    temp_fixed.clear();
                    counter_seek = 0; //Restart
                }
            }
        }
      else{
          return false;
      }
      
      //What  if too broken down and really no big segments formed
      if (segments.size() == 0){
          return false;
      }
      
      clusterFloes = segments;
      clusterFixed = fixed_segments;
      return true;
  }
  
  //Clustering of Floes depending on Bond Connectivity (Critical Function) (newer function)
  bool bond_clustering_v2(vector<vector<size_t>> & clusterFloes, vector<vector<size_t>> & clusterFixed, vector<size_t> & fixed_floes){
      vector<vector<size_t>> segments;
      vector<vector<size_t>> fixed_segments;
      vector<size_t> temp_segment;
      vector<size_t> temp_fixed;
      
      //Generate Bond Connectivity list using grain Bond Data
      //Data to get bond info for all grains
        vector<Vector3i>    bondInfoIJ;  //Place holder to go all over
        vector<Vector3i>    bondV;       //List of all 0: on/off 1: leader index 2: follower index Bond Connectivity Data
        size_t floeL; //For follower
        size_t floeF; //For leader
        size_t min_size = 4; //Minimum cluster size to even consider for damping
         
        for (size_t i = 0; i < _ngrains; i++){
            bondInfoIJ = _grains[i].getBondInfo();
            for (size_t j = 0; j < bondInfoIJ.size(); ++j) {
                if (bondInfoIJ[j](0) == 1) { //Only save active bonds
                    if (size_t(bondInfoIJ[j](2)) < (_ngrains+1)) { //Additional restriction, please check
                        bondV.push_back(bondInfoIJ[j]);
                    }    
                }
            }
        }
        //Bond included to speed up process
        vector<size_t> bond_included(bondV.size()); //
        for (size_t i = 0; i < bondV.size(); i++){
            bond_included[i] = 0; //Initialized all zeros or not included. 1 means already included in cluster and must be skipped.
        }
        
        if (bondV.size() > 0){  //Else do nothing
            while( sum_included(bond_included) < bondV.size() ){ //Only leave if all bonds are added to clusters
                for (size_t i = 0; i < bondV.size(); i++){   //Loop over all bonds to cluster
                    if (bond_included[i] == 0){  //Only analyze if the bonds is not added to the cluster already
                        floeL = size_t(bondV[i](1));
                        floeF = size_t(bondV[i](2));
                        
                        if ( exists_seg(segments, floeL) == false && exists_seg(segments, floeF) == false ){ 
                            temp_segment = exhaust_options(bondV, floeL, bond_included);
                            if (temp_segment.size() >= min_size){
                                segments.push_back(temp_segment);
                                temp_fixed = fixed_detect(temp_segment, fixed_floes); //Make the same for fixed floe cluster list for convenience
                                fixed_segments.push_back(temp_fixed);
                            }
                            temp_segment.clear(); //Clear to start again
                            temp_fixed.clear();
                        }
                        else{
                            cout << "HUH!!!! No more???: Clustering error!!!!" << endl;
                            exit(1);
                        }
                            
                   }
               }
            }
      }
      else{
          return false;
      }
      
      //What  if too broken down and really no big segments formed
      if (segments.size() == 0){
          return false;
      }
      
      clusterFloes = segments;
      clusterFixed = fixed_segments;
      return true;
  }
  
  //Needed for Under Ice Detection and Distance From Edge (Critical Function)
  void clusterStats(vector<vector<size_t>> & clusterFloes, vector<Vector2d> & clusterCentroidOut, vector<double> & clusterBradiusOut){
      size_t nClusters = clusterFloes.size();
      size_t lenList;
      vector<Vector2d> clusterCentroid(nClusters); //Vector of cluster centroids
      vector<double> clusterBradius(nClusters);     //Vector of cluster bounding circles
      size_t jdx;
      
      double avex, avey, max_dist;
      for (size_t i = 0; i < nClusters; ++i) {
          lenList = clusterFloes[i].size();
          avex = 0;
          avey = 0;
          //Get average values over all positions to get centroid in x and y. 
          for (size_t j = 0; j < lenList; ++j) {
              jdx = clusterFloes[i][j];
              avex += _grains[jdx].getPosition()(0);
              avey += _grains[jdx].getPosition()(1);
          }
          if (lenList > 0){
              clusterCentroid[i](0) = avex/lenList;
              clusterCentroid[i](1) = avey/lenList;
          }
          else{
              clusterCentroid[i](0) = 0.000;
              clusterCentroid[i](1) = 0.000;
          }
          //Get max distance among all floes of cluster
          max_dist = max_distFx(clusterFloes[i]);
          clusterBradius[i] = 0.575*max_dist; //0.7 is a bit better than 0.5 to avoid being too close to edge just in case, 0.575 is to be faster
      }
      clusterCentroidOut = clusterCentroid;
      clusterBradiusOut = clusterBradius;
  }
      
  //Distance of Grid Point from Ice Edge No fixed floes (Critical Function)
    void under_iceDEM_bond(double & dist_find, Vector2d & ptGrid, vector<vector<size_t>> & clusterFloes, vector<Vector2d> & clusterCentroid, vector<double> & clusterBradius){   
      size_t nClusters = clusterFloes.size();
      
      bool close_ice =  false; //Assume none, prove proximity as true
      for (size_t i = 0; i < nClusters; ++i) {
          if ( (ptGrid - clusterCentroid[i]).norm() < clusterBradius[i] ){
              close_ice = true;
              break;
          }
      }
        
        dist_find = 0.0; //For reference
        size_t edge_index;     
      double d_grid, d_edge, d_p, cosG, sinG, Rc; //Distance between grid_Point-to-centroid, edge_point-to-centroid, grid_Point-to-edge_point, cos and sine angle for d_grid respect to horizontal
      Vector2d cPoint; //Point in projected circle
      Vector2d edgePoint;
      if (close_ice){
          for (size_t i = 0; i < nClusters; ++i) {
              if ( (ptGrid - clusterCentroid[i]).norm() < clusterBradius[i] ){
                  //We are within action circle so we need to get all critical distances
                  d_grid = (ptGrid - clusterCentroid[i]).norm();
                  cosG = (ptGrid(0) - clusterCentroid[i](0)) / d_grid;
                  sinG = (ptGrid(1) - clusterCentroid[i](1)) / d_grid;
                  cPoint(0) = clusterCentroid[i](0) + clusterBradius[i] * cosG;
                  cPoint(1) = clusterCentroid[i](1) + clusterBradius[i] * sinG;
                  edge_index = min_distFx(clusterFloes[i], cPoint);
                  Rc = _grains[edge_index].getRadius();
                  edgePoint(0) = _grains[edge_index].getPosition()(0);
                  edgePoint(1) = _grains[edge_index].getPosition()(1);
                  d_edge = (edgePoint - clusterCentroid[i]).norm();
                  d_p = (ptGrid - edgePoint).norm();
                  if (d_grid > d_edge + Rc){
                     dist_find = 0.0;
                  }
                  else if (d_grid > d_edge && d_grid < d_edge + Rc){
                     dist_find = max(d_edge + Rc - d_grid, dist_find); 
                  }
                  else{
                     dist_find = max(d_p + Rc, dist_find); 
                  }
              }
          }
      }
      else{
          dist_find = 0.0;
      }
      
  } 
  
  //Distance of Grid Point from Ice Edge with fixed floes (Critical Function)
  void under_iceDEM_bond_fixed(double & dist_find, Vector2d & ptGrid, vector<vector<size_t>> & clusterFloes, vector<vector<size_t>> & clusterFixed, vector<Vector2d> & clusterCentroid, vector<double> & clusterBradius){   
      size_t nClusters = clusterFloes.size();
      size_t intersec_index;
      
      bool close_ice =  false; //Assume none, prove proximity as true
      for (size_t i = 0; i < nClusters; ++i) {
          if  ( cluster_intersects(clusterFloes[i], ptGrid) ){  //All floes but first intersect at small floe level
          //if ( (ptGrid - clusterCentroid[i]).norm() < clusterBradius[i] ){ //All floes based on Bradius
              close_ice = true;
              intersec_index = i;
              break;
          }
      }
      
      bool printC = false;
      if ( (ptGrid(0) == 100.0 && ptGrid(1) == 230.0) || (ptGrid(0) == 300.0 && ptGrid(1) == 230.0) || (ptGrid(0) == 200.0 && ptGrid(1) == 50.0) ){
          cout <<"Print samples" << endl;
          printC = false;
      }
        
      double start_angle, angle_val;
      double dist_min_fixed;
      bool fixed_floe_q;
      dist_find = 0.0; //For reference
      size_t edge_index, edge_indexT;    
      double d_grid, d_edge, d_p, cosG, sinG, Rc; //Distance between grid_Point-to-centroid, edge_point-to-centroid, grid_Point-to-edge_point, cos and sine angle for d_grid respect to horizontal
      Vector2d cPoint; //Point in projected circle
      Vector2d edgePoint;
      if (close_ice){
          //for (size_t i = 0; i < nClusters; ++i) {
              //if ( (ptGrid - clusterCentroid[i]).norm() < clusterBradius[i] ){
                  //We are within action circle so we need to get all critical distances
                  d_grid = (ptGrid - clusterCentroid[intersec_index]).norm();
                  cosG = (ptGrid(0) - clusterCentroid[intersec_index](0)) / d_grid;
                  sinG = (ptGrid(1) - clusterCentroid[intersec_index](1)) / d_grid;
                  cPoint(0) = clusterCentroid[intersec_index](0) + clusterBradius[intersec_index] * cosG;
                  cPoint(1) = clusterCentroid[intersec_index](1) + clusterBradius[intersec_index] * sinG;
                  //edge_index = min_distFx(clusterFloes[i], cPoint);
                  //fixed_floe_q = fixed_find(clusterFloes[i], clusterFixed[i], edge_index); //Assume fixed and then improve
                  fixed_floe_q = true;
                  //Hard way (Assume bottom fixed) (Modify when needed)
                  dist_min_fixed = 10000000000000;
                  if (fixed_floe_q){
                      //Fix rotational constant for bottom fix (sweep from 180 to 360 degrees)  //Overkill would be 0 ot 360. (More expensive)
                      start_angle = 0;  //360
                      for (size_t adx = 0; adx < 721; ++adx) {  //Inspect for 180 degrees   //for (size_t adx = 0; adx < 361; ++adx) { 
                          angle_val = start_angle + 0.5*double(adx); //In degrees
                          angle_val *= M_PI/180; //Change to radians
                          cPoint(0) = clusterCentroid[intersec_index](0) + clusterBradius[intersec_index] * cos(angle_val); //New iterative projected point
                          cPoint(1) = clusterCentroid[intersec_index](1) + clusterBradius[intersec_index] * sin(angle_val);
                          edge_indexT = min_distFx(clusterFloes[intersec_index], cPoint);
                          fixed_floe_q = fixed_find(clusterFloes[intersec_index], clusterFixed[intersec_index], edge_indexT);
                          if (fixed_floe_q == false && (ptGrid - _grains[edge_indexT].getPosition()).norm() < dist_min_fixed){
                              dist_min_fixed = (ptGrid - _grains[edge_indexT].getPosition()).norm();
                              edge_index = edge_indexT;
                          }
                      }
                  }
                  //Easy way
                  Rc = _grains[edge_index].getRadius();
                  edgePoint(0) = _grains[edge_index].getPosition()(0);
                  edgePoint(1) = _grains[edge_index].getPosition()(1);
                  d_edge = (edgePoint - clusterCentroid[intersec_index]).norm();
                  d_p = (ptGrid - edgePoint).norm();
                  if (d_grid > d_edge + Rc){
                     dist_find = 0.0;
                  }
                  else if (d_grid > d_edge && d_grid < d_edge + Rc){
                     dist_find = max(d_edge + Rc - d_grid, dist_find); 
                  }
                  else{
                     dist_find = max(d_p + Rc, dist_find); 
                  }
                  
                  if (printC){
                      cout << "dist_find: " << dist_find << endl;
                      cout << "dist_min_fixed: " << dist_min_fixed << endl;
                      cout << "Bradius: " << clusterBradius[intersec_index] << endl;
                      cout << "ptGrid: " << ptGrid(0) << " , " << ptGrid(1) <<endl;
                      cout << "Centroid: " << clusterCentroid[intersec_index](0) << " , " << clusterCentroid[intersec_index](1) <<endl;
                      cout << "cPoint: " << cPoint(0) << " , " << cPoint(1) <<endl;
                      cout << "edgePoint: " << edgePoint(0) << " , " << edgePoint(1) <<endl;
                      cout << "d_grid: " << d_grid << endl;
                      cout << "d_edge: " << d_edge << endl;
                      cout << "d_p: " << d_p << endl;
                      cout << "Rc: " << Rc << endl;
                      cout << "Angles: " << cosG << " , " << sinG <<endl;
                  }
              //}
          //}
      }
      else{
          dist_find = 0.0;
      }
      
  } 
  
  //YFRONT functions
  double vec_ave_y(vector<size_t> & clusterV){
      size_t idx;
      double avesum = 0;
      for (size_t i = 0; i < clusterV.size(); ++i) {
          idx = clusterV[i];
          avesum += _grains[idx].getPosition()(1);
      }
      if (clusterV.size() < 1) {
          return 0;
      }
      else{
          return avesum / clusterV.size();
      }
  }
  
  double maxY_and_cut(vector<size_t> & clusterV){
      double maxY = 0.00;
      double posY;
      size_t idx, max_ind;
      for (size_t i = 0; i < clusterV.size(); ++i) {
          idx = clusterV[i];
          posY = _grains[idx].getPosition()(1);
          if (posY > maxY){
              maxY = posY;
              max_ind = i;
          }
      }
      clusterV.erase(clusterV.begin()+max_ind);
      return maxY;
  }
  
  double maxY_and_buffer(vector<size_t> & clusterV){
      double maxY = 0.00;
      double posY;
      double buffer_dist = 30.00; //Orig 30.00
      size_t idx;
      for (size_t i = 0; i < clusterV.size(); ++i) {
          idx = clusterV[i];
          posY = _grains[idx].getPosition()(1);
          if (posY > maxY){
              maxY = posY;
          }
      }
      return maxY - buffer_dist;
  }
  
  double maxY_ave(vector<size_t> & clusterV, size_t & nave){
      double aveY = 0.00;
      vector<size_t> clusterVcopy; //To avoid problems
      clusterVcopy = clusterV;
      bool buffer_use = true; //Else is average of max points
      //Just do the average if too few points
      if ( nave >= clusterVcopy.size() ){
          return vec_ave_y(clusterVcopy); //Defines average y position of cluster
      }
      //Filter for real top nave or top - buffer
      else{
          if (buffer_use){
              aveY = maxY_and_buffer(clusterVcopy);  //Defines max y position of cluster minus a small buffer
              return aveY;
          }
          else{
              for (size_t i = 0; i < nave; ++i) {
                  aveY += maxY_and_cut(clusterVcopy); //Find max Y value and also cut vector to find next best and so on.
              }
              return aveY/nave;
          }
      }
  }
  
  double maxY_ave_fix(vector<size_t> & clusterV, size_t & nave){
      double aveY = 0.00;
      vector<size_t> clusterVcopy; //To avoid problems
      clusterVcopy = clusterV;
      bool buffer_use = true; //Else is average of max points
      //Return zero if no fixed floes
      if (clusterV.size()<1){
          return 0.0;
      }
      
      //Just do the average if too few points
      if ( nave >= clusterVcopy.size() ){      
          return vec_ave_y(clusterVcopy);       //Defines average y position of cluster
      }
      //Filter for real top nave or top - buffer
      else{
          if (buffer_use){
              //aveY = maxY_and_buffer(clusterVcopy);  //Defines max y position of fixed cluster minus a small buffer
              aveY = vec_ave_y(clusterVcopy) + 20.00;  //Defines mean y position of fixed cluster plus a small buffer
              return aveY;
          }
          else{
              for (size_t i = 0; i < nave; ++i) {
                  aveY += maxY_and_cut(clusterVcopy); //Find max Y value and also cut vector to find next best and so on.
              }
              return aveY/nave;
          }
      }
  }
  
  bool have_fixed(vector<size_t> & vFixed){
      bool fixed_exists = false;
      for (size_t i = 0; i < vFixed.size(); ++i) {
          if (vFixed[i] == 1){
              fixed_exists = true;
              break;
          }
      }
      return fixed_exists;
  }
  
  double find_yfront(vector<vector<size_t>> & clusterFloes, vector<vector<size_t>> & clusterFixed){ 
      
      double y_front;
      double y_front_all = 300.00; //Let's say that the y front is the maximum point of the floe with fixed floes, based on how it fails, only bottom floe will have fixed floes. (Example value)
      y_front_all = 100.00; //Assume is close to original ice geometry
      double y_limit = 11.00 + 40.00; //Cannot go below this to simulate cantilever //54.00 good for Fram Up, 12 for lower fram rec
      size_t nave = 0;  //20
      size_t naveF = 0; //1
      double y_front_fix = 10.00 +  40.00;
      double y_limit_max = 90.00 + 150.00; //Damp cannot go higher than that
      //Remove damping for now:
      //y_front = 0.0;
      //return y_front;
    
      //RETURN AGAIN TO THIS      
      size_t nClusters = clusterFloes.size();
      size_t nClustersF = clusterFixed.size();
      
    //   //For moving floes //Not applies for damp since they are loose.
    //   for (size_t i = 0; i < nClusters; ++i) {
    //           y_front_all = maxY_ave(clusterFloes[i], nave); //Average of top X floes (say 20 since it will not be reduced that much)
    //   }
      size_t idx_max = 0;
      size_t nfloe_w_fix = 0;
      double clusterymax = 0.0;
      //For fixed floes
      for (size_t i = 0; i < nClusters; ++i) {
          if ( have_fixed(clusterFixed[i]) ){ //Now assume more than 1 clusters are fixed. Use max value of fixed but no more than initial ref value //Assume only 1 cluster is fixed, rest aren't
              //Max value of only fixed floes, keep accumulating the max value of fixed floe clusters
              clusterymax = maxY_ave_fix(clusterFloes[i], naveF);
              if ( clusterymax >  y_front_fix  ){
                  y_front_fix =  clusterymax;
                  idx_max = i;
              }
              nfloe_w_fix++;
          }
      }
      
      //Use the smallest of the two options (fix limit should win)
      //y_front = min(y_front_all, y_front_fix);
      y_front = min(y_limit_max, y_front_fix);
      cout << "Y_front found from clusters: " << y_front_fix << endl;
      cout << "Y_front found before limit y = 11: " << y_front << endl;
      cout << "Y_front index found: " << idx_max << endl;
      cout << "Number of floes with fixed: " << nfloe_w_fix << endl;
      //RE-DO this part
      //   for (size_t i = 0; i < nClusters; ++i) {
      //       if ( have_fixed(clusterFixed[i]) ){ //Assume only 1 cluster is fixed, rest aren't
      //           y_front_all = maxY_ave(clusterFloes[i], nave); //Average of top X floes (say 20 since it will not be reduced that much)
      //           //y_front_fix = maxY_ave_fix(clusterFixed[i], naveF); //Average of only fixed floes
      //           y_front = min(y_front_all, y_front_fix);
      //           break;
      //       }
      //   }
      return max(y_front, y_limit);
  }  

  //Distance of Grid Point from a y front Ice Edge in 1D
  double under_iceDEM_bond_fixed_yfront(Vector2d & ptGrid, double & y_front){   
      double dist_find;
      bool below_front = false;
      
      if (ptGrid(1) <= y_front){
          below_front = true;
      }
     
      if (below_front){
          dist_find = y_front - ptGrid(1);
      }
      else{
          dist_find = 0.0;
      }
      return dist_find;
  } 
  
  //Non Fixed Floe Version and Fixed version together (WARNING DOES NOT YET ACCOUNT FOR PBC!!!!!!!). If fixed_floes is empty we always we 0 for fixedCluster so it's okay.
  double updateDamp(vector<size_t> & fixed_floes){
      double Ds =  2.5; //Limit for damping //2.5 //20
      vector<double> newdampMat(_x_cells * _y_cells);
      vector<double> dist_edge(_x_cells * _y_cells); //Distance used to find damping matrix
      
      //1. Do floe segmentation function
      //INPUT ALL FLOE BONDS FROM ALL GRAINS
      //OUTPUT: Vector of floe bond indices and Vector of grain indices joined into distinct floes, at least 1 per simulation, unless all is broken
      //IDENTIFY BEST DATA STRUCTURE FOR BONDS AND FLOES (CREATE!!!) //Also best outputs like centroids and x-y dists
      vector<vector<size_t>> clusterFloes; 
      vector<vector<size_t>> clusterFixed; 
      //Add special one for fixed floes in the other function
      //bond_clustering function that goes over all bonds and functions (CREATE!!!)
      bool bond_exist = bond_clustering_v2(clusterFloes, clusterFixed, fixed_floes);  //bond_clustering_fixed(clusterFloes, clusterFixed, fixed_floes);
      cout << "Number of clusters: " << clusterFloes.size() << endl; 
      //Accelerate things if no bonds exist
      if (bond_exist == false){
          for (size_t i = 0; i < _y_cells; i++) {
                for (size_t j = 0; j < _x_cells; j++) {
                    newdampMat[j+i*_x_cells] = 1.0; //By defect it will be 1 if D = 0 or no floes are bonded to damp
                }
            }
            _dampMat = newdampMat;
            return 0.0;
      }
      
      //Get useful stats for distance
      vector<Vector2d> clusterCentroid; //Vector of cluster centroids
      vector<double> clusterBradius;     //Vector of cluster bounding circles
      clusterStats(clusterFloes, clusterCentroid, clusterBradius); //Get stats for distance
      
      bool simple_y = true; 
      double y_front = 0.0; //YFRONT
      if (simple_y){
        cout << "Use simple y damp front!!!" << endl;
        y_front = find_yfront(clusterFloes, clusterFixed); //YFRONT
        cout << "y damp front: " << y_front << endl;
      }
     
      //2. Update dist_edge based on NO or contact with ice
      //INPUT dist_edge empty
      //OUTPUT dist_edge with zero values, which will be skipped, only work on negative values position relative to edges.
      Vector2d ptGrid;
      double dist_find = 0.0;
      double dist_Change;
      for (size_t i = 0; i < _y_cells; i++) {
            for (size_t j = 0; j < _x_cells; j++) {
                ptGrid = _fluid_coord[j+i*_x_cells];
                // under_iceDEM_bond function (CREATE!!!)
                //under_iceDEM_bond(dist_find, ptGrid, clusterFloes, clusterCentroid, clusterBradius);  //Assume no contact with ice by default, then disprove. (NO FIXED floes) (Left for reference)
                if (simple_y == false){
                    under_iceDEM_bond_fixed(dist_find, ptGrid, clusterFloes, clusterFixed, clusterCentroid, clusterBradius); //FOR fixed floes //radial in 2D
                }
                else{
                    dist_find = under_iceDEM_bond_fixed_yfront(ptGrid, y_front); //FOR fixed floes //linear in y direction, damping underneath front //YFRONT
                }
                dist_Change = dist_find;
                dist_edge[j+i*_x_cells] = dist_Change; 
            }
      }
      
      //Resmooth to avoid holes
      for (size_t i = 0; i < _y_cells; i++) {
          for (size_t j = 0; j < _x_cells; j++) {
              if (  (i > 0 && i <_y_cells-1) && (j > 0 && i <_x_cells-1)  ){
                  if (dist_edge[j+i*_x_cells] == 0 && dist_edge[j+(i+1)*_x_cells] > 0 && dist_edge[j+(i-1)*_x_cells] > 0 && dist_edge[(j+1)+i*_x_cells] > 0 && dist_edge[(j-1)+i*_x_cells] > 0 ){
                      dist_edge[j+i*_x_cells] = 0.25 * ( dist_edge[j+(i+1)*_x_cells] + dist_edge[j+(i-1)*_x_cells] + dist_edge[(j+1)+i*_x_cells] + dist_edge[(j-1)+i*_x_cells] );
                  }
                  if (dist_edge[j+i*_x_cells] == 0 && dist_edge[j+(i+1)*_x_cells] > 0 && dist_edge[j+(i-1)*_x_cells] > 0){
                      dist_edge[j+i*_x_cells] = 0.5 * ( dist_edge[j+(i+1)*_x_cells] + dist_edge[j+(i-1)*_x_cells] );
                  }
                  if (dist_edge[j+i*_x_cells] == 0 && dist_edge[(j+1)+i*_x_cells] > 0 && dist_edge[(j-1)+i*_x_cells] > 0 ){
                      dist_edge[j+i*_x_cells] = 0.5 * ( dist_edge[(j+1)+i*_x_cells] + dist_edge[(j-1)+i*_x_cells] );
                  }
              }
          }
      }

      //3. Find Damping Matrix using Distance Info
      for (size_t i = 0; i < _y_cells; i++) {
            for (size_t j = 0; j < _x_cells; j++) {
                newdampMat[j+i*_x_cells] = 1.0 * exp(-dist_edge[j+i*_x_cells]/Ds); //By defect it will be 1 if D = 0
                ////Debugging print damp matrix
                ////cout << "Damp Matrix value at pos : " << i << " , " << j << " --> " << newdampMat[j+i*_x_cells] << endl;
                ////cout << "Dist edge value at pos : " << i << " , " << j << " --> " << dist_edge[j+i*_x_cells] << endl;
            }
        }
      
      _dampMat = newdampMat;
      
     // cout << "Segments #: " << clusterFloes.size() << endl;
     // //Debugging print cluster floes:
     // for (size_t i = 0; i < clusterFloes.size(); i++){
     //     cout << "Cluster segment number: " << i << endl;
     //     cout << "Included floes #: " << clusterFloes[i].size() <<endl;
     //     for (size_t j = 0; j < clusterFloes[i].size(); j++){
     //         cout << j << " , " << clusterFloes[i][j] << endl;
     //     }
     // }
      //cout << "End Debugging print" << endl;
      //exit(1); //Leave for debugging
      return y_front;
  }
  

  void changeDamp(bool & valueIn){
      _ice_damp = valueIn;
  }
  
    void changedampMat (vector<double> & newdampMat)
    {
        _dampMat = newdampMat;
    }
    
    vector<double> getDampMat(){
        return _dampMat;
    }
    
    //END DAMPING FUNCTIONS
    
    void initBondNumbers(){
        _ntension = 0;
        _nshear = 0;
    }
    
    size_t getntension(){
        return _ntension;
    }
    
    size_t getnshear(){
        return _nshear;
    }
    
    void changeCurrFactor(double & newval){
        _curr_factor = newval;
    }
  
  
  
  
private:
    vector<Grain2d>         _grains;        // vector of grain objects

    vector<Grain2d>         _grainsWall;                // vector of grain objects that simulate land BCs, unmovable, unmutable only constrain grains

  vector<Grain2d>       _wallBottomGrains;    // vector of wall bottom grains
    vector<Wall2d>          _walls;
    Vector2d            _offset;        // offset (longit. dim) of the periodic cell
    double          _dt;          // time increment
    double          _gDamping;        // global damping
    size_t                  _stepup;  //Get time step for time dependent processes
    double                  _slopedir; //Slope for sinusoidal variations (simple)
    double                  _flowangle; //Angle for moving grains
    double                  _flowspeed; //Just speed
    double                  _flowforce;  //Force for moving grains
  GrainState2d      _globalGrainState;    // grain state of entire assembly
  size_t           _ngrains;        // number of grains 
    size_t                  _ngrainsWall;           // number of LAND grains   
    size_t                  _nwalls;
    FracProps2D             _fracProps;
    size_t                  _maxId;  
    vector<Vector2d>        _fluid_coord;           //Constant fluid grid coords
    vector<Vector2d>        _Uwg;                   //Constant fluid grid speed
    size_t                  _x_cells;                //Number of columns for fluid grid 
    size_t                  _y_cells;                //Number of rows for fluid grid
    vector<double>          _oceanTemp;              //GlobalOceanTemp grid

    size_t                  _START_TEMP; //Can also be used as Step if melt Step is Constant
    double                  _THERM_COEFF; //0.200 //0.005 for 100 steps //High speed makes unstable grains
    double                  _MELT_MASS; //0.02 eg for mass //Divided Break and Melt //10000 for point size  //SHOULD IT BE LARGER FOR FASTER SLOPE????
    size_t                  _BREAK_STEP;
    int                     _BREAK_PROB; // 10 for 100 steps
    
    size_t                  _cell_sizex;  //fluid cell size
    size_t                  _cell_sizey;  //fluid cell size
    size_t                  _fluid_mode;  //For different fluid functions
    
    vector<size_t>          _movable_floes;  //For different fixing floes
    Vector2d                _loadvel;  //Load velocity for external load app.
    Vector2d                _loadvelD;  //Load velocity for external load app. //Down velocity control
    
    bool                    _only_dem; //Flip for only DEM
    
    bool                    _ice_damp; //Flip for ice damping
    vector<double>          _dampMat;  //Ocean Damping grid
    
    size_t                  _ntension; //Number of bonds broken by tension
    size_t                  _nshear; //Number of bonds broken by shear
    
    double                  _curr_factor; //Adjust more easily ocean velocity
};

#endif /* World2d_H_ */
