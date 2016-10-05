#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <random>
#include <stdlib.h>
#include <ctime>
#include <cstring>
#include <queue>
#include <signal.h>
#include <unistd.h>

#ifdef USE_GUI
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#endif

#include "tests.h"

#include "data.h"
#include "branchedChain.h"
#include "onlineVariance.h"


//use cmake . in base folder to generate Makefile

void writeState(branchedChain *mol1, std::fstream *fs);
void writeState(branchedChain *mol1, std::fstream *fs){
	*fs << "D:,"; //Write header specify data
	for(int i = 0; i < mol1->numMonomers; i++){
		*fs << i <<","<<mol1->monomers[i].r.vector[0] << "," << mol1->monomers[i].r.vector[1] <<"," << mol1->monomers[i].r.vector[2];
		if(i < mol1->numMonomers -1) {
			*fs << ",";
		}
	}
	*fs << std::endl;
}

void ThermodynamicIntegrate(branchedChain *mol1, const double dnu, const double emax, stateGen *generator, const char *logfile);
void ThermodynamicIntegrate(branchedChain *mol1, const double dnu, const double emax, stateGen *generator, const char *logfile) {
	//Do thermodynamic integration
	double nu = 0.0;  //thermodynamic intergration parameter
	double FreeEnergy = 0.0; //result of the integral
	int nucount = 0;
	double DS = 0.001; //How to take derivative
	std::fstream ls;
	ls.open(logfile, std::ofstream::out);

	std::cout << "CHAINLENGTH: " << CHAINLENGTH << std::endl;
	while(nu <= 1.0+dnu) {
		double pair_energy = 0.0; //Pair energy
		double rog2 = 0.0; //Radius of gyration squared
		//derivative of pair potential
		double dexpect_potential = 0.0;
		double num_potential_points = 0.0;

		mol1->epsilon = nu*emax;
		mol1->runMC(EQUILIBRATION,generator); //run the simulation to EQUILIBRATE
		ls << "PARAMS," << "nu," << nu << ", DS," << DS << ",EPSILON," << mol1->epsilon << ",SPACING," << \
			mol1->spacing << ",INTERACTION DISTANCE," << mol1->sigma << ",SNAPSHOT," << SNAPSHOT << "Energy "<< mol1->energy <<  std::endl;

		ls.flush();
		//Run the simulation and take periodic snapshots
		for(int i = 0; i < 1000; i++ ){
			mol1->runMC(SNAPSHOT,generator);
			double trg2 = mol1->findRg();
			rog2 += trg2; //find radius of gyration
	        //calculate derivative
	        mol1->epsilon = nu*emax+DS;
			double denergy = mol1->totalLJ();
			dexpect_potential += (denergy - mol1->energy)/DS;
			mol1->epsilon = (nu*emax);

			writeState(mol1, &ls); //write the location of the molecules
			ls << "C:," << trg2 << "," << mol1->energy << "," <<  (denergy - mol1->energy)/DS << std::endl;
			ls.flush();

			pair_energy += mol1->energy;
			num_potential_points += 1.0;
    	}

		//Integrate Trapezoid Rule
		double dul = (dexpect_potential/num_potential_points);
		double displayFreeEnergy = 0.0;
		if(nucount == 0) {
			FreeEnergy += dul;
			displayFreeEnergy = 0.0;
		}
		else {
			FreeEnergy += 2*dul;
			displayFreeEnergy = nu*(FreeEnergy-dul)/(2.0*(double)nucount);
		}
		std::cout << nu << ","<< dul << "," << pair_energy/num_potential_points << "," << displayFreeEnergy << "," << rog2/num_potential_points << std::endl;
		nucount += 1;
		nu += dnu; //increment integration parameter
	}
	ls.close();
}


void ThermodynamicIntegrateSig(branchedChain *mol1, const double dnu, const double smax, stateGen *generator, const char *logfile);
void ThermodynamicIntegrateSig(branchedChain *mol1, const double dnu, const double smax, stateGen *generator, const char *logfile) {
	//Do thermodynamic integration
	double nu = 0.0;  //thermodynamic intergration parameter
	double FreeEnergy = 0.0; //result of the integral
	int nucount = 0;
	std::fstream ls;
	ls.open(logfile, std::fstream::out);

	std::cout << "CHAINLENGTH: " << CHAINLENGTH << std::endl;
	while(nu <= 1.0+dnu) {
		double pair_energy = 0.0; //Pair energy
		double rog2 = 0.0; //Radius of gyration squared
		//derivative of pair potential
		double dexpect_potential = 0.0;
		double num_potential_points = 0.0;

		double DS = 0.0005; //How to take derivative

		mol1->sigma = nu*smax;
		mol1->runMC(EQUILIBRATION,generator); //run the simulation to EQUILIBRATE
		ls << "PARAMS," << "nu," << nu << ", DS," << DS << ",EPSILON," << mol1->epsilon << ",SPACING," << \
			mol1->spacing << ",INTERACTION DISTANCE," << mol1->sigma << ",SNAPSHOT," << SNAPSHOT << std::endl;
		ls.flush();
		//Run the simulation and take periodic snapshots
		for(int i = 0; i < 1000; i++ ){
			mol1->runMC(SNAPSHOT,generator);
			double trg2 = mol1->findRg();
			rog2 += trg2; //find radius of gyration
	        //calculate derivative
	        mol1->sigma = nu*smax + DS;
			double denergy = mol1->totalLJ();
			mol1->sigma = nu*smax;
			dexpect_potential += (denergy - mol1->energy)/DS;

			writeState(mol1, &ls); //write the location of the molecules
			ls << "C:," << trg2 << "," << mol1->energy << "," <<  (denergy - mol1->energy)/DS << std::endl;
			ls.flush();

			pair_energy += mol1->energy;
			num_potential_points += 1.0;
    	}

		//Integrate Trapezoid Rule
		double dul = (dexpect_potential/num_potential_points);
		double displayFreeEnergy = 0.0;
		if(nucount == 0) {
			FreeEnergy += dul;
			displayFreeEnergy = 0.0;
		}
		else {
			FreeEnergy += 2*dul;
			displayFreeEnergy = nu*(FreeEnergy-dul)/(2.0*(double)nucount);
		}
		std::cout << nu << ","<< dul << "," << pair_energy/num_potential_points << "," << displayFreeEnergy << "," << rog2/num_potential_points << std::endl;
		nucount += 1;
		nu += dnu; //increment integration parameter
	}
	ls.close();
}

void EnergyRelaxation(branchedChain *mol1, double dnu, const double maxeps, stateGen *generator);
void EnergyRelaxation(branchedChain *mol1, double dnu, const double maxeps, stateGen *generator) {
	double nu = 0.0; //controls the number of iterations in this context
	double data[EQUILIBRATION];
	double num_points = 0.0;

	for(int i = 0; i < EQUILIBRATION; i++) { //initialize
		data[i] = 0.0;
	}

	std::cout << "CHAINLENGTH: " << mol1->numMonomers - mol1->numPhantoms << std::endl;
	mol1->printParams(std::cout);
	while(nu <= 1.0+dnu) {
		std::cerr << "Progress is " << nu << std::endl;
		//reset
		mol1->epsilon = maxeps*0.0; //turn off interactions
		mol1->runMC(CHAINLENGTH*10,generator);

		mol1->epsilon = maxeps*1.0; //turn on interactions
		//Run the simulation and take periodic snapshots
		for(int i = 0; i < EQUILIBRATION; i++ ){
			mol1->runMC(10, generator);
			data[i] += mol1->energy;
    	}
    	num_points += 1.0;
		nu += dnu; //increment integration parameter
	}

	for(int i = 0; i < EQUILIBRATION; i++) { //normalize
		data[i] = data[i]/num_points;
		std::cout << data[i] << std::endl;
	}
}

void ETECorr(branchedChain *mol1, int frames, int numframes, int mi, int mj, const char *fname, stateGen *generator);
void ETECorr(branchedChain *mol1, int frames, int numframes, int mi, int mj, const char *fname, stateGen *generator) {
	mol1->runMC(EQUILIBRATION,generator);  //equiibrate

	std::cout << "CHAINLENGTH: " << CHAINLENGTH << std::endl;

	 std::fstream fs;
	 fs.open (fname, std::fstream::out);

	for(int i = 0; i < numframes; i++){
		mol1->runMC(frames,generator);
		double d2 = diffSquared(&mol1->monomers[mi].r, &mol1->monomers[mj].r);
		jVector res = {{0,0,0}};
		vectorSubP(&res,&mol1->monomers[mj].r, &mol1->monomers[mi].r);
		fs << d2 << "," << mol1->findRg() << "," << res.vector[0] << "," << res.vector[1] << "," << res.vector[2] << std::endl;
	}

	fs.close();
}

//This is not multithreaded.  DO NOT WRITE TO THE SAME fname with simultaneous processes!!!
void JustRun(branchedChain *mol1, const int batches, const int frames, const int snapshot, stateGen *generator, const char *fname);
void JustRun(branchedChain *mol1, const int batches, const int frames, const int snapshot, stateGen *generator, const char *fname) {
	dataFP fp; dataFP pd;
	char bname[500];
	char pdname[500]; char ssname[500];
	char command[1000];

	int b = 0;

	OnlineVar org2;  //calculating variance online
	OnlineVar odesum;
	OnlineVar aveCos; //used for calculation average cos angle

	int edgeCount = mol1->numMonomers - mol1->numPhantoms - 1;
	OnlineVar *ete_dists = new OnlineVar[edgeCount];
	OnlineVar *avg_dist = new OnlineVar[edgeCount];

	mol1->runMC(EQUILIBRATION,generator); //equilibrate the run

	sprintf(pdname,"Avedat.%llu.dat",generator->Seed); //write processed data
	sprintf(ssname,"Avedat.%llu.save",generator->Seed); //write processed data

	//check if save state exists, and load the save file
	if (access( ssname, F_OK ) != -1 ) {
		mol1->loadState(ssname, generator, &b);
		std::cerr << "Loaded resume file, batch " << b << std::endl;
		if (b == batches) {
			std::cerr << "This program needs to do no more work!" << std::endl;
			exit(0);
		}
	}
	//begin running and taking data
	for(b; b < batches; b++){
		std::cerr << "Starting Batch #" << b << std::endl;
		sprintf(bname,"%s_%d_%llu.tf",FRAME_PREFIX,b,generator->Seed);
		fp = fopen(bname,"wb"); //open a block data file
		if(ferror(fp)){
			std::cerr << "Coult not open file for writing: " << bname << std::endl;
			exit(0);
		}
		mol1->writeHeader(fp,generator);
		mol1->writeParameterLine(fp);
		#ifndef LOG_STATE
		mol1->writeDataLine(fp,generator);
		#endif 

		LogLine ll; memset((void *)&ll, 0, sizeof(ll));
		org2.setState(0,0,0);  //reset the rg counters
		odesum.setState(0,0,0);
		aveCos.setState(0,0,0);
		for(int i = 0; i < edgeCount; i++){
			ete_dists[i].setState(0,0,0);
			avg_dist[i].setState(0,0,0);
		}

		for(int f = 0; f < frames; f++) { //run this many times
			mol1->runMC(snapshot,generator);
			#ifdef LOG_STATE
			if(frames%LOG_STATE == 0) { //Log the state every so often
				mol1->writeDataLine(fp,generator);
			}
			#endif
			double _trg2 = mol1->findRg();
			ll.rog2_sum += _trg2; org2.addValue(_trg2);
			double _tde = mol1->energyDerivative(DE);
			ll.de_sum += _tde; odesum.addValue(_tde);
			mol1->edgeLength(ete_dists, avg_dist); //find edge length and position

			jVector l1,l2; double tbe; //find average angle
			vectorSubP(&l1, &mol1->monomers[0].r, &mol1->monomers[2].r);
			vectorSubP(&l2, &mol1->monomers[1].r, &mol1->monomers[2].r);
			tbe = std::abs(dotP(&l1, &l2));
			aveCos.addValue(tbe);

			ll.rmax = DENSITY_MAX_DOMAIN;
			mol1->findDensity(ll.density_sum, DENSITY_BINS,DENSITY_MAX_DOMAIN);
			
#ifdef FREE_ENERGY_DENSITY
			double _tde2 = mol1->binParticles(DE);
			if( abs(_tde - _tde2 ) > 10e-6) {
				std::cerr <<"Warning:  There is error in the derivative calculation of size: " << abs(_tde - _tde2) << std::endl;
				exit(0);
			}
			//Now average to online variables
#endif

			ll.N++;
		}
		fclose(fp);
		std::cerr << "Rg2 is " << ll.rog2_sum/(double)ll.N << " +- " << sqrt(org2.Variance()/(double)org2.n)  << " <cos> " << aveCos.mean << std::endl;
		//std::cerr << "End to end distance is " << ete_dist.mean << std::endl;
		std::cerr << "Segment averages are ";
		for(int i = 0; i < edgeCount; i++){
			std::cerr << ete_dists[i].mean << ",";
		}

		ll.epsilon = mol1->epsilon;
		ll.spacing = mol1->spacing;
		ll.sigma = mol1->sigma;
		//online variance data.  Save so we can calculate full fluctuations with parallel algorithm
		ll.rg2_n = org2.n; ll.rg2_mean = org2.mean; ll.rg2_m2 = org2.m2;
		ll.de_n = odesum.n; ll.de_mean = odesum.mean; ll.de_m2 = odesum.m2;

		pd = fopen(pdname,"ab");
		if(ferror(pd)){
			std::cerr << "Could not open " << pdname << " for writing!" << std::endl;
			exit(0);
		}

		fwrite((void *)&ll.N,  sizeof(ll.N), (size_t)1, pd); //write Log Line
		fwrite((void *)&ll.rog2_sum,  sizeof(ll.rog2_sum), (size_t)1, pd); //write Log Line
		fwrite((void *)&ll.de_sum,  sizeof(ll.de_sum), (size_t)1, pd); //write Log Line
		fwrite((void *)&ll.rmax,  sizeof(ll.rmax), (size_t)1, pd); //write Log Line
		fwrite((void *)ll.density_sum,  sizeof(ll.density_sum[0]), (size_t)DENSITY_BINS, pd); //write Log Line
		fwrite((void *)&ll.spacing,  sizeof(ll.spacing), (size_t)1, pd); //write Log Line
		fwrite((void *)&ll.sigma,  sizeof(ll.sigma), (size_t)1, pd); //write Log Line
		fwrite((void *)&ll.epsilon,  sizeof(ll.epsilon), (size_t)1, pd); //write Log Line
		fwrite((void *)&ll.rg2_n,  sizeof(ll.rg2_n), (size_t)1, pd); //write Log Line
		fwrite((void *)&ll.rg2_mean,  sizeof(ll.rg2_mean), (size_t)1, pd); //write Log Line
		fwrite((void *)&ll.rg2_m2,  sizeof(ll.rg2_m2), (size_t)1, pd); //write Log Line
		fwrite((void *)&ll.de_n,  sizeof(ll.de_n), (size_t)1, pd); //write Log Line
		fwrite((void *)&ll.de_mean,  sizeof(ll.de_mean), (size_t)1, pd); //write Log Line
		fwrite((void *)&ll.de_m2,  sizeof(ll.de_m2), (size_t)1, pd); //write Log Line
		fwrite((void *)&edgeCount, sizeof(edgeCount), (size_t)1, pd); //write number of non-phantom edges
		for(int i = 0; i < edgeCount; i++) { //write edge lengths
			fwrite((void *)&ete_dists[i].n,  sizeof(ll.de_n), (size_t)1, pd); //write Log Line
			fwrite((void *)&ete_dists[i].mean,  sizeof(ll.de_mean), (size_t)1, pd); //write Log Line
			fwrite((void *)&ete_dists[i].m2,  sizeof(ll.de_m2), (size_t)1, pd); //write Log Line
		}
		for(int i = 0; i < edgeCount; i++) { //write edge distances
			fwrite((void *)&avg_dist[i].n,  sizeof(ll.de_n), (size_t)1, pd); //write Log Line
			fwrite((void *)&avg_dist[i].mean,  sizeof(ll.de_mean), (size_t)1, pd); //write Log Line
			fwrite((void *)&avg_dist[i].m2,  sizeof(ll.de_m2), (size_t)1, pd); //write Log Line
		}

		fflush(pd);
		fclose(pd);

		//saveState here
		int _b = b + 1;
		mol1->saveState(ssname,generator,&_b);  //this is for resuming incomplete batches (must run from beginning)

		sprintf(command,"zip -rv %s.zip %s",fname,bname); 		system(command); //save block to zip
		sprintf(command,"rm %s",bname);		system(command);  //remove file
		std::cout << command << std::endl;
	}
	delete[] ete_dists;
	delete[] avg_dist;
}

//global timer molecule
branchedChain *timerMol;

void sigSpeed(int sig);
void sigSpeed(int sig){
	timerMol->speedFlag = 1;
}

int main(){

#ifdef RUN_TESTS
	testCenterOfMass();
	exit(0);
#endif
	//initialize randomizer
	//std::default_random_engine generator;
	stateGen generator; generator.seed((uint64_t)INIT_SEED);
	std::uniform_real_distribution<double> azDistribution(0,M_PI*2.0); //azimuthal angle
	std::uniform_real_distribution<double> thetaDistribution(-1.0,1.0); //cos theta
	std::uniform_int_distribution<int> 	intDistribution(0, CHAINLENGTH - 1);
    std::uniform_real_distribution<double> accdist(0.0,1.0); //metropolis
    std::uniform_int_distribution<int> latticeRot(0,3);

    //Write Header info
    std::cerr << "Initial seed " << (uint64_t)INIT_SEED << std::endl;
    std::cerr << "Zip Name " << ZIPNAME << std::endl;
    std::cerr << "Batches " << BATCHES << std::endl;
    std::cerr << "Lines " << LINES << std::endl;
    std::cerr << "Epsilon " << MOL_EPSILON << std::endl;
    std::cerr << "Sigma " << MOL_SIGMA << std::endl;
    std::cerr << "Spacing " << MOL_SPACING << std::endl;

	//Build chain
	branchedChain mol1 = branchedChain(CHAINLENGTH);
	//set parameters
	mol1.sigma = MOL_SIGMA;
	mol1.epsilon = MOL_EPSILON;
	mol1.spacing = MOL_SPACING;
	mol1.beta = MOL_BETA;
	timerMol = &mol1;
	signal(SIGALRM, &sigSpeed ); alarm(INDICATOR_PERIOD);
#ifdef LINEAR
	mol1.createLinear(CHAINLENGTH - PHANTOMS);
	mol1.setLinearPositions(MOL_SPACING);
	if(PHANTOMS > 0){	mol1.insertPhantoms(CHAINLENGTH-PHANTOMS); }
	std::cerr << "Linear: Num PHANTOMS " << mol1.numPhantoms << std::endl;
#endif
#ifdef ARBITRARY
	std::cerr <<" Loading " << ADJ_FILE << std::endl;
	mol1.loadAdjacency(CHAINLENGTH - PHANTOMS, ADJ_FILE, &generator);
	if(PHANTOMS > 0){	mol1.insertPhantoms(CHAINLENGTH-PHANTOMS); }
#endif
#ifdef DENDRIMER
	int mygen = 1;
	std::cout << "Dendrimer Mass " << mol1.dendrimerMass(3,mygen) <<std::endl;
	mol1.createDendrimerR(CHAINLENGTH-PHANTOMS, 3,mygen,&generator);
	if(PHANTOMS > 0){	mol1.insertPhantoms(CHAINLENGTH-PHANTOMS); }
#endif
	mol1.saveTopofile(TOPOLOGY_FILE);
	mol1.enableDerivativeDensity(100.0, 20.0);
	clock_t b,e;
	b = clock();

	//ThermodynamicIntegrateSig(&mol1,0.05,0.9,&generator,TI_FILE);  //38kt
	//ThermodynamicIntegrate(&mol1,0.05,1.0,&generator,TI_FILE);
#ifdef MODE_RELAX
	EnergyRelaxation(&mol1,1.0/20.0,5.0,&generator);
#endif
	//mol1.nu = 1.0;
	//ETECorr(&mol1,10000,10,0,CHAINLENGTH-1,"corr94.csv",&generator);
#ifdef MODE_RUN
	JustRun(&mol1, BATCHES,LINES, DECORR_TIME, &generator, ZIPNAME);
#endif

	e = clock();
	//std::cout << "Clock ticks was " << e - b << std::endl;
	std::cout << "Arc Length " << mol1.arcLength(0,0) << std::endl;

#ifdef USE_GUI

	// Create a polydata object
	vtkSmartPointer<vtkPolyData> point = showMonomers(&mol1);//randomPointsOnSphere(1000);//pointsOnSphere(10,100);
 
  	// Visualize
  	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	#if VTK_MAJOR_VERSION <= 5
  		mapper->SetInput(point);
	#else
  		mapper->SetInputData(point);
	#endif
 
  	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
 	actor->SetMapper(mapper);
  	actor->GetProperty()->SetPointSize(15);
 
  	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  	renderWindow->AddRenderer(renderer);
  	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = 
    			vtkSmartPointer<vtkRenderWindowInteractor>::New();
  	renderWindowInteractor->SetRenderWindow(renderWindow);
 
  	renderer->AddActor(actor);
 
  	renderWindow->Render();
  	renderWindowInteractor->Start();

 #endif
 
  	return EXIT_SUCCESS;
}
