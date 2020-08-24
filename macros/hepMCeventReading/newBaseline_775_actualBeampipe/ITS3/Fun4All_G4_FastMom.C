#pragma once
#include <phgenfit/Track.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <g4detectors/PHG4DetectorSubsystem.h>
#include <g4detectors/PHG4CylinderSubsystem.h>

#include <g4histos/G4HitNtuple.h>
#include <gdmlimporter/GdmlImportDetectorSubsystem.h>
#include <g4main/PHG4ParticleGenerator.h>
#include <g4main/PHG4ParticleGeneratorBase.h>
#include <g4main/PHG4ParticleGun.h>
#include <g4main/PHG4Reco.h>
#include <g4main/PHG4SimpleEventGenerator.h>
#include <g4main/PHG4TruthSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4trackfastsim/PHG4TrackFastSimEval.h>
#include <phool/recoConsts.h>
#include <g4detectors/PHG4GDMLSubsystem.h>

#include <g4main/HepMCNodeReader.h>
#include <phhepmc/Fun4AllHepMCInputManager.h>

#include <gdmlimporter/GdmlImportDetectorSubsystem.h>
#include <gdmlimporter/SimpleNtuple.h>
#include <gdmlimporter/TrackFastSimEval.h>

#include <charmanalysis/CharmAnalysis.h>


#include "detector_setup.h"


R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libgdmlimportdetector.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)
R__LOAD_LIBRARY(libcharmanalysis.so)

void Fun4All_G4_FastMom(
			int nEvents = -1,
			const char *outputFile = "out_allSi",
			const char *genpar = "e-")
{
	// ======================================================================================================
	// Input from the user
	const int particle_gen = 1;	// 1 = particle generator, 2 = particle gun
	const int magnetic_field = 1;	// 1 = uniform 1.5T, 2 = uniform 3.0T, 3 = sPHENIX 1.4T map, 4 = Beast 3.0T map
	// ======================================================================================================
	// Parameters for black holes
	string radialBH   = "DIRC";		// Cylindrical surface object name
	double radialBHradius = 80.0;		// [cm] 
	double radialBHlength     = 400.0;		// [cm]
	// ---
	double thickness    = 0.1;		// black hole thickness, needs to be taken into account for the z positions
	// ---
	string forwardBH   = "FOR"; 		// Forward plane object name
	double forwardBHz   = 130.0+thickness/2.0;	// [cm]
	double forwardBHradius = 80.0;		// [cm]
	// ---
	string backwardBH   = "BACK";		// Backward plane object name
	double backwardBHz   = -(130.0+thickness/2.0);// [cm]
	double backwardBHradius = 80.0;		// [cm]
	// ======================================================================================================
	// Make the Server
	Fun4AllServer *se = Fun4AllServer::instance();
	se->Verbosity(1);
	// If you want to fix the random seed for reproducibility
	recoConsts *rc = recoConsts::instance();
	rc->set_IntFlag("RANDOMSEED", 12345);
	// ======================================================================================================
	// Particle Generation
	cout << "Particle that will be generated: " << std::string(genpar) << endl;
	// --------------------------------------------------------------------------------------
	// Particle Generator Setup
	PHG4ParticleGenerator *gen = new PHG4ParticleGenerator();
	gen->set_name(std::string(genpar));	// geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ... (currently passed as an input)
	gen->set_vtx(0,0,0);			// Vertex generation range
	gen->set_mom_range(0.0, 30.0);		// Momentum generation range in GeV/c
	gen->set_z_range(0.0, 0.0); //THIS NEEDS TO BE SET, TO NOT DEFAULT TO -10 TO 10
	gen->set_eta_range(0.0, 1.0);		// Detector coverage corresponds to |η|< 4
	gen->set_phi_range(0.0, 2.0*TMath::Pi());
	
	// --------------------------------------------------------------------------------------
	// Particle Gun Setup
	PHG4ParticleGun *gun = new PHG4ParticleGun();
	gun->set_name(std::string(genpar));	// geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ...
	gun->set_vtx(0,0,0);
	gun->set_mom(0,1,0);
	// --------------------------------------------------------------------------------------
	//if     (particle_gen==1){se->registerSubsystem(gen); cout << "Using particle generator" << endl;}
	//else if(particle_gen==2){se->registerSubsystem(gun); cout << "Using particle gun"       << endl;}
	

	// ======================================================================================================
	PHG4Reco *g4Reco = new PHG4Reco();
	// ======================================================================================================
	// Magnetic field setting
	TString B_label;
	if(magnetic_field==1){		// uniform 1.5T
		B_label = "_B_1.5T";
		g4Reco->set_field(1.5);
	}
	else if(magnetic_field==2){	// uniform 3.0T
		B_label = "_B_3.0T";
		g4Reco->set_field(3.0);
	}
	else if(magnetic_field==3){	// sPHENIX 1.4T map
		B_label = "_sPHENIX";
		g4Reco->set_field_map(string(getenv("CALIBRATIONROOT")) + string("/Field/Map/sPHENIX.2d.root"), PHFieldConfig::kField2D);
		g4Reco->set_field_rescale(-1.4/1.5);
	}
	else if(magnetic_field==4){	// Beast 3.0T map
		B_label = "_Beast";
		g4Reco->set_field_map(string(getenv("CALIBRATIONROOT")) + string("/Field/Map/mfield.4col.dat"), PHFieldConfig::kFieldBeast);
	}
	else{				// The user did not provide a valid B field setting
		cout << "User did not provide a valid magnetic field setting. Set 'magnetic_field'. Bailing out!" << endl;
	}
	// ======================================================================================================
	// Physics list (default list is "QGSP_BERT")
	//g4Reco->SetPhysicsList("FTFP_BERT_HP"); // This list is slower and only useful for hadronic showers.
	// ======================================================================================================
	
	#ifdef _BEAMPIPE_
	// Central pipe dimension
	// Extracted via mechanical model: Detector chamber 3-20-20
	// directly implimenting the central Be section in G4 cylinder for max speed simulation in the detector region.
	// The jointer lip structure of the pipe R = 3.2cm x L=5mm is ignored here
	const double be_pipe_radius = 3.1000;
	const double be_pipe_thickness = 3.1762 - be_pipe_radius;  // 760 um for sPHENIX
	const double be_pipe_length_plus = 66.8;                   // +z beam pipe extend.
	const double be_pipe_length_neg = -79.8;                   // -z beam pipe extend.
	bool overlapcheck = false;
	bool absorberactive = false;
	
	// process pipe extentions?
	bool use_forward_pipes = true;
	const bool do_pipe_hadron_forward_extension = use_forward_pipes && true;
	const bool do_pipe_electron_forward_extension = use_forward_pipes && true;

	const double be_pipe_length = be_pipe_length_plus - be_pipe_length_neg;  // pipe length
	const double be_pipe_center = 0.5 * (be_pipe_length_plus + be_pipe_length_neg);
	
	// mid-rapidity beryillium pipe
	PHG4CylinderSubsystem* cyl = new PHG4CylinderSubsystem("VAC_BE_PIPE", 0);
	cyl->set_double_param("radius", 0.0);
	cyl->set_int_param("lengthviarapidity", 0);
	cyl->set_double_param("length", be_pipe_length);
	cyl->set_double_param("place_z", be_pipe_center);
	cyl->set_string_param("material", "G4_Galactic");
	cyl->set_double_param("thickness", be_pipe_radius);
	cyl->SuperDetector("PIPE");
	cyl->OverlapCheck(overlapcheck);
	if (absorberactive) cyl->SetActive();
	g4Reco->registerSubsystem(cyl);

	cyl = new PHG4CylinderSubsystem("BE_PIPE", 1);
	cyl->set_double_param("radius", be_pipe_radius);
	cyl->set_int_param("lengthviarapidity", 0);
	cyl->set_double_param("length", be_pipe_length);
	cyl->set_double_param("place_z", be_pipe_center);
	cyl->set_string_param("material", "G4_Be");
	cyl->set_double_param("thickness", be_pipe_thickness);
	cyl->SuperDetector("PIPE");
	cyl->OverlapCheck(overlapcheck);
	if (absorberactive) cyl->SetActive();
	g4Reco->registerSubsystem(cyl);
	
	if (do_pipe_electron_forward_extension)
	{
	PHG4GDMLSubsystem* gdml = new PHG4GDMLSubsystem("ElectronForwardEnvelope");
	gdml->set_string_param("GDMPath", string(getenv("CALIBRATIONROOT")) + "/Beam/Detector chamber 3-20-20.G4Import.gdml");
	gdml->set_string_param("TopVolName", "ElectronForwardEnvelope");
	gdml->set_int_param("skip_DST_geometry_export", 1);  // do not export extended beam pipe as it is not supported by TGeo and outside Kalman filter acceptance
	gdml->OverlapCheck(1);
	g4Reco->registerSubsystem(gdml);
	}

	if (do_pipe_hadron_forward_extension)
	{
	PHG4GDMLSubsystem* gdml = new PHG4GDMLSubsystem("HadronForwardEnvelope");
	gdml->set_string_param("GDMPath", string(getenv("CALIBRATIONROOT")) + "/Beam/Detector chamber 3-20-20.G4Import.gdml");
	gdml->set_string_param("TopVolName", "HadronForwardEnvelope");
	gdml->set_int_param("skip_DST_geometry_export", 1);  // do not export extended beam pipe as it is not supported by TGeo and outside Kalman filter acceptance
	gdml->OverlapCheck(1);
	g4Reco->registerSubsystem(gdml);
	}
	#endif
	
	
	#ifdef _BARREL_
	// Loading barrel from gdml file
	GdmlImportDetectorSubsystem * svtPart = new GdmlImportDetectorSubsystem();
	svtPart->set_string_param("GDMPath","Vst_GDML_.gdml"); 
	//svtPart->OverlapCheck(); //Doesn't do anything, it seems.
	svtPart->AddAssemblyVolume("VST");	// Barrel. ********** THIS IS THE NAME IN THE GDML FILE, FROM EICROOT ********
	svtPart->SuperDetector("SVT");
	svtPart->SetActive();          // this saves hits in the MimosaCore volumes
	svtPart->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)
	g4Reco->registerSubsystem(svtPart);
	#endif
	
	#ifdef _TIMING_
	//Time-stamping layer
	GdmlImportDetectorSubsystem * timePart = new GdmlImportDetectorSubsystem();
	timePart->set_string_param("GDMPath","Timing_GDML_.gdml"); 
	//timePart->OverlapCheck(); //Doesn't do anything, it seems.
	timePart->AddAssemblyVolume("TIMING");	// Barrel. ********** THIS IS THE NAME IN THE GDML FILE, FROM EICROOT ********
	timePart->SuperDetector("TimeStamping");
	timePart->SetActive();          // this saves hits in the MimosaCore volumes
	timePart->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)
	g4Reco->registerSubsystem(timePart);
	#endif
	
	#ifdef _SITPC_
	//Silicon TPC replacement
	GdmlImportDetectorSubsystem * siTPCPart = new GdmlImportDetectorSubsystem();
	siTPCPart->set_string_param("GDMPath","Sitpc_GDML_775.gdml"); 
	//siTPCPart->OverlapCheck(); //Doesn't do anything, it seems.
	siTPCPart->AddAssemblyVolume("SITPC");	// Barrel. ********** THIS IS THE NAME IN THE GDML FILE, FROM EICROOT ********
	siTPCPart->SuperDetector("SiTPCReplacement");
	siTPCPart->SetActive();          // this saves hits in the MimosaCore volumes
	siTPCPart->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)
	g4Reco->registerSubsystem(siTPCPart);
	#endif
	
	#ifdef _FORWARD_SILICON_DISKS_
	//Forward disks
	GdmlImportDetectorSubsystem * fstPart = new GdmlImportDetectorSubsystem();
	fstPart->set_string_param("GDMPath","Fst_GDML_775.gdml"); 
	//fstPart->OverlapCheck(); //Doesn't do anything, it seems.
	fstPart->AddAssemblyVolume("FST");	// Disks. ********** THIS IS THE NAME IN THE GDML FILE, FROM EICROOT ********
	fstPart->SuperDetector("FstDisks");
	fstPart->SetActive();          // this saves hits in the MimosaCore volumes
	fstPart->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)
	g4Reco->registerSubsystem(fstPart);
	#endif
	
	#ifdef _FORWARD_SILICON_DISKS_
	//Backward disks
	GdmlImportDetectorSubsystem * bstPart = new GdmlImportDetectorSubsystem();
	bstPart->set_string_param("GDMPath","Bst_GDML_775.gdml"); 
	//bstPart->OverlapCheck(); //Doesn't do anything, it seems.
	bstPart->AddAssemblyVolume("BST");	// Disks. ********** THIS IS THE NAME IN THE GDML FILE, FROM EICROOT ********
	bstPart->SuperDetector("BstDisks");
	bstPart->SetActive();          // this saves hits in the MimosaCore volumes
	bstPart->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)
	g4Reco->registerSubsystem(bstPart);
	#endif
	
	//g4Reco->SetWorldMaterial("G4_Galactic"); //For material scans as in EICROOT

	#ifdef _BLACKHOLES_
	// ======================================================================================================	
	//Setting up black holes in the forward, backward, and radial directions
	PHG4CylinderSubsystem * bhCyl;
	bhCyl = new PHG4CylinderSubsystem(radialBH,0);
	bhCyl->set_double_param("length", radialBHlength);
	bhCyl->set_double_param("radius", radialBHradius); // dirc radius
	bhCyl->set_double_param("thickness", 0.1); // needs some thickness
	bhCyl->set_string_param("material", "G4_AIR");
	bhCyl->SetActive(1);
	bhCyl->SuperDetector(radialBH);
	bhCyl->BlackHole();
	bhCyl->set_color(0,1,1,0.3); //blue
	g4Reco->registerSubsystem(bhCyl);
	
	bhCyl = new PHG4CylinderSubsystem(forwardBH,0);
	bhCyl->set_double_param("length", thickness);
	bhCyl->set_double_param("radius", 0); // 
	bhCyl->set_double_param("thickness", forwardBHradius); // 
	bhCyl->set_string_param("material", "G4_AIR");
	bhCyl->set_double_param("place_z", forwardBHz);
	bhCyl->SetActive(1);
	bhCyl->SuperDetector(forwardBH);
	bhCyl->BlackHole();
	bhCyl->set_color(0,1,1,0.3); //blue
	g4Reco->registerSubsystem(bhCyl);
	
	bhCyl = new PHG4CylinderSubsystem(backwardBH,0);
	bhCyl->set_double_param("length", thickness);
	bhCyl->set_double_param("radius", 0); // 
	bhCyl->set_double_param("thickness", backwardBHradius); // 
	bhCyl->set_string_param("material", "G4_AIR");
	bhCyl->set_double_param("place_z", backwardBHz);
	bhCyl->SetActive(1);
	bhCyl->SuperDetector(backwardBH);
	bhCyl->BlackHole();
	bhCyl->set_color(0,1,1,0.3); //blue
	g4Reco->registerSubsystem(bhCyl);
	#endif
	// ======================================================================================================
	
	//For HepMC
	HepMCNodeReader * hr = new HepMCNodeReader();
	se->registerSubsystem(hr);
	
	//Register the truth hits
	PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
	g4Reco->registerSubsystem(truth);

	g4Reco->Verbosity(0);
	se->registerSubsystem(g4Reco);


	// ======================================================================================================
	//Reconstruction
	// Fast pattern recognition and full Kalman filter
	// output evaluation file for truth track and reco tracks are PHG4TruthInfoContainer
	PHG4TrackFastSim *kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
	kalman->set_use_vertex_in_fitting(false);
	kalman->set_sub_top_node_name("SVTX");
	kalman->set_trackmap_out_name("SvtxTrackMap");
	
	//10: SVT
	//20: Time-stamping layer(s)
	//30: Silicon TPC
	//40: Forward disks
	//50: Backward disks
	double pixelResolution = 10e-4/sqrt(12);
	std::cout << pixelResolution << std::endl;
	#ifdef _BARREL_
	//SVT
	//This adds the hit nodes for EACH LAYER OF THE DETECTOR.
	for (int i = 10; i < 10+_NO_OF_BARREL_LAYERS_; i++) {
		std::string nodeName = "G4HIT_SVT_" + std::to_string(i); // ****** THIS NEEDS TO MATCH THE NAME FROM GetName() in GdmlImportDetector.cc ******
		kalman->add_phg4hits(
				nodeName,				// const std::string& phg4hitsNames
				PHG4TrackFastSim::Cylinder,		// const DETECTOR_TYPE phg4dettype
				999.,					// radial-resolution [cm] (this number is not used in cylindrical geometry)
				pixelResolution,					// azimuthal (arc-length) resolution [cm]
				pixelResolution,					// longitudinal (z) resolution [cm]
				1,					// efficiency (fraction)
				0					// hit noise
		);
	}
	#endif
	
	#ifdef _TIMING_
	//Time-stamping
	for (int i = 20; i < 20+_NO_OF_TIMING_LAYERS_; i++) {
		std::string nodeName = "G4HIT_TimeStamping_" + std::to_string(i); // ****** THIS NEEDS TO MATCH THE NAME FROM GetName() in GdmlImportDetector.cc ******
		kalman->add_phg4hits(
				nodeName,				// const std::string& phg4hitsNames
				PHG4TrackFastSim::Cylinder,		// const DETECTOR_TYPE phg4dettype
				999.,					// radial-resolution [cm] (this number is not used in cylindrical geometry)
				pixelResolution,					// azimuthal (arc-length) resolution [cm]
				pixelResolution,					// longitudinal (z) resolution [cm]
				1,					// efficiency (fraction)
				0					// hit noise
		);
	}
	#endif
	
	#ifdef _SITPC_
	//Silicon TPC replacement
	for (int i = 30; i < 30+_NO_OF_SITPC_LAYERS_; i++) {
		std::string nodeName = "G4HIT_SiTPCReplacement_" + std::to_string(i); // ****** THIS NEEDS TO MATCH THE NAME FROM GetName() in GdmlImportDetector.cc ******
																			  //Comes from SuperDetector name when adding it to Geant4.
		kalman->add_phg4hits(
				nodeName,				// const std::string& phg4hitsNames
				PHG4TrackFastSim::Cylinder,		// const DETECTOR_TYPE phg4dettype
				999.,					// radial-resolution [cm] (this number is not used in cylindrical geometry)
				pixelResolution,					// azimuthal (arc-length) resolution [cm]
				pixelResolution,					// longitudinal (z) resolution [cm]
				1,					// efficiency (fraction)
				0					// hit noise
		);
	}
	#endif
	
	#ifdef _FORWARD_SILICON_DISKS_
	//Forward silicon disks
	for (int i = 40; i < 40+_NO_OF_FORWARD_DISKS_ ; i++) {
		std::string nodeName = "G4HIT_FstDisks_" + std::to_string(i); // ****** THIS NEEDS TO MATCH THE NAME FROM GetName() in GdmlImportDetector.cc ******
																			  //Comes from SuperDetector name when adding it to Geant4.
		kalman->add_phg4hits(
				nodeName,				// const std::string& phg4hitsNames
				PHG4TrackFastSim::Vertical_Plane,		// const DETECTOR_TYPE phg4dettype
				pixelResolution,					// radial-resolution [cm] 
				pixelResolution,					// azimuthal (arc-length) resolution [cm]
				999.0,					// longitudinal (z) resolution [cm](not used for vertical plane)
				1,					// efficiency (fraction)
				0					// hit noise
		);
	}
	#endif
	#ifdef _BACKWARD_SILICON_DISKS_
	//Backward silicon disks
	for (int i = 50; i < 50+_NO_OF_BACKWARD_DISKS_ ; i++) {
		std::string nodeName = "G4HIT_BstDisks_" + std::to_string(i); // ****** THIS NEEDS TO MATCH THE NAME FROM GetName() in GdmlImportDetector.cc ******
																			  //Comes from SuperDetector name when adding it to Geant4.
		kalman->add_phg4hits(
				nodeName,				// const std::string& phg4hitsNames
				PHG4TrackFastSim::Vertical_Plane,		// const DETECTOR_TYPE phg4dettype
				pixelResolution,					// radial-resolution [cm] 
				pixelResolution,					// azimuthal (arc-length) resolution [cm]
				999.0,					// longitudinal (z) resolution [cm] (not used for vertical plane)
				1,					// efficiency (fraction)
				0					// hit noise
		);
	}
	#endif
	
	
	
	se->registerSubsystem(kalman);
	
	
	

	// ======================================================================================================
	TrackFastSimEval *fast_sim_eval = new TrackFastSimEval("FastTrackingEval");
	fast_sim_eval->set_filename(TString(outputFile)+B_label+"_FastTrackingEval.root");
	//fast_sim_eval->AddProjection(radialBH);
	//fast_sim_eval->AddProjection(forwardBH);
	//fast_sim_eval->AddProjection(backwardBH);
	se->registerSubsystem(fast_sim_eval);




	SimpleNtuple * hits = new SimpleNtuple("Hits");
	#ifdef _BARREL_
	//SVT starts at 10, to allow for absorbers before
	//hits->AddNode("ABSORBER_SVT",0); // hits in the passive volumes
	for (int i = 10; i < 10+_NO_OF_BARREL_LAYERS_; i++) { // hits in the  MimosaCore volumes
		std::string nodeName = "SVT_" + std::to_string(i); // ****** THIS NEEDS TO MATCH WHAT'S INPUT IN GdmlImportDetector.cc, stripped of the G4HIT, because currently SimpleNtuple adds that.
														   // See SimpleNtuple.cc ******
		hits->AddNode(nodeName, i);
	}
	#endif
	#ifdef _TIMING_
	//Time-stamping
	//hits->AddNode("ABSORBER_TimeStamping",1); // hits in the passive volumes
	for (int i = 20; i < 20+_NO_OF_TIMING_LAYERS_; i++) { // hits in the  MimosaCore volumes
		std::string nodeName = "TimeStamping_" + std::to_string(i); // ****** THIS NEEDS TO MATCH WHAT'S INPUT IN GdmlImportDetector.cc, stripped of the G4HIT, because currently SimpleNtuple adds that.
														   // See SimpleNtuple.cc ******
		hits->AddNode(nodeName, i);
	}
	#endif
	#ifdef _SITPC_
	//Silicon TPC replacement
	//hits->AddNode("ABSORBER_SiTPCReplacement",2); // hits in the passive volumes
	for (int i = 30; i < 30+_NO_OF_SITPC_LAYERS_; i++) { // hits in the  MimosaCore volumes
		std::string nodeName = "SiTPCReplacement_" + std::to_string(i); // ****** THIS NEEDS TO MATCH WHAT'S INPUT IN GdmlImportDetector.cc, stripped of the G4HIT, because currently SimpleNtuple adds that.
														   // See SimpleNtuple.cc ******
		hits->AddNode(nodeName, i);
	}
	#endif
	#ifdef _FORWARD_SILICON_DISKS_
	//Forward silicon disks
	//hits->AddNode("ABSORBER_FstDisks",3); // hits in the passive volumes
	for (int i = 40; i < 40+_NO_OF_FORWARD_DISKS_; i++) { // hits in the  MimosaCore volumes
		std::string nodeName = "FstDisks_" + std::to_string(i); // ****** THIS NEEDS TO MATCH WHAT'S INPUT IN GdmlImportDetector.cc, stripped of the G4HIT, because currently SimpleNtuple adds that.
														   // See SimpleNtuple.cc ******
		hits->AddNode(nodeName, i);
	}
	#endif
	#ifdef _BACKWARD_SILICON_DISKS_
	//Backward silicon disks
	//hits->AddNode("ABSORBER_BstDisks",4); // hits in the passive volumes
	for (int i = 50; i < 50+_NO_OF_FORWARD_DISKS_; i++) { // hits in the  MimosaCore volumes
		std::string nodeName = "BstDisks_" + std::to_string(i); // ****** THIS NEEDS TO MATCH WHAT'S INPUT IN GdmlImportDetector.cc, stripped of the G4HIT, because currently SimpleNtuple adds that.
														   // See SimpleNtuple.cc ******
		hits->AddNode(nodeName, i);
	}
	#endif
	
	se->registerSubsystem(hits);


	//The charm analysis
	CharmAnalysis * charmAnalysis = new CharmAnalysis("charmAnalysis", string(outputFile) + "_charmAnalysis.root");
	charmAnalysis->setMinJetPt(10.);
	charmAnalysis->Verbosity(100);
	charmAnalysis->analyzeTracks(false);
	charmAnalysis->analyzeClusters(false);
	charmAnalysis->analyzeJets(false);
	charmAnalysis->analyzeTruth(true);
	//Let's skip the fancier stuff for now
	//se->registerSubsystem(charmAnalysis);

	///////////////////////////////////////////
	// IOManagers...
	///////////////////////////////////////////
	const std::string dst_name = std::string(outputFile)+"_gdmlimporter.root";
	//Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT",TString(outputFile)+"_gdmlimporter.root");
	Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT",dst_name);
	out->Verbosity(10);
	//Let's not output right now
	//se->registerOutputManager(out);
	
	Fun4AllInputManager *in2 = new Fun4AllHepMCInputManager("DSTIN");
	in2->Verbosity(10);
	se->registerInputManager(in2);
	//in2->fileopen(in2->Name().c_str(), inputFile); //Need to add both file name and the name of the input manager to bind it to!
	se->fileopen(in2->Name().c_str(), "../../eventFiles/ggm2ccbar_28_1_100.dat"); //Need to add both file name and the name of the input manager to bind it to!

	//Fun4AllInputManager *in = new Fun4AllDummyInputManager("JADE");
	//se->registerInputManager(in);
	if (nEvents <= 0)
	{
		return;
	}
	se->run(nEvents);
	se->End();
	delete se;
	gSystem->Exit(0);
}
