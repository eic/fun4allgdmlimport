# Fun4AllGdmlImport

Importing GDML files for the silicon created in EICROOT into the Fun4All framework.
Using the Fun4All structure of things.

Recompile by setting up the Fun4All setup;
`singularity shell -B /home  -B /cvmfs /cvmfs/eic.opensciencegrid.org/singularity/rhic_sl7_ext.simg`

`source /cvmfs/eic.opensciencegrid.org/x8664_sl7/opt/fun4all/core/bin/eic_setup.sh -n`

`source /cvmfs/eic.opensciencegrid.org/x8664_sl7/opt/fun4all/core/bin/setup_local.sh $HOME/myinstall`

then going into source/build and running

`../autogen.sh --prefix=$HOME/myinstall`

and

`make install`

Then `source /cvmfs/eic.opensciencegrid.org/x8664_sl7/opt/fun4all/core/bin/setup_local.sh $HOME/myinstall`.

Then we can go into macros and run `root -l Fun4All_G4_FastMom.C\(-1\)`.

The most recent version of this file is located in **macros/hepMCeventReading/newBaseline_775_actualBeampipe/ITS3**

The active detector parts can be changed via the file detector_setup.h.

# To add a detector:

In https://gitlab.com/hwennlof/fun4allgdmlimport/-/blob/master/source/GdmlImportDetector.cc: Add a relevant string to detprefix. This is for the detector name creation, and assignment of detector indices.
There is a potential difficulty here, in extracting the active part of the detector. For the silicon it is "MimosaCore", but the detector prefix varies, and comes from EICROOT.
As it is now, it pulls the VstStave detectors out, selects "MimosaCore" as the active part, puts the rest of the geometry as passive parts, and gives the layers numbers starting from 10.
The same is then done for the time-stamping layer(s), with indices starting from 20.

In the macro you want to run:
`   

    // Loading All-Si Tracker from gdml file
	GdmlImportDetectorSubsystem * svtPart = new GdmlImportDetectorSubsystem();
	svtPart->set_string_param("GDMPath","Vst_GDML_.gdml"); 
	svtPart->AddAssemblyVolume("VST");	// Barrel. ********** THIS IS THE NAME IN THE GDML FILE, FROM EICROOT ********
	svtPart->SuperDetector("SVT");
	svtPart->SetActive();          // this saves hits in the MimosaCore volumes
	svtPart->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)
	g4Reco->registerSubsystem(svtPart);`
	
Make sure the AssemblyVolume name matches the one in the GDML file.


Also in the macro to run, if reconstruction is desired, each detector part needs to be added to the (for example) Kalman filter;
`

    //SVT
	//This adds the hit nodes for EACH LAYER OF THE DETECTOR.
	for (int i = 10; i < 10+_NO_OF_BARREL_LAYERS_; i++) {
		std::string nodeName = "G4HIT_SVT_" + std::to_string(i); // ****** THIS NEEDS TO MATCH THE NAME FROM GetName() in GdmlImportDetector.cc ******
		kalman->add_phg4hits(
				nodeName,				// const std::string& phg4hitsNames
				PHG4TrackFastSim::Cylinder,		// const DETECTOR_TYPE phg4dettype
				999.,					// radial-resolution [cm] (this number is not used in cylindrical geometry)
				5.8e-4,					// azimuthal (arc-length) resolution [cm]
				5.8e-4,					// longitudinal (z) resolution [cm]
				1,					// efficiency (fraction)
				0					// hit noise
		);
	}`

This is where the detector number comes in. It is linked to the G4Hits, so it's important to get right.