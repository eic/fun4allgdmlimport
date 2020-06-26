# Fun4AllGdmlImport

Importing GDML files for the silicon created in EICROOT into the Fun4All framework.
Using the Fun4All structure of things.

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