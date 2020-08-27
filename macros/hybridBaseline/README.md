# Hybrid detector silicon baseline

This folder contains the GDML files for the silicon parts of the hybrid detector baseline.

It has 5 silicon barrel layers, consisting of 3 vertexing layers and 2 outer layers. The vertexing layers have a material budget of 0.05% X/X0 per layer, and the outer layers have 0.55% X/X0 per layer.

The disks are in a tapered configuration for the innermost two, hitting full available radius at the third disk. The material thickness for the disks is 0.24% X/X0.

The detector parts present are changed in detector_setup.h.

To get a display of this configuration:
First, make sure you have done the steps in the basic readme (https://gitlab.com/hwennlof/fun4allgdmlimport/-/tree/master/), to get Singularity going with the GdmlImporter installed and sourced. Then do

```
root -l Fun4All_G4_HybridBaseline.C\(-1\)
.L DisplayOn()
PHG4Reco * g4 = DisplayOn()
```

Do not worry if the three vertexing layers look odd in the display; this is an artefact of their thinness and how it is handled via EICROOT. They function correctly, and have the correct material.