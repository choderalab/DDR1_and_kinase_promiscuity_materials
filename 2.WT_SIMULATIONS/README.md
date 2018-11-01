# WT DDR1 simulations

Code used to run WT DDR1 simulations and generate parts of Figures 2, 3, and S4 of the paper 'What makes a kinase promiscuous for inhibitors?'.  

#### Data
* Complete Folding@home simulation data analyzed for this publication is available via the [Open Science Framework](https://osf.io/4r8x2/).
#### Python Simulation Setup Scripts
* `models` - scripts used make DDR1 models and setup Folding@home simulations
  * Here we have included the full setup scripts, but in the end we only present the analysis of the protonated DDR1 models on already existing DDR1 structures
#### Python Analysis Scripts
* `dot_plot` - scripts and PDBs used to make Figure 2D
* `MSM_analysis` - scripts used to perform HMM and MSM analysis of Folding@home trajectories and produce Figure 3A and Figure S4
