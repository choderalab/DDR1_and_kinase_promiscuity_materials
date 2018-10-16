# Reseeding WT and mutant DDR1 simulations

Code used to run reseeded WT and mutant DDR1 simulations and generate parts of Figures 3 and S4 of the paper 'What makes a kinase promiscuous for inhibitors?'.  

#### Data
* Complete Folding@home simulation data analyzed for this publication is available via the [Open Science Framework](https://osf.io/).
#### Python Simulation Setup Scripts
* `models` - scripts used make DDR1 models and setup Folding@home simulations
 * This also includes the folder `multi_model_DDR1_Dasatinib` with scripts used to make the model of DDR1 bound to Dasatinib based on multiple templates.
* `picking_new_starting` - scripts and plots used to pick new starting configurations to reseed Folding@home simulations
#### Python Analysis Scripts
* `MSM_analysis` - scripts used to perform HMM and MSM analysis of Folding@home trajectories and produce Figure 3C and Figure S4
