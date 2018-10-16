# This script selects structures from simulations that have sampled the DFG-flip to restart simulations.

# Imports

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import seaborn
import mdtraj as md
import re
from glob import glob
import os

# First we will define our trajectories from which we are sampling

base_path = '/cbio/jclab/projects/fah/fah-data/munged3/no-solvent/11403'

filenames = [
    '%s/run1-clone23.h5'%base_path,
    '%s/run1-clone47.h5'%base_path,
    '%s/run3-clone5.h5'%base_path,
    '%s/run3-clone9.h5'%base_path]

# Then we will define the DFG flip coordinate, which we will use to adequately sample this transition.

  # In this case we are trying the distance between the
  # CZ of the F in the DFG motif and the CA of a stable glycine/alanine.
  # QIAS_G_MAY in Src and QISS_A_MEY in Abl
  # and QIAS_G_MRY in DDR1

DDR1_DFG_distance = ['resid 149 and resname GLY and name CA','resid 181 and resname PHE and name CZ']

# We also want to define the KER distance for our 2D plot

KER_hbond_coords = [[51,68],[68,185]]

# This will be our csv file that stores how we've selected these structures
writer = open('selected-structures.csv','w')
writer.write('%s,%s,%s,%s\n'%('trajectory','DFG distance','frame','new pdb file'))

# We will make a 1D plot of this for our initial trajectory.

for filename in filenames:

    short_filename = filename.split('/')[-1]
    clone_name = short_filename.split('.')[0]

    print 'working on %s' %short_filename

    traj = md.load(filename)
    topology = traj.topology

    def_DFG_atom_1 = topology.select(DDR1_DFG_distance[0])
    def_DFG_atom_2 = topology.select(DDR1_DFG_distance[1])   

    #print 'Atom distances computed between %s and %s' %(topology.atom(def_DFG_atom_1),topology.atom(def_DFG_atom_2))       
    def_DFG_atoms = [def_DFG_atom_1[0], def_DFG_atom_2[0]]
    #print 'These correspond to atom numbers %s.' %def_DFG_atoms

    distances = md.compute_distances(traj,[def_DFG_atoms])

    flattened_distances = [val for sublist in distances for val in sublist]

    plt.figure(figsize=(5,3))

    plt.plot(distances,'o',ms='3',label='%s'%short_filename)
    plt.axhline(y=1.2,color='red',label='cut-off')

    plt.ylim((0.6,2.2))
    plt.title('%s: Structure Selection'%clone_name)
    plt.ylabel('Distance(nm)',fontsize=16)
    plt.yticks(fontsize=13)
    plt.xlabel('Frame',fontsize=16)
    plt.xticks(fontsize=13)
    plt.legend(loc=0)
    plt.savefig('%s-distances-1D.png'%clone_name,bbox_inches='tight')
    plt.close()

# We will then sample along this coordinate to choose 10 structures (plotting where these structures fall).
 
    sorted_flattened_distances = np.sort(flattened_distances)

    #ten_distances = sorted_flattened_distances[100::len(sorted_flattened_distances)/10]

    #can we pick ten evenly spaced distances from these sorted_distances

    my_min = sorted_flattened_distances.min()
    my_max = sorted_flattened_distances.max()

    ideal_distances = np.arange(my_min+.1,my_max-.1,(my_max-my_min)/10.)
    print ideal_distances

    ten_distances = []
    for distance in ideal_distances:
        closest_distance = min(sorted_flattened_distances, key = lambda x:abs(x-distance))
        ten_distances.append(closest_distance)
    print ten_distances

    ten_frames = []
    for distance in ten_distances:
        ten_frames.append(flattened_distances.index(distance))
    print ten_frames

    #Write out these structures.

    for i,frame in enumerate(ten_frames):
        traj[frame].save_pdb('%s_%s.pdb'%(clone_name,frame))

        #and save in csv file
        writer.write('%s,%s,%s,%s_%s.pdb\n'%(short_filename,ten_distances[i],frame,clone_name,frame))

    plt.figure(figsize=(5,3))

    plt.plot(distances, 'o', ms='3', label='%s'%short_filename)
    plt.plot(ten_frames,ten_distances, 'o', color='m', label = 'selected structures')
    plt.axhline(y=1.2,color='red',label='cut-off')

    plt.ylim((0.6,2.2))
    plt.title('%s: Structure Selection'%clone_name)
    plt.ylabel('Distance(nm)',fontsize=16)
    plt.yticks(fontsize=13)
    plt.xlabel('Frame',fontsize=16)
    plt.xticks(fontsize=13)
    plt.legend(loc=0)
    plt.savefig('%s-selected-1D.png'%clone_name,bbox_inches='tight')
    plt.close()


# Let's also plot where these structures fall on a 2D plot using a coordinate representative of the active/inactive transition.

    k295e310 = md.compute_contacts(traj, [KER_hbond_coords[0]])
    e310r409 = md.compute_contacts(traj, [KER_hbond_coords[1]])
    KER = 10*(e310r409[0] - k295e310[0]) # 10x because mdtraj is naturally in nm

    plt.figure(figsize=(5,4))

    plt.plot(distances,KER,'.', alpha=0.25, ms='3', label='all frames')
    plt.plot(distances[ten_frames[0]], KER[ten_frames[0]],'o',color='m', ms='3',label = 'selected structures')
    for frame in ten_frames:
        plt.plot(distances[frame],KER[frame],'o',color='m',ms='3')
    plt.axvline(x=1.2,color='0.7',linestyle='--')
    plt.text(0.6,17.5,'DFG-in',fontweight='bold')
    plt.text(2.1,17.5,'DFG-out',fontweight='bold')

    plt.xlabel('DFG Distance (nm)')
    plt.ylabel('d(E635-R752) - d(618-E635) ($\AA$)')
    plt.xlim(0.5,2.5)
    plt.ylim(-20,20)
    plt.title('%s: Structure Selection' %clone_name)
    plt.legend(loc=4)

    plt.savefig('%s-selected-2D.png'%clone_name,bbox_inches='tight',dpi=700)
    plt.close()





