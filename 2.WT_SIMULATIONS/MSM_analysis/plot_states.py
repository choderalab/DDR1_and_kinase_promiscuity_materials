import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from msmbuilder import dataset
import numpy as np
import mdtraj as md
from glob import glob

files = glob('hmm*.xtc')
 
top = 'hmm1_100samples.pdb'

protein = 'DDR1'

DFG = dict()
DFG['DDR1'] = ['resid 149 and resname GLY and name CA','resid 181 and resname PHE and name CZ']

KER_hbond = {}
KER_hbond['DDR1'] =  ['resid 51 and resname LYS','resid 68 and resname GLU','resid 185 and resname ARG']

def convert_atom_list_to_resid(atom_list,topology):
    resid_set = set()
    for atom_index in atom_list:
        atom = topology.atom(atom_index)
        resid_set.add(atom.residue.index)
    resid = list(resid_set)
    if len(resid) == 1:
        resid = resid[0]
    else:
        print 'your selection does not select a single residue'
    return resid

def DFG_KER_byrun(files,KER,def_DFG):

    difference = []
    DFG = []

    difference_combinetrajs = []
    DFG_combinetrajs = []

    for file in files:
        
        print 'working on %s' %file

        trajectories = dataset.MDTrajDataset(file,topology=top)

        for traj in trajectories:

            topology = traj.topology

           # append difference
            KER_K_atoms = topology.select(KER[0])
            KER_E_atoms = topology.select(KER[1])
            KER_R_atoms = topology.select(KER[2])

            KER_K = convert_atom_list_to_resid(KER_K_atoms,topology)
            KER_E = convert_atom_list_to_resid(KER_E_atoms,topology)
            KER_R = convert_atom_list_to_resid(KER_R_atoms,topology)

            #print 'Atom distances computed between %s, %s, and %s' %(topology.residue(KER_K),topology.residue(KER_E),topology.residue(KER_R))

            # note the default for compute_contacts is 'closest-heavy'
            k295e310 = md.compute_contacts(traj, [[KER_K,KER_E]])
            e310r409 = md.compute_contacts(traj, [[KER_E,KER_R]])

            difference_combinetrajs.append(10*(e310r409[0] - k295e310[0])) # 10x because mdtraj is naturally in nm

            # append DFG
            def_DFG_atom_1 = topology.select(def_DFG[0])
            def_DFG_atom_2 = topology.select(def_DFG[1])

            #print 'Atom distances computed between %s and %s' %(topology.atom(def_DFG_atom_1),topology.atom(def_DFG_atom_2))       
            def_DFG_atoms = [def_DFG_atom_1[0], def_DFG_atom_2[0]]
            #print 'These correspond to atom numbers %s.' %def_DFG_atoms

            DFG_combinetrajs.append(md.compute_distances(traj,[def_DFG_atoms]))

        # flatten list of arrays
        difference_combinetrajs = np.asarray([val for sublist in difference_combinetrajs for val in sublist])
        DFG_combinetrajs = np.asarray([val for sublist in DFG_combinetrajs for val in sublist])

        difference.append(difference_combinetrajs)
        difference_combinetrajs = []

        DFG.append(DFG_combinetrajs)
        DFG_combinetrajs = []

    return [DFG, difference]

[DFG_separate,difference_separate] = DFG_KER_byrun(files,KER_hbond[protein],DFG[protein])

#save DFG and difference data
np.save('DFG_hmm_states.npy',DFG_separate)
np.save('difference_hmm_states.npy',difference_separate)

import matplotlib
colors = ['red','cyan','blue','purple','magenta','brown','black','white','orange','pink','yellowgreen']

for i in range(len(DFG_separate)):
    plt.plot(DFG_separate[i],difference_separate[i],'o',color=colors[i])

plt.ylabel('d(E635-R752)-d(K618-E635) $\AA$')
plt.xlabel('G716 CA - F748 CZ $\AA$')
plt.xlim(0.5,2.5)
plt.ylim(-20,20)

plt.savefig('hmm-states-on-DFG.png',dpi=700)
