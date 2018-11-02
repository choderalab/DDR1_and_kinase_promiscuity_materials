import matplotlib
matplotlib.use('Agg')

#import general libraries
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt

from glob import glob

D671N_path_to_trajs = 'D671N-pro-trajectories/*.h5'
Y755A_path_to_trajs = 'Y755A-pro-trajectories/*.h5'
Y759A_path_to_trajs = 'Y759A-pro-trajectories/*.h5'

filenames_D671N = sorted(glob(D671N_path_to_trajs))
ref_traj_D671N = md.load(filenames_D671N[0], stride=10)
top_D671N = ref_traj_D671N.topology

filenames_Y755A = sorted(glob(Y755A_path_to_trajs))
ref_traj_Y755A = md.load(filenames_Y755A[0], stride=10)
top_Y755A = ref_traj_Y755A.topology

filenames_Y759A = sorted(glob(Y759A_path_to_trajs))
ref_traj_Y759A = md.load(filenames_Y759A[0], stride=10)
top_Y759A = ref_traj_Y759A.topology

def convert_atom_list_to_resid(atom_list,topology):
    resid_set = set()
    for atom_index in atom_list:
        atom = topology.atom(atom_index)
        resid_set.add(atom.residue.index)
    resid = list(resid_set)
    if len(resid) == 1:
        resid = resid[0]
    else:
        print('your selection does not select a single residue')
    return resid

def add_kinase_coords_featurizer(feat,top):
   
    # Roux pseudo-dihedral
    Roux = ['resid 179 and resname ALA and name CB',
        'resid 179 and resname ALA and name CA',
        'resid 180 and resname ASP and name CA',
        'resid 180 and resname ASP and name CG']

    Roux_atoms = np.array([],dtype=np.int64)
    for atom in Roux:
        atom_select = int(top.select(atom))
        Roux_atoms = np.append(Roux_atoms, atom_select )

    feat.add_dihedrals([Roux_atoms])

    # DFG distance
    DFG_distance = ['resid 149 and resname GLY and name CA',
                'resid 181 and resname PHE and name CZ']

    def_DFG_atom_1 = top.select(DFG_distance[0])
    def_DFG_atom_2 = top.select(DFG_distance[1])

    feat.add_distances([def_DFG_atom_1[0],def_DFG_atom_2[0]])
    
    # Two KER distances
    KER =  ['resid 51 and resname LYS',
        'resid 68 and resname GLU',
        'resid 185 and resname ARG']

    KER_res = np.array([],dtype=np.int64)
    for res in KER:
        atom_select = top.select(res)
        res_select = convert_atom_list_to_resid(atom_select,top)
        KER_res = np.append(KER_res, res_select )

    feat.add_residue_mindist([[KER_res[0],KER_res[1]],[KER_res[1],KER_res[2]]])
    
    #Let's also add our D671 - R752 distance
    DR = ['resid 104 and (resname ASP or resname ASN)',
      'resid 185 and resname ARG']

    DR_res = np.array([],dtype=np.int64)
    for res in DR:
        atom_select = top.select(res)
        res_select = convert_atom_list_to_resid(atom_select,top)
        DR_res = np.append(DR_res, res_select )

    feat.add_residue_mindist([[ DR_res[0], DR_res[1] ]])

    # DFG phi, psi, chi angles
      #Note that because of the way dihedrals are defined, we're also including
      # one residue before and after
    DFG_residues = ['resid 179 and resname ALA',
                'resid 180 and resname ASP',
                'resid 181 and resname PHE',
                'resid 182 and resname GLY',
                'resid 183 and resname MET']
    for residue in DFG_residues:
        feat.add_backbone_torsions(residue)

    for residue in DFG_residues:
        try:
            feat.add_chi1_torsions(residue)
        except:
            pass


    # other pseudo dihedrals
    # pseudo dihedrals from Roux Phys Chem B 2015 2D Umbrella Sampling paper
    US2D_Abl = ['resid 179 and resname ALA and name CB',
        'resid 179 and resname ALA and name CA',
        'resid 181 and resname PHE and name CA',
        'resid 181 and resname PHE and name CG']

    US2D_Abl_atoms = np.array([],dtype=np.int64)
    for atom in US2D_Abl:
        atom_select = int(top.select(atom))
        US2D_Abl_atoms = np.append(US2D_Abl_atoms, atom_select )

    feat.add_dihedrals([US2D_Abl_atoms])

    US2D_Src = ['resid 183 and resname MET and name CB',
        'resid 183 and resname MET and name CA',
        'resid 181 and resname PHE and name CA',
        'resid 181 and resname PHE and name CG']

    US2D_Src_atoms = np.array([],dtype=np.int64)
    for atom in US2D_Src:
        atom_select = int(top.select(atom))
        US2D_Src_atoms = np.append(US2D_Src_atoms, atom_select )

    feat.add_dihedrals([US2D_Src_atoms])
	
  # pseudo dihedrals from Mobitz are defined by four consecutive CA's
    Mobitz1 = ['resid 179 and resname ALA and name CA',
           'resid 180 and resname ASP and name CA',
           'resid 181 and resname PHE and name CA',
           'resid 182 and resname GLY and name CA']

    Mobitz1_atoms = np.array([],dtype=np.int64)
    for atom in Mobitz1:
        atom_select = int(top.select(atom))
        Mobitz1_atoms = np.append(Mobitz1_atoms, atom_select )

    feat.add_dihedrals([Mobitz1_atoms])

    Mobitz2 = ['resid 180 and resname ASP and name CA',
           'resid 181 and resname PHE and name CA',
           'resid 182 and resname GLY and name CA',
           'resid 183 and resname MET and name CA']

    Mobitz2_atoms = np.array([],dtype=np.int64)
    for atom in Mobitz2:
        atom_select = int(top.select(atom))
        Mobitz2_atoms = np.append(Mobitz2_atoms, atom_select )

    feat.add_dihedrals([Mobitz2_atoms])

    # minimal distances to include R-spine
    Rspine =  ['resid 72 and resname MET',
           'resid 83 and resname LEU',
           'resid 160 and resname HIS',
           'resid 181 and resname PHE']

    Rspine_res = np.array([],dtype=np.int64)
    for res in Rspine:
        atom_select = top.select(res)
        res_select = convert_atom_list_to_resid(atom_select,top)
        Rspine_res = np.append(Rspine_res, res_select )

    feat.add_residue_mindist([[Rspine_res[0],Rspine_res[1]],[Rspine_res[1],Rspine_res[2]],[Rspine_res[2],Rspine_res[3]]])

    print('Final Features Dimensions: %s '%feat.dimension())

    return feat

## D671N ##
# Make our featurizers
feat_D671N = coor.featurizer(top_D671N)
feat_D671N = add_kinase_coords_featurizer(feat_D671N,top_D671N)

# Write out files for these features for our D671N trajs
src_D671N  = coor.source(filenames_D671N, features=feat_D671N)
calculated_features_D671N = src_D671N.get_output()

print('len(calculated_features_D671N): %s' %len(calculated_features_D671N))

for i, traj in enumerate(calculated_features_D671N):
    np.save('D671N-pro/calculated_features_D671N_%s.npy'%i, traj)
    
## Y755A ##
# Make our featurizers
feat_Y755A = coor.featurizer(top_Y755A)
feat_Y755A = add_kinase_coords_featurizer(feat_Y755A,top_Y755A)

# Write out files for these features for our Y755A trajs
src_Y755A  = coor.source(filenames_Y755A, features=feat_Y755A
calculated_features_Y755A = src_Y755A.get_output()

print('len(calculated_features_Y755A): %s' %len(calculated_features_Y755A)

for i, traj in enumerate(calculated_features_Y755A):
    np.save('Y755A-pro/calculated_features_Y755A_%s.npy'%i, traj)
    
## Y759A ##
# Make our featurizers
feat_Y759A = coor.featurizer(top_Y759A)
feat_Y759A = add_kinase_coords_featurizer(feat_Y759A,top_Y759A)

# Write out files for these features for our Y759A trajs
src_Y759A = coor.source(filenames_Y759A, features=feat_Y759A
calculated_features_Y759A = src_Y759A.get_output()

print('len(calculated_features_Y759A): %s' %len(calculated_features_Y759A)

for i, traj in enumerate(calculated_features_Y759A):
    np.save('Y759A-pro/calculated_features_Y759A_%s.npy'%i, traj)