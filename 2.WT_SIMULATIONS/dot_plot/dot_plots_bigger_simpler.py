import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set_style("whitegrid")

# Load in PDBs
DDR1_3ZOS = md.load('3ZOS.pdb')
DDR1_4BKJ = md.load('4BKJ.pdb')
DDR1_Dasatinib = md.load('DDR1_Dasatinib.pdb')
DDR1_VX680= md.load('DDR1_VX680.pdb')

DDR1_4CKR = md.load('4CKR.pdb')
DDR1_5FDP = md.load('5fdp.pdb')
DDR1_5FDX = md.load('5fdx.pdb')

DDR1_5BVK = md.load('5bvk.pdb')
DDR1_5BVN = md.load('5bvn.pdb')
DDR1_5BVO = md.load('5bvo.pdb')
DDR1_5BVW = md.load('5bvw.pdb')

ABL_2F4J = md.load('2f4j.pdb')
ABL_2GQG = md.load('2gqg.pdb')
ABL_1IEP = md.load('1iep.pdb')
ABL_3OXZ = md.load('3oxz.pdb')

ABL_1OPJ = md.load('1OPJ.pdb')

DFG = dict()
DFG['DDR1_B'] = ['resSeq 753 and resname GLY and name CA and chainid 0','resSeq 785 and resname PHE and name CZ and chainid 0']
DFG['DDR1_A'] = ['resSeq 716 and resname GLY and name CA and chainid 0','resSeq 748 and resname PHE and name CZ and chainid 0']
DFG['ABL'] = ['resSeq 350 and resname ALA and name CA and chainid 0','resSeq 382 and resname PHE and name CZ and chainid 0']
DFG['ABL_1OPJ'] = ['resSeq 369 and resname ALA and name CA and chainid 0','resSeq 401 and resname PHE and name CZ and chainid 0']

DFG_dih = dict()
DFG_dih['DDR1_B'] = ['resSeq 783 and resname ALA and name CB and chainid 0','resSeq 783 and resname ALA and name CA and chainid 0','resSeq 784 and resname ASP and name CA and chainid 0','resSeq 784 and resname ASP and name CG and chainid 0']
DFG_dih['DDR1_A'] = ['resSeq 746 and resname ALA and name CB and chainid 0','resSeq 746 and resname ALA and name CA and chainid 0','resSeq 747 and resname ASP and name CA and chainid 0','resSeq 747 and resname ASP and name CG and chainid 0']
DFG_dih['ABL'] = ['resSeq 380 and resname ALA and name CB and chainid 0','resSeq 380 and resname ALA and name CA and chainid 0','resSeq 381 and resname ASP and name CA and chainid 0','resSeq 381 and resname ASP and name CG and chainid 0']
DFG_dih['ABL_1OPJ'] = ['resSeq 399 and resname ALA and name CB and chainid 0','resSeq 399 and resname ALA and name CA and chainid 0','resSeq 400 and resname ASP and name CA and chainid 0','resSeq 400 and resname ASP and name CG and chainid 0']

ASP_phi = dict()
ASP_phi['DDR1_B'] = ['resSeq 784 and resname ASP and name C and chainid 0','resSeq 784 and resname ASP and name CA and chainid 0','resSeq 784 and resname ASP and name N and chainid 0','resSeq 783 and resname ALA and name C and chainid 0']
ASP_phi['DDR1_A'] = ['resSeq 747 and resname ASP and name C and chainid 0','resSeq 747 and resname ASP and name CA and chainid 0','resSeq 747 and resname ASP and name N and chainid 0','resSeq 746 and resname ALA and name C and chainid 0']
ASP_phi['ABL'] = ['resSeq 381 and resname ASP and name C and chainid 0','resSeq 381 and resname ASP and name CA and chainid 0','resSeq 381 and resname ASP and name N and chainid 0','resSeq 380 and resname ALA and name C and chainid 0']
ASP_phi['ABL_1OPJ'] = ['resSeq 400 and resname ASP and name C and chainid 0','resSeq 400 and resname ASP and name CA and chainid 0','resSeq 400 and resname ASP and name N and chainid 0','resSeq 399 and resname ALA and name C and chainid 0']

def DFG_distance(trajectories,def_DFG):

    distance = []

    for traj in trajectories:

        topology = traj.topology

        def_DFG_atom_1 = topology.select(def_DFG[0])
        def_DFG_atom_2 = topology.select(def_DFG[1])
   
        def_DFG_atoms = [def_DFG_atom_1[0], def_DFG_atom_2[0]]

        distance.append(md.compute_distances(traj,[def_DFG_atoms]))

    flattened_distance = np.asarray([val for sublist in distance for val in sublist])

    return [flattened_distance]

def DFG_dihedral(trajectories,def_DFG):

    dihedral = []

    for traj in trajectories:

        topology = traj.topology

        def_DFG_atom_1 = topology.select(def_DFG[0])
        def_DFG_atom_2 = topology.select(def_DFG[1])        
        def_DFG_atom_3 = topology.select(def_DFG[2])
        def_DFG_atom_4 = topology.select(def_DFG[3])  

        def_DFG_atoms = [def_DFG_atom_1[0], def_DFG_atom_2[0],def_DFG_atom_3[0], def_DFG_atom_4[0]]
   
        dihedral.append(md.compute_dihedrals(traj,[def_DFG_atoms]))

    flattened_dihedral = np.asarray([val for sublist in dihedral for val in sublist])

    return [flattened_dihedral]
	
[DDR1_3ZOS_dih] = DFG_dihedral(DDR1_3ZOS, DFG_dih['DDR1_B'])
[DDR1_4BKJ_dih] = DFG_dihedral(DDR1_4BKJ, DFG_dih['DDR1_B'])
[DDR1_Dasatinib_dih] = DFG_dihedral(DDR1_Dasatinib, DFG_dih['DDR1_A'])
[DDR1_VX680_dih] = DFG_dihedral(DDR1_VX680, DFG_dih['DDR1_A'])

[DDR1_4CKR_dih] = DFG_dihedral(DDR1_4CKR, DFG_dih['DDR1_B'])
[DDR1_5FDP_dih] = DFG_dihedral(DDR1_5FDP, DFG_dih['DDR1_B'])
[DDR1_5FDX_dih] = DFG_dihedral(DDR1_5FDX, DFG_dih['DDR1_B'])

[DDR1_5BVK_dih] = DFG_dihedral(DDR1_5BVK, DFG_dih['DDR1_B'])
[DDR1_5BVN_dih] = DFG_dihedral(DDR1_5BVN, DFG_dih['DDR1_B'])
[DDR1_5BVO_dih] = DFG_dihedral(DDR1_5BVO, DFG_dih['DDR1_B'])
[DDR1_5BVW_dih] = DFG_dihedral(DDR1_5BVW, DFG_dih['DDR1_B'])

[ABL_2F4J_dih] = DFG_dihedral(ABL_2F4J, DFG_dih['ABL'])
[ABL_2GQG_dih] = DFG_dihedral(ABL_2GQG, DFG_dih['ABL'])
[ABL_1IEP_dih] = DFG_dihedral(ABL_1IEP, DFG_dih['ABL'])
[ABL_3OXZ_dih] = DFG_dihedral(ABL_3OXZ, DFG_dih['ABL'])

[ABL_1OPJ_dih] = DFG_dihedral(ABL_1OPJ, DFG_dih['ABL_1OPJ'])

plt.figure(figsize=(5,2.5))

plt.scatter(DDR1_VX680_dih,0.5, edgecolors="slateblue", marker='o', linewidth='3', s=80, facecolors='none',label='DDR1_VX680')
plt.scatter(DDR1_Dasatinib_dih,0.25, edgecolors="slateblue", marker='o', linewidth='3', s=80, facecolors='none',label='DDR1_Dasatinib')
plt.scatter(DDR1_4BKJ_dih,0, edgecolors="slateblue", marker='o', linewidth='3', s=80, facecolors='none',label='4BKJ (DDR1_Imatinib)')
plt.scatter(DDR1_3ZOS_dih,-0.25, edgecolors="slateblue", marker='o', linewidth='3', s=80, facecolors='none',label='3ZOS (DDR1_Ponatinib)')

plt.scatter(ABL_2F4J_dih,-0.75, edgecolors="green", marker='*', linewidth='3', s=80, facecolors='none',label='2F4J (ABL_VX680)')
plt.scatter(ABL_2GQG_dih,-0.5, edgecolors="green", marker='*', linewidth='3', s=80, facecolors='none',label='2GQG (ABL_Dasatinib)')
plt.scatter(ABL_1IEP_dih,-1, edgecolors="green", marker='*', linewidth='3', s=80, facecolors='none',label='1IEP (ABL_Imatinib)')
plt.scatter(ABL_3OXZ_dih,-1.25, edgecolors="green", marker='*', linewidth='3', s=80, facecolors='none',label='3OXZ (ABL_Ponatinib)')

plt.xlabel('ALA_CB ALA_CA ASP_CA ASP_CG pseudo-dihedral (radians)')
plt.ylim((-1.4, 0.65))
plt.xlim((-4,4))
plt.yticks([])
#plt.legend()

ax = plt.gca()
ax.invert_xaxis()

plt.tight_layout()
plt.savefig('DDR1_ABL_DFG_pseudo_dihedral_1D-simpler.png',dpi=300)

plt.clf()

[DDR1_3ZOS_distance] = DFG_distance(DDR1_3ZOS, DFG['DDR1_B'])
[DDR1_4BKJ_distance] = DFG_distance(DDR1_4BKJ, DFG['DDR1_B'])
[DDR1_Dasatinib_distance] = DFG_distance(DDR1_Dasatinib, DFG['DDR1_A'])
[DDR1_VX680_distance] = DFG_distance(DDR1_VX680, DFG['DDR1_A'])

[DDR1_4CKR_distance] = DFG_distance(DDR1_4CKR, DFG['DDR1_B'])
[DDR1_5FDP_distance] = DFG_distance(DDR1_5FDP, DFG['DDR1_B'])
[DDR1_5FDX_distance] = DFG_distance(DDR1_5FDX, DFG['DDR1_B'])

[DDR1_5BVK_distance] = DFG_distance(DDR1_5BVK, DFG['DDR1_B'])
[DDR1_5BVN_distance] = DFG_distance(DDR1_5BVN, DFG['DDR1_B'])
[DDR1_5BVO_distance] = DFG_distance(DDR1_5BVO, DFG['DDR1_B'])
[DDR1_5BVW_distance] = DFG_distance(DDR1_5BVW, DFG['DDR1_B'])

[ABL_2F4J_distance] = DFG_distance(ABL_2F4J, DFG['ABL'])
[ABL_2GQG_distance] = DFG_distance(ABL_2GQG, DFG['ABL'])
[ABL_1IEP_distance] = DFG_distance(ABL_1IEP, DFG['ABL'])
[ABL_3OXZ_distance] = DFG_distance(ABL_3OXZ, DFG['ABL'])

[ABL_1OPJ_distance] = DFG_distance(ABL_1OPJ, DFG['ABL_1OPJ'])

plt.figure(figsize=(5,2.5))

plt.scatter(DDR1_VX680_distance,0.5, edgecolors="slateblue", marker='o', linewidth='3', s=80, facecolors='none',label='DDR1_VX680')
plt.scatter(DDR1_Dasatinib_distance,0.25, edgecolors="slateblue", marker='o', linewidth='3', s=80, facecolors='none',label='DDR1_Dasatinib')
plt.scatter(DDR1_4BKJ_distance,0, edgecolors="slateblue", marker='o', linewidth='3', s=80, facecolors='none',label='4BKJ (DDR1_Imatinib)')
plt.scatter(DDR1_3ZOS_distance,-0.25, edgecolors="slateblue", marker='o', linewidth='3', s=80, facecolors='none',label='3ZOS (DDR1_Ponatinib)')

plt.scatter(ABL_2F4J_distance,-0.75, edgecolors="green", marker='*', linewidth='3', s=80, facecolors='none',label='2F4J (ABL_VX680)')
plt.scatter(ABL_2GQG_distance,-0.5, edgecolors="green", marker='*', linewidth='3', s=80, facecolors='none',label='2GQG (ABL_Dasatinib)')
plt.scatter(ABL_1IEP_distance,-1, edgecolors="green", marker='*', linewidth='3', s=80, facecolors='none',label='1IEP (ABL_Imatinib)')
plt.scatter(ABL_3OXZ_distance,-1.25, edgecolors="green", marker='*', linewidth='3', s=80, facecolors='none',label='3OXZ (ABL_Ponatinib)')

plt.xlabel('F748 CZ - G716 CG2 distance (nm)')
plt.ylim((-1.4, 0.65))
plt.xlim((0, 2.0))
plt.yticks([])
#plt.legend()

plt.tight_layout()
plt.savefig('DDR1_ABL_DFG_distance_1D-simpler.png',dpi=300)

plt.clf()

[DDR1_3ZOS_phi] = DFG_dihedral(DDR1_3ZOS, ASP_phi['DDR1_B'])
[DDR1_4BKJ_phi] = DFG_dihedral(DDR1_4BKJ, ASP_phi['DDR1_B'])
[DDR1_Dasatinib_phi] = DFG_dihedral(DDR1_Dasatinib, ASP_phi['DDR1_A'])
[DDR1_VX680_phi] = DFG_dihedral(DDR1_VX680, ASP_phi['DDR1_A'])

[DDR1_4CKR_phi] = DFG_dihedral(DDR1_4CKR, ASP_phi['DDR1_B'])
[DDR1_5FDP_phi] = DFG_dihedral(DDR1_5FDP, ASP_phi['DDR1_B'])
[DDR1_5FDX_phi] = DFG_dihedral(DDR1_5FDX, ASP_phi['DDR1_B'])

[DDR1_5BVK_phi] = DFG_dihedral(DDR1_5BVK, ASP_phi['DDR1_B'])
[DDR1_5BVN_phi] = DFG_dihedral(DDR1_5BVN, ASP_phi['DDR1_B'])
[DDR1_5BVO_phi] = DFG_dihedral(DDR1_5BVO, ASP_phi['DDR1_B'])
[DDR1_5BVW_phi] = DFG_dihedral(DDR1_5BVW, ASP_phi['DDR1_B'])

[ABL_2F4J_phi] = DFG_dihedral(ABL_2F4J, ASP_phi['ABL'])
[ABL_2GQG_phi] = DFG_dihedral(ABL_2GQG, ASP_phi['ABL'])
[ABL_1IEP_phi] = DFG_dihedral(ABL_1IEP, ASP_phi['ABL'])
[ABL_3OXZ_phi] = DFG_dihedral(ABL_3OXZ, ASP_phi['ABL'])

[ABL_1OPJ_phi] = DFG_dihedral(ABL_1OPJ, ASP_phi['ABL_1OPJ'])

plt.figure(figsize=(5,2.5))

plt.scatter(DDR1_VX680_phi,0.5, edgecolors="slateblue", marker='o', linewidth='3', s=80, facecolors='none',label='DDR1_VX680')
plt.scatter(DDR1_Dasatinib_phi,0.25, edgecolors="slateblue", marker='o', linewidth='3', s=80, facecolors='none',label='DDR1:Dasatinib')
plt.scatter(DDR1_4BKJ_phi,0, edgecolors="slateblue", marker='o', linewidth='3', s=80, facecolors='none',label='4BKJ (DDR1_Imatinib)')
plt.scatter(DDR1_3ZOS_phi,-0.25, edgecolors="slateblue", marker='o', linewidth='3', s=80, facecolors='none',label='3ZOS (DDR1_Ponatinib)')

plt.scatter(ABL_2F4J_phi,-0.75, edgecolors="green", marker='o', linewidth='3', s=80, facecolors='none',label='2F4J (ABL_VX680)')
plt.scatter(ABL_2GQG_phi,-0.5, edgecolors="green", marker='o', linewidth='3', s=80, facecolors='none',label='2GQG (ABL_Dasatinib)')
plt.scatter(ABL_1IEP_phi,-1, edgecolors="green", marker='o', linewidth='3', s=80, facecolors='none',label='1IEP (ABL_Imatinib)')
plt.scatter(ABL_3OXZ_phi,-1.25, edgecolors="green", marker='o', linewidth='3', s=80, facecolors='none',label='3OXZ (ABL_Ponatinib)')

plt.xlabel('ASP phi (radians)',fontsize=14, fontweight='bold')
plt.xticks(fontsize=14,fontweight='bold')
plt.ylim((-1.4, 0.65))
plt.xlim((-4,4))
plt.yticks([])
#plt.legend()

ax = plt.gca()
ax.invert_xaxis()

plt.tight_layout()
plt.savefig('DDR1_ABL_DFG_ASP_phi_1D-simpler.png',dpi=300)

plt.figure(figsize=(5,2.5))

plt.scatter(DDR1_VX680_phi* (180.0 / np.pi),0.5, edgecolors="slateblue", marker='o', linewidth='3', s=80, facecolors='none',label='DDR1_VX680')
plt.scatter(DDR1_Dasatinib_phi* (180.0 / np.pi),0.25, edgecolors="slateblue", marker='o', linewidth='3', s=80, facecolors='none',label='DDR1:Dasatinib')
plt.scatter(DDR1_4BKJ_phi* (180.0 / np.pi),0, edgecolors="slateblue", marker='o', linewidth='3', s=80, facecolors='none',label='4BKJ (DDR1_Imatinib)')
plt.scatter(DDR1_3ZOS_phi* (180.0 / np.pi),-0.25, edgecolors="slateblue", marker='o', linewidth='3', s=80, facecolors='none',label='3ZOS (DDR1_Ponatinib)')

plt.scatter(ABL_2F4J_phi* (180.0 / np.pi),-0.75, edgecolors="green", marker='o', linewidth='3', s=80, facecolors='none',label='2F4J (ABL_VX680)')
plt.scatter(ABL_2GQG_phi* (180.0 / np.pi),-0.5, edgecolors="green", marker='o', linewidth='3', s=80, facecolors='none',label='2GQG (ABL_Dasatinib)')
plt.scatter(ABL_1IEP_phi* (180.0 / np.pi),-1, edgecolors="green", marker='o', linewidth='3', s=80, facecolors='none',label='1IEP (ABL_Imatinib)')
plt.scatter(ABL_3OXZ_phi* (180.0 / np.pi),-1.25, edgecolors="green", marker='o', linewidth='3', s=80, facecolors='none',label='3OXZ (ABL_Ponatinib)')

plt.xlabel('DFG ASP $\phi$',fontsize=14,fontweight='bold')
plt.xticks([-180,-90,0,90,180],fontsize=14,fontweight='bold')
plt.ylim((-1.4, 0.65))
plt.xlim((-180,180))
plt.yticks([])
#plt.legend()

#ax = plt.gca()
#ax.invert_xaxis()

plt.tight_layout()
plt.savefig('DDR1_ABL_DFG_ASP_phi_1D-simpler-degrees.png',dpi=300)

