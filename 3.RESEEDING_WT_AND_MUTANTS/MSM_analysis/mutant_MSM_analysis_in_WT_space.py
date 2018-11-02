import matplotlib
matplotlib.use('Agg')

# import general libraries
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt

import corner

print('loading trajectories')

# load trajectories
path_to_trajs = '../DDR1_WT_pro_plus_OG_trajectories/*.h5'
from glob import glob
filenames = sorted(glob(path_to_trajs))

from tables.exceptions import NoSuchNodeError

trajs = list()
filenames_checked=list()
for filename in filenames:
    try:
        traj = md.load(filename, stride=10)
        good_file = filename
        trajs.append(traj)
        filenames_checked.append(good_file)
    except NoSuchNodeError as e:
        print('Trajectory file %s is corrupted: %s' % (filename, str(e)))

np.save('filenames_checked.npy',filenames_checked)

top = trajs[0].topology

print('featurizing')

# Make our featurizer
feat = coor.featurizer(top)

# add our known kinase coordinates
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

print('Final number of features')
print(feat.dimension())

print('running tica')

src  = coor.source(filenames_checked, features=feat)
tica = coor.tica(src, stride=10,lag=100,kinetic_map=False,commute_map=True)
Y = tica.get_output()

Y1 = [y[:,0] for y in Y]
Y2 = [y[:,1] for y in Y]

np.save('tica_projection-otherpro_repeat.npy',Y)

plt.clf()
plt.figure(figsize=(8,5))
mplt.plot_free_energy(np.hstack(Y1),np.hstack(Y2))
plt.xlabel('tic 1')
plt.ylabel('tic 2')

plt.savefig('tic1-tic2-otherpro_repeat.png')

colors = ['red','cyan','blue','purple','magenta','brown','black','white','orange','pink','yellowgreen']

print('Finished with WT-pro tica analysis')

print('Begin with other pro analysis projected on WT-pro tics')

#define our other three features

other_pro = {}
other_pro['D671N-pro'] = glob('D671N-pro/calculated_features_D671N_*.npy')
other_pro['Y755A-pro'] = glob('Y755A-pro/calculated_features_Y755A_*.npy')
other_pro['Y759A-pro'] = glob('Y759A-pro/calculated_features_Y759A_*.npy')

for mutant in other_pro:

    calculated_features = other_pro[mutant]

    otherpro_features = []
    for traj in calculated_features:
        otherpro_features.append(np.load(traj))
    
    Y_otherpro = tica.transform(otherpro_features)

    np.save('%s/tica_projection.npy'%mutant,Y_otherpro)

    Y1_otherpro = [y[:,0] for y in Y_otherpro]
    Y2_otherpro = [y[:,1] for y in Y_otherpro]
    Y3_otherpro = [y[:,2] for y in Y_otherpro]
    Y4_otherpro = [y[:,3] for y in Y_otherpro]
    Y5_otherpro = [y[:,4] for y in Y_otherpro]
    Y6_otherpro = [y[:,5] for y in Y_otherpro]

    sixtics=np.vstack((np.hstack(Y1_otherpro),np.hstack(Y2_otherpro),np.hstack(Y3_otherpro),np.hstack(Y4_otherpro),np.hstack(Y5_otherpro),np.hstack(Y6_otherpro))).T

    corner.corner(sixtics, labels=[r"$tic1$", r"$tic2$", r"$tic3$",r"$tic4$", r"$tic5$", r"$tic6$", r"$\Gamma \, [\mathrm{parsec}]$"], quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12});

    plt.savefig('%s/corner.png'%mutant)

    plt.clf()
    plt.figure(figsize=(8,5))
    mplt.plot_free_energy(np.hstack(Y1_otherpro),np.hstack(Y2_otherpro))
    plt.xlabel('tic 1')
    plt.ylabel('tic 2')

    plt.savefig('%s/tic1-tic2.png'%mutant)

    print('running %s kmeans'%mutant)

    clkmeans = coor.cluster_kmeans(Y_otherpro,300,max_iter=300)

    plt.clf()
    plt.figure(figsize=(8,5))
    plt.plot(clkmeans.clustercenters[:,0], clkmeans.clustercenters[:,1],' ok')
    mplt.plot_free_energy(np.hstack(Y1),np.hstack(Y2))
    plt.xlabel('tic 1')
    plt.ylabel('tic 2')

    plt.savefig('%s/kmeans_cluster-on_tic1tic2.png'%mutant)

    np.save('%s/clkmeans_dtrajs.npy'%mutant,clkmeans.dtrajs)
    np.save('%s/clkmeans_clustercenters.npy'%mutant,clkmeans.clustercenters)

    print('running %s MSM'%mutant)

    MSM = msm.estimate_markov_model(clkmeans.dtrajs, 50)

    time_scale_sep = MSM.timescales()[:-1]/MSM.timescales()[1:]
    ##slow_indices = [index for index, value in enumerate(time_scale_sep) if value >= 1.5]
    ##last_slow_index = slow_indices[-1]
    ##n_macrostates_timescales = last_slow_index + 2
    #We are actually going to over ride all this machinery
    n_macrostates = 2

    plt.clf()
    plt.figure(figsize=(5,3))
    plt.plot(time_scale_sep, linewidth=0,marker='o')
    #plt.axvline(x=last_slow_index+0.5,color='r')
    plt.axvline(x=0.5,color='g')
    plt.xlabel('index'); plt.ylabel('timescale separation');
    plt.xlim(0,30)

    plt.savefig('%s/timescale_separation.png'%mutant)

    #print('%s macrostates chosen from timescale separation' %n_macrostates_timescales)
    print('%s macrostates chosen because thats what we want' %n_macrostates)

    plt.clf()
    lags = [1,2,5,10,20,50,100,200,400]
    its = msm.its(clkmeans.dtrajs, lags=lags)
    mplt.plot_implied_timescales(its)

    plt.savefig('%s/implied_timescale_plot.png'%mutant)
    plt.clf()

    print('fraction of states used = ', MSM.active_state_fraction)
    print('fraction of counts used = ', MSM.active_count_fraction)

    mplt.plot_cktest(MSM.cktest(3));

    plt.savefig('%s/cktest_msm.png'%mutant)

    plt.clf()
    plt.figure(figsize=(8,5))
    mplt.plot_free_energy(np.hstack(Y1_otherpro),np.hstack(Y2_otherpro),weights=np.hstack(MSM.trajectory_weights()))
    plt.xlabel('tic 1')
    plt.ylabel('tic 2')

    plt.savefig('%s/reweighted_MSM_tic1_tic2.png'%mutant)

    print('running %s hmm'%mutant)
    HMM = msm.bayesian_hidden_markov_model(clkmeans.dtrajs,n_macrostates,lag=100,nsamples=1000)
    hmm_dist = HMM.metastable_distributions
    hmm_sets = HMM.metastable_sets

    plt.clf()
    plt.figure(figsize=(8,5))
    mplt.plot_free_energy(np.hstack(Y1_otherpro),np.hstack(Y2_otherpro),weights=np.hstack(HMM.trajectory_weights()))
    plt.xlabel('tic 1')
    plt.ylabel('tic 2')

    plt.savefig('%s/reweighted_HMM_tic1_tic2.png'%mutant)

    plt.clf()
    plt.figure(figsize=(8,5))
    mplt.plot_free_energy(np.hstack(Y1_otherpro),np.hstack(Y2_otherpro),weights=np.hstack(HMM.trajectory_weights()))
    colors = ['red','cyan','blue','purple','magenta','brown','black','white','orange','pink','yellowgreen']
    for i in range(n_macrostates):
        plt.scatter(clkmeans.clustercenters[hmm_sets[i],0], clkmeans.clustercenters[hmm_sets[i],1], color=colors[i],s=50)
    plt.xlabel('tic 1')
    plt.ylabel('tic 2')

    plt.savefig('%s/hmm_clusters_on_tic1_tic2.png'%mutant)

    state_1_colors = HMM.metastable_distributions[0]/max(HMM.metastable_distributions[0])
    state_2_colors = HMM.metastable_distributions[1]/max(HMM.metastable_distributions[1])

    plt.clf()
    plt.figure(figsize=(8,5))
    #NOTE this needs to be changed if more than two macrostates are of interest
    mplt.plot_free_energy(np.hstack(Y1_otherpro),np.hstack(Y2_otherpro),weights=np.hstack(HMM.trajectory_weights()))
    for i in range(len(clkmeans.clustercenters[:,0])):
        plt.scatter(clkmeans.clustercenters[i,0], clkmeans.clustercenters[i,1], color='red', s=100,alpha=(state_1_colors[i]))
        plt.scatter(clkmeans.clustercenters[i,0], clkmeans.clustercenters[i,1], color='cyan', s=100,alpha=(state_2_colors[i]))
    plt.xlabel('tic 1')
    plt.ylabel('tic 2')

    plt.savefig('%s/hmm_clusters_on_tic1_tic2_fade.png'%mutant)

    hmm_samples = HMM.sample_by_observation_probabilities(100)

    #for i in range(n_macrostates):
    #    coor.save_traj(src, hmm_samples[i], './hmm%s_100samples.xtc'%i)
    #    t = md.load( './hmm%s_100samples.xtc'%i,top=top)
    #    t[0].save_pdb('./hmm%s_100samples.pdb'%i)

    statdist = MSM.stationary_distribution
    relative_counts = MSM.count_matrix_active.sum(0) / np.sum(MSM.count_matrix_active)

    plt.clf()
    plt.figure()
    plt.scatter(statdist, relative_counts)
    plt.xlabel('MSM stationary distribution')
    plt.ylabel('Relative counts')

    plt.savefig('%s/msm_sanity_check.png'%mutant)

    statdist = HMM.stationary_distribution
    relative_counts = HMM.count_matrix_EM.sum(0) / np.sum(HMM.count_matrix_EM)

    plt.clf()
    plt.figure()
    plt.scatter(statdist, relative_counts)
    plt.xlabel('HMMstationary distribution')
    plt.ylabel('Relative counts')

    plt.savefig('%s/hmm_sanity_check.png'%mutant)
    plt.clf()

    mplt.plot_cktest(HMM.cktest(3));

    plt.savefig('%s/cktest_hmm.png'%mutant)
    plt.clf()

    # macrostate free energies
    f_i = -np.log(sorted(HMM.stationary_distribution))[::-1]
    f_i -= f_i.min()
    plt.figure()
    plt.plot(f_i, '.')
    plt.xlabel('Macrostate')
    plt.ylabel(r'$\Delta G$ $(k_B T)$')
    plt.title('Macrostate free energies')

    plt.savefig('%s/macrostate_free_energies.png'%mutant)
    plt.clf()

    delG_distribution = -np.log(np.transpose(HMM.sample_f('stationary_distribution'))[1]/
                            np.transpose(HMM.sample_f('stationary_distribution'))[0])

    delG_interval = np.percentile(a=delG_distribution, q=[2.5, 50.0, 97.5])

    plt.clf()
    plt.figure()
    plt.bar(0,delG_interval[1],yerr=[[delG_interval[1]-delG_interval[0]],[delG_interval[2]-delG_interval[1]]],color=colors,error_kw=dict(ecolor='gray',lw=2,capsize=5,capthick=2))
    plt.ylabel(r'$\Delta G$ $(k_B T)$')

    plt.savefig('%s/delG_states.png'%mutant)
    plt.clf()

    np.save('%s/delG_distribution.npy'%mutant,delG_distribution)

    print('HMM stationary distribution')
    print(HMM.stationary_distribution)

    #plot stationary distribution with errors
    pi = HMM.stationary_distribution
    piL,piR = HMM.sample_conf('stationary_distribution')
    plt.xlim((-0.5,1.5))
    plt.ylabel('stationary distribution')
    plt.xlabel('state')
    for i in range(2):
        plt.errorbar(i, pi[i], yerr=[[piL[i]], [piR[i]]],fmt='o',color=colors[i])

    plt.savefig('%s/stationary_distribution.png'%mutant)
    plt.clf()




