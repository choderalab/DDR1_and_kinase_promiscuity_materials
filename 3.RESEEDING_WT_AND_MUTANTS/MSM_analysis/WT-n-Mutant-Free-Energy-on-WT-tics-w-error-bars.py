import matplotlib
matplotlib.use('Agg')

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt

systems = ['WT-pro',
           'D671N-pro',
           'Y755A-pro',
           'Y759A-pro']

# These file paths are for the mutant data that has been analyzed on our WT-pro tics
file_paths = {}
file_paths['WT-pro'] = '../MSM-analysis/transform_all_eight/transform_all_eight/'
file_paths['D671N-pro'] = '../MSM-analysis/transform_all_eight/transform_all_eight/D671N-pro/'
file_paths['Y755A-pro'] = '../MSM-analysis/transform_all_eight/transform_all_eight/Y755A-pro/'
file_paths['Y759A-pro'] = '../MSM-analysis/transform_all_eight/transform_all_eight/Y759A-pro/'

#First we are analyzing the WT data to get the line separating our states

print('calculating line from WT-pro HMM')

system = 'WT-pro'

Y = np.load('%s/tica_projection.npy'%file_paths[system])
Y1 = [y[:,0] for y in Y]
Y2 = [y[:,1] for y in Y]
    
clkmeans_clustercenters = np.load('%s/clkmeans_clustercenters.npy'%file_paths[system])
clkmeans_dtrajs = np.load('%s/clkmeans_dtrajs.npy'%file_paths[system])
    
clkmeans_dtrajs = clkmeans_dtrajs.tolist()
    
HMM = msm.bayesian_hidden_markov_model(clkmeans_dtrajs,nstates=2,lag=100,nsamples=1000,mincount_connectivity=100)

hmm_sets = HMM.metastable_sets

state0_samples = [clkmeans_clustercenters[hmm_sets[0],0],clkmeans_clustercenters[hmm_sets[0],1]]
state1_samples = [clkmeans_clustercenters[hmm_sets[1],0],clkmeans_clustercenters[hmm_sets[1],1]]

X_HMM = np.concatenate((state0_samples,state1_samples), axis = 1)
Y_HMM = np.array([0]*len(state0_samples[0]) + [1]*len(state1_samples[0]))

import scipy
from sklearn import svm

C = 1.0  # SVM regularization parameter
clf = svm.SVC(kernel = 'linear',  gamma=0.7, C=C )
clf.fit(X_HMM.transpose(), Y_HMM)

plt.figure(figsize=(4,3))
plt.scatter(clkmeans_clustercenters[hmm_sets[0],0], clkmeans_clustercenters[hmm_sets[0],1], color='violet')
plt.scatter(clkmeans_clustercenters[hmm_sets[1],0], clkmeans_clustercenters[hmm_sets[1],1], color='red')
plt.xlabel('tic 1')
plt.ylabel('tic 2')
w = clf.coef_[0]
a = -w[0] / w[1]
xx = np.linspace(-60, 60)
yy = a * xx - (clf.intercept_[0]) / w[1]

plt.plot(xx, yy, 'k-')
plt.xlabel('TIC 1',fontsize=12)
plt.xticks(fontsize=12)
plt.xlim((-75,87.5))
plt.ylabel('TIC 2',fontsize=12)
plt.yticks(fontsize=12)
plt.ylim((-60,60))
plt.gca().invert_xaxis()

plt.savefig('EB-defining-line-WT-pro.png',dpi=300,bbox_inches='tight')

import seaborn as sns
sns.set_style("ticks")
import matplotlib as mpl

# Now that we've found our line we can just plot our free energy surfaces, reweighted by their MSM, for each mutant.

#Create a dictionary in which to collect free energy differences:
delG = {}
delG_nstates = {}
delG_nstates_prob = {}
delG_interval = {}

for system in systems:
    print('calculating free energy landscape for %s'%system)

    delG_nstates_prob[system] = []

    Y = np.load('%s/tica_projection.npy'%file_paths[system])
    Y1 = [y[:,0] for y in Y]
    Y2 = [y[:,1] for y in Y]

    clkmeans_dtrajs = np.load('%s/clkmeans_dtrajs.npy'%file_paths[system])
    clkmeans_dtrajs = clkmeans_dtrajs.tolist()

    MSM = msm.bayesian_markov_model(clkmeans_dtrajs, lag=200, nsamples=1000)

    #Calculate delG
    Y1_hstack = np.hstack(Y1)
    Y2_hstack = np.hstack(Y2)
    MSM_weights_hstack = np.hstack(MSM.trajectory_weights())    

    Y1_DFG_in = []
    Y2_DFG_in = []
    Y1_DFG_out = []
    Y2_DFG_out = []
    weights_DFG_in = []
    weights_DFG_out = []

    for i,y1 in enumerate(Y1_hstack):
        this_y1 = y1
        this_y2 = Y2_hstack[i]
        this_weight = MSM_weights_hstack[i]
        if this_y2 > a * y1 - (clf.intercept_[0]) / w[1]:
            Y1_DFG_out.append(this_y1)
            Y2_DFG_out.append(this_y2)
            weights_DFG_out.append(this_weight)
        else:
            Y1_DFG_in.append(this_y1)
            Y2_DFG_in.append(this_y2)
            weights_DFG_in.append(this_weight)

    delG[system] = -np.log(np.sum(weights_DFG_in)/np.sum(weights_DFG_out))

    plt.figure(figsize=(5,3))
    mplt.plot_free_energy(np.hstack(Y1),np.hstack(Y2),weights=np.hstack(MSM.trajectory_weights()),
                      cmap='pink',ncountours=11,vmax=9.6,cbar=False)
    plt.xlabel('TIC 1',fontsize=12)
    plt.xticks(fontsize=12)
    plt.xlim((-75,87.5))
    plt.xticks(np.arange(-75,76,25),['','','','','','',''])
    plt.ylabel('TIC 2',fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylim((-60,60))
    plt.yticks(np.arange(-60,61,20),['','','','','','',''])
    plt.plot(xx, yy, '--',color='0.25')
    plt.gca().invert_xaxis()
    ax, _ = mpl.colorbar.make_axes(plt.gca())
    cbar = mpl.colorbar.ColorbarBase(ax,cmap='pink',
                       norm=mpl.colors.Normalize(vmin=0, vmax=9.6))
    cbar.set_label('Free energy (kT)')
    cbar.set_clim(0,9.6)

    plt.savefig('EB-countour-MSM-line-%s.png'%system,dpi=300,bbox_inches='tight')

# Okay, so that was calculating the free energy difference using trajectory weights.
# Do we get a very different answer if instead we use the stationary distribution of the nstates?

    #In this case, instead of getting the tica space components of each frame in the trajectories
    # we only care about the tica space components of our nstates.

    clkmeans_clustercenters = np.load('%s/clkmeans_clustercenters.npy'%file_paths[system])

    Y1 = clkmeans_clustercenters[:,0]
    Y2 = clkmeans_clustercenters[:,1]

    #Calculate delG
    Y1_hstack = np.hstack(Y1)
    Y2_hstack = np.hstack(Y2)
    MSM_statdist_hstack = np.hstack(MSM.stationary_distribution)

    Y1_DFG_in = []
    Y2_DFG_in = []
    Y1_DFG_out = []
    Y2_DFG_out = []
    weights_DFG_in = []
    weights_DFG_out = []

    for i,y1 in enumerate(Y1_hstack):
        this_y1 = y1
        this_y2 = Y2_hstack[i]
        this_weight = MSM_statdist_hstack[i]
        if this_y2 > a * y1 - (clf.intercept_[0]) / w[1]:
            Y1_DFG_out.append(this_y1)
            Y2_DFG_out.append(this_y2)
            weights_DFG_out.append(this_weight)
        else:
            Y1_DFG_in.append(this_y1)
            Y2_DFG_in.append(this_y2)
            weights_DFG_in.append(this_weight)

    delG_nstates[system] = -np.log(np.sum(weights_DFG_in)/np.sum(weights_DFG_out))

# Okay, now we want to see if we can do this with the errors from the Bayesian sampling

    #Calculate delG
    Y1_hstack = np.hstack(Y1)
    Y2_hstack = np.hstack(Y2)
    #MSM_statdist_hstack = np.hstack(MSM.stationary_distribution)
    MSM_statdist_hstack = np.transpose(MSM.sample_f('stationary_distribution'))

    Y1_DFG_in = []
    Y2_DFG_in = []
    Y1_DFG_out = []
    Y2_DFG_out = []
    weights_DFG_in = []
    weights_DFG_out = []

    for i,y1 in enumerate(Y1_hstack):
        this_y1 = y1
        this_y2 = Y2_hstack[i]
        this_weight = MSM_statdist_hstack[i]
        if this_y2 > a * y1 - (clf.intercept_[0]) / w[1]:
            Y1_DFG_out.append(this_y1)
            Y2_DFG_out.append(this_y2)
            weights_DFG_out.append(this_weight)
        else:
            Y1_DFG_in.append(this_y1)
            Y2_DFG_in.append(this_y2)
            weights_DFG_in.append(this_weight)
    
    for i in range(len(weights_DFG_in[0])):
        delG_nstates_prob[system].append(-np.log(np.sum(np.transpose(weights_DFG_in)[i])/np.sum(np.transpose(weights_DFG_out)[i])))

    delG_interval[system] = np.percentile(a= delG_nstates_prob[system], q =[2.5,50.0,92.5])

# From this we can now calculate the free energy differences
print('plotting free energy bar graph')

index = np.arange(4)
width = 0.35

fig, ax = plt.subplots(figsize=(6,3.5))
rects1 = ax.barh(index, [delG_interval['WT-pro'][1],delG_interval['D671N-pro'][1],delG_interval['Y755A-pro'][1],delG_interval['Y759A-pro'][1]], 
               width,
               xerr=[[delG_interval['WT-pro'][1]-delG_interval['WT-pro'][0],delG_interval['D671N-pro'][1]-delG_interval['D671N-pro'][0],delG_interval['Y755A-pro'][1]-delG_interval['Y755A-pro'][0],delG_interval['Y759A-pro'][1]-delG_interval['Y759A-pro'][0]],
                     [delG_interval['WT-pro'][2]-delG_interval['WT-pro'][1],delG_interval['D671N-pro'][2]-delG_interval['D671N-pro'][1],delG_interval['Y755A-pro'][2]-delG_interval['Y755A-pro'][1],delG_interval['Y759A-pro'][2]-delG_interval['Y759A-pro'][1]]],
               color='C4', error_kw=dict(ecolor='gray',lw=2,capsize=5,capthick=2))

ax.set_xlabel(r'$\Delta G$ $(k_B T)$')
ax.set_xlim((-0.6,1.2))
ax.set_yticks(index )
ax.set_yticklabels(('WT-pro', 'D671N-pro','Y755A-pro','Y759A-pro'))

ax.invert_yaxis()
plt.savefig('EB-bargraph-MSM-line.png',dpi=300,bbox_inches='tight')


