import matplotlib
matplotlib.use('Agg')

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt

Y = np.load('../MSM-analysis/WT-pro-2macro-27-commute-lag100-stride-working/tica_projection.npy')

Y1 = [y[:,0] for y in Y]
Y2 = [y[:,1] for y in Y]

clkmeans_clustercenters = np.load('../MSM-analysis/WT-pro-2macro-27-commute-lag100-stride-working/clkmeans_clustercenters.npy')

clkmeans_dtrajs = np.load('../MSM-analysis/WT-pro-2macro-27-commute-lag100-stride-working/clkmeans_dtrajs.npy')

clkmeans_dtrajs = clkmeans_dtrajs.tolist()

HMM = msm.bayesian_hidden_markov_model(clkmeans_dtrajs,2,lag=100)

state_1_colors = HMM.metastable_distributions[0]/max(HMM.metastable_distributions[0])
state_2_colors = HMM.metastable_distributions[1]/max(HMM.metastable_distributions[1])

path_to_trajs = '../MSM-analysis/DDR1_WT_pro_plus_OG_trajectories/*.h5'
from glob import glob
filenames = sorted(glob(path_to_trajs))

print('These are our trajectories that flip.')
print(filenames[16],filenames[42],filenames[145],filenames[149])

import seaborn as sns
sns.set_style("ticks")

#We pick one of these trajectories and the frames of interest, as described in the directory 'picking_new_starting'

plt.figure(figsize=(5,3))
mplt.plot_free_energy(np.hstack(Y1),np.hstack(Y2),weights=np.hstack(HMM.trajectory_weights()),cmap='YlGnBu_r')
for i in range(len(clkmeans_clustercenters[:,0])):
    plt.scatter(clkmeans_clustercenters[i,0], clkmeans_clustercenters[i,1], color='red', s=50,alpha=(state_1_colors[i]))
    plt.scatter(clkmeans_clustercenters[i,0], clkmeans_clustercenters[i,1], color='violet', s=50,alpha=(state_2_colors[i]))
plt.plot([Y1[16][30],Y1[16][520]],[Y2[16][30],Y2[16][520]],color='yellow',linestyle='--')
plt.plot(Y1[16][30],Y2[16][30],'o',color='yellow',markersize=15,markeredgecolor='black',markeredgewidth=2)
plt.plot([Y1[16][520],Y1[16][3114]],[Y2[16][520],Y2[16][3114]],color='yellow',linestyle='--')
plt.plot(Y1[16][520],Y2[16][520],'o',color='yellow',markersize=15,markeredgecolor='black',markeredgewidth=2)
plt.plot([Y1[16][3114],Y1[16][3130]],[Y2[16][3114],Y2[16][3130]],color='yellow',linestyle='--')
plt.plot([Y1[16][3130],Y1[16][3845]],[Y2[16][3130],Y2[16][3845]],color='yellow',linestyle='--')
plt.plot([Y1[16][3845],Y1[16][5610]],[Y2[16][3845],Y2[16][5610]],color='yellow',linestyle='--')
plt.plot([Y1[16][5610],Y1[16][9906]],[Y2[16][5610],Y2[16][9906]],color='yellow',linestyle='--')
plt.plot([Y1[16][9906],Y1[16][10887]],[Y2[16][9906],Y2[16][10887]],color='yellow',linestyle='--')
plt.plot([Y1[16][10887],Y1[16][11113]],[Y2[16][10887],Y2[16][11113]],color='yellow',linestyle='--')
plt.plot(Y1[16][3114],Y2[16][3114],'o',color='yellow',markersize=15,markeredgecolor='black',markeredgewidth=2)
plt.plot(Y1[16][3130],Y2[16][3130],'o',color='0.8',markersize=15,alpha=0.8)
plt.plot(Y1[16][3845],Y2[16][3845],'o',color='0.8',markersize=15,alpha=0.8)
plt.plot(Y1[16][5610],Y2[16][5610],'o',color='0.8',markersize=15,alpha=0.8)
plt.plot(Y1[16][9906],Y2[16][9906],'o',color='0.8',markersize=15,alpha=0.8)
plt.plot(Y1[16][11113],Y2[16][11113],'o',color='0.8',markersize=15,alpha=0.8)
plt.plot(Y1[16][10887],Y2[16][10887],'o',color='yellow',markersize=15,markeredgecolor='black',markeredgewidth=2)
plt.xlabel('TIC 1',fontsize=12)
plt.xticks([])
plt.xlim((-75,87.5))
plt.ylabel('TIC 2',fontsize=12)
plt.yticks([])
plt.ylim((-60,60))
plt.gca().invert_xaxis()

plt.savefig('WT-Free-Energy-HMM.png',dpi=300,bbox_inches='tight')



