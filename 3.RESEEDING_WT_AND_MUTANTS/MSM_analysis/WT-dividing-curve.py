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
file_paths['WT-pro'] = '../MSM-analysis/transform_all_eight/'
file_paths['D671N-pro'] = '../MSM-analysis/transform_all_eight/D671N-pro/'
file_paths['Y755A-pro'] = '../MSM-analysis/transform_all_eight/Y755A-pro/'
file_paths['Y759A-pro'] = '../MSM-analysis/transform_all_eight/Y759A-pro/'

#First we are analyzing the WT data to get the line separating our states

print('calculating curve from WT-pro HMM')

system = 'WT-pro'

Y = np.load('%s/tica_projection.npy'%file_paths[system])
Y1 = [y[:,0] for y in Y]
Y2 = [y[:,1] for y in Y]
    
clkmeans_clustercenters = np.load('%s/clkmeans_clustercenters.npy'%file_paths[system])
clkmeans_dtrajs = np.load('%s/clkmeans_dtrajs.npy'%file_paths[system])
    
clkmeans_dtrajs = clkmeans_dtrajs.tolist()
    
HMM = msm.bayesian_hidden_markov_model(clkmeans_dtrajs,nstates=2,lag=100,nsamples=1000,mincount_connectivity=100)

print('calculating probability of each frame to be in each state')

state1_weights = lambda x: HMM.observation_probabilities[0][x]
state2_weights = lambda x: HMM.observation_probabilities[1][x]

discrete_traj_state1_weights = list(map(state1_weights, HMM.discrete_trajectories_obs))
discrete_traj_state2_weights = list(map(state2_weights, HMM.discrete_trajectories_obs))

plt.figure(figsize=(8,5))
mplt.plot_free_energy(np.hstack(Y1),np.hstack(Y2),
                      weights=np.hstack(discrete_traj_state1_weights),cmap='Reds_r',
                      vmax=40, cbar_label=None)
plt.title('%s - DFG out'%system)
plt.xlabel('TIC 1',fontsize=12)
plt.xticks(fontsize=12)
plt.xlim((-75,87.5))
plt.ylabel('TIC 2',fontsize=12)
plt.yticks(fontsize=12)
plt.ylim((-60,60))
plt.gca().invert_xaxis()

plt.savefig('DFG-OUT-WT-pro.png',dpi=300,bbox_inches='tight')

plt.clf()
plt.figure(figsize=(8,5))
mplt.plot_free_energy(np.hstack(Y1),np.hstack(Y2),
                      weights=np.hstack(discrete_traj_state2_weights),cmap='Blues_r',
                      vmax=40, cbar_label=None)
plt.title('%s - DFG in'%system)
plt.xlabel('TIC 1',fontsize=12)
plt.xticks(fontsize=12)
plt.xlim((-75,87.5))
plt.ylabel('TIC 2',fontsize=12)
plt.yticks(fontsize=12)
plt.ylim((-60,60))
plt.gca().invert_xaxis()

plt.savefig('DFG-IN-WT-pro.png',dpi=300,bbox_inches='tight')

print('calculating which state is most likely at each frame')

z_1, xedge_1, yedge_1 = np.histogram2d(np.hstack(Y1),np.hstack(Y2),bins=100, weights=np.hstack(discrete_traj_state1_weights))
x_1 = 0.5*(xedge_1[:-1] + xedge_1[1:])
y_1 = 0.5*(yedge_1[:-1] + yedge_1[1:])
kT=1.0
F_1 = -kT * np.log(z_1)
# minener zero
F_1 -= np.min(F_1)
# avoid zero count
zmin_nonzero = np.min(z_1[np.where(z_1 > 0)])
z_1 = np.maximum(z_1, zmin_nonzero)
extent_1 = [yedge_1[0], yedge_1[-1], xedge_1[0], xedge_1[-1]]

z_2, xedge_2, yedge_2 = np.histogram2d(np.hstack(Y1),np.hstack(Y2),bins=100, weights=np.hstack(discrete_traj_state2_weights))
x_2 = 0.5*(xedge_2[:-1] + xedge_2[1:])
y_2 = 0.5*(yedge_2[:-1] + yedge_2[1:])
F_2 = -kT * np.log(z_2)
# minener zero
F_2 -= np.min(F_2)
# avoid zero count
zmin_nonzero = np.min(z_2[np.where(z_2 > 0)])
z_2 = np.maximum(z_2, zmin_nonzero)
extent_2 = [yedge_2[0], yedge_2[-1], xedge_2[0], xedge_2[-1]]

# Note that x_1 = x_2 and y_1 = y_2, double check if you are skeptical

probability_1 = F_1.T
probability_2 = F_2.T

binary_probability_state_1 = np.zeros([100,100])
binary_probability_state_2 = np.zeros([100,100])

for i,entry in enumerate(probability_2):
    for j,value in enumerate(entry):
        if value > probability_1[i][j]:
            binary_probability_state_1[i][j] = probability_1[i][j]
        else:
            binary_probability_state_1[i][j] = np.inf

for i,entry in enumerate(probability_2):
    for j,value in enumerate(entry):
        if value < probability_1[i][j]:
            binary_probability_state_2[i][j] = probability_2[i][j]
        else:
            binary_probability_state_2[i][j] = np.inf

plt.clf()

plt.figure(figsize=(4,3))

plt.contourf(x_1, y_1, binary_probability_state_1, 100, extent=extent_1, cmap='Reds_r',offset=-1)
plt.contourf(x_1, y_1, binary_probability_state_2, 100, extent=extent_1, cmap='Blues_r',offset=-1)

plt.xlabel('TIC 1',fontsize=12)
plt.xticks(fontsize=12)
plt.xlim((-75,87.5))
plt.ylabel('TIC 2',fontsize=12)
plt.yticks(fontsize=12)
plt.ylim((-60,60))
plt.gca().invert_xaxis()

plt.savefig('Boundary-WT-pro.png',dpi=300,bbox_inches='tight')

print('calculating fit line to the boundary between these two states')

boundary_100 = np.zeros([100,100])
for i,entry in enumerate(probability_2):
    for j,value in enumerate(entry):
        if probability_1[i][j]-0.7 < value < probability_1[i][j]+0.7:
            boundary_100[i][j] = probability_2[i][j]
        else:
            boundary_100[i][j] = np.inf

x_boundary = []
y_boundary = []
for i,entry in enumerate(boundary_100):
    for j,value in enumerate(entry):
        if boundary_100[i][j] == np.inf:
            continue
        else:
            x_boundary.append(x_1[j])
            y_boundary.append(y_1[i])

import scipy.interpolate as inter
s2 = inter.UnivariateSpline(x_boundary[::-5], y_boundary[::-5],s=5e5)

linear_x = np.linspace(-75,75, 101)

plt.clf()

plt.figure(figsize=(4,3))

plt.contourf(x_1, y_1, binary_probability_state_1, 100, extent=extent_1, cmap='Reds_r',offset=-1)
plt.contourf(x_1, y_1, binary_probability_state_2 , 100, extent=extent_2, cmap='Blues_r',offset=-1)
plt.scatter(x_boundary, y_boundary,color ='green',marker='.')
plt.plot (linear_x, s2(linear_x), 'k--', label='Spline, fit')

plt.xlabel('TIC 1',fontsize=12)
plt.xticks(fontsize=12)
plt.xlim((-75,87.5))
plt.ylabel('TIC 2',fontsize=12)
plt.yticks(fontsize=12)
plt.ylim((-60,60))
plt.gca().invert_xaxis()

plt.savefig('Boundary-fitspline-WT-pro.png',dpi=300,bbox_inches='tight')

