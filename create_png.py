#!/usr/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from scipy.stats import kde
from scipy import integrate

import gizmo_analysis as ga

import matplotlib
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True
#matplotlib.rcParams.update({'font.size': 10})
matplotlib.rcParams['font.family'] = 'serif'

import matplotlib as mpl
mpl.use('Agg')

from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
from matplotlib import ticker, cm, colors

from sys import argv
import time

NAME = argv[1]
I = int(argv[2])
OUTPUT = argv[3]
fn = NAME % I

vals = ga.par_path(os.getcwd())
M = vals[0]
R = vals[1]
Res = int(vals[5])

rc, G, Mc = R/1e3, 4.3e4, M/1e10
tff = np.pi/2 *np.sqrt(rc**3/G/Mc/2)
tunit = 206265*1000*1.5e8/(86400*365)
lim = R/1e3*2
ext = R/1e3*1.5
# Res = 64
nbase = Res**3
ids_toshow = np.arange(0, nbase-1, nbase/64**3, dtype=int) 
periodic = False

selected_BH = True
if selected_BH:
    sim = ga.simulation('./output/')
    interesting_BHs = sim.find_interesting_BHs()
    #print(interesting_BHs)

sp = ga.snapshot(fn)
with h5py.File(fn, 'r') as f:
    # print(i, end=' ')
    xyz=f['PartType0']['Coordinates'][()]
    if (xyz>0).all():
        periodic = True
        xyz[:,0] -= R/1e3
        xyz[:,1] -= R/1e3
    pids = f['PartType0']['ParticleIDs'][()]
    rho = f['PartType0']['Density'][()]
    pids_b = np.where(np.in1d(pids, ids_toshow))[0]
    
    #print(len(pids_b))
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(xyz[:,0][pids_b], xyz[:,1][pids_b], s=0.01, \
                c=np.log10(rho)[pids_b]+xyz[:,2][pids_b]/lim, alpha=0.2, cmap='Blues'\
                )
    
    if 'PartType4' in list(f.keys()):
        xyz=f['PartType4']['Coordinates'][()]
        pids = f['PartType4']['ParticleIDs'][()]
        pids_b = np.where(np.in1d(pids, ids_toshow))[0]
        #print(pids_b)
        ax.scatter(xyz[:,0][pids_b], xyz[:,1][pids_b], s=0.01,\
                    #c='orange', \
                    c=xyz[:,2][pids_b]/lim, alpha=0.2, cmap='Wistia_r', vmax=1, vmin=-1\
                    )
    
    if 'PartType5' in list(f.keys()):
        if selected_BH:
            xyz = []
            pids = interesting_BHs + nbase
            for idd in interesting_BHs:
                xyz.append(sp.single_bh(idd, 'Coordinates'))
            xyz = np.array(xyz)
        else:
            xyz=f['PartType5']['Coordinates'][()]
            pids = f['PartType5']['ParticleIDs'][()]
        if periodic:
            xyz[:,0] -= R/1e3
            xyz[:,1] -= R/1e3
        ids = pids - nbase
        bh_sink = f['PartType5']['SinkRadius'][()][0]/2.8
        ax.scatter(xyz[:,0], xyz[:,1], s=5, c=xyz[:,2]/lim, alpha=0.95, cmap='bone_r', vmax=1, vmin=-1)

    for j in range(len(xyz)):
        if np.abs(xyz[j,0]) < lim*0.95 and np.abs(xyz[j,1]) < lim*0.95:
            ax.text(xyz[j,0], xyz[j,1], str(ids[j]), fontsize=8)
    
    ax.set_aspect('equal')
    ax.text(-ext, +ext, r'\rm Res=$%d^3$'%Res)
    ax.text(+0.0, +ext, r'\rm $r_{\rm sink}$=%.2e pc'%(bh_sink*1000))
    ax.text(-ext, -ext, r'\rm $t$=%.2e yr'%(f['Header'].attrs['Time']*tunit))
    ax.text(+0.0, -ext, r'\rm $t$=%.2f $t_{\rm ff}$'%(f['Header'].attrs['Time']/tff))
    if not periodic:
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
    else:
        ax.set_xlim(-lim/2, lim/2)
        ax.set_ylim(-lim/2, lim/2)
    ax.set_xlabel(r'\rm $x$/kpc')
    ax.set_ylabel(r'\rm $y$/kpc')
    plt.savefig('%s%03d.png'%(OUTPUT, I), bbox_inches='tight')
    plt.close()

print("Finished processing #%d" % I)
