
import numpy as np
import matplotlib.pyplot as plt
import mplhep as mpl
mpl.style.use("LHCb2")

import ROOT
from ROOT import TFile, gSystem, gInterpreter
from ROOT import TH1D, TH2D, TCanvas, TChain
import os,sys
from math import * 
from os import path,listdir
import random

basedir=path.dirname(path.realpath(__file__))

#from Selections import load_event_library
#load_event_library()

sys.path.insert(0,basedir)
from MCTools import * 
gInterpreter.AddIncludePath( f'{basedir}/../include')
gSystem.Load( f'{basedir}/../build/libEvent.so')

from ROOT import uParticle

def read_Ntuple(file_name, tuple_name, charge):
    '''reads the data from an NTuple of the shape (data,isTrue)
        isTrue - is the reconstructed particle actually that particle (using MC PID)
    file_name - name of the file to be read
    tuple_name - name of the tuple to be read
    charge - e+ or e-
    '''
    dataFile = ROOT.TFile.Open(file_name)
    dataNtuple = dataFile.Get(tuple_name)

    pT = [ np.sqrt(particle.px**2 + particle.py**2) for particle in dataNtuple if particle.charge == charge ]
    p = [ np.sqrt(particle.px**2 + particle.py**2 + particle.pz**2) for particle in dataNtuple if particle.charge == charge ]
    eta = [ np.arctanh(particle.pz / np.sqrt(particle.px**2 + particle.py**2 + particle.pz**2))
            for particle in dataNtuple if particle.charge == charge]
    theta = [ 2 * np.arctan( np.exp( -x ) ) for x in eta ]
    phi = [ np.arctan( particle.py / particle.px ) for particle in dataNtuple if particle.charge == charge ]

    return pT, p, eta, theta, phi

def plot_2d_hist(file_name, tuple_name, title, N, bin_range, fig):
    ''' plots several bar charts on the same axis
    file_name - name of file to be opened
    tuple_name - array of names of tuples to be extracted
    title - array of [title, xaxislabel, yaxislabel]
    N - number of bins
    bin_range - size of bins
    fig - matplotlib figure
    '''
    ax = fig.add_subplot(1,1,1)

    pT1, p1, eta1, theta1, phi1 = read_Ntuple(file_name,tuple_name, 1) #e+
    pT2, p2, eta2, theta2, phi2 = read_Ntuple(file_name,tuple_name, -1) #e-

    H, xedges, yedges = np.histogram2d(eta1, eta2, bins=N)#, range=[bin_range, bin_range])
    mpl.hist2dplot(H, xedges, yedges, ax=ax, cbar=True)

    ax.set_title(title[0])
    ax.set_xlabel(title[1])
    ax.set_ylabel(title[2])



file_name = f'{basedir}/../graphing/test_results_electron_pT.root'
tuple_name = "electron_momentum"
#labels = [ r'$\mu$+ pT', r'$\mu$- pT' ]
#title = [ r'$\mu$+ $\eta$ vs $\mu$- $\eta$', r'$\mu$+ $\eta$', r'$\mu$- $\eta$']
title = [ r'e+ $\eta$ vs e- $\eta$', r'e+ $\eta$', r'e- $\eta$']
N = 50
bin_range = [ 0, 0.25 ]
fig = plt.figure(figsize=(15,10))

plot_2d_hist(file_name, tuple_name, title, N, bin_range, fig)

fig.savefig(f'{basedir}/../graphing/e+eta_vs_e-eta_30_10_25.pdf', bbox_inches="tight")
#plt.show()