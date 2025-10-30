
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

def read_Ntuple(file_name, tuple_name, isTrue):
    '''reads the data from an NTuple of the shape (data,isTrue)
        isTrue - is the reconstructed particle actually that particle (using MC PID)
    file_name - name of the file to be read
    tuple_name - name of the tuple to be read
    isTrue - boolean, using only true Bs decays or all decays?
    '''
    dataFile = ROOT.TFile.Open(file_name)
    dataNtuple = dataFile.Get(tuple_name)

    mass = [ particle.mass for particle in dataNtuple if not isTrue or (isTrue and particle.true)]

    return mass

def plot_many_bar_graphs(file_name, tuple_names, hist_labels, title, N, bin_range, fig, isTrue):
    ''' plots several bar charts on the same axis
    file_name - name of file to be opened
    tuple_names - array of names of tuples to be extracted
    hist_labels - array of labels for the data sets, should be same length as many_data
    title - array of [title, xaxislabel, yaxislabel]
    N - number of bins
    bin_range - size of bins
    fig - matplotlib figure
    isTrue - list of booleans, using only true Bs decays or all decays?
    '''
    ax = fig.add_subplot(1,1,1)

    for i, tuple_name in enumerate(tuple_names):

        data = read_Ntuple(file_name,tuple_name,isTrue[i])

        #filtered_data = [ x for x in data if 1<x<1.05] # for if the peaks haven't been filtered already

        count = np.size(data)
        mean = np.sum(data)/count
        std_dev = np.sqrt( np.sum([ (x - mean)**2 for x in data ]) / (count - 1) )

        hist_labels[i] += f'\n  Count = {count} \n  Mean = {mean:.3f} \n  Std Dev = {std_dev:.3f} \n'
            
        mpl.histplot(*np.histogram(data, bins=N, range=bin_range), ax=ax, label=hist_labels[i])

    ax.set_title(title[0])
    ax.legend()
    ax.set_xlabel(title[1])
    ax.set_ylabel(title[2])


file_name = f'{basedir}/../graphing/test_results.root'
#tuple_names = [ "Ds_mass", "phi_mass", "K+pi_mass", "K-pi_mass", ]
tuple_names = [ "phi_mass", "phi_mass" ]
#hist_labels = [ r'$\text{D}_{\text{s}}$ mass', r'phi mass', r'$\text{K}^{\text{+}} \pi$ mass', r'$\text{K}^{\text{-}} \pi$ mass' ]
#hist_labels = [ r'J/$\psi$ mass', r'J/$\psi$ mass true' ]
hist_labels = [ r'Phi mass', r'Phi mass true' ]
#title = [ r'Reconstructed J/$\psi$ mass for muons no selection', r'J/$\psi$ mass [GeV]', r'Count']
title = [ r'Reconstructed phi mass for muons no selection', r'Phi mass [GeV]', r'Count']
N = 100
bin_range = [ 0, 4 ]
isTrue = [ False, True]
fig = plt.figure(figsize=(15,10))

plot_many_bar_graphs(file_name, tuple_names, hist_labels, title, N, bin_range, fig, isTrue)

fig.savefig(f'{basedir}/../graphing/histogram_phi_muons_no_selections_29_10_25.pdf', bbox_inches="tight")
#plt.show()