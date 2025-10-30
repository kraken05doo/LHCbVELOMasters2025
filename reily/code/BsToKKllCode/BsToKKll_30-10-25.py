import ROOT
from ROOT import TFile, gSystem, gInterpreter
from ROOT import TH1D, TH2D, TCanvas, TChain
import os,sys
from math import * 
from os import path,listdir
import random
import array

basedir=path.dirname(path.realpath(__file__))

#from Selections import load_event_library
#load_event_library()

sys.path.insert(0,basedir)
from MCTools import * 
gInterpreter.AddIncludePath( f'{basedir}/../include')
gSystem.Load( f'{basedir}/../build/libEvent.so')

from ROOT import uParticle

def dira_bpv( particle, vertices, max_dt):
  vertex = vertices[0]
  ipw = lambda particle, vertex : particle.ip(vertex) if particle.ip_t(vertex) < max_dt else 9999

  ip = ipw(particle, vertex ) 
  for index in range(1, len(vertices ) ) :
    if ipw(particle, vertices[index]) < ip : 
      vertex = vertices[index]
      ip = particle.ip(vertices[index])
  dx = particle.firstState.x0 - vertex.x 
  dy = particle.firstState.y0 - vertex.y 
  dz = particle.firstState.z  - vertex.z
  p = particle.p4()
  return (dx * p.x() + dy * p.y() + dz * p.z() ) / sqrt( (dx**2 + dy**2 + dz**2 )*p.P2() )

events = TChain("Events")
#events.AddFile(("/disk/moose/general/djdt/lhcbUII_masters_2025/scoping_phijpsiee_1k/u2Selec_tuple_Long.root"))
events.AddFile(("/disk/moose/general/djdt/lhcbUII_masters_2025/scoping_phijpsimm_1k/u2Selec_tuple_Long.root")) #for muons

entry = 0

# create an empty file to write results to
hfile = TFile( 'h_test_results.root', 'RECREATE',)

phi_tuple = ROOT.TNtuple("phi_mass", "phi_mass", "mass:true")
Jpsi_tuple = ROOT.TNtuple("Jpsi_mass", "Jpsi_mass", "mass:true")

electron_momentum_tuple = ROOT.TNtuple("electron_momentum", "electron_momentum", "charge:px:py:pz:E")

for event in events:
    displaced_tracks = ROOT.select( event.Particles, event.Vertices, 250, 1500, 6 )

    #look at all PIDs
    #all_particles_ID = [ track.trueID for track in displaced_tracks]
    #print(all_particles_ID)

    #get electrons and kaons
    #electrons = [ track for track in displaced_tracks if abs(track.trueID) == 11 ]
    electrons = [ track for track in displaced_tracks if abs(track.trueID) == 13 ] #for muons
    kaons = [ track for track in displaced_tracks if abs(track.trueID) == 321 ]

    #get photons
    #photons = [ track for track in displaced_tracks if abs(track.trueID) == 22 ]
    #photons_from_electrons = [ photon for photon in photons if is_from(photon, event, 11) ]
    #print(f"photons: {len(photons)}, photons from electrons: {len(photons_from_electrons)}")

    #separate +ve/-ve kaons and electrons
    ep = [track for track in electrons if track.charge() > 0]
    em = [track for track in electrons if track.charge() < 0]
    kp = [track for track in kaons if track.charge() > 0]
    km = [track for track in kaons if track.charge() < 0]      

    entry += 1
    doca_cut = 0.10
    nPVs = npvs(event) 
    found_signal = False

    #reconstruct phi and Jpsi
    phi_candidates = ROOT.combine(kp, km, doca_cut, 15, 0)
    Jpsi_candidates = ROOT.combine(ep, em, doca_cut, 15, 0)

    phi_list = []
    phi_vtx_list = []
    Jpsi_list = []
    Jpsi_vtx_list = []

    print( f"Entry # = {entry}, Number of PVs Found =  {nPVs}, Number of Electrons Found = {len(electrons)} ({len(ep)} positive, {len(em)} negative), Number of Kaons found = {len(kaons)} ({len(kp)} positive, {len(km)} negative), Number of Phis created = {len(phi_candidates)}, Number of J/psis created = {len(Jpsi_candidates)}" )

    for i,phi_candidate in enumerate(phi_candidates):
        k1, k2, phi, phi_vtx = phi_candidate

        phi_list.append(phi)
        phi_vtx_list.append(phi_vtx)

        is_phi_signal = is_from(k1, event, 333) and is_from(k2, event, 333) #are the K's from a phi using MC PID

        #selections for phi
        pv_phi = phi.bpv_4d(event.Vertices)
        if phi_vtx.chi2 / phi_vtx.ndof > 5 : continue
        if k1.pt() + k2.pt() < 900 : continue
        if phi.mass < 900 or phi.mass > 1100 : continue
        if phi_vtx.chi2_distance(pv_phi) > 50 : continue
        if dira_bpv(phi,event.Vertices,0.050) < 0.9 : continue

        if is_phi_signal:
          phi_tuple.Fill( array.array("f", [ phi.mass * 0.001, True ]))
        else:
          phi_tuple.Fill( array.array("f", [ phi.mass * 0.001, False ]))
  
    for j,Jpsi_candidate in enumerate(Jpsi_candidates):
        e1, e2, Jpsi, Jpsi_vtx = Jpsi_candidate
        is_Jpsi_signal = is_from(e1, event, 443) and is_from(e2, event, 443) #are the e's from a J/psi using MC PID

        Jpsi_list.append(Jpsi)
        Jpsi_vtx_list.append(Jpsi_vtx)

        electron_momentum_tuple.Fill( array.array("f", [ e1.charge(), e1.p4().x() * 0.001, e1.p4().y() * 0.001, e1.p4().z() * 0.001, e1.p4().E() * 0.001 ]) )
        electron_momentum_tuple.Fill( array.array("f", [ e2.charge(), e2.p4().x() * 0.001, e2.p4().y() * 0.001, e2.p4().z() * 0.001, e2.p4().E() * 0.001 ]) )

        #selections for J/psi
        pv_Jpsi = Jpsi.bpv_4d(event.Vertices)
        if Jpsi_vtx.chi2 / Jpsi_vtx.ndof > 5 : continue
        if e1.pt() + e2.pt() < 2500 : continue
        if Jpsi.mass < 2500 or Jpsi.mass > 3500 : continue
        if Jpsi_vtx.chi2_distance(pv_Jpsi) > 50 : continue
        if dira_bpv(Jpsi,event.Vertices,0.050) < 0.9 : continue
        
        if is_Jpsi_signal:
          Jpsi_tuple.Fill( array.array("f", [ Jpsi.mass * 0.001, True ]))
        else:
          Jpsi_tuple.Fill( array.array("f", [ Jpsi.mass * 0.001, False ]))

hfile.cd()
phi_tuple.Write()
Jpsi_tuple.Write()
electron_momentum_tuple.Write()
hfile.Close()
