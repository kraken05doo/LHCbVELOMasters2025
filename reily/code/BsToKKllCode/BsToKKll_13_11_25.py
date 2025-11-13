
import ROOT
from ROOT import TFile, gSystem, gInterpreter
from ROOT import TH1D, TH2D, TCanvas, TChain
import os,sys
from math import * 
from os import path,listdir
import random
import array
import numpy as np
#import argparser

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

def MCselect( input, vertices, minP, maxZ):
  selected_particles = []
  selected_vtx = []
  for i,particle in enumerate(input):
    p = particle.p4()
    vtx_i = particle.vertexIndex
    if vtx_i == -1 : return False
    z = vertices[vtx_i].pos.z()
    pAbs = np.sqrt( p.x() ** 2 + p.y() ** 2 + p.z() ** 2 )

    if z < maxZ and pAbs > minP:
      selected_particles.append(particle)
      selected_vtx.append(vertices[vtx_i])
  return selected_particles, selected_vtx

def electron_photon_reconstruction(electrons, photons):
  photon_motherIndices = [ photon.motherMCParticleIndex for photon in photons if photon.p4().E() > 100 ]
  electron_indices = [ electron.mcParticleIndex for electron in electrons ]
  #print(photon_motherIndices)
  #print(electron_indices)

  used_photon_indices = []

  for i,photon_motherIndex in enumerate(photon_motherIndices):
    if photon_motherIndex in electron_indices:
      j = np.where( np.array(electron_indices) == photon_motherIndex )[0][0]
      #print(electrons[j].p4())
      e_charge = electrons[j].charge()

      p_new = electrons[j].p4() + photons[i].p4()
      electrons[j].firstState.tx = p_new.x() / p_new.z()
      electrons[j].firstState.ty = p_new.y() / p_new.z()
      electrons[j].firstState.qop = abs(1 / p_new.P()) * e_charge
      used_photon_indices.append(i)

      #print(electrons[j].p4())
      #print("recomb!")
  return electrons, used_photon_indices

#parser = argparse.ArgumentParser()
#parser.add_argument("-ev_i", metavar="eventIndex", required=True,
#                    help="event index to be processed")
#args = parser.parse_args()

#particle_type = "ee"
#event_index = args.ev_i
events = TChain("Events")
#event_file = f"/disk/moose/general/djdt/lhcbUII_masters_2025/U2SelecTuples/Bs2PhiJps{particle_type}/{event_index}/u2Selec_tuple_Long.root"
events.AddFile(("/disk/moose/general/djdt/lhcbUII_masters_2025/scoping_phijpsiee_1k/u2Selec_tuple_Long.root"))
#events.AddFile((event_file))

entry = 0
checked_entry = 0


# debug variable
DEBUG = True


# create an empty file to write results to
#hfile = TFile( f'test_results{event_index}.root', 'RECREATE',)
hfile = TFile( f'h_test_results.root', 'RECREATE',)

phi_tuple = ROOT.TNtuple("phi_mass", "phi_mass", "mass:true")
Jpsi_tuple = ROOT.TNtuple("Jpsi_mass", "Jpsi_mass", "mass:true")
Bs_tuple = ROOT.TNtuple("Bs_mass", "Bs_mass", "mass:true")

electron_momentum_tuple = ROOT.TNtuple("electron_momentum", "electron_momentum", "charge:px:py:pz:E")
kaon_momentum_tuple = ROOT.TNtuple("kaon_momentum", "kaon_momentum", "charge:px:py:pz:E")

photon_tuple = ROOT.TNtuple("photons", "photons", "x:y:z:E:P:PT")
photon_hist = ROOT.TH3D("photon pos", "photon pos", 100, -3000, 3000, 100, -3000, 3000, 100, 0, 15000)

mc_electron_momentum_tuple = ROOT.TNtuple("mcElectron_momentum", "mcElectron_momentum", "charge:px:py:pz:E")
mc_kaon_momentum_tuple = ROOT.TNtuple("mcKaon_momentum", "mcKaon_momentum", "charge:px:py:pz:E")

electron_count = 0
kaon_count = 0

#test_list = [7, 31, 86, 115, 124, 147, 169, 214, 238, 256, 258, 265, 291, 292, 299, 332, 339, 340, 363, 385, 418, 436, 448, 453, 466, 472, 488, 490, 515, 516, 520, 528, 570, 593, 600, 604, 608, 611, 633, 641, 666, 669, 677, 692, 699, 720, 721, 744, 757, 770, 777, 795, 914, 918, 927, 949, 952, 956, 971, 975]

for index,event in enumerate(events):
    #if not index in test_list : continue #testing Jpsi candidates

    particles = event.Particles
    vertices = event.Vertices
    MCparticles = event.MCParticles
    MCvertices = event.MCVertices

    displaced_tracks = ROOT.select( particles, vertices, 200, 350, 6 )
    MC_tracks, MC_vertices = MCselect( MCparticles, MCvertices, 100, 6000 )
    #print([ tracks.ID for tracks in MC_tracks ])

    entry += 1
    doca_cut = 0.10
    chi2_ndof = 15
    charge = 0
    nPVs = npvs(event) 
    found_signal = False

    #look at all PIDs
    if DEBUG and False:
      all_particles_ID = [ track.trueID for track in displaced_tracks]
      print(all_particles_ID)

    #get electrons and kaons
    electrons_not_reconstructed = [ track for track in displaced_tracks if abs(track.trueID) == 11 ]
    #electrons_not_reconstructed = [ track for track in displaced_tracks if abs(track.trueID) == 13 ] #for muons
    kaons = [ track for track in displaced_tracks if abs(track.trueID) == 321 and track.p4().P() > 1600 and np.sqrt(track.p4().x()**2 + track.p4().y()**2) > 350 ]

    if len(electrons_not_reconstructed) < 2 or len(kaons) < 2:
      if DEBUG and False:
        print(f"removing entry: {entry} \n")
      continue #remove all data samples where there isn't at least set of KKee
    
    print( f"electron charges: {[ f'{track.mcParticleIndex}  {track.charge()}' for track in electrons_not_reconstructed ]}")

    checked_entry += 1

    electron_count += np.size(electrons_not_reconstructed)
    kaon_count += np.size(kaons)

    #get photons
    photons = [ track for i,track in enumerate(MC_tracks) if abs(track.ID) == 22 and MC_vertices[i].type == 101 ]
    photon_vtxs = [ MC_vertices[i] for i,track in enumerate(MC_tracks) if abs(track.ID) == 22 and MC_vertices[i].type == 101 ]
    print(f"photons: {len(photons)}") #, photons from electrons: {len(photons_from_electrons)}")

    #get mc electrons and kaons
    if DEBUG and False:
      mc_electrons = [ track for i,track in enumerate(MC_tracks) if abs(track.ID) == 11 and abs(track.motherID) == 443 ]#or abs(track.GDmotherID) == 443 or abs(track.GGDmotherID) == 443 or abs(track.GGGDmotherID) == 443 ]
      mc_kaons = [ track for i,track in enumerate(MC_tracks) if abs(track.ID) == 321 and abs(track.motherID) == 333 ]#or abs(track.GDmotherID) == 333 or abs(track.GGDmotherID) == 333 or abs(track.GGGDmotherID) == 333 ]
    

    if DEBUG and False:
      for electron in mc_electrons:
        mc_electron_momentum_tuple.Fill(array.array("f", [np.sign(electron.ID), electron.p4().x(), electron.p4().y(), electron.p4().z(), electron.p4().E()]) )
      for kaon in mc_kaons:
        mc_kaon_momentum_tuple.Fill(array.array("f", [np.sign(kaon.ID), kaon.p4().x(), kaon.p4().y(), kaon.p4().z(), kaon.p4().E()]) )
      

    #brem recon
    electrons, used_photon_indices = electron_photon_reconstruction(electrons_not_reconstructed, photons)
    

    if DEBUG and False:
      for i in used_photon_indices:
        photon_tuple.Fill(array.array("f", [photon_vtxs[i].pos.x(), photon_vtxs[i].pos.y(), photon_vtxs[i].pos.z(), photons[i].p4().E(), photons[i].p4().P(), np.sqrt(photons[i].p4().x()**2 + photons[i].p4().y()**2)]) )
        photon_hist.Fill(photon_vtxs[i].pos.x(), photon_vtxs[i].pos.y(), photon_vtxs[i].pos.z())

      electrons_from_Jpsi = [ track for track in electrons if is_from(track, MCparticles, 443) ]
      kaons_from_phi = [ track for track in kaons if is_from(track, MCparticles, 333) ]

      for electron in electrons_from_Jpsi:
        electron_momentum_tuple.Fill( array.array("f", [electron.charge(), electron.p4().x(), electron.p4().y(), electron.p4().z(), electron.p4().E()]) )
      for kaon in kaons_from_phi:
        kaon_momentum_tuple.Fill( array.array("f", [kaon.charge(), kaon.p4().x(), kaon.p4().y(), kaon.p4().z(), kaon.p4().E()]) )
      

      #print( f"modified electron charges: {[ f'{track.mcParticleIndex}  {track.charge()}' for track in electrons ]}")


    #separate +ve/-ve kaons and electrons
    ep = [track for track in electrons if track.charge() > 0]
    em = [track for track in electrons if track.charge() < 0]
    kp = [track for track in kaons if track.charge() > 0]
    km = [track for track in kaons if track.charge() < 0]      

    #reconstruct phi and Jpsi
    phi_candidates = ROOT.combine(kp, km, doca_cut, chi2_ndof, charge)
    Jpsi_candidates = ROOT.combine(ep, em, doca_cut, chi2_ndof, charge)

    phi_list = []
    phi_vtx_list = []
    Jpsi_list = []
    Jpsi_vtx_list = []

    print( f"entry # = {index}, PVs =  {nPVs}, electrons = {len(electrons)} ({len(ep)} positive, {len(em)} negative), kaons = {len(kaons)} ({len(kp)} positive, {len(km)} negative), phis = {len(phi_candidates)}, J/psis = {len(Jpsi_candidates)}" )
    
    if DEBUG and False:
      print( f"electron energy: {[ track.p4().E() for track in em ]}, positron energy: {[ track.p4().E() for track in ep ]} ")
      print( f"electrons from Jpsi: {np.sum([ 1 if is_from(track, MCparticles, 443) else 0 for track in em ])}. positrons from Jpsi: {np.sum([ 1 if is_from(track, MCparticles, 443) else 0 for track in ep ])}")
      print( f"Jpsi ID for electrons: {[ MCparticles[track.mcParticleIndex].motherMCParticleIndex for track in em if is_from(track, MCparticles, 443) ]}, Jpsi ID for positrons: {[ MCparticles[track.mcParticleIndex].motherMCParticleIndex for track in ep if is_from(track, MCparticles, 443) ]}")

    for i,phi_candidate in enumerate(phi_candidates):
        k1, k2, phi, phi_vtx = phi_candidate

        phi_list.append(phi)
        phi_vtx_list.append(phi_vtx)

        is_phi_signal = is_from(k1, MCparticles, 333) and is_from(k2, MCparticles, 333) #are the K's from a phi using MC PID
        is_phi_signal *= are_from(k1, k2, MCparticles) #are the K's from the same phi using MCParticleIndex

        #print(is_from(k1, MCparticles, 333) and is_from(k2, MCparticles, 333), are_from(k1, k2, MCparticles))

        #selections for phi
        #pv_phi = phi.bpv_4d(vertices)
        #if phi_vtx.chi2 / phi_vtx.ndof > 5 : continue
        #if k1.pt() + k2.pt() < 900 : continue
        #if phi.mass < 900 or phi.mass > 1100 : continue
        #if phi_vtx.chi2_distance(pv_phi) > 50 : continue
        #if dira_bpv(phi,vertices,0.050) < 0.9 : continue

        if is_phi_signal:
          phi_tuple.Fill( array.array("f", [ phi.mass * 0.001, True ]))
          #print(f"Phi charge: {phi.charge()}")
        else:
          phi_tuple.Fill( array.array("f", [ phi.mass * 0.001, False ]))
  
    for j,Jpsi_candidate in enumerate(Jpsi_candidates):
        e1, e2, Jpsi, Jpsi_vtx = Jpsi_candidate

        Jpsi_list.append(Jpsi)
        Jpsi_vtx_list.append(Jpsi_vtx)

        is_Jpsi_signal = is_from(e1, MCparticles, 443) and is_from(e2, MCparticles, 443) #are the e's from a J/psi using MC PID
        is_Jpsi_signal *= are_from(e1, e2, MCparticles) #are the J/psi's from the same phi using MCParticleIndex

        #print(is_from(e1, MCparticles, 443) and is_from(e2, MCparticles, 443), are_from(e1, e2, MCparticles))

        #electron_momentum_tuple.Fill( array.array("f", [ e1.charge(), e1.p4().x() * 0.001, e1.p4().y() * 0.001, e1.p4().z() * 0.001, e1.p4().E() * 0.001 ]) )
        #electron_momentum_tuple.Fill( array.array("f", [ e2.charge(), e2.p4().x() * 0.001, e2.p4().y() * 0.001, e2.p4().z() * 0.001, e2.p4().E() * 0.001 ]) )

        #selections for J/psi
        #pv_Jpsi = Jpsi.bpv_4d(vertices)
        #if Jpsi_vtx.chi2 / Jpsi_vtx.ndof > 5 : continue
        #if e1.pt() + e2.pt() < 2500 : continue
        #if Jpsi.mass < 2500 or Jpsi.mass > 3500 : continue
        #if Jpsi_vtx.chi2_distance(pv_Jpsi) > 50 : continue
        #if dira_bpv(Jpsi,vertices,0.050) < 0.9 : continue
        
        if is_Jpsi_signal:
          Jpsi_tuple.Fill( array.array("f", [ Jpsi.mass * 0.001, True ]))
          #print(f"Jpsi charge: {Jpsi.charge()}")
        else:
          Jpsi_tuple.Fill( array.array("f", [ Jpsi.mass * 0.001, False ]))
    
    if DEBUG and False:
      print(f"actual Jpsi: {is_Jpsi_signal} , actual phi: {is_phi_signal}")
      print(Jpsi_list, phi_list)
      print([ Jpsi.mass * 0.001 for Jpsi in Jpsi_list ])

    #Bs_fake_charge = 2 # combine function incorrectly handles charge 0 particles, meaning Jpsi and phi are both given charge 1
    Bs_list = ROOT.combine(phi_list, Jpsi_list, doca_cut, chi2_ndof, charge)

    for i, Bs_candidate in enumerate(Bs_list):
      phi, Jpsi, Bs, Bs_vtx = Bs_candidate
      print("Bs candidate:")

      is_Bs_signal = is_from( phi, MCparticles, 531 ) and is_from( Jpsi, MCparticles, 531 )
      print("phi: ", MCparticles[phi.mcParticleIndex].motherID, MCparticles[phi.mcParticleIndex].GDmotherID, MCparticles[phi.mcParticleIndex].GGDmotherID, MCparticles[phi.mcParticleIndex].GGGDmotherID )
      print("Jpsi: ", MCparticles[Jpsi.mcParticleIndex].motherID, MCparticles[Jpsi.mcParticleIndex].GDmotherID, MCparticles[Jpsi.mcParticleIndex].GGDmotherID, MCparticles[Jpsi.mcParticleIndex].GGGDmotherID)
      print(is_Bs_signal)
      is_Bs_signal *= are_from( phi, Jpsi, MCparticles )
      print(is_Bs_signal, are_from( phi, Jpsi, MCparticles ))

      #selections for J/psi
      #pv_Bs = Bs.bpv_4d(vertices)
      #if Bs_vtx.chi2 / Bs_vtx.ndof > 5 : continue
      #if phi.pt() + Jpsi.pt() < 5200 : continue
      #if phi.mass < 980 or phi.mass > 1050 : continue
      #if Bs.mass < 5200 or Bs.mass > 5500 : continue
      #if Bs_vtx.chi2_distance(pv_Bs) > 50 : continue
      #if dira_bpv(Bs,vertices,0.050) < 0.9 : continue

      if is_Bs_signal:
        Bs_tuple.Fill( array.array("f", [ Bs.mass * 0.001, True ]))
      else:
        Bs_tuple.Fill( array.array("f", [ Bs.mass * 0.001, False ]))
    
    if DEBUG:
      print("\n\n")


print(f"total no. of electrons: {electron_count} \ntotal no. of kaons: {kaon_count} \nchecked entried: {checked_entry}")

hfile.cd()

phi_tuple.Write()
Jpsi_tuple.Write()
Bs_tuple.Write()
electron_momentum_tuple.Write()
kaon_momentum_tuple.Write()
mc_electron_momentum_tuple.Write()
mc_kaon_momentum_tuple.Write()
photon_tuple.Write()

photon_hist.SetTitle("Photon Brem pos")
photon_hist.GetXaxis().SetTitle("x")
photon_hist.GetYaxis().SetTitle("y")
photon_hist.GetZaxis().SetTitle("z")
photon_hist.Draw()
photon_hist.Write()

hfile.Close()
