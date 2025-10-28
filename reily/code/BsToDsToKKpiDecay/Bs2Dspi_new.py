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

#define a function that returns "dira" with the best primary vertex (see the LHCb Glossary https://lhcb.github.io/glossary/glossary/D.html)
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

def kaon_eff( pT ):
  #returns kaon misidentification as pion efficiency (NOT CURRENTLY wrt momentum)
  k_to_pi_eff = 0.01
  return k_to_pi_eff

events = TChain("Events") # create a TChain for our event files

for i in range(666000,666001,1): #just one for now
    # add events to the tChain
    events.AddFile(("/disk/moose/general/djdt/lhcbUII_masters/dataStore/"
                    "Beam7000GeV-md100-nu38-VerExtAngle_vpOnly/13264021/"
                    "VP_U2_ParamModel-SX/SX_10um50s_75umcylindr3p5_nu38_Bs2Dspi_2111/"
                    f"moore/U2Tuple_u2_250um_4d-{i}-SX_10um50s_75umcylindr3p5_nu38_Bs2Dspi_2111.root"))
    
entry=0 # entry counter
plot = ROOT.TH1D("m_ds","",100,1.8,2.1) # create two empty plots
plot_TrueDs = ROOT.TH1D("m_ds_true","",100,1.8,2.1)
# you could do the above with a python histogram instead!

# plot phi mass
plot_phi = ROOT.TH1D("m_phi", "", 100,0.95,1.5)
plot_phi_true = ROOT.TH1D("m_phi_true", "", 100,0.95,1.5)

#plot K+pi mass
plot_Kplus_pi = ROOT.TH1D("m_Kplus_pi", "", 100,0.6,1.6)
plot_Kplus_pi_true = ROOT.TH1D("m_Kplus_pi_true", "", 100,0.6,1.6)

#plot K-pi mass
plot_Kminus_pi = ROOT.TH1D("m_Kminus_pi", "", 100,0.6,1.6)
plot_Kminus_pi_true = ROOT.TH1D("m_Kminus_pi_true", "", 100,0.6,1.6)

#plot Bs mass
plot_Bs = ROOT.TH1D("m_Bs", "", 100,5,5.5)
plot_Bs_true = ROOT.TH1D("m_Bs_true", "", 100,5,5.5)

# create an empty file to write results to
hfile = TFile( 'h_test_results.root', 'RECREATE',)

# create NTuples for Bs, Ds, phi, K-pi, K+pi
Bs_tuple = ROOT.TNtuple("Bs_mass", "Bs_mass", "mass:true")
Ds_tuple = ROOT.TNtuple("Ds_mass", "Ds_mass", "mass:true")
phi_tuple = ROOT.TNtuple("phi_mass", "phi_mass", "mass:true")
Kplus_pi_tuple = ROOT.TNtuple("K+pi_mass", "K+pi_mass", "mass:true")
Kminus_pi_tuple = ROOT.TNtuple("K-pi_mass", "K-pi_mass", "mass:true")

n_signal=0 # number of signal events counter

for event in events: # loop through each event
  # select some displaced tracks, see definition in selections/src/uParticle.cpp
  displaced_tracks = ROOT.select( event.Particles, event.Vertices, 250, 1500, 6 )
 
  # from within these displaced tracks, select the kaons and pions using their true ID
  good_pions = [ track for track in displaced_tracks if abs( track.trueID ) == 211 ]
  good_kaons = [ track for track in displaced_tracks if abs( track.trueID ) == 321 ]
  # think about how this would change if we implemented "realistic Particle identification", what happens to the result? 


  # particle misIDing simulation
  # creates a seed number for each xxon (pion/kaon)
  seed_pion = [ random.random() for track in displaced_tracks if abs( track.trueID ) == 211 ]
  seed_kaon = [ random.random() for track in displaced_tracks if abs( track.trueID ) == 321 ]

  # removes all xxons where the seed value is less than the efficiency
  realistic_pions = [ track for i,track in enumerate(good_pions) if seed_pion[i] > kaon_eff( track.pt() ) ]
  realistic_kaons = [ track for i,track in enumerate(good_kaons) if seed_kaon[i] > kaon_eff( track.pt() ) ]

  # gathers all of the removed xxons
  mis_IDed_pions = [ track for i,track in enumerate(good_pions) if seed_pion[i] < kaon_eff( track.pt() ) ]
  mis_IDed_kaons = [ track for i,track in enumerate(good_kaons) if seed_kaon[i] < kaon_eff( track.pt() ) ]

  for i,new_pion in enumerate(mis_IDed_kaons):
    mis_IDed_kaons[i].mass = 139.57
  for i,new_kaon in enumerate(mis_IDed_pions):
    mis_IDed_pions[i].mass = 493.68

  # inserts removed xxons into the other list (i.e. one has been misIDed as the other)
  realistic_pions.extend(mis_IDed_kaons)
  realistic_kaons.extend(mis_IDed_pions)
  # the rest of the code has been amended to use realistic_xxons instead of good_xxons


  #if len(realistic_kaons) > 0 : print(realistic_kaons[0].mass)
  #if len(realistic_pions) > 0 : print(realistic_pions[0].mass)
  # change kaon particles to have kaon mass
  #for i,kaon in enumerate(realistic_kaons):
  #  realistic_kaons[i].mass = 493.7
  #if len(realistic_kaons) > 0 : print(realistic_kaons[0].mass)


  # separate the kaons into positive and negative
  # kp = [track for track in good_kaons if track.charge() > 0 ]
  # km = [track for track in good_kaons if track.charge() < 0 ]
  kp = [track for track in realistic_kaons if track.charge() > 0 ]
  km = [track for track in realistic_kaons if track.charge() < 0 ]


  doca_cut = 0.10 # cut value for DOCA cut (see https://lhcb.github.io/glossary/glossary/D.html)
  entry = entry + 1
  nPVs = npvs( event ) 
  found_signal = False
#   print( f"Entry # = {entry}, Number of PVs Found =  {nPVs}, Number of Pions Found = {len(good_pions)}, Number of Kaons found = {len(good_kaons)}")
  
  phi_candidates = ROOT.combine( kp, km, doca_cut, 15, 0) # build phi -> KK candidates using combine,  see definition in selections/src/uParticle.cpp
  print( f"Entry # = {entry}, Number of PVs Found =  {nPVs}, True Number of Pions Found = {len(good_pions)}, True Number of Kaons found = {len(good_kaons)} ({len(kp)} positive, {len(km)} negative), Number of Phis created = {len(phi_candidates)}")
  print( f"Entry # = {entry}, Number of PVs Found =  {nPVs}, Number of Pions Found = {len(realistic_pions)}, Number of Kaons found = {len(realistic_kaons)} ({len(kp)} positive, {len(km)} negative), Number of Phis created = {len(phi_candidates)}")

  # for pion in good_pions : # loop through our pions
  for i,pion in enumerate(realistic_pions) : # loop through our pions

    for my_phi in phi_candidates: # loop through our phi candidates for each pion
    #   print(help(my_phi))
      k1,k2,phi,phi_vtx = my_phi
      ds_vtx = ROOT.uVertex( [k1,k2,pion] ) # build the Ds common vertex
      ds     = ROOT.uParticle( [k1,k2,pion] ) # build the Ds particle

      bs_vtx = [ ROOT.uVertex([ds,pion2]) for j,pion2 in enumerate(realistic_pions) if j != i ] # and (ds.charge() * pion2.charge() == -1) ]
      bs = [ ROOT.uParticle([ds,pion2]) for j,pion2 in enumerate(realistic_pions) if j != i ] # and (ds.charge() * pion2.charge() == -1) ]    

      #phi = ROOT.uParticle( [k1,k2] ) # build the phi particle

      if k1.charge() > 0: # build k1 + pion particle
        Kplus_pi = ROOT.uParticle( [k1,pion] )
      elif k1.charge() < 0:  
        Kminus_pi = ROOT.uParticle( [k1,pion] )

      if k2.charge() > 0: # build k2 + pion particle
        Kplus_pi = ROOT.uParticle( [k2,pion] )
      elif k2.charge() < 0:  
        Kminus_pi = ROOT.uParticle( [k2,pion] )


      # print(k1.mass, k2.mass)

      # build a boolean statement instructing us whether these three particles did come from a Ds
      is_signal = is_from(k1, event, 431) and is_from(k2, event, 431) and is_from(pion, event, 431)
      is_signal_bs = [ is_signal  # k1,k2,pi1 from Ds
                        and (is_from(k1, event, 531) and is_from(k2, event, 531) and is_from(pion, event, 531)) #k1,k2,pi1 froom Bs
                        and is_from(pion2, event, 531) for j,pion2 in enumerate(realistic_pions,i+1) ]  # pi2 from Bs
          # investigate potential Bs being filtered out by this

      ## Selection cuts on the quality of the vertex and the tracks we select to construct it
      if ds_vtx.chi2 / ds_vtx.ndof > 5 : continue
      if k1.pt() + k2.pt() + pion.pt() < 1800 : continue
      if ds.mass < 1800 or ds.mass  > 2100 : continue 

      pv  = ds.bpv_4d( event.Vertices ) # find definition in selections/src/uParticle.cpp

      # More selections, look up what the chi2_distance is! definition in selections/src/uVertex.cpp
      if ds_vtx.chi2_distance(pv) < 50 : continue
      if dira_bpv(ds,event.Vertices,0.050)  < 0.9 : continue 

      plot.Fill( ds.mass * 0.001 ) # if all of these signal selections are passed, fill the histograms
      plot_phi.Fill( phi.mass * 0.001 )
      plot_Kplus_pi.Fill( Kplus_pi.mass * 0.001 )
      plot_Kminus_pi.Fill( Kminus_pi.mass * 0.001 )
      for single_bs in bs:
        plot_Bs.Fill( single_bs.mass * 0.001 )
      # plotting in units of GeV!

      found_signal |= is_signal 

      # if it passes the is_signal requirments we set, also plot it in the true plot
      if is_signal: 
        print("Adding a Ds")
        plot_TrueDs.Fill( ds.mass * 0.001)
        plot_phi_true.Fill( phi.mass * 0.001 )
        plot_Kplus_pi_true.Fill( Kplus_pi.mass * 0.001 )
        plot_Kminus_pi_true.Fill( Kminus_pi.mass * 0.001)

        Ds_tuple.Fill( array.array("f", [ ds.mass * 0.001, True ]))
        phi_tuple.Fill( array.array("f", [ phi.mass * 0.001, True ]))
        Kplus_pi_tuple.Fill( array.array("f", [ Kplus_pi.mass * 0.001, True ]))
        Kminus_pi_tuple.Fill( array.array("f", [ Kminus_pi.mass * 0.001, True ]))

      else: 
        # otherwise print info about why this failed!
        print( "Background")
        print_mc_particle( k1, event.MCParticles) 
        print_mc_particle( k2, event.MCParticles) 
        print_mc_particle( pion, event.MCParticles) 

        Ds_tuple.Fill( array.array("f", [ ds.mass * 0.001, False ]))
        phi_tuple.Fill( array.array("f", [ phi.mass * 0.001, False ]))
        Kplus_pi_tuple.Fill( array.array("f", [ Kplus_pi.mass * 0.001, False ]))
        Kminus_pi_tuple.Fill( array.array("f", [ Kminus_pi.mass * 0.001, False ]))
      
      #signal test for bs
      for k,single_bs in enumerate(bs):
        ## Selection for Bs
        if bs_vtx[k].chi2 / bs_vtx[k].ndof > 5 : continue
        if k1.pt() + k2.pt() + pion.pt() + realistic_pions[k].pt() < 5250 : continue
        if single_bs.mass < 5250 or single_bs.mass  > 5500 : continue 
        if ds.mass < 1940 or ds.mass > 2000 : continue

        pv_bs  = single_bs.bpv_4d( event.Vertices )

        if bs_vtx[k].chi2_distance(pv_bs) < 50 : continue
        if dira_bpv(single_bs,event.Vertices,0.050)  < 0.9 : continue

        if is_signal_bs[k]:
          plot_Bs_true.Fill( single_bs.mass * 0.001 )
          Bs_tuple.Fill( array.array("f", [ single_bs.mass * 0.001, True ]))
        else:
          Bs_tuple.Fill( array.array("f", [ single_bs.mass * 0.001, False ]))


  n_signal = n_signal + found_signal 

hfile.cd()

#titles
plot.SetTitle("D_{s} mass")
plot_TrueDs.SetTitle("True D_{s} mass")
plot.GetXaxis().SetTitle("D_{s} mass [GeV]")
plot_TrueDs.GetXaxis().SetTitle("True D_{s} mass [GeV]")
plot.GetYaxis().SetTitle("Count")
plot_TrueDs.GetYaxis().SetTitle("Count")

plot_phi.SetTitle("Phi mass")
plot_phi_true.SetTitle("True Phi mass")
plot_phi.GetXaxis().SetTitle("Phi mass [GeV]")
plot_phi_true.GetXaxis().SetTitle("True phi mass [GeV]")
plot_phi.GetYaxis().SetTitle("Count")
plot_phi_true.GetYaxis().SetTitle("Count")

plot_Kplus_pi.SetTitle("K+pi mass")
plot_Kplus_pi_true.SetTitle("True K+pi mass")
plot_Kplus_pi.GetXaxis().SetTitle("K+pi mass [GeV]")
plot_Kplus_pi_true.GetXaxis().SetTitle("True K+pi mass [GeV]")
plot_Kplus_pi.GetYaxis().SetTitle("Count")
plot_Kplus_pi_true.GetYaxis().SetTitle("Count")

plot_Kminus_pi.SetTitle("K-pi mass")
plot_Kminus_pi_true.SetTitle("True K-pi mass")
plot_Kminus_pi.GetXaxis().SetTitle("K-pi mass [GeV]")
plot_Kminus_pi_true.GetXaxis().SetTitle("True K-pi mass [GeV]")
plot_Kminus_pi.GetYaxis().SetTitle("Count")
plot_Kminus_pi_true.GetYaxis().SetTitle("Count")

plot_Bs.SetTitle("B_{s} mass")
plot_Bs_true.SetTitle("True B_{s} mass")
plot_Bs.GetXaxis().SetTitle("B_{s} mass [GeV]")
plot_Bs_true.GetXaxis().SetTitle("True B_{s} mass [GeV]")
plot_Bs.GetYaxis().SetTitle("Count")
plot_Bs_true.GetYaxis().SetTitle("Count")

plot.Draw()
plot_TrueDs.Draw()
plot_phi.Draw()
plot_phi_true.Draw()
plot_Kplus_pi.Draw()
plot_Kplus_pi_true.Draw()
plot_Kminus_pi.Draw()
plot_Kminus_pi_true.Draw()
plot_Bs.Draw()
plot_Bs_true.Draw()

plot.SetDirectory(0)
plot_TrueDs.SetDirectory(0)
plot_phi.SetDirectory(0)
plot_phi_true.SetDirectory(0)
plot_Kplus_pi.SetDirectory(0)
plot_Kplus_pi_true.SetDirectory(0)
plot_Kminus_pi.SetDirectory(0)
plot_Kminus_pi_true.SetDirectory(0)
plot_Bs.SetDirectory(0)
plot_Bs_true.SetDirectory(0)

plot.Write()
plot_TrueDs.Write()
plot_phi.Write()
plot_phi_true.Write()
plot_Kplus_pi.Write()
plot_Kplus_pi_true.Write()
plot_Kminus_pi.Write()
plot_Kminus_pi_true.Write()
plot_Bs.Write()
plot_Bs_true.Write()

Bs_tuple.Write()
Ds_tuple.Write()
phi_tuple.Write()
Kplus_pi_tuple.Write()
Kminus_pi_tuple.Write()

hfile.Close()
# more plotting options and tutorials here: https://root.cern/doc/v636/hsimple_8py.html

