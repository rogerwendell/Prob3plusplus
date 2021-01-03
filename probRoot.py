#!/usr/local/bin/python2.7

#
# Use python wrapper to make 2D oscillagrams
# 

import numpy
from   math  import *

from BargerPropagator import *

from ROOT import TFile
from ROOT import TH2D
from ROOT import TCanvas
from ROOT import ROOT
import ROOT

#ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

# Create the method to compute probability
bNu = BargerPropagator()

# Binning	
NBinsEnergy    = 200
Dm12ZenithNbin = 200 

# Zenith Angle Range
Dm12ZenithEdge = numpy.zeros(Dm12ZenithNbin+1)
cz_start       = -1.001
cz_end         =  1.101
cz_step        = ( cz_end - cz_start)/float(Dm12ZenithNbin)

# Energy Range
EnergyBins     = numpy.zeros(NBinsEnergy+1)
e_start        = 0.110000001
e_end          = 300.0
e_step         = log10(e_end/e_start)/float(NBinsEnergy)
   
# Oscillation Parameters
kSquared = True  # are we using sin^2(x) variables?
kNuBar  =  1
DM2     = 2.5e-3
Theta23 = 0.5
Theta13 = 0.0238
dm2     = 7.9e-5
Theta12 = 0.302
#dcp     = 4.7123
dcp     = 0

# Print useful information
print "Using\tDM2\t", DM2, \
    "\n\tTheta23\t", Theta23, \
    "\n\tTheta13\t", Theta13, \
    "\n\tdm2\t", dm2, \
    "\n\tTheta12\t", Theta12, \
    "\n\tdcp\t", dcp 

print "From [", e_start, "-", e_end, "] GeV"

# Setup the axis for the histograms

for i in range(0, NBinsEnergy):
    Entry = e_start*pow( 10.0 , float(i)*e_step )
    EnergyBins[i] =Entry

EnergyBins[NBinsEnergy] = EnergyBins[NBinsEnergy-1]*1.001

Dm12ZenithEdge[0]= cz_start*0.9999
for i in range(1, Dm12ZenithNbin):
    Dm12ZenithEdge[i] = Dm12ZenithEdge[0] + float(i)*cz_step

Dm12ZenithEdge[Dm12ZenithNbin] = Dm12ZenithEdge[Dm12ZenithNbin-1]*1.001

# NuE Histograms
NuEToNuE3f     = TH2D("NuEToNuE3f","3 Flavor P_{#nu_{e}#rightarrow#nu_{e}}",
                      NBinsEnergy -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge)
NuEToNuMu3f    = TH2D("NuEToNuMu3f","3 Flavor P_{#nu_{e}#rightarrow#nu_{#mu}}",
                      NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge)
NuEToNuTau3f   = TH2D("NuEToNuTau3f","3 Flavor P_{#nu_{e}#rightarrow#nu_{#tau}}",
                      NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge)
NuEToNuX3f     = TH2D("NuEToNuX3f","3 Flavor P_{#nu_{e}#rightarrow#nu_{x}}",
                      NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge)
   
#NuMu Histograms
NuMuToNuE3f   = TH2D("NuMuToNuE3f","3 Flavor P_{#nu_{#mu}#rightarrow#nu_{e}}",
                     NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge)
NuMuToNuMu3f  = TH2D("NuMuToNuMu3f","3 Flavor P_{#nu_{#mu}#rightarrow#nu_{#mu}}",
                     NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge)
NuMuToNuTau3f = TH2D("NuMuToNuTau3f","3 Flavor P_{#nu_{#mu}#rightarrow#nu_{#tau}}",
                     NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge)
NuMuToNuX3f   = TH2D("NuMuToNuX3f","3 Flavor P_{#nu_{#mu}#rightarrow#nu_{x}}",
                     NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge)

#NuTau Histograms   
NuTauToNuE3f  = TH2D("NuTauToNuE3f","3 Flavor P_{#nu_{#tau}#rightarrow#nu_{e}}",
                     NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge)
NuTauToNuMu3f = TH2D("NuTauToNuMu3f","3 Flavor P_{#nu_{#tau}#rightarrow#nu_{#mu}}",
                     NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge)
NuTauToNuTau3f= TH2D("NuTauToNuTau3f","3 Flavor P_{#nu_{#tau}#rightarrow#nu_{#tau}}",
                     NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge)
NuTauToNuX3f  = TH2D("NuTauToNuX3f","3 Flavor P_{#nu_{#tau}#rightarrow#nu_{x}}",
                     NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge)

# Two Flavor Histogram                                
NuMuToNuTau2f = TH2D("NuMuToNuTau2f","2 Flavor P_{#nu_{#mu}#rightarrow#nu_{#tau}}",
                     NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge)

# Fill the Histograms
for i in range(0, NBinsEnergy):
    energy = e_start*pow(10.0, float(i)*e_step)
    for j in range(0, Dm12ZenithNbin):
        
        cosineZ = cz_start + float(j)*cz_step
  
        bNu.SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, dcp , energy, kSquared, kNuBar ) 
        bNu.DefinePath( cosineZ, 25.00, True  )
        bNu.propagate( 1*kNuBar )
  
        NuEToNuE3f.Fill( energy, cosineZ, bNu.GetProb(1,1) ) 
        NuEToNuMu3f.Fill( energy, cosineZ, bNu.GetProb(1,2) ) 
        NuEToNuTau3f.Fill( energy, cosineZ, bNu.GetProb(1,3) ) 
        NuEToNuX3f.Fill( energy, cosineZ, 1.0 - bNu.GetProb(1,1) ) 

        NuMuToNuE3f.Fill( energy, cosineZ, bNu.GetProb(2,1) ) 
        NuMuToNuMu3f.Fill( energy, cosineZ, bNu.GetProb(2,2) ) 
        NuMuToNuTau3f.Fill( energy, cosineZ, bNu.GetProb(2,3) ) 
        NuMuToNuX3f.Fill( energy, cosineZ, 1.0 - bNu.GetProb(2,2) ) 

        NuTauToNuE3f.Fill( energy, cosineZ, bNu.GetProb(3,1) ) 
        NuTauToNuMu3f.Fill( energy, cosineZ, bNu.GetProb(3,2) ) 
        NuTauToNuTau3f.Fill( energy, cosineZ, bNu.GetProb(3,3) ) 
        NuTauToNuX3f.Fill( energy, cosineZ, 1.0 - bNu.GetProb(3,3) ) 

tmp = TFile("RawProb.root", "recreate")

NuEToNuE3f.Write()
NuEToNuMu3f.Write()
NuEToNuTau3f.Write()
NuEToNuX3f.Write()
   
NuMuToNuE3f.Write()
NuMuToNuMu3f.Write()
NuMuToNuTau3f.Write()
NuMuToNuX3f.Write()

NuTauToNuE3f.Write()
NuTauToNuMu3f.Write()
NuTauToNuTau3f.Write()
NuTauToNuX3f.Write()

NuMuToNuTau2f.Write()

tmp.Close()
   	
print "Happy oscillating! Job is finished."
