#!/usr/local/bin/python2.7

#
# Simple Example of using python  
# wrapper to BargerPropagator
# 

from BargerPropagator import * 

from ROOT import TH1D
from ROOT import TCanvas
from ROOT import ROOT

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

# Get the Propagator
b = BargerPropagator()

# Define the oscillation parameters 
x12    = 0.825 
x13    = 0.10
x23    = 1.0
m21    = 7.9e-5
mAtm   = 2.5e-3
delta  = 0.
energy = 0.650

# ks: 0 - sin2(2q) variables
# ks: 1 - sin2( q) variables
ks     = 0

# nutype:  1 - neutrino
# nutype: -1 - antineutrino
nutype = 1


h = TH1D('h' , '#nu_{#mu} #rightarrow #nu_{e} ' , 10000, 0.1, 20 )

# Compute some oscillation probabilities
for i in range(1, h.GetNbinsX() ):

  energy = h.GetBinCenter(i)

  b.SetMNS( x12, x13, x23, m21, mAtm, delta, energy, ks, nutype )
  b.propagateLinear( 1 , 295.0 , 2.6 )
  h.SetBinContent( i , b.GetProb(2,1) )


h.GetXaxis().SetTitle( ' E_{#nu} ' )
h.GetYaxis().SetTitle( ' Osc. Probability ' )

c = TCanvas()
c.SetFillColor(0)
c.SetLogx(1)

h.Draw()
c.SaveAs('t2k.beamline.test.png')


# end of script
##################


  





