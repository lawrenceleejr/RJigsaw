#!/usr/bin/env python -b 

# from ROOT import *

from rootpy import *
from ROOT import *

gSystem.Load("./RJigsawTool");


#stuff for getting gamma distribution from PDF's
type=0; #0-q qbar 1-gg

eta_1 = 0.27871;
eta_2 = 3.3627;
eps_u = 4.4343;
g_u = 38.599;

del_S = -0.11912;
eta_S = 9.4189;
eps_S = -2.6287;
g_S = 18.065;

A_g = 3.4055;
del_g = -0.12178;
eta_g = 2.9278;
eps_g = -2.3210;
g_g = 1.9233;
A_g1 = -1.6189;
del_g1 = -0.23999;
eta_g1 = 24.792;

P = TLorentzVector()


#give transverse momenta to CM system in lab frame?
PT = 0.1; #In units of sqrt{shat}

#gamma factor associated with 'off-threshold-ness' of tops
gamma = 1.2;

#Now, we also have the option to take gamma, event-by-event, 
#from a more realistic distribution
#to do this, set 'b_gamma' to true and the rest
#of the variables below appropriately
b_gamma = true;
rootS = 13.;
type = 1; #0 quark-antiquark  1 gluon-gluon
M = 175./1000.; #TeV units

#Number of toy events to throw
N = 100;

RJTool = RJigsawTool.RJigsawTool()

RJTool.initialize()
RJTool.setHemisphereMode(1)
RJTool.getHemisphereMode()


# #here, we set up the stuff for dynamic gamma
# h_gamma = TH1D("newgamma","newgamma",500,1.0,6.0);
# if b_gamma:
# 	print "generating gamma distribution"
# 	for ibin in xrange(500):
# 		g = h_gamma->GetBinCenter(ibin);
# 		entry = Calc_dsigma_dgamma(g);

# 		if entry > 0.:
# 			h_gamma->SetBinContent(ibin, entry);


# for i in xrange(N): // event generation loop
# 	if b_gamma:
#       gamma = h_gamma->GetRandom();


#     # Event generation

#     Mt1 = 175.;
# 	MW1 = 80.;
# 	Mt2 = 175.;
# 	MW2 = 80.;
# 	Mnu1 = 0.;
# 	Mnu2 = 0.;




