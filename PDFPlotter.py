#import sys, os
import numpy as np
#import scipy.stats as st
from scipy import interpolate
#import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.path import Path
#from matplotlib.ticker import MultipleLocator

from ROOT import gROOT, TCanvas, TF1, TFile, TH2D
import root_numpy

#sys.path.insert(2,'/global/homes/x/xxiang/gitlab/LZStats/Tools/PDFVisualiser/pyModules/')
#import pyVisualiser as pyv
from MyTools import *

'''
For every ROOT histograms
 1. Load root files and get histograms
 2. Convert histogram to numpy array
 3. Get bin cenerts from historam
 4. make the histograms finer (optional)
'''

exposure=5.6e6

# B8 Background
fBkg = TFile.Open("background.root")
hS1S2 = fBkg.Get("h-B8_fine") # per kg d
hS1S2.Scale(exposure)
bkg_Z = np.transpose(root_numpy.hist2array(hS1S2))
bkg_x, bkg_y = get_TH2_bin_centers(hS1S2)
bkg_X, bkg_Y = np.meshgrid(bkg_x, bkg_y)
bkg_x_fine, bkg_y_fine, bkg_Z_fine = make_finer_pdf_2d(bkg_x, bkg_y, bkg_Z)
bkg_X_fine, bkg_Y_fine = np.meshgrid(bkg_x_fine, bkg_y_fine)


# Backgrounds others than B8
fBkger = TFile.Open("/global/project/projectdirs/lz/data/LZStats_data/LowE_ER_group/root_pdfs/Background/AllBkg_MDC2.root")
er_str=["hatm", "hhep", "hDSN", "hPP", "h2vBB", "hKr", "hRn", "hRn220", "hDetE", "hDetN"]
for i in range(0, len(er_str)):
    if i==0:
        hS1S2er = fBkger.Get(er_str[i])
    else:
        hS1S2er.Add(fBkger.Get(er_str[i]))
hS1S2er.Scale(exposure)
bkger_Z = np.transpose(root_numpy.hist2array(hS1S2er))
bkger_x, bkger_y = get_TH2_bin_centers(hS1S2er)
bkger_ex, bkger_ey = get_TH2_bin_edges(hS1S2er)
bkger_X, bkger_Y = np.meshgrid(bkger_x, bkger_y)
plt.figure()
plt.pcolormesh(bkger_ex,bkger_ey, bkger_Z)
plt.savefig('figs/tmp_bkger.pdf')


# WIMP signal
fSig = TFile.Open("signal.root")
hSig = []
mass_list=[6,8,10,20]
for i in range(0, len(mass_list)):
    hist_name = str("h-wimp-m%d_fine" % mass_list[i])
    hSig.append(fSig.Get(hist_name))
    hSig[i].Scale(1./hSig[i].Integral())
    # hSig[i].Scale(exposure)

sig_x=[]; sig_y=[]
sig_X=[]; sig_Y=[]; sig_Z=[]
sig_X_fine=[]; sig_Y_fine=[]; sig_Z_fine=[];
for i in range(0, len(mass_list)):
    sig_Z.append(np.transpose(root_numpy.hist2array(hSig[i])))
    tmp_x, tmp_y = get_TH2_bin_centers(hSig[i])
    sig_x.append(tmp_x); sig_y.append(tmp_y)
    a,b = np.meshgrid(tmp_x, tmp_y)
    sig_X.append(a); sig_Y.append(b)
    tmp_x, tmp_y, tmp_Z = make_finer_pdf_2d(tmp_x, tmp_y, sig_Z[i])
    sig_Z_fine.append(tmp_Z)
    tmp_X, tmp_Y = np.meshgrid(tmp_x, tmp_y)
    sig_X_fine.append(tmp_X); sig_Y_fine.append(tmp_Y)


'''
    Plotting Contours
'''
fig = plt.figure()
sig_xvals=[];
sig_yvals=[];
sig_conts=[]
for i in range(0, len(mass_list)):
    a, b = sample_pdf_2d(sig_X_fine[i], sig_Y_fine[i], sig_Z_fine[i], nsamples=500000)
    sig_xvals.append(a);
    sig_yvals.append(b);
    sig_conts.append(plt.contour(sig_X_fine[i], sig_Y_fine[i], sig_Z_fine[i], 150))
#print(sig_xvals, sig_yvals)
#plt.hist2d(sig_xvals[2], sig_yvals[2], bins=[80, 50])
fig.savefig('figs/to_be_removed.pdf')
bkg_xvals, bkg_yvals = sample_pdf_2d(bkg_X_fine, bkg_Y_fine, bkg_Z_fine, nsamples=500000)
bkg_conts = plt.contour(bkg_X_fine, bkg_Y_fine, bkg_Z_fine, 150)


#plt.figure()
#plt.pcolormesh(bkger_x, bkger_y, bkger_Z)
#fig.savefig('figs/tmp_pcolormesh.pdf')

#Levs = sig_conts.levels
##for i in range(len(Levs)):
##    print(i, ": ", Levs[i])
colors=['blue','green', 'orange', 'purple']
labels=['6 GeV', '8 GeV', '10 GeV', '20 GeV'] #note: the output always outside->inside
fig = plt.figure()
all_sig_segs =[]
for m in range(0, len(sig_conts)):
    sig_segs = find_percentile_contour(sig_conts[m], sig_xvals[m], sig_yvals[m], tolerance=0.01)
    all_sig_segs.append(sig_segs)
    aColor = colors[m]
    aLabel = labels[m]
    if sig_segs is None:
        print("ERROR: problem with finding all the contours")
        continue;
    for i in range(0, len(sig_segs)):
        for j in range(0, len(sig_segs[i])):
            aGraph = sig_segs[i][j].T
            if j==0 and i==0:
                plt.plot(aGraph[0],aGraph[1], '-', color=aColor, label=aLabel)
            else:
                plt.plot(aGraph[0],aGraph[1], '-', color=aColor)


bkg_segs = find_percentile_contour(bkg_conts, bkg_xvals, bkg_yvals, tolerance=0.01)
for i in range(0, len(bkg_segs)):
    for j in range(0, len(bkg_segs[i])):
        aGraph = bkg_segs[i][j].T
        if (i==0 and j==0):
            plt.plot(aGraph[0],aGraph[1], '-', color='red', label='B8 Bkg')
        else:
            plt.plot(aGraph[0],aGraph[1], '-', color='red')

plt.legend()
plt.xlim(0, 20)
plt.ylim(2.5, 4.5)
fig.savefig("figs/temp_contours.pdf")


'''
    Print Numbers to Screen
'''
target_pct = ['99.7%','95%','68%']
#WIMP number (fraction) in different contours
print "WIMP fractions:"
print "{} {} {} {}".format("mass","pct", "num1", "num2")
for m in range(0, len(all_sig_segs)):
    sig_segs=all_sig_segs[m]
    for i in range(0, len(sig_segs)):
        num1=0; num2=0;
        for j in range(0, len(sig_segs[i])):
            num1+= hist2d_sum(sig_x[m], sig_y[m], sig_Z[m], sig_segs[i][j])
            num2+= hist2d_sum(sig_x[m], sig_y[m], sig_Z[m], bkg_segs[i][j])
        print "{} {} {:.2f} {:.2f}".format(labels[m], target_pct[i], num1, num2)

# ER background numbers in various contours
print "---------------------"
print "ER numbers:"
print "{} {} {} {}".format("mass","pct", "num1", "num2")
for m in range(0, len(all_sig_segs)):
    sig_segs=all_sig_segs[m]
    for i in range(0, len(sig_segs)):
        num1=0; num2=0;
        for j in range(0, len(sig_segs[i])):
            num1 +=hist2d_sum(bkger_x, bkger_y, bkger_Z, sig_segs[i][j]) 
            num2 +=hist2d_sum(bkger_x, bkger_y, bkger_Z, bkg_segs[i][j])
        print "{} {} {:.2f} {:.2f}".format(labels[m], target_pct[i], num1, num2)
