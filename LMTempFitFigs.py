import numpy as np
import ROOT
import scipy
from scipy import interpolate

import sys;
sys.path.append("JPyPlotRatio");
import JPyPlotRatio

fload = ROOT.TFile("outputs/CorrFit_2022.root","read"); #Opens figs

dataTypePlotParams = 	[
							{'plotType':'data','color':'k','fmt':'o','markersize':2.0},
							{'plotType':'theory','facecolor':'b','edgecolor':'b','alpha':0.3,'linestyle':'solid','linecolor':'b'},
							{'plotType':'data','color':'r','fmt':'s','markersize':1.0},
							{'plotType':'theory','facecolor':'y','edgecolor':'y','alpha':0.3,'linestyle':'solid','linecolor':'y'},
							{'plotType':'theory','facecolor':'g','edgecolor':'g','alpha':0.3,'linestyle':'solid','linecolor':'g'},
							{'plotType':'theory','facecolor':'r','edgecolor':'r','alpha':0.3,'linestyle':'solid','linecolor':'r'}
						];

xtitle = ["$\\Delta\\varphi (\\mathrm{rad})$"];
ytitle = ["$\\frac{1}{N_{\\mathrm{trig}}}\\frac{\\mathrm{d}N^{\\mathrm{pair}}}{\\mathrm{d}\\Delta\\varphi}$"];
nrow = 1;
ncol = 1;

'''
# Limits the panel 
xlimits = {0:(-1.4,4.5)};
ylimits = {0:(0.903,0.99)};
'''
rlimits = [(0.8,1.2)];

plot = JPyPlotRatio.JPyPlotRatio(panels=(1,1),
	panelsize=(13,13), # change the size
	#panelLabel="",
	#rowBounds=ylimits,  # for nrow
	#colBounds=xlimits,  # for ncol
	ratioBounds=rlimits,# for nrow
	panelLabelLoc=(0.2,1.1),panelLabelSize=13,panelLabelAlign="left",
	legendPanel=0,
	legendLoc=(0.2,1.12),legendSize=11,xlabel=xtitle[0],ylabel=ytitle[0]); # x- and y-coordinate labels

fig = fload.Get("CMSeta_projection_1_3;49");
fit = fload.Get("fFit;49");
fit_v2 = fload.Get("fit_v2");
fit_v3 = fload.Get("fit_v3");
#fit_G = fload.Get("F*Y+G;49");

data = plot.Add(0, fig, **dataTypePlotParams[0], label='signal');
data_fit = plot.Add(0, fit, **dataTypePlotParams[1], label='fit');
data_fit_v2 = plot.Add(0, fit_v2, **dataTypePlotParams[3], label='v2');
data_fit_v3 = plot.Add(0, fit_v3, **dataTypePlotParams[4], label='v3');
#data_fit_G = plot.Add(0, fit_G, **dataTypePlotParams[5], label='F*Y+G');

fload.Close();

plot.GetPlot().text(0.35,0.98,"Peripheral 2-particle correlation",fontsize=12);
plot.GetPlot().text(0.1,0.90,"pp $\\sqrt{s}$ = ",fontsize=11);
plot.GetPlot().text(0.1, 0.92,"$ < p_\\mathrm{T,trig(assoc)} < \\,\\mathrm{GeV}/c$",fontsize=10);
plot.GetPlot().text(0.5, 0.925,"$-4.0 < \\eta < 4.0 $");
plot.GetRatioAxes(0).xaxis.set_ticks_position('both');
plot.GetRatioAxes(0).yaxis.set_ticks_position('both');
plot.Ratio(data, data_fit); # Plots theory vs data ratio

plot.Plot();


# plot.Save("figs/Fig14_FlowExt.pdf");
# plot.Save("figs/Fig14_FlowExt.png");
plot.Show();