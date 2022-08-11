import numpy as np
import ROOT
import scipy
from scipy import interpolate

import sys;
sys.path.append("JPyPlotRatio");
import JPyPlotRatio

fload = ROOT.TFile("CorrFit_2022.root","read"); #Opens figs

dataTypePlotParams = 	[
							{'plotType':'data','color':'k','fmt':'o','markersize':2.0},
							{'plotType':'theory','facecolor':'b','edgecolor':'b','alpha':0.5,'linestyle':'solid','linecolor':'b'},
							{'plotType':'theory','facecolor':'r','edgecolor':'r','alpha':0.5,'linestyle':'solid','linecolor':'r'}
						];

xtitle = ["$\\Delta\\varphi (\\mathrm{rad})$"];
ytitle = ["$\\frac{1}{N_{\\mathrm{trig}}}\\frac{\\mathrm{d}N^{\\mathrm{pair}}}{\\mathrm{d}\\Delta\\varphi}$"];


plot = JPyPlotRatio.JPyPlotRatio(panels=(1,1),
	panelsize=(8,8), # change the size
	#panelLabel="",
	panelLabelLoc=(0.2,1.1),panelLabelSize=10,panelLabelAlign="left",
	legendPanel=0,
	legendLoc=(0.2,1.1),legendSize=11,xlabel=xtitle[0],ylabel=ytitle[0]); # x- and y-coordinate labels

fig = fload.Get("CMSeta_projection_0_1;1");
fit = fload.Get("fit;1");
fit2 = fload.Get("CMSeta_projection_1_1;1");
#if (!fit_2) { print("No fit found!" ); }

data = plot.Add(0, fig, **dataTypePlotParams[0], label='signal');
data_fit = plot.Add(0, fit, **dataTypePlotParams[1], label='fit');
data_fit2 = plot.Add(0, fit2, **dataTypePlotParams[2], label='Y_HM-F*Y_MB');

Signal_entries = fig.GetEntries();
LM_entries = fit2.GetEntries();

fload.Close();

plot.GetPlot().text(0.35,0.98,"Peripheral 2-particle correlation",fontsize=12);
plot.GetPlot().text(0.2,0.83,"pp $\\sqrt{s}$ = ",fontsize=11);
plot.GetPlot().text(0.5, 0.90,"$ < p_\\mathrm{T,trig(assoc)} < \\,\\mathrm{GeV}/c$",fontsize=10);
plot.GetPlot().text(0.5, 0.925,"$-4.0 < \\eta < 4.0 $");
plot.GetPlot().text(0.05, 0.05, "Signal entries");
plot.GetPlot().text(0.05, 0.03, "Y_HM-F*Y_MB entries");
plot.GetPlot().text(0.2, 0.05, Signal_entries);
plot.GetPlot().text(0.2, 0.03, LM_entries);
plot.GetRatioAxes(0).xaxis.set_ticks_position('both');
plot.GetRatioAxes(0).yaxis.set_ticks_position('both');
#plot.Ratio(data_fit, data_fit2); # Plots theory vs data ratio

plot.Plot();


# plot.Save("figs/Fig14_FlowExt.pdf");
# plot.Save("figs/Fig14_FlowExt.png");
plot.Show();