#include "TFile.h"

void testing() 
{

	TFile* fin = new TFile ("Corr_1_3_GeV.root", "read");

	TH1D* hY = (TH1D*) fin->Get("CMSeta_projection_0_0");

	hY->Draw();

	for (int i = 0; i < hY->GetNbinsX(); i++) cout << hY->GetBinContent(i) << endl;
}