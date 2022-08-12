#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"


/*
	Program loads two data sets and finds the best fit for them using Chi2.
	Outputs Chi2 statistical value and also parameters from the best fit.
*/
void LMTempFit_2022() {

	// Options
	bool showChi = false; // When "true" shows chi2 for every histogram

	// Initializings 
	Int_t numbOfFVar = 100; // Number of F values
	Double_t factorF[numbOfFVar];
	Double_t stepsize = (3-1)/(double) 100;
 	Double_t chi2_best;
 	Double_t factorF_best;

 	/*
 	Double_t G_par;
 	Double_t V1_par;
 	Double_t V2_par;
 	Double_t V3_par;
 	Double_t V4_par;
 	Double_t V5_par;
 	*/
 	
 	Int_t indexVal;
	Int_t NH = 5;
 	TH1D*  hY_a[numbOfFVar];
	TF1 *fitvn[NH];
	Double_t vn[NH];
	Double_t vnError[NH];

	TString paramNames[NH+1] = {'const', 'v1', 'v2', 'v3', 'v4', 'v5'};
	Double_t params[NH+1];

 	// Opens data 
	TFile* fIn = new TFile ("Corr_1_3_GeV.root", "read");
 	
	// Loads MB and HM (CMS) data
	TH1D* hY;
	fIn->GetObject("CMSeta_projection_1_0", hY); // HM data - appropriate values: 1_0, 1_1, 1_3
	if (!hY) { std::cout << "Histogram NOT found!" << std::endl; }		

	TH1D* hY_MB;
	fIn->GetObject("CMSeta_projection_0_0", hY_MB); // // MB data - appropriate values: 1_0, 1_1, 1_3
	if (!hY_MB) { std::cout << "Histogram NOT found!" << std::endl; }

 	// Initializing Chi2 function
 	Double_t Chi2(TH1D *hY_a, TF1 *fFit, bool showChi);


	//	Fit function creation -> G(fourier) or the Y_ridge
 	string cosine = "[0]*(1";
	for (int i=1; i<=NH; i++) {
		ostringstream app;
		app << "+2*[" << i << "]*TMath::Cos(" << i << "*(x-[" << i + NH << "]))"; 
		string append = app.str();
		cosine = cosine + append;
	}
	cosine = cosine + ")";
	cout << cosine << endl;
	const char* cos = cosine.c_str();


	TF1* fFit = new TF1("fFit", cos, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);

	fFit->SetParName(0, "G_param");
	fFit->SetParLimits(0, -10, 10);
	fFit->SetParameter(0, 1);

	for (int i = 1; i <= NH; i++) 
	{
		fFit->SetParName(i, Form("V%d,%d", i, i)); // Param 1: V2,2 and Param 2: V3,3 ...
		fFit->SetParLimits(i, -1, 1);
	}
	
	// Initial values for parameters for fitting
	fFit->SetParameter(1, 0.1); 
	fFit->SetParameter(2, 0.06); 
	fFit->SetParameter(3, 0.03); 
	fFit->SetParameter(4, 0.03); 
	fFit->SetParameter(5, 0.03); 



	// Creating factor F values
 	for (int i = 0; i <= numbOfFVar; i++) 
 	{
 		factorF[i] = 1 + (i*stepsize);
 	}	

	TFile* fOut = new TFile ("CorrFit_2022.root", "recreate");


 	// 	Multiplying, subtracting, fitting and Chi2 testing
 	for (int j = 0; j < numbOfFVar; j++) 
 	{
 		hY_a[j] = (TH1D*) hY->Clone(); 
 		hY_a[j]->Add(hY_MB, -factorF[j]);
 		hY_a[j]->Fit("fFit", "W", "", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
	

 		//Double_t min_val = fFit->GetChisquare();
 		Double_t min_val = Chi2(hY_a[j], fFit, showChi);	

 		if (j == 0) chi2_best = min_val;

 		// Saving
 		if (min_val < chi2_best) 
 		{
 			chi2_best = min_val;
 			indexVal = j; // Index for the histo (hY - F*Y_MB) closest to theory G(Fourier)
 			factorF_best = factorF[j];
			G_par = fFit->GetParameter("G_param");
			V1_par = TMath::Sqrt(fFit->GetParameter("V1,1"));
			V2_par = TMath::Sqrt(fFit->GetParameter("V2,2")); // Squaring Vn,n to get Vn
			V3_par = TMath::Sqrt(fFit->GetParameter("V3,3"));
			V4_par = TMath::Sqrt(fFit->GetParameter("V4,4"));
			V5_par = TMath::Sqrt(fFit->GetParameter("V5,5"));
			hY_a[j]->Write(); // Saves fitted histogram
			fFit->Write();
 		}	
 	}


 	// Saving harmonics
	for (Int_t n=0; n<NH; n++)
	{
		TString formula = Form("[0]*(1 + 2*[1]*TMath::Cos(%d*x))",n);									
		fitvn[n]= new TF1(Form("fit_v%d", n),formula, -TMath::Pi()/2.0, 3.0/2.0*TMath::Pi());			
		vn[n] = fFit->GetParameter(n+1);																	
		vnError[n] = fFit->GetParError(n+1);																
		fitvn[n]->SetParameter(1,vn[n]);																
		fitvn[n]->SetParameter(0, G_par);
		fitvn[n]->Write();
	}
	

 	// Outputs
 	cout << "\n\n" << "Lowest Chi2: " << chi2_best << "\n" << endl;
 	cout << "PARAMETERS \n" << endl; 
 	cout << "F value: " << factorF_best << "\n" << endl;
 	cout << "G parameter: " << G_par << "\n" << endl;
 	cout << "V1,1: " << V2_par << "\n" << endl;
 	cout << "V2,2: " << V2_par << "\n" << endl;
 	cout << "V3,3: " << V3_par << "\n\n" << endl;
 	cout << "V4,4: " << V4_par << "\n\n" << endl;
 	cout << "V5,5: " << V5_par << "\n\n" << endl;
	cout << "Index: " << indexVal << "\n\n" << endl;

 	fIn->Close();
 	fOut->Close();
} 


/*	
	Chi2 Test

	Parameters: hY' , fFit -> "Yield and fit function" 
	Returns: Double_t -> "Chi2 statistic value"  

*/
Double_t Chi2(TH1D *hY_a, TF1 *fFit, bool chiOp) 
{
	Double_t chi2 = 0.0;

	for (int i = 1; i < hY_a->GetNbinsX(); i++) 
	{
		Double_t bincent = hY_a->GetXaxis()->GetBinCenter(i); // x-value for bin center
		Double_t obs = hY_a->GetBinContent(i); // bin value
		Double_t exp = fFit->Eval(bincent); // fit value for each bin center

		// Calculating errors 
		Double_t fit_G_Err = TMath::Power(fFit->GetParError(0), 2);
		Double_t fit_V1_err = TMath::Power(fFit->GetParError(1), 2);
		Double_t fit_V2_err = TMath::Power(fFit->GetParError(2), 2);
		Double_t fit_V3_err = TMath::Power(fFit->GetParError(3), 2);
		Double_t fit_V4_err = TMath::Power(fFit->GetParError(4), 2);
		Double_t fit_V5_err = TMath::Power(fFit->GetParError(5), 2);
		Double_t hist_err = TMath::Power(hY_a->GetBinError(i), 2);

		Double_t total_err = fit_G_Err + fit_V1_err + fit_V2_err + fit_V3_err + fit_V4_err + fit_V5_err + hist_err;
		Double_t val = obs - exp;
		Double_t chi2_temp = 0.0;

		if (total_err != 0) // Exclude zero denominator
		{
			chi2_temp = TMath::Power(val, 2) / total_err;
			chi2 = chi2 + chi2_temp;
		}
	}

	if ( chiOp ) cout << "Chi2 statistical value: " << chi2 << endl;
	return chi2;
}


/* 

STEPS IN THE ALGORITHM:

1. 	Loads two input histos
2. 	Creates G(fourier) fit function
3. 	Gives fit some initial values
4. 	Multiplies hY_MB histo with F_i value
5. 	Substracts hY_MB from hY histo to create hY' histo
6. 	Fits G(fourier) to hY'
7. 	Calculates Chi2 value
8. 	Compares Chi2 value to the previous best Chi2 value
9. 	Saves histos and fits to .root file
10.	Outputs best chi2 value and associated parameters on screen

*/