#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"


void LMTempFit() {

	Int_t numbOfFVar = 100;
	Double_t factorF[numbOfFVar];
	Double_t stepsize = (1-3)/100;
	
	// Input histos 
	TH1D *hY; 
	TH1D *hY_LM;
 	TH1D *hY_a; // hY'

	// Harmonics
	Int_t NH = 2;

	// Functions
	void LoadData_hY(infile_hY);
	void LoadData_hY_LM(infile_hY_LM);
	void FitFNC();
	Double_t Chi2();


	// Loading input data from hY (high multiplicity yield)
	void LoadData_hY(TString filename_hY)
	{
			TFile *fIn = TFile::Open(filename_hY,"read");
			hY = (TH1D*)fIn->Get("hY");
	}

	// Loading input data from hY_LM (low multiplicity yield)
	void LoadData_hY_LM(TString filename_hY_LM)
	{
			TFile *fIn_LM = TFile::Open(filename_hY_LM,"read");
			hY_LM = (TH1D*)fIn_LM->Get("hY_LM");
	}



	//	Fit function G(fourier)
 	string cosine = "[0]*(1";
	for (int i=1; i<=NH; i++) {
		ostringstream app;
		app << "+2*[" << i << "]*TMath::Cos(" << i+1 << "*x)"; // G*{1+2*Vn,n*Cos(n*[delta]Phi)}
		string append = app.str();
		cosine = cosine + append;
	}

	cosine = cosine + ")";
	cout << cosine << endl;
	const char* cos = cosine.c_str();

	TF1* fFit = new TF1("fFit", cos, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
	fFit->SetParName(0, "G");

	for (int i = 1; i <= NH; i++) 
	{
		fFit->SetParName(i, Form("V%d,%d", i)); // V2,2 and V3,3
	}

	// Initial values for parameters for fitting
	fFit->SetParameter(0, 1000);
	fFit->SetParameter(1, 0.1);
	fFit->SetParameter(2, 0.01); 


 	// Dividing interval [1, 3] into 100 step values
 	factorF_samp = 1;
 	for (int i = 0; i < numbOfFVar; i++) 
 	{
 		factorF[i] = factorF_samp + stepsize; // F factor variations
 	}	

 	// Initializing for Chi2 test
 	Double_t chi2_best = 0;
 	Double_t factorF_best = 0;

 	/* 	Multiplying, subtracting, fitting and Chi2 testing
	/
	/	Result: 
	*/
 	for (int j = 0; j < numbOfFVar; j++) 
 	{
 		hY_LM->Scale(factorF[j]); 							// F * hY_LM
 		hY_a = hY->Add(hY_LM, -1); 							// hY' 
 		hY_a->Fit("fFit")									// Fitting j'th fit to hY'
 		Double_t min_val = Chi2(hY_a);						// Estimating fit 

 		if (min_val < chi2_best) 
 		{ 
 			chi2_best = min_val; 
 			factorF_best = factorF[j] 
 		}
 	}

 	/*	Output
 	/	
 	/	Parameters: G , V2,2 , V3,3 , F 
 	*/
 	cout << Form("Best fit was hY_a%d \n", min_fit_id) << endl; 

 	cout << "Parameters are: \n" << endl;

 	for (int i = 1; i <= 3; i++) 
 	{
 		cout << Form("Param%d: ", i) << hY_a[min_fit_id]->GetParameter(i); << "\n" << endl;	
 		cout << "F value: " << factorF[min_fit_id] << endl;
 	}
 	
} 



/*	Chi2 Calc
/		
/	Parameters: hY' -> "For fitting"
/	Returns: Double_t -> "Chi2 statistic value"  
*/
Double_t Chi2(TH1D *hY_a) 
{
	Double_t chi2_val = 0.0;

	// Goes through every bin in the current hY'
	for (int i = 0; i < hY_a->GetNbinsX(); i++) 
	{
		// Calculates x-value for current bin to get fit function value at x
		Double_t bincent = hY_a->GetXaxis()->GetBinCenter(i);

		// Calculates (hY'-fv2)^2 / Sigma^2 value per bin and adds them up for chi2_val
		// Sigma^2 = (G_err + V2,2_err + V3,3_err)^2 + (hY_a_err)^2
		Double_t chi2_val = chi2_val + (TMath::Pow((hY_a->GetBinContent(i) - fFit->Eval(bincent)), 2) / 
		TMath::Power(fFit->GetParError(0) + fFit->GetParError(1) + fFit->GetParError(2), 2) + TMath:Power(hY_a->GetBinError(i), 2));
	}
			
	return chi2_val;
}


/* 

NOTES:

- Chi2 statistical value should be checked


