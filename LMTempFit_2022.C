#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"


void LMTempFit() {

	// Initializings for F value calc
	Int_t numbOfFVar = 100;
	Double_t factorF[numbOfFVar];
	Double_t stepsize = (1-3)/100;

	// Initializings for Chi2 test
 	Double_t chi2_best = 0;
 	Double_t factorF_best = 0;
 	Double_t G_par = 0;
 	Double_t V2_par = 0;
 	Double_t V3_par = 0;
	
	// Input histos 
	TH1D *hY; 
	TH1D *hY_LM;
 	TH1D *hY_a; // hY'

	// Harmonics in fit function
	Int_t NH = 2;

	// Function introductions with set parameters
	void LoadData_hY(infile_hY);
	void LoadData_hY_LM(infile_hY_LM);


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


	//	Fit function creation -> G(fourier) or the Y_ridge
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

 	// 	Multiplying, subtracting, fitting and Chi2 testing
 	for (int j = 0; j < numbOfFVar; j++) 
 	{
 		hY_LM->Scale(factorF[j]); 				// F * hY_LM
 		hY_a = hY->Add(hY_LM, -1); 				// hY' 
 		hY_a->Fit("fFit");						// Fitting j'th fit to hY'
 		Double_t min_val = Chi2(hY_a, fFit);	// Calculating Chi2 value 

 		// Saving best values
 		if (min_val < chi2_best) 
 		{ 
 			chi2_best = min_val; 
 			factorF_best = factorF[j];
			G_par = hY_a->GetParameter(0);
			V2_par = hY_a->GetParameter(1);
			V3_par = hY_a->GetParameter(2);	
 		}
 	}

 	// Outputs
 	cout << "Lowest Chi2 value: " << chi2_best << "\n\n" << endl;
 	cout << "PARAMETERS \n" << endl; 
 	cout << "F value: " << factorF_best << "\n" << endl;
 	cout << "G parameter: " << G_par << "\n" << endl;
 	cout << "V2,2: " << V2_par << "\n" << endl;
 	cout << "V3,3: " << V3_par << endl;
} 



/*	Chi2 Test
/		
/	Parameters: hY' , fFit -> "Yield and fit function" 
/	Returns: Double_t -> "Chi2 statistic value"  
*/
Double_t Chi2(TH1D *hY_a, TF1 *fFit) 
{
	Double_t chi2 = 0.0;

	for (int i = 0; i < h->GetNbinsX(); i++) 
	{
		Double_t bincent = hY_a->GetXaxis()->GetBinCenter(i);
		Double_t obs = hY_a->GetBinContent(i);
		Double_t exp = fFit->Eval(bincent);
		Double_t err = TMath::Power(hY_a->GetBinError(i), 2) + TMath::Power(fFit->GetParError(0), 2);

		if (err != 0) // Exclude sigma^2 = 0
		{
			chi2 = chi2 + (TMath::Power(obs - exp, 2) / err);
		}
	}

	return chi2;
}


/* 


NOTES:

- Chi2 statistical value should be checked


STEPS IN THE ALGORITHM:

1. Loads two input histos
2. Creates G(fourier) fit function
3. Gives fit some initial values
4. Multiplies hY_LM histo with F_i value
5. Substracts hY_LM from hY histo to create hY' histo
6. Fits G(fourier) to hY'
7. Calculates Chi2 value
8. Compares Chi2 value to the previous best Chi2 value
9. Prints fit parameters on screen
 
