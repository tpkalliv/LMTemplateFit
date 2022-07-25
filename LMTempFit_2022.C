#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"

	Int_t variationsF = 100;
	Double_t factorF[variationsF];
	Double_t stepsize = (1-3)/100;
	
	// Input histos 
	TH1D *hY[variationsF]; 
	TH1D *hY_LM[variationsF];

	// Harmonics
	Int_t NH = 2;

	// Functions
	void LoadData_hY(infile_hY);
	void LoadData_hY_LM(infile_hYLM);
	void FitFNC();
	void Chi2();


void LMTempFit() {

	// Loading input data from hY 
	void LoadData_hY(TString filename_hY)
	{
		for (int i = 0; i < variationsF; i++) 
		{	
			// Loading given high multiplicity yield histo into hY histo
			TFile *fIn = TFile::Open(filename_hY,"read");
			hY[i] = (TH1D*)fIn->Get("hY");
		}
	}

	// Loading input data from hY_LM 
	void LoadData_hY_LM(TString filename_hY_LM)
	{
		for (int i = 0; i < variationsF; i++)
		{	
			// Loading given low multiplicity yield histo into hY_LM histo
			TFile *fIn_LM = TFile::Open(filename_hY_LM,"read");
			hY_LM[i] = (TH1D*)fIn_LM->Get("hY_LM");
		}
	}

	//	Fitting function
 	string cosine = "[0]*(1";
	for (int i=1; i<=NH; i++) {
		ostringstream app;
		app << "+2*[" << i << "]*TMath::Cos(" << i+1 << "*x"; // G*{1+2*Vn,n*Cos(n*[delta]Phi)}
		string append = app.str();
		cosine = cosine + append;
	}

	TF1* fFit = new TF1("fFit", cosine, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
	fFit->SetParName(0, "G");
	
	for (int i = 1; i <= NH; i++) 
	{
		fFit->SetParName(i, Form("V%d,%d", i)); // V2,2 and V3,3
	}

 	// Dividing interval [1, 3] into 100 step values
 	factorF_samp = 1;
 	for (int i = 0; i < variationsF; i++) 
 	{
 		factorF[i] = factorF_samp + stepsize; // F factor variations
 	}	

 	// Multiplying hY_LM with F 
 	for (int j = 0; j < variationsF; j++) 
 	{
 		hY_LM[j]->Scale(factorF[j]);
 	}
 	
 	// Subtracting hY_LM from hY to get hY'
 	TH1D *hY_a = hY->Clone("hY_a");

 	for (int j = 0; j < variationsF; j++) 
 	{
		hY_a[j] = hY[j]->Add(hY_LM[j], -1); // hY[j] + ( (-1) * hY_LM[j] ) 
 	}

 	FitFNC(); // Fittings

 	
 	int min_fit_id = Chi2(); // Chi2 estimation

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


/* 	G(fourier) fitting
/
/	Fits same fit function for different hY' histos
*/
void FitFNC()
{
	for (int i = 0; i < variationsF; i++) 
	{
		hY_a[i]->Fit("fFit"); // Fits G*(1 + 2*V2,2*Cos(2x) + 3*V3,3*Cos(3x)) 
	}
}


/*	Chi2 Estimation
/	
/	Returns: Int - Histogram index number with minimum Chi2 value  
*/
int Chi2() 
{

	Double_t min = 999; // chi2 value
	int fit_id = -999; // Best fit ID

	// For every hY' with with different F value, we calculate the chi2 value
	for (int j = 0; j < variationsF; j++) 
	{
		Double_t chi2_val = 0;

		// Goes through every bin in the current hY'
		for (int i = 0; i < hY_a[j]->GetNbinsx(); i++) 
		{
			// Calculates x-value for current bin (so to get fit function value at x)
			Double_t bincent = hY_a[j]->GetXAxis()->GetBinCenter(i);

			// Calculates (hY'-fv2)^2 / sigma^2 value for every bin 
			Double_t chi2_val += TMath::Pow((hY_a[j]->GetBinContent(i) - fv2->Eval(bincent)), 2) / 
			TMath::Pow(fv2->GetParError(0) + fv2->GetParError(1) + fv2->GetParError(2), 2) + TMath:Pow(hY_a[j]->GetBinError(i), 2);
			// Sigma^2 = (G_err + V2,2_err + V3,3_err)^2 + (hY_a[j]_err)^2
		}

		if (chi2_val < min) { min = chi2_val; fit_id = j;  } // Memorizes fit index
	}

	cout << "Minimum chi2 value: " << min << endl;

	return fit_id;

}


/* 

NOTES:

-