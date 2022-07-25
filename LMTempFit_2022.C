#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"



/*----------------------INITIALIZINGS------------------------*/

	// For F value variations
	Int_t variationsF = 100;
	Double_t factorF[variationsF];
	Double_t stepsize = (1-3)/100;
	
	// Input histos (per-trigger-particle yields), initializing 200 histos for finding F 
	TH1D *hY[variationsF]; 
	TH1D *hY_LM[variationsF];

	Int_t NH = 2;

	// Functions for loading the histogram data
	void LoadData_hY(infile_hY);
	void LoadData_hY_LM(infile_hYLM);

	// Function for fitting G(fourier)
	void FitFNC();

	// Chi-squared function for estimating the best parameters 
	void Chi2();


// PROGRAM STARTS HERE
void LMTempFit() {


/* -------------------LOADING INPUTS---------------------------*/

	// Loading input data from hY 
	void LoadData_hY(TString filename_hY)
	{
		// Get one hY histo for every F value
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
		// Get one hY_LM histo for every F value
		for (int i = 0; i < variationsF; i++)
		{	
			// Loading given low multiplicity yield histo into hY_LM histo
			TFile *fIn_LM = TFile::Open(filename_hY_LM,"read");
			hY_LM[i] = (TH1D*)fIn_LM->Get("hY_LM");
		}
	}

/*--------------------PRODUCING FIT FUNCTION------------------*/

	// Creating G(fourier) function for fitting

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
		fFit->SetParName(0, "G");
		fFit->SetParameter(i, Form("V%d,%d", i)); // V2,2 and V3,3
	}



/*-------------------CREATING F FACTOR VARIATIONS-------------*/

 	// Dividing [1, 3] interval into 100 step values
 	factorF_samp = 1;
 	for (int i = 0; i < variationsF; i++) 
 	{
 		factorF[i] = factorF_samp + stepsize;
 	}	


/*-------------------MULTIPLYING HISTOS WITH F FACTOR---------*/

 	// Goes through every factor F value and creates 100 F*hY_LM histos with different bin values
 	for (int j = 0; j < variationsF; j++) 
 	{
 		for (int i = 0; i < hY_LM->GetNbinsx(i); i++) 
 		{
			hY_LM_binVal = hY_LM->GetBinContent(i); // fetches the i'th bin center value from x-axis
			hY_LM->SetBinContent(i, hY_LM_binVal * factorF[j]); // multiplies i'th bin value with factorF value (F*hY_LM)
		}

 	}


/*-------------------SUBTRACTING F*hY_LM from hY--------------*/

 	// Cloning hY to create hY'
 	TH1D *hY_a = hY->Clone("hY_a");

 	for (int i = 1; i <= hY->GetBinSize(); i++) 
 	{
 		hY_a->SetBinContent(i, hY->GetBinContent(i) - hY_LM->GetBinContent(i)); // Subtracting bin contents and adding them into hY'
 	}

/*-------------------FITTING-----------------------------------*/

 	// Fits
 	FitFNC();


/*-------------------CHI2 ESTIMATION---------------------------*/

 	// Chi2 estimation, returns ID (index number) for best fit
 	double min_fit_id = Chi2();


/*-------------------OUTPUTS-----------------------------------*/

 	cout << Form("Best fit was hY%d \n", min_fit_id) << endl; 

 	cout << "Parameters are: \n" << endl;

 	for (int i = 1; i <= 3; i++) 
 	{
 		
 		cout << Form("Param%d: ", i) << hY_a[min_fit_id]->GetParameter(i); << "\n" << endl;	
 	}
 	
 


} // PROGRAM ENDS HERE




/*-------------------CALLED FUNCTIONS---------------------------*/

// Fitting:
void FitFNC()
{

	// Goes through all different hY' histos 
	for (int i = 0; i < variationsF; i++) 
	{
		hY_a[i]->Fit("fFit");
	}

}


// Chi2 test: 
double Chi2() 
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

			// Calculates (hY'-fv2)^2 / fv2 value for every bin and adds them up (this might be wrong way to do chi2)
			Double_t chi2_val += TMath::Pow((hY_a[j]->GetBinContent(i) - fv2->Eval(bincent)), 2) / fv2->Eval(bincent);
		}

		if (chi2_val < min) { min = chi2_val; fit_id = j;  } // Memorizes fit index
	}

	cout << "Minimum chi2 value: " << chi2_val << endl;

	return fit_id;

}