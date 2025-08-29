#include <iostream>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>


struct PenalizedLikelihood {
    PenalizedLikelihood(TH1D* hist_, TF1* func_, double true_sigma_, double sigma_error_) 
        : hist(hist_), func(func_), true_sigma(true_sigma_), sigma_error(sigma_error_) {}

    double operator() (const double *par) const {
        func->SetParameters(par);
        double NLL = 0;
        
        for(int i = 1; i <= hist->GetNbinsX(); i++){
            double mu = func->Eval(hist->GetBinCenter(i));
            mu = TMath::Max(mu, 1e-10);
            double n = hist->GetBinContent(i);
            if (n > 0) {
                NLL += (mu - n + n*(TMath::Log(n) - TMath::Log(mu)));
            } else {
                NLL += mu;
            }
        }
        double sigma = par[3];
        double penalty = 0.5 * TMath::Power((sigma - true_sigma) / sigma_error, 2);
        
        return NLL + penalty;
    }
    
    TH1D* hist; 
    TF1* func; 
    double true_sigma;
    double sigma_error;
};

void ex3_1() {
    TH1D *hist = new TH1D("hist", ";X;Y", 100, -50, 50);

    TRandom3 rand(0);
    for (int i = 0; i < 100000; ++i) {
        hist->Fill(rand.Uniform(-50, 50));
    }

    for (int i = 0; i < 700; ++i) {
        hist->Fill(rand.Gaus(0, 5));
    }
    
    TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*exp(-0.5*((x-[2])/[3])**2)", -50, 50);

    fitFunc->SetParameter(0, 1000);
    fitFunc->SetParameter(1, 10);
    fitFunc->SetParameter(2, 0);
    fitFunc->SetParameter(3, 5);
    
    std::cout << "=== No penalty function ===" << std::endl;
    TFitResultPtr fitResult1 = hist->Fit("fitFunc", "LS"); 
    
    double amplitude1 = fitFunc->GetParameter(1);
    double amplitudeErr1 = fitFunc->GetParError(1);
    double sigma1 = fitFunc->GetParameter(3);
    double sigmaErr1 = fitFunc->GetParError(3);
    
    std::cout << amplitude1 << " ± " << amplitudeErr1 << std::endl;
    std::cout <<  sigma1 << " ± " << sigmaErr1 << std::endl;
    
    std::cout << "\n=== With penalty function (error_s = 1) ===" << std::endl;
    
    PenalizedLikelihood penalizedLL(hist, fitFunc, 5.0, 1.0);
    
    ROOT::Fit::Fitter fitter;
    double par[4] = {1000, 10, 0, 5};
    
    fitter.FitFCN(4, penalizedLL, par, hist->GetNbinsX(), false);
    ROOT::Fit::FitResult result = fitter.Result();
    
    double amplitude2 = result.Parameter(1);
    double amplitudeErr2 = result.ParError(1);
    double sigma2 = result.Parameter(3);
    double sigmaErr2 = result.ParError(3);
    
    std::cout <<  amplitude2 << " ± " << amplitudeErr2 << std::endl;
    std::cout << sigma2 << " ± " << sigmaErr2 << std::endl;
    
    std::cout << sigmaErr1 << std::endl;
    std::cout << " (error_s=1): " << sigmaErr2 << std::endl;
    
    

    std::cout << "\n=== With penalty function  (error_s = 0.1) ===" << std::endl;
    
    PenalizedLikelihood penalizedLL2(hist, fitFunc, 5.0, 0.1);
    
    ROOT::Fit::Fitter fitter2;
    fitter2.FitFCN(4, penalizedLL2, par, hist->GetNbinsX(), false);
    ROOT::Fit::FitResult result2 = fitter2.Result();
    
    double amplitude3 = result2.Parameter(1);
    double amplitudeErr3 = result2.ParError(1);
    double sigma3 = result2.Parameter(3);
    double sigmaErr3 = result2.ParError(3);
    
    std::cout << amplitude3 << " ± " << amplitudeErr3 << std::endl;
    std::cout  << sigma3 << " ± " << sigmaErr3 << std::endl;
    
    std::cout << "\nNo penalty: " << sigmaErr1 << std::endl;
    std::cout << "error_s=1: " << sigmaErr2 << std::endl;
    std::cout << "error_s=0.1: " << sigmaErr3 << std::endl;
    

    TCanvas *c1 = new TCanvas("c1", "Fit Results", 1200, 800);
    hist->Draw();
    fitFunc->Draw("same");
    c1->Update();
}