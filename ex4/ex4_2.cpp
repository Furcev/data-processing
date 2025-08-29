#include <iostream>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TStyle.h>

void ex4_2() {
    TFitResultPtr r, rold;
    const int N = 5000;
    const int Ntoys = 10000; 


    TH1D *hist_data = new TH1D("hist_data", "Data;X;Y", 100, -2.0, 2.0);
    TH1D *hist_toy = new TH1D("hist_toy", "Toy Data;X;Y", 100, -2.0, 2.0);
    TH1D *hfcn_true = new TH1D("hfcn_true", "-2logL (true model toys);#chi^{2}", 100, 30, 130);
    TH1D *hfcn_false = new TH1D("hfcn_false", "-2logL (false model toys);#chi^{2}", 100, 30, 130);


    TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x + [2]*x*x", -2.0, 2.0);
    gROOT->GetListOfFunctions()->Add(fitFunc);

    TRandom3 rand(0);
    double fcn;
    int nfcn = 0;
    int j;

    TCanvas *c1 = new TCanvas("c1", "Fit and ToyMC", 600, 600);

    for (int i = 0; i < N; ++i) {
        hist_data->Fill(rand.Gaus(0, 1));
    }
    hist_data->Draw("E");

    rold = hist_data->Fit("fitFunc", "LS");
    fcn = rold->MinFcnValue();

    std::cout << "\nReal data fit:" << std::endl;
    std::cout << "chi2 = " << rold->Chi2() << ", ndf = " << rold->Ndf()
              << ", p = " << rold->Prob() << ", MinFcnValue = " << fcn << std::endl;

    hfcn_false->Reset();
    for (int i = 0; i < Ntoys; ++i) {
        hist_toy->Reset();
        hist_toy->FillRandom("gaus", rand.Poisson(N));
        r = hist_toy->Fit("pol2", "LSQN");              
        hfcn_false->Fill(r->MinFcnValue());
    }

    hfcn_true->Reset();
    for (int i = 0; i < Ntoys; ++i) {
        hist_toy->Reset();
        hist_toy->FillRandom("fitFunc", rand.Poisson(N));
        r = hist_toy->Fit("pol2", "LSQN");
        double fcn_toy = r->MinFcnValue();
        hfcn_true->Fill(fcn_toy);
        if (fcn_toy > fcn) nfcn++;
    }


    for (j = 100; j > 0; --j) {
        if (hfcn_true->Integral(j, 100) > 500) break;
    }

    gStyle->SetOptStat(1111);
    hfcn_false->SetLineColor(kRed);
    hfcn_false->SetFillColor(kRed);
    hfcn_false->SetFillStyle(3004);
    hfcn_true->SetFillColor(kBlue);
    hfcn_true->SetFillStyle(3005);

    hfcn_true->Draw("HIST");
    hfcn_false->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(hfcn_true, "Toys: pol2", "F");
    leg->AddEntry(hfcn_false, "Toys: gaus", "F");
    leg->Draw();

    c1->Update();

    double p_value = static_cast<double>(nfcn) / Ntoys;
    double fcn_max = hfcn_true->GetBinCenter(j);
    double p_false = hfcn_false->Integral(1, j) / Ntoys;

    std::cout << "\n--- ToyMC Results ---" << std::endl;
    std::cout << "p = " << p_value
              << ", fcn_data = " << fcn
              << ", fcn_95% = " << fcn_max << std::endl;
    std::cout << "p_false (power proxy) = " << p_false << std::endl;
}