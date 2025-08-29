#include <iostream>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TMath.h>

void ex4_3() {
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(111);

    TFitResultPtr r;
    const int N = 1000;
    TH1D *hist = new TH1D("hist", "Gaussian Data;X;Y", 100, -2, 2);
    TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x + [2]*x*x", -2, 2);

    TRandom3 rand(0);

    for (int i = 0; i < N; ++i) {
        double x = rand.Gaus(0, 1);
        if (x >= -2 && x <= 2) {  
            hist->Fill(x);
        }
    }


    r = hist->Fit("fitFunc", "SQ"); 

    std::cout << "\n=== Результат подгонки ===" << std::endl;
    std::cout << "chi2 = " << r->Chi2() << ", ndf = " << r->Ndf()
              << ", p-value (chi2) = " << r->Prob() << std::endl;

    TH1D *hist_fit = new TH1D("hist_fit", "Fit Model (Pol2);X;Y", 100, -2, 2);

    TF1 *func_clone = (TF1*)fitFunc->Clone("func_clone");
    func_clone->SetRange(-2, 2);
    hist_fit->FillRandom("func_clone", N);

    // Критерий Колмогорова-Смирнова
    double ks_prob = hist->KolmogorovTest(hist_fit, "M");
    std::cout << "\n=== Критерий Колмогорова-Смирнова ===" << std::endl;
    std::cout << "K-S p-value = " << ks_prob << std::endl;

    // Критерий Андерсона-Дарлинга
    double ad_prob = hist->AndersonDarlingTest(hist_fit);
    std::cout << "\n=== Критерий Андерсона-Дарлинга ===" << std::endl;
    std::cout << "A-D test p-value = " << ad_prob << std::endl;

    TCanvas *c1 = new TCanvas("c1", "Fit and Goodness-of-Fit Tests", 900, 600);
    hist->SetLineColor(kBlack);
    hist->SetLineWidth(2);
    hist->Draw("E");

    fitFunc->SetLineColor(kBlue);
    fitFunc->SetLineWidth(2);
    fitFunc->Draw("same");


    TLegend *legend = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend->AddEntry(hist, "Data (Gaussian)", "lep");
    legend->AddEntry(fitFunc, "Fit: Pol2", "l");
    legend->Draw();

    TPaveText *pt = new TPaveText(0.15, 0.7, 0.45, 0.88, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->AddText(Form("K-S p = %5.3f", ks_prob));
    pt->AddText(Form("A-D p = %5.3f", ad_prob));
    pt->AddText(Form("#chi^{2}/ndf = %4.2f/%d", r->Chi2(), r->Ndf()));
    pt->AddText(Form("P(#chi^{2}) = %5.3f", r->Prob()));
    pt->Draw();

    c1->Update();
}