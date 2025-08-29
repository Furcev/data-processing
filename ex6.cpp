#include <iostream>
#include <TMath.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TStyle.h>

void ex6() {

    gStyle->SetOptStat(111111);

    double a = 2.0;
    double b = 1.5;
    double c = 1.0;

    double V_box = 8.0 * a * b * c;


    const long N = 1000000;


    TRandom3 rand(0);

    double sum_f = 0.0;
    double sum_f2 = 0.0;


    for (long i = 0; i < N; ++i) {

        double x = rand.Uniform(-a, a);
        double y = rand.Uniform(-b, b);
        double z = rand.Uniform(-c, c);

        double value = (x*x)/(a*a) + (y*y)/(b*b) + (z*z)/(c*c);
        double f = (value <= 1.0) ? 1.0 : 0.0;


        sum_f += f;
        sum_f2 += f * f;
    }


    double f_avg = sum_f / N;
    double f2_avg = sum_f2 / N;


    double var_f = f2_avg - f_avg * f_avg;

    double V = V_box * f_avg;


    double sigma_V = V_box * TMath::Sqrt(var_f / N);


    double V_true = (4.0/3.0) * TMath::Pi() * a * b * c;


    double rel_error = TMath::Abs(V - V_true) / V_true;


    std::cout << "\n=== Метод Монте-Карло: объём эллипсоида ===" << std::endl;
    std::cout << "Полуоси: a = " << a << ", b = " << b << ", c = " << c << std::endl;
    std::cout << "Число точек N = " << N << std::endl;
    std::cout << "Объём (Монте-Карло) = " << V << " ± " << sigma_V << std::endl;
    std::cout << "Аналитический объём = " << V_true << std::endl;
    std::cout << "Разница = " << V - V_true << std::endl;
    std::cout << "Относительная ошибка = " << rel_error * 100.0 << " %" << std::endl;
    std::cout << "Статистическая ошибка (sigma_V) = " << sigma_V << std::endl;


    TH1D *h_f = new TH1D("h_f", "Indicator Function f_i;Inside (1) / Outside (0);Count", 2, -0.5, 1.5);
    h_f->SetBinContent(1, N * (1 - f_avg));
    h_f->SetBinContent(2, N * f_avg);
    h_f->SetFillColor(kBlue);
    h_f->SetFillStyle(3001);


    TCanvas *c1 = new TCanvas("c1", "Monte Carlo: Ellipsoid Volume", 800, 600);
    h_f->Draw();

    TPaveText *pt = new TPaveText(0.5, 0.7, 0.88, 0.88, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->AddText(Form("V_{MC} = %.4f #pm %.4f", V, sigma_V));
    pt->AddText(Form("V_{true} = %.4f", V_true));
    pt->AddText(Form("#epsilon_{rel} = %.3f%%", rel_error * 100));
    pt->Draw();

    c1->Update();

}