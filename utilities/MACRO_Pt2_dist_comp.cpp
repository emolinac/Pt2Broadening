// Code that outputs a plot comparing two Pt2 distributions
// Specifically, I used it to compare distributions with different corrections configurations

void MACRO_Pt2_dist_comp(Int_t Q2_bin, Int_t Nu_bin, Int_t Zh_bin, TString material)
{
  TFile* f1 = new TFile("corr_data_Pt2.root","READ");
  TFile* f2 = new TFile("./ACC/corr_data_Pt2.root","READ");

  TH1F* h1 = (TH1F*) f1->Get(Form("corr_data_Pt2_"+material+"_%i%i%i",Q2_bin,Nu_bin,Zh_bin));
  TH1F* h2 = (TH1F*) f2->Get(Form("corr_data_Pt2_"+material+"_%i%i%i",Q2_bin,Nu_bin,Zh_bin));

  h2->SetMarkerColor(kRed);
  h2->SetLineColor(kRed);

  Int_t min_bin = h1->GetMinimumBin();
  Int_t xlim = h1->GetBinCenter(min_bin+1);
  
  THStack* h = new THStack("h","Pt2 comp; Pt2; counts");
  h->Add(h1);
  h->Add(h2);

  h->Draw("NOSTACK");
  h->GetXaxis()->SetRangeUser(0,xlim);

  TLegend* legend = new TLegend(0.7,0.8,0.9,0.92,"","NDC");
  legend->AddEntry(h1,"ACC+RC","lpf");
  legend->AddEntry(h2,"ACC","lpf");
  legend->SetFillStyle(0);                                                                                                                                    
  legend->SetBorderSize(0);                                                                                                                                   
  legend->SetTextFont(62);                                                                                                                                    
  legend->SetTextSize(0.04);

  h->Draw("NOSTACKSAME");
  legend->Draw("SAME");
  gPad->SetLogy(1);
}
