Int_t Q2_bin_selected, Nu_bin_selected;

Int_t N_Q2 = 3;
Int_t N_Nu = 3;
const Int_t N_Zh = 8;
Int_t N_Pt2 = 90;
Int_t N_Phi = 12;

Double_t Zh_min = 0.;
Double_t Zh_max = 1.;

Float_t Zh_limits[N_Zh+1] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1};

void set_error_X_null(TGraphErrors* g);
void set_X_shift(TGraphErrors* g, Double_t shift);

void MACRO_Ptbroadening_Zh_Integrated()
{
  TCanvas* c = new TCanvas("c","",800,600);
  c->Draw();

  TFile* fsource  = new TFile(Form("meanPt2_Zh_%i.root",N_Zh),"READ");

  TH1F* h_sol[3];  TH1F* h_liq[3];  TH1F* h_broad[3];
  
  h_sol[0] = (TH1F*) fsource->Get("meanPt2_C_CLEAN_INTERPOLATED");
  h_sol[1] = (TH1F*) fsource->Get("meanPt2_Fe_CLEAN_INTERPOLATED");
  h_sol[2] = (TH1F*) fsource->Get("meanPt2_Pb_CLEAN_INTERPOLATED");

  h_liq[0] = (TH1F*) fsource->Get("meanPt2_DC_CLEAN_INTERPOLATED");
  h_liq[1] = (TH1F*) fsource->Get("meanPt2_DFe_CLEAN_INTERPOLATED");
  h_liq[2] = (TH1F*) fsource->Get("meanPt2_DPb_CLEAN_INTERPOLATED");

  for(Int_t i = 0 ; i < 3 ; i++){
    h_broad[i] = new TH1F(Form("h_broad_%i",i),"",N_Zh,Zh_limits);
    h_broad[i]->Add(h_sol[i],h_liq[i],1,-1);
  }

  TGraphErrors* g[3];
  for(Int_t i = 0 ; i < 3 ; i++){
    g[i] = (TGraphErrors*) TH1TOTGraph(h_broad[i]);
    set_error_X_null(g[i]);
  }

  g[0]->SetMarkerColor(kRed);
  g[0]->SetLineColor(kRed);
  g[1]->SetMarkerColor(kBlue);
  g[1]->SetLineColor(kBlue);
  g[2]->SetMarkerColor(kBlack);
  g[2]->SetLineColor(kBlack);

  set_X_shift(g[0],-0.015);
  set_X_shift(g[2],0.015);
  
  TLegend* legend = new TLegend(0.15,0.75,0.3,0.92,"","NDC");
  legend->AddEntry(g[0],"C","lpf");
  legend->AddEntry(g[1],"Fe","lpf");
  legend->AddEntry(g[2],"Pb","lpf");
  legend->SetFillStyle(0);                                                                                                                                    
  legend->SetBorderSize(0);                                                                                                                                   
  legend->SetTextFont(62);                                                                                                                                    
  legend->SetTextSize(0.04);

  TLatex* kinematics2 = new TLatex();
  kinematics2->SetTextFont(52);
  kinematics2->SetTextSize(0.04);

  TPad* p = new TPad("p","p",0,0,1,1);
  p->SetLeftMargin(0.13);
  p->SetTopMargin(0.06);
  p->SetGridx(1);
  p->SetGridy(1);
  p->Draw();
  p->cd();

  TMultiGraph* mg = new TMultiGraph();
  mg->Add(g[0]);
  mg->Add(g[1]);
  mg->Add(g[2]);

  mg->Draw("APE0");

  gStyle->SetTitleFont(62,"XY");
  gStyle->SetTitleSize(0.04,"XY");
  
  mg->GetYaxis()->SetRangeUser(0.,0.059);
  mg->GetXaxis()->SetRangeUser(Zh_min,Zh_max);
  mg->GetXaxis()->SetTitle("z_{h}");
  mg->GetXaxis()->CenterTitle();
  mg->GetXaxis()->SetTitleOffset(1.1);

  mg->GetYaxis()->SetTitle("#Delta P^{2}_{T} [GeV^{2}]");
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetTitleOffset(1.3);

  mg->Draw("APE0");

  legend->Draw();
  //  kinematics2->DrawLatexNDC(0.3,0.9,"CLAS PRELIMINARY");

  // c->Print(Form("./PLOTS/Pt_broad_Zh.pdf",Q2_bin_selected,Nu_bin_selected));
  // c->Print(Form("./PLOTS/Pt_broad_Zh.png",Q2_bin_selected,Nu_bin_selected));

  return c;
}

void set_error_X_null(TGraphErrors* g)
{
  Double_t* errors_Y = g->GetEY();
  
  for(Int_t point = 0 ; point < N_Zh ; point++)
    {
      g->SetPointError(point,0,errors_Y[point]);  
    }
}

void set_X_shift(TGraphErrors* g, Double_t shift)
{
  Double_t* content_Y = g->GetY();
  Double_t* content_X = g->GetX();

  for(Int_t point = 0; point < N_Zh ; point++)
    {
      g->SetPoint(point, content_X[point] + shift , content_Y[point]);
    }
}

