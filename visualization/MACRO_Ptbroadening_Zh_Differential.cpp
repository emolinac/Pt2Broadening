Int_t Q2_bin_selected, Nu_bin_selected;

Int_t N_Q2 = 3;
Int_t N_Nu = 3;
const Int_t N_Zh = 8;
Int_t N_Pt2 = 90;
Int_t N_Phi = 12;

Double_t Zh_min = 0.;
Double_t Zh_max = 1.;

const Double_t Zh_limits[N_Zh+1] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1};

void set_error_X_null(TGraphErrors* g);
void set_X_shift(TGraphErrors* g, Double_t shift);
void set_kinematic_limits(TFile* fbinning, Double_t& Q2_min, Double_t& Q2_max, Double_t& Nu_min, Double_t& Nu_max);

void MACRO_Ptbroadening_Zh_Differential(Int_t Q2_bin_selected_input, Int_t Nu_bin_selected_input)
{
  Q2_bin_selected = Q2_bin_selected_input;
  Nu_bin_selected = Nu_bin_selected_input;
  
  TCanvas* c = new TCanvas("c","",800,600);
  c->Draw();

  TFile* fbinning = new TFile(Form("binning_%i%i%i%i%i_Zh_fancy.root",N_Q2,N_Nu,N_Zh,N_Pt2,N_Phi),"READ");
  TFile* fsource  = new TFile("corr_data_Pt2_processed.root","READ");

  Double_t Q2_min, Q2_max, Nu_min, Nu_max;
  set_kinematic_limits(fbinning, Q2_min, Q2_max, Nu_min, Nu_max);

  TH1F* h_sol[3];   TH1F* h_liq[3];  TH1F* h_broad[3];
  
  h_sol[0] = (TH1F*) fsource->Get(Form("meanPt2_C_%i%i_CLEAN_INTERPOLATED",Q2_bin_selected,Nu_bin_selected));
  h_sol[1] = (TH1F*) fsource->Get(Form("meanPt2_Fe_%i%i_CLEAN_INTERPOLATED",Q2_bin_selected,Nu_bin_selected));
  h_sol[2] = (TH1F*) fsource->Get(Form("meanPt2_Pb_%i%i_CLEAN_INTERPOLATED",Q2_bin_selected,Nu_bin_selected));

  h_liq[0] = (TH1F*) fsource->Get(Form("meanPt2_DC_%i%i_CLEAN_INTERPOLATED",Q2_bin_selected,Nu_bin_selected));
  h_liq[1] = (TH1F*) fsource->Get(Form("meanPt2_DFe_%i%i_CLEAN_INTERPOLATED",Q2_bin_selected,Nu_bin_selected));
  h_liq[2] = (TH1F*) fsource->Get(Form("meanPt2_DPb_%i%i_CLEAN_INTERPOLATED",Q2_bin_selected,Nu_bin_selected));

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


  TLatex* kinematics = new TLatex();
  kinematics->SetTextFont(62);
  kinematics->SetTextSize(0.04);

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
  
  mg->GetYaxis()->SetRangeUser(0.,0.085);
  mg->GetXaxis()->SetRangeUser(Zh_min,Zh_max);
  mg->GetXaxis()->SetTitle("z_{h}");
  mg->GetXaxis()->CenterTitle();
  mg->GetXaxis()->SetTitleOffset(1.1);

  mg->GetYaxis()->SetTitle("#Delta P^{2}_{T} [GeV^{2}]");
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetTitleOffset(1.3);

  //  mg->SetTitle(";z_{h};#Delta P^{2}_{T} [GeV^{2}]");
  mg->Draw("APE0");

  legend->Draw();
  kinematics->DrawLatexNDC(2.4*p->GetLeftMargin(),1.024*(1-p->GetTopMargin()),Form("%.2f<Q^{2}[GeV^{2}]<%.2f     %.2f<#nu[GeV]<%.2f",Q2_min,Q2_max,Nu_min,Nu_max));
  //  kinematics2->DrawLatexNDC(0.3,0.9,"CLAS PRELIMINARY");

  // c->Print(Form("./PLOTS/Pt_broad_Zh_%i%i.pdf",Q2_bin_selected,Nu_bin_selected));
  // c->Print(Form("./PLOTS/Pt_broad_Zh_%i%i.png",Q2_bin_selected,Nu_bin_selected));

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

void set_kinematic_limits(TFile* fbinning, Double_t& Q2_min, Double_t& Q2_max, Double_t& Nu_min, Double_t& Nu_max)
{
  TNtuple* limits_tuple = (TNtuple*) fbinning->Get("limits_tuple");
  Float_t Q2_min_local, Q2_max_local, Nu_min_local, Nu_max_local, Q2_bin_local, Nu_bin_local;
  limits_tuple->SetBranchAddress("Q2_min",&Q2_min_local);
  limits_tuple->SetBranchAddress("Q2_max",&Q2_max_local);
  limits_tuple->SetBranchAddress("Nu_min",&Nu_min_local);
  limits_tuple->SetBranchAddress("Nu_max",&Nu_max_local);
  limits_tuple->SetBranchAddress("Q2_bin",&Q2_bin_local);
  limits_tuple->SetBranchAddress("Nu_bin",&Nu_bin_local);
  
  for(Int_t tuple_entry = 0 ; tuple_entry < limits_tuple->GetEntries() ; tuple_entry++)
  {
    limits_tuple->GetEntry(tuple_entry);
    if(Q2_bin_selected==Q2_bin_local&&Nu_bin_selected==Nu_bin_local)
    {
      Q2_min = Q2_min_local;
      Q2_max = Q2_max_local;
      Nu_min = Nu_min_local;
      Nu_max = Nu_max_local;

      break;
    }
  }
}
