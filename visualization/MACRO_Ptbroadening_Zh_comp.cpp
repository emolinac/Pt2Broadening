const Int_t N_comparison = 2;

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

void MACRO_Ptbroadening_Zh_comp(Int_t Q2_bin_selected_input, Int_t Nu_bin_selected_input)
{
  Q2_bin_selected = Q2_bin_selected_input;
  Nu_bin_selected = Nu_bin_selected_input;
  
  TCanvas* c = new TCanvas("c","",800,600);
  c->Draw();

  TFile* fsource[N_comparison];
  fsource[0]  = new TFile("corr_data_Pt2_processed.root","READ");
  fsource[1]  = new TFile("./ACC/corr_data_Pt2_processed.root","READ");
  
  TFile* fbinning  = new TFile(Form("binning_%i%i%i%i%i_Zh_fancy.root",N_Q2,N_Nu,N_Zh,N_Pt2,N_Phi),"READ");

  Double_t Q2_min, Q2_max, Nu_min, Nu_max;
  set_kinematic_limits(fbinning, Q2_min, Q2_max, Nu_min, Nu_max);

  TH1F* h_sol[N_comparison][3];   TH1F* h_liq[N_comparison][3];  TH1F* h_broad[N_comparison][3];
  TGraphErrors* g[N_comparison][3];

  for(Int_t comp_index = 0 ; comp_index < N_comparison ; comp_index++){
    h_sol[comp_index][0] = (TH1F*) fsource[comp_index]->Get(Form("meanPt2_C_%i%i_CLEAN_INTERPOLATED",Q2_bin_selected,Nu_bin_selected));
    h_sol[comp_index][1] = (TH1F*) fsource[comp_index]->Get(Form("meanPt2_Fe_%i%i_CLEAN_INTERPOLATED",Q2_bin_selected,Nu_bin_selected));
    h_sol[comp_index][2] = (TH1F*) fsource[comp_index]->Get(Form("meanPt2_Pb_%i%i_CLEAN_INTERPOLATED",Q2_bin_selected,Nu_bin_selected));

    h_liq[comp_index][0] = (TH1F*) fsource[comp_index]->Get(Form("meanPt2_DC_%i%i_CLEAN_INTERPOLATED",Q2_bin_selected,Nu_bin_selected));
    h_liq[comp_index][1] = (TH1F*) fsource[comp_index]->Get(Form("meanPt2_DFe_%i%i_CLEAN_INTERPOLATED",Q2_bin_selected,Nu_bin_selected));
    h_liq[comp_index][2] = (TH1F*) fsource[comp_index]->Get(Form("meanPt2_DPb_%i%i_CLEAN_INTERPOLATED",Q2_bin_selected,Nu_bin_selected));

    for(Int_t i = 0 ; i < 3 ; i++){
      h_broad[comp_index][i] = new TH1F(Form("h_broad_%i_%i",i,comp_index),"",N_Zh,Zh_limits);
      h_broad[comp_index][i]->Add(h_sol[comp_index][i],h_liq[comp_index][i],1,-1);
    }

    for(Int_t i = 0 ; i < 3 ; i++){
      g[comp_index][i] = (TGraphErrors*) TH1TOTGraph(h_broad[comp_index][i]);
      set_error_X_null(g[comp_index][i]);
    }
    g[comp_index][0]->SetMarkerColor(kRed);
    g[comp_index][0]->SetLineColor(kRed);
    g[comp_index][1]->SetMarkerColor(kBlue);
    g[comp_index][1]->SetLineColor(kBlue);
    g[comp_index][2]->SetMarkerColor(kBlack);
    g[comp_index][2]->SetLineColor(kBlack);

    g[comp_index][0]->SetMarkerStyle(20+comp_index*4);
    g[comp_index][1]->SetMarkerStyle(20+comp_index*4);
    g[comp_index][2]->SetMarkerStyle(20+comp_index*4);
   
    
    set_X_shift(g[comp_index][0],-0.015);
    set_X_shift(g[comp_index][2],0.015);

  }
  
  TLegend* legend = new TLegend(0.15,0.75-0.17*(N_comparison-1),0.3,0.92,"","NDC");
  legend->AddEntry(g[0][0],"C ACC+RC","lpf");
  legend->AddEntry(g[0][1],"Fe ACC+RC","lpf");
  legend->AddEntry(g[0][2],"Pb ACC+RC","lpf");
  legend->AddEntry(g[1][0],"C ACC","lpf");
  legend->AddEntry(g[1][1],"Fe ACC","lpf");
  legend->AddEntry(g[1][2],"Pb ACC","lpf");

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

  TLine* line = new TLine();
  
  TPad* p = new TPad("p","p",0,0,1,1);
  p->SetLeftMargin(0.13);
  p->SetTopMargin(0.06);
  p->SetGridy(1);
  p->Draw();
  p->cd();

  TMultiGraph* mg = new TMultiGraph();
  for(Int_t comp_index = 0 ; comp_index < N_comparison ; comp_index++){
    mg->Add(g[comp_index][0]);
    mg->Add(g[comp_index][1]);
    mg->Add(g[comp_index][2]);
  }
  
  mg->Draw("APE0");

  gStyle->SetTitleFont(62,"XY");
  gStyle->SetTitleSize(0.04,"XY");
  
  mg->GetYaxis()->SetRangeUser(0.,0.085);
  mg->GetXaxis()->SetRangeUser(Zh_min,Zh_max);
  mg->GetXaxis()->SetTitle("z_{h}");
  mg->GetXaxis()->SetTitleOffset(1.1);

  mg->GetYaxis()->SetTitle("#Delta P^{2}_{T} [GeV^{2}]");
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetTitleOffset(1.3);

  mg->SetTitle(";z_{h};#Delta P^{2}_{T} [GeV^{2}]");
  mg->Draw("APE0");

  legend->Draw();
  kinematics->DrawLatexNDC(1.5*p->GetLeftMargin(),1.02*(1-p->GetTopMargin()),Form("%.2f<Q^{2}[GeV^{2}]<%.2f     %.2f<#nu[GeV]<%.2f",Q2_min,Q2_max,Nu_min,Nu_max));
  kinematics2->DrawLatexNDC(0.3,0.9,"CLAS PRELIMINARY");
  line->DrawLine(0.1,0,1,0);

  c->Print(Form("./PLOTS/Pt_broad_Zh_%i%i_ACCRCCOMP.pdf",Q2_bin_selected,Nu_bin_selected));
  c->Print(Form("./PLOTS/Pt_broad_Zh_%i%i_ACCRCCOMP.png",Q2_bin_selected,Nu_bin_selected));

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
