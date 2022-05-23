const Int_t N_Q2 = 3;
const Int_t N_Nu = 3;
const Int_t N_Zh = 8;
const Int_t N_Pt2 = 90;
const Int_t N_Phi = 12;

Float_t Nu_limits[N_Nu+1] = {2.2,3.2,3.7,4.26};
Float_t Nu_centroids[3][3] = {{2.87,3.46,3.96},
			      {2.87,3.46,3.96},
			      {2.87,3.46,3.97}};

void set_error_X_null(TGraphErrors* g);
void set_X_shift_and_centroid(TGraphErrors* g, Double_t shift, Int_t TargIndex);

void MACRO_Ptbroadening_Nu_Integrated_BinCentroid()
{
  TCanvas* c = new TCanvas("c","",800,600);
  c->Draw();

  TFile* fsource  = new TFile(Form("meanPt2_Nu_%i.root",N_Nu),"READ");

  TH1F* h_meanPt2_sol[3];   TH1F* h_meanPt2_liq[3];

  //C
  h_meanPt2_sol[0] = (TH1F*) fsource->Get("meanPt2_C_Nu_CLEAN_INTERPOLATED");
  h_meanPt2_liq[0] = (TH1F*) fsource->Get("meanPt2_DC_Nu_CLEAN_INTERPOLATED");
  //Fe
  h_meanPt2_sol[1] = (TH1F*) fsource->Get("meanPt2_Fe_Nu_CLEAN_INTERPOLATED");
  h_meanPt2_liq[1] = (TH1F*) fsource->Get("meanPt2_DFe_Nu_CLEAN_INTERPOLATED");
  //Pb
  h_meanPt2_sol[2] = (TH1F*) fsource->Get("meanPt2_Pb_Nu_CLEAN_INTERPOLATED");
  h_meanPt2_liq[2] = (TH1F*) fsource->Get("meanPt2_DPb_Nu_CLEAN_INTERPOLATED");
  
  TH1F* h_broad[3];  
  for(Int_t i = 0 ; i < 3 ; i++){
    h_broad[i] = new TH1F(Form("h_broad_%i",i),"",N_Nu,Nu_limits);
    h_broad[i]->Add(h_meanPt2_sol[i],h_meanPt2_liq[i],1,-1);
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

  set_X_shift_and_centroid(g[0],0,0);
  set_X_shift_and_centroid(g[1],0,1);
  set_X_shift_and_centroid(g[2],0,2);

  //FIT TEST
  std::cout<<"C FIT"<<std::endl;
  g[0]->Fit("pol1","N");
  std::cout<<"Fe FIT"<<std::endl;
  g[1]->Fit("pol1","N");
  std::cout<<"Pb FIT"<<std::endl;
  g[2]->Fit("pol1","N");

  
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

  TLine* line = new TLine();
  
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
  mg->GetXaxis()->SetRangeUser(2.1,4.3);
  mg->GetXaxis()->SetTitle("#nu[GeV]");
  mg->GetXaxis()->CenterTitle();
  mg->GetXaxis()->SetTitleOffset(1.1);

  mg->GetYaxis()->SetTitle("#Delta P^{2}_{T} [GeV^{2}]");
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetTitleOffset(1.3);

  //  mg->SetTitle(";z_{h};#Delta P^{2}_{T} [GeV^{2}]");
  mg->Draw("APE0");

  legend->Draw();
  //  kinematics2->DrawLatexNDC(0.3,0.9,"CLAS PRELIMINARY");
  //  line->DrawLine(0.1,0,1,0);

  //c->Print("./PLOTS/Pt_broad_Nu.pdf");
  //c->Print("./PLOTS/Pt_broad_Nu.png");

  return c;
}

void set_error_X_null(TGraphErrors* g)
{
  Double_t* errors_Y = g->GetEY();
  
  for(Int_t point = 0 ; point < N_Nu ; point++)
    {
      g->SetPointError(point,0,errors_Y[point]);  
    }
}

void set_X_shift_and_centroid(TGraphErrors* g, Double_t shift, Int_t TargIndex)
{
  Double_t* content_Y = g->GetY();
  Double_t* content_X = g->GetX();

  for(Int_t point = 0; point < N_Nu ; point++)
    {
      g->SetPoint(point, Nu_centroids[TargIndex][point] + shift , content_Y[point]);
    }
}
