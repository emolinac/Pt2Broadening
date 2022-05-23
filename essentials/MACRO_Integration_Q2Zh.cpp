const Int_t N_Q2 = 3;
const Int_t N_Nu = 3;
const Int_t N_Zh = 8;
const Int_t N_Pt2 = 90;
const Int_t N_Phi = 12;

const Double_t Q2_min = 1.;
const Double_t Q2_max = 4.;
const Double_t Zh_min = 0.;
const Double_t Zh_max = 1;
const Double_t Pt2_min = 0.;
const Double_t Pt2_max = 3.;

Float_t Nu_limits[N_Nu+1] = {2.2,3.2,3.7,4.26};

Int_t empty_histo(TH1F* h);
void integration(TFile* ftarget, TFile* fsource, TString material, Float_t* Nu_limits);

void MACRO_Integration_Q2Zh()
{
  
  TFile* fsource  = new TFile("corr_data_Pt2_processed.root","READ");
  TFile* ftarget  = new TFile(Form("meanPt2_Nu_%i.root",N_Nu),"RECREATE");

  integration(ftarget, fsource, "DC", Nu_limits);
  integration(ftarget, fsource, "DFe", Nu_limits);
  integration(ftarget, fsource, "DPb", Nu_limits);
  
  integration(ftarget, fsource, "C", Nu_limits);
  integration(ftarget, fsource, "Fe", Nu_limits);
  integration(ftarget, fsource, "Pb", Nu_limits);
}

Int_t empty_histo(TH1F* h)
{
  Int_t empty = 0;
  for(Int_t i = 1 ; i <= h->GetNbinsX() ; i++)
    {
      if(h->GetBinContent(i)==0){
	empty++;}
    }
  if(empty==h->GetNbinsX())
    {return 1;}
  else
    {return 0;}
}

void integration(TFile* ftarget, TFile* fsource, TString material, Float_t* Nu_limits)
{
  std::cout<<"Material "<<material<<std::endl;
  
  TH1F* h_Nu = new TH1F("corr_data_"+material+"_Nu","",N_Nu,Nu_limits);
  for(Int_t Nu_bin = 0 ; Nu_bin < N_Nu ; Nu_bin++)
    {
      TH1F* h_Pt2_Integrated = new TH1F(Form("corr_data_Pt2_"+material+"_%i",Nu_bin),"",N_Pt2,Pt2_min,Pt2_max);
      for(Int_t Zh_bin = 2 ; Zh_bin < N_Zh ; Zh_bin++)
	{
	  for(Int_t Q2_bin = 0 ; Q2_bin < N_Q2 ; Q2_bin++)
	    {
	      TH1F* h_Pt2 =  (TH1F*) fsource->Get(Form("corr_data_Pt2_"+material+"_%i%i%i_CLEAN_INTERPOLATED",Q2_bin,Nu_bin,Zh_bin));

	      h_Pt2_Integrated->Add(h_Pt2);

	      delete h_Pt2;
	    }// End Q2 loop
	}// End Zh loop

      h_Nu->SetBinContent(Nu_bin+1, h_Pt2_Integrated->GetMean());
      h_Nu->SetBinError(Nu_bin+1, h_Pt2_Integrated->GetMeanError());

      delete h_Pt2_Integrated;
    }// End Nu loop

  ftarget->cd();
  h_Nu->Write("meanPt2_"+material+"_Nu_CLEAN_INTERPOLATED");

  delete h_Nu;

}
