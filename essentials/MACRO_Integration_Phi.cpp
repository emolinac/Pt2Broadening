Int_t N_Q2 = 3;
Int_t N_Nu = 3;
Int_t N_Zh = 8;
Int_t N_Pt2 = 90;
Int_t N_Phi = 12;

Double_t Pt2_min = 0.;
Double_t Pt2_max = 3.;

Double_t delta_Pt2 = (Pt2_max - Pt2_min)/N_Pt2;

Int_t empty_histo(TH1F* h);
void Phi_integration(TString material, TFile* ftarget, TFile* fsource);

void MACRO_Integration_Phi()
{
  TFile* fsource = new TFile("corr_data_Phi.root","READ");

  TFile* ftarget = new TFile("corr_data_Pt2.root","RECREATE");

  Phi_integration("DC", ftarget, fsource);
  Phi_integration("DFe", ftarget, fsource);
  Phi_integration("DPb", ftarget, fsource);
  
  Phi_integration("C", ftarget, fsource);
  Phi_integration("Fe", ftarget, fsource);
  Phi_integration("Pb", ftarget, fsource);
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

void Phi_integration(TString material, TFile* ftarget, TFile* fsource)
{
  for(Int_t Q2_bin = 0 ; Q2_bin < N_Q2 ; Q2_bin++)
  {
    for(Int_t Nu_bin = 0 ; Nu_bin < N_Nu ; Nu_bin++)
    {
      for(Int_t Zh_bin = 0 ; Zh_bin < N_Zh ; Zh_bin++)
      {
	
	TH1F* h_Pt2 = new TH1F(Form("corr_data_Pt2_"+material+"_%i%i%i",Q2_bin,Nu_bin,Zh_bin),"",N_Pt2,Pt2_min,Pt2_max);
	
        for(Int_t Pt2_bin = 0 ; Pt2_bin < N_Pt2 ; Pt2_bin++)
        {
	  TH1F* h_Phi = (TH1F*) fsource->Get(Form("corr_data_"+material+"_%i%i%i%i",Q2_bin,Nu_bin,Zh_bin,Pt2_bin));
	  
	  if(h_Phi==NULL){continue;}
	  if(empty_histo(h_Phi)==1){continue;}
	  Double_t errors_array[1] = {};
	  Double_t integral = h_Phi->IntegralAndError(1,12,*errors_array/*,"width"*/);
	  //Integral
	  h_Pt2->SetBinContent(Pt2_bin+1, integral);
	  h_Pt2->SetBinError(Pt2_bin+1, errors_array[0]);
	  
	}// End Pt2 loop

	if(empty_histo(h_Pt2)==0)
	  {
	    ftarget->cd();
	    h_Pt2->Write();
	  }
      }// End Zh loop
    }// End Nu loop
  }// End Q2 loop
  
}
