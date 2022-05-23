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

Float_t Zh_limits[N_Zh+1] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1};

Int_t empty_histo(TH1F* h);
void integration(TFile* ftarget, TFile* fsource, TString material);

void MACRO_Integration_Q2NuZh()
{
  TFile* fsource = new TFile("corr_data_Pt2_processed.root","READ");
  TFile* ftarget = new TFile("meanPt2.root","RECREATE");


  integration(ftarget, fsource, "DC");
  integration(ftarget, fsource, "DFe");
  integration(ftarget, fsource, "DPb");

  integration(ftarget, fsource, "C");
  integration(ftarget, fsource, "Fe");
  integration(ftarget, fsource, "Pb");

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

void integration(TFile* ftarget, TFile* fsource, TString material)
{
  TH1F* h         = new TH1F("corr_data_"+material,"",1,0,1);
  TH1F* h_Pt2_INT = new TH1F("corr_data_"+material+"_Pt2","",N_Pt2,Pt2_min,Pt2_max);
  
  for(Int_t Zh_bin = 1 ; Zh_bin < N_Zh ; Zh_bin++)
    {	
      for(Int_t Nu_bin = 0 ; Nu_bin < N_Nu ; Nu_bin++)
	{
	  for(Int_t Q2_bin = 0 ; Q2_bin < N_Q2 ; Q2_bin++)
	    {
	      TH1F* h_Pt2 = (TH1F*) fsource->Get(Form("corr_data_Pt2_"+material+"_%i%i%i_CLEAN_INTERPOLATED",Q2_bin,Nu_bin,Zh_bin));	  

	      h_Pt2_INT->Add(h_Pt2);
	      
	      delete h_Pt2;
	    }// End Q2 loop
	}// End Nu loop
    }// End Zh loop

  h->SetBinContent(1, h_Pt2_INT->GetMean());
  h->SetBinError(1, h_Pt2_INT->GetMeanError());  
  
  ftarget->cd();
  h->Write("meanPt2_"+material+"_CLEAN_INTERPOLATED");
  
  delete h;
  delete h_Pt2_INT;
}
