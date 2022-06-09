// Name convention:
// X                    = Distribution without treatment
// X_CLEAN              = Distribution with Pt2 cutoff applied
// X_CLEAN_INTERPOLATED = Distribution with Pt2 cutoff applied and interpolation


// Constants

const Int_t N_Q2  = 3;
const Int_t N_Nu  = 3;
const Int_t N_Zh  = 8;
const Int_t N_Pt2 = 90;
const Int_t N_Phi = 12;

const Double_t Zh_min            = 0;
const Double_t Zh_max            = 1;
const Double_t Zh_limits[N_Zh+1] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1};

const Double_t Pt2_min   = 0.;
const Double_t Pt2_max   = 3.;
const Double_t delta_Pt2 = (Pt2_max-Pt2_min)/N_Pt2; 

const Double_t ChiSQndf_weight = 0.95;
const Double_t ndf_weight      = 0.05;

const Int_t number_of_fits  = 15;
const Double_t max_ChiSQndf = 10;  

// Functions declaration

void     fit_and_cutoff(TString material, TFile* fsource, TFile* ftarget);
void     GetCleanPt2Distribution(TH1F* h, Double_t cutoff);
void     GetInterpolatedPt2Distribution(TH1F* h);

Int_t    GetFirstEmptyBin(TH1F* h); 

Double_t GetWeightedPt2Cutoff(Double_t* ChiSQndf_array,Double_t* ndf_array,Double_t* Pt2_cutoff, Int_t array_size);
Double_t GetInterpolatedPt2Point(Double_t x, Double_t prev_bin, Double_t next_bin, Double_t prev_bin_center, Double_t next_bin_center);
Double_t GetInterpolatedPt2PointError(Double_t x, Double_t prev_bin_error, Double_t next_bin_error, Double_t prev_bin_center, Double_t next_bin_center);

// MAIN
void MACRO_Pt2_processing(){
  TFile* fsource = new TFile("corr_data_Pt2.root","READ");
  TFile* ftarget = new TFile("corr_data_Pt2_processed.root","RECREATE");

  fit_and_cutoff("DC",fsource,ftarget);
  fit_and_cutoff("DFe",fsource,ftarget);
  fit_and_cutoff("DPb",fsource,ftarget);

  fit_and_cutoff("C",fsource,ftarget);
  fit_and_cutoff("Fe",fsource,ftarget);
  fit_and_cutoff("Pb",fsource,ftarget);
}

void fit_and_cutoff(TString material, TFile* fsource, TFile* ftarget)
{
  for(Int_t Q2_bin = 0 ; Q2_bin<N_Q2 ; Q2_bin++){
    for(Int_t Nu_bin = 0 ; Nu_bin<N_Nu ; Nu_bin++){

      TH1F* meanPt2                    = new TH1F("meanPt2","",N_Zh,Zh_limits);
      TH1F* meanPt2_CLEAN              = new TH1F("meanPt2_CLEAN","",N_Zh,Zh_limits);
      TH1F* meanPt2_CLEAN_INTERPOLATED = new TH1F("meanPt2_CLEAN_INTERPOLATED","",N_Zh,Zh_limits);
      
      for(Int_t Zh_bin = 0 ; Zh_bin<N_Zh ; Zh_bin++){
	std::cout<<"_________________________Zh_bin="<<Zh_bin<<"__________________________"<<std::endl;
	TH1F* h        = (TH1F*)fsource->Get(Form("corr_data_Pt2_"+material+"_%i%i%i",Q2_bin,Nu_bin,Zh_bin));
	TH1F* h_clone  = (TH1F*)h->Clone();

	// The tail treatment in zh<0.3 is futile since Pt2 distributions fall quickly
	if(Zh_bin<4){
	  ftarget->cd();
	  h->Write();
	  h_clone->Write(Form("corr_data_Pt2_"+material+"_%i%i%i_CLEAN",Q2_bin,Nu_bin,Zh_bin));
	  GetInterpolatedPt2Distribution(h_clone);
	  h_clone->Write(Form("corr_data_Pt2_"+material+"_%i%i%i_CLEAN_INTERPOLATED",Q2_bin,Nu_bin,Zh_bin));

	  meanPt2->SetBinContent(Zh_bin+1, h->GetMean());
	  meanPt2->SetBinError(Zh_bin+1, h->GetMeanError());

	  meanPt2_CLEAN->SetBinContent(Zh_bin+1, h->GetMean());
	  meanPt2_CLEAN->SetBinError(Zh_bin+1, h->GetMeanError());

	  meanPt2_CLEAN_INTERPOLATED->SetBinContent(Zh_bin+1, h_clone->GetMean());
	  meanPt2_CLEAN_INTERPOLATED->SetBinError(Zh_bin+1, h_clone->GetMeanError());

	  continue;
	}

	Double_t ChiSQndf_array[number_of_fits];
	Double_t ndf_array[number_of_fits];
	Double_t Pt2_cutoff[number_of_fits];
	
	//START FIT LOOP
	Int_t fits_completed = 0;
	Double_t Pt2_fit_min, Pt2_fit_max;

	for(Int_t i = 0 ; i < number_of_fits ; i++){
	  Pt2_fit_min = i*delta_Pt2;
	  Pt2_fit_max = 3.;

	  TF1* fit_func = new TF1("fit_func","[0]*TMath::Exp(-x/[1])",Pt2_fit_min,Pt2_fit_max);
	  fit_func->SetParameter(0,h_clone->GetBinContent(1));
	  fit_func->SetParameter(1,0.2);

	  h_clone->Fit(fit_func,"SRQ");

	  //REGION VALIDITY TEST FUNCTION
	  TF1* validity_func = new TF1("validity_func","[0]*TMath::Exp(-x/[1])",Pt2_fit_min,1000);
	  validity_func->SetParameter(0,fit_func->GetParameter(0));
	  validity_func->SetParameter(1,fit_func->GetParameter(1));

	  // Check if ChiSQndf is lower than the maximum allowed
	  // Check if width of exponential is reasonable
	  if((fit_func->GetChisquare()/fit_func->GetNDF()<max_ChiSQndf)&&(fit_func->GetParameter(1)<0.3)){

	    // Check whether the validity function gives a cutoff value bigger or smaller than the maximum experimental value
	    // and assign the cutoff 
	    if(validity_func->GetX(1)>=Pt2_fit_max){
	      Pt2_cutoff[i]   = Pt2_fit_max;
	    }
	    else{
	      Pt2_cutoff[i]   = fit_func->GetX(1);
	    }

	    ChiSQndf_array[i] = fit_func->GetChisquare()/fit_func->GetNDF();
	    ndf_array[i]      = fit_func->GetNDF();
	    ftarget->cd();
	    h_clone->Write(Form("corr_data_Pt2_"+material+"_%i%i%i_fit%i",Q2_bin,Nu_bin,Zh_bin,i));

	    fits_completed++;
	  }//END FIT QUALLITY CHECK
	  else{
	    Pt2_cutoff[i]     = 0;
	    ChiSQndf_array[i] = 0;
	    ndf_array[i]      = 0;
	  }
	}//END FIT LOOP

	// Obtaining final cutoff
	Double_t weighted_Pt2_cutoff = 0;

	// If no fit passed the quality test write the standard distribution with the same names
	if(fits_completed==0){
	  ftarget->cd();
	  h->Write();
	  h_clone->Write(Form("corr_data_Pt2_"+material+"_%i%i%i_CLEAN",Q2_bin,Nu_bin,Zh_bin));
	  GetInterpolatedPt2Distribution(h_clone);
	  h_clone->Write(Form("corr_data_Pt2_"+material+"_%i%i%i_CLEAN_INTERPOLATED",Q2_bin,Nu_bin,Zh_bin));

	  meanPt2->SetBinContent(Zh_bin+1, h->GetMean());
	  meanPt2->SetBinError(Zh_bin+1, h->GetMeanError());
	  meanPt2_CLEAN->SetBinContent(Zh_bin+1, h->GetMean());
	  meanPt2_CLEAN->SetBinError(Zh_bin+1, h->GetMeanError());
	  meanPt2_CLEAN_INTERPOLATED->SetBinContent(Zh_bin+1, h_clone->GetMean());
	  meanPt2_CLEAN_INTERPOLATED->SetBinError(Zh_bin+1, h_clone->GetMeanError());

	  continue;
	}
	// If one fit passed quality test return the only non-null value in the Pt2_cutoff array
	else if(fits_completed==1){
	  weighted_Pt2_cutoff = TMath::MaxElement(number_of_fits,Pt2_cutoff);
	  std::cout<<"Single fit completed. Cutoff="<<weighted_Pt2_cutoff<<std::endl;
	}
	// If more than one fit passes the quality tesst, obtain the weighted average 
	else{
	  weighted_Pt2_cutoff = GetWeightedPt2Cutoff(ChiSQndf_array, ndf_array, Pt2_cutoff, number_of_fits);
	  std::cout<<"Multiple fit completed. Cutoff="<<weighted_Pt2_cutoff<<std::endl;
	}
	//END GETTING CUTOFF

	// Removal of events passed the Pt2 cutoff
	GetCleanPt2Distribution(h_clone, weighted_Pt2_cutoff);
	TH1F* h_clone2 = (TH1F*)h_clone->Clone();

	// Interpolation of the resulting distribution
	GetInterpolatedPt2Distribution(h_clone2);

	ftarget->cd();
	h->Write();
	h_clone->Write(Form("corr_data_Pt2_"+material+"_%i%i%i_CLEAN",Q2_bin,Nu_bin,Zh_bin));	
	h_clone2->Write(Form("corr_data_Pt2_"+material+"_%i%i%i_CLEAN_INTERPOLATED",Q2_bin,Nu_bin,Zh_bin));

	meanPt2->SetBinContent(Zh_bin+1, h->GetMean());
	meanPt2->SetBinError(Zh_bin+1, h->GetMeanError());
	meanPt2_CLEAN->SetBinContent(Zh_bin+1, h_clone->GetMean());
	meanPt2_CLEAN->SetBinError(Zh_bin+1, h_clone->GetMeanError());
	meanPt2_CLEAN_INTERPOLATED->SetBinContent(Zh_bin+1, h_clone2->GetMean());
	meanPt2_CLEAN_INTERPOLATED->SetBinError(Zh_bin+1, h_clone2->GetMeanError());

	gROOT->cd();
	delete h;
	delete h_clone;
	delete h_clone2;
      } //END ZH LOOP

      ftarget->cd();
      
      meanPt2->Write(Form("meanPt2_"+material+"_%i%i",Q2_bin,Nu_bin));
      meanPt2_CLEAN->Write(Form("meanPt2_"+material+"_%i%i_CLEAN",Q2_bin,Nu_bin));
      meanPt2_CLEAN_INTERPOLATED->Write(Form("meanPt2_"+material+"_%i%i_CLEAN_INTERPOLATED",Q2_bin,Nu_bin));

      gROOT->cd();
      delete meanPt2;
      delete meanPt2_CLEAN;
      delete meanPt2_CLEAN_INTERPOLATED;
      
    } //END NU LOOP
  } //END Q2 LOOP
}

Int_t GetFirstEmptyBin(TH1F* h)
{
  for(Int_t bin = 1 ; bin <= h->GetNbinsX() ; bin++){
    if(h->GetBinContent(bin)==0){
	  return bin;
	}
    }
}

Double_t GetWeightedPt2Cutoff(Double_t* ChiSQndf_array,Double_t* ndf_array,Double_t* Pt2_cutoff,Int_t array_size){
  Double_t total_ChiSQndf	= 0;
  Double_t total_ndf		= 0;
  Double_t total_weight		= 0;
  Double_t weight		= 0;
  Double_t sum_weightPt2	= 0;
  
  //VARIABLES TOTAL WEIGHT LOOP
  for(Int_t index = 0 ; index < array_size ; index++){
    // This loop only considers usable cutoff values
    if(TMath::IsNaN(Pt2_cutoff[index])||Pt2_cutoff[index]==0){continue;}
    total_ChiSQndf += ChiSQndf_array[index];
    total_ndf += ndf_array[index];
  }
  //END VARIABLES TOTAL WEIGHT LOOP

  //SUM WEIGHT*PT2 LOOP
  for(Int_t index = 0 ; index <array_size ; index++){
    // This loop only considers usable cutoff values
    if(TMath::IsNaN(Pt2_cutoff[index])||Pt2_cutoff[index]==0){continue;}
    weight         = (TMath::Gaus(ChiSQndf_array[index],1,0.2)*ChiSQndf_weight + (ndf_array[index]/total_ndf)*ndf_weight)/(ChiSQndf_weight+ndf_weight);
    
    //    std::cout<<"ChiSQndf="<<ChiSQndf_array[index]<<"  ndf="<<ndf_array[index]<<"   weight="<<weight<<std::endl;
    
    total_weight  += weight;
    sum_weightPt2 += weight*Pt2_cutoff[index];
  }
  //END SUM WEIGHT*PT2 LOOP

  std::cout<<"weighted_Pt2_cutoff="<<sum_weightPt2/total_weight<<std::endl;
  return sum_weightPt2/total_weight;
}

void GetCleanPt2Distribution(TH1F* h, Double_t cutoff)
{
  for(Int_t Pt2_bin = 1 ; Pt2_bin <= h->GetNbinsX() ; Pt2_bin++){
    if(h->GetBinCenter(Pt2_bin)>cutoff){
      for(Int_t Pt2_sub_bin = Pt2_bin ; Pt2_sub_bin <= h->GetNbinsX() ; Pt2_sub_bin++){
	h->SetBinContent(Pt2_sub_bin, 0.);
	h->SetBinError(Pt2_sub_bin, 0.);
      }
      break;
    }
  }
}

void GetInterpolatedPt2Distribution(TH1F* h)
{
  for(Int_t bin = 1 ; bin <= h->GetNbinsX() ; bin++)
    {
      if(h->GetBinContent(bin)==0)
	{
	  Double_t prev_bin = h->GetBinContent(bin-1);
	  Double_t next_bin = h->GetBinContent(bin+1);
	      
	  Int_t counter = 0;
	  Int_t exit_counter = 0;
	  while(next_bin==0.){
	    counter++;
	    next_bin = h->GetBinContent(bin+1+counter);
	    if(bin+1+counter == N_Pt2){
	      exit_counter++;
	      break;
	    }
	  }

	  if(exit_counter!=0){break;}

	  Double_t prev_bin_error  = h->GetBinError(bin-1);
	  Double_t next_bin_error  = h->GetBinError(bin+1+counter);
	  Double_t prev_bin_center = h->GetBinCenter(bin-1);
	  Double_t next_bin_center = h->GetBinCenter(bin+1+counter);
	  Double_t missing_bin_center = h->GetBinCenter(bin);
	      
	  Double_t interpolated_Pt2 = GetInterpolatedPt2Point(missing_bin_center,prev_bin,next_bin,prev_bin_center,next_bin_center);
	  h->SetBinContent(bin, interpolated_Pt2);
	      
	  Double_t interpolated_Pt2error = GetInterpolatedPt2PointError(missing_bin_center,prev_bin,next_bin,prev_bin_center,next_bin_center);
	  h->SetBinError(bin, interpolated_Pt2error);
	}
    }
}

Double_t GetInterpolatedPt2Point(Double_t x,Double_t prev_bin, Double_t next_bin, Double_t prev_bin_center, Double_t next_bin_center)
{
  Double_t m     = (next_bin - prev_bin)/(next_bin_center - prev_bin_center);
  Double_t y_ref = prev_bin;
  Double_t x_ref = prev_bin_center;

  //returns straigh line interpolation
  return (x - x_ref)*m + y_ref;
}

Double_t GetInterpolatedPt2PointError(Double_t x, Double_t prev_bin_error, Double_t next_bin_error, Double_t prev_bin_center, Double_t next_bin_center)
{
  Double_t m     = (next_bin_error - prev_bin_error)/(next_bin_center - prev_bin_center);
  Double_t y_ref = prev_bin_error;
  Double_t x_ref = prev_bin_center;

  //returns straigh line interpolation
  return (x - x_ref)*m + y_ref;
}
