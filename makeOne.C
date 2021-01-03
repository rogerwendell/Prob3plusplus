{

  gStyle->SetLineStyleString(5,"[]");
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTitleSize(0.04);
  gStyle->SetLabelSize(0.04);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetPadBorderMode(0);



  gStyle->SetPalette(1); 

	TFile * Bob = _file0;
	Bob->cd();


// Raw Probability Plots	

	TH2F *_NuEToNuE3f 	=  (TH2F*) Bob->Get("NuEToNuE3f");
	TH2F *_NuEToNuMu3f 	=  (TH2F*) Bob->Get("NuEToNuMu3f");
	TH2F *_NuEToNuTau3f 	=  (TH2F*) Bob->Get("NuEToNuTau3f");
	TH2F *_NuEToNuX3f 	=  (TH2F*) Bob->Get("NuEToNuX3f");
	
	
	TH2F *_NuMuToNuE3f 	=  (TH2F*) Bob->Get("NuMuToNuE3f");
	TH2F *_NuMuToNuMu3f 	=  (TH2F*) Bob->Get("NuMuToNuMu3f");
	TH2F *_NuMuToNuTau3f 	=  (TH2F*) Bob->Get("NuMuToNuTau3f");
	TH2F *_NuMuToNuX3f 	=  (TH2F*) Bob->Get("NuMuToNuX3f");

	TH2F *_NuTauToNuE3f 	=  (TH2F*) Bob->Get("NuTauToNuE3f");
	TH2F *_NuTauToNuMu3f 	=  (TH2F*) Bob->Get("NuTauToNuMu3f");
	TH2F *_NuTauToNuTau3f 	=  (TH2F*) Bob->Get("NuTauToNuTau3f");
	TH2F *_NuTauToNuX3f 	=  (TH2F*) Bob->Get("NuTauToNuX3f");
//------------------- End of Raw Probability Plots	


  TCanvas *c4 = new TCanvas("c4","",5,5,800,800);
  c4->SetLeftMargin(0.153226);
  c4->SetRightMargin(0.153226);
  c4->SetBottomMargin(0.153226);
  c4->SetTopMargin(0.153226);
  c4->SetRightMargin(0.153226);
  c4->SetFillColor(kWhite);
  c4->SetLogx(1);

  gStyle->SetPalette(1); 
  gStyle->SetPalette(1); 

	
	_NuMuToNuE3f->GetXaxis()->SetTitleSize(0.03);
	_NuMuToNuE3f->GetXaxis()->SetTitleOffset(2.0);
	_NuMuToNuE3f->GetYaxis()->SetTitleSize(0.03);
	_NuMuToNuE3f->GetYaxis()->SetRangeUser(-1.00, 1.00);
	//_NuMuToNuE3f->GetZaxis()->SetRangeUser(0.00, 0.55);
	_NuMuToNuE3f->GetYaxis()->SetTitleOffset(2.5);
	//_NuMuToNuMu3f->SetTitle->("Prob 3 Flavor NuMu->NuE");
	_NuMuToNuE3f->GetXaxis()->SetTitle("Energy [GeV]");
	_NuMuToNuE3f->GetYaxis()->SetTitle("Cosine Zenith Angle");
	_NuMuToNuE3f->Draw("cont4z");

	
  TCanvas *c5 = new TCanvas("c5","",5,5,800,800);
  c5->SetLeftMargin(0.153226);
  c5->SetRightMargin(0.153226);
  c5->SetBottomMargin(0.153226);
  c5->SetTopMargin(0.153226);
  c5->SetRightMargin(0.153226);
  c5->SetFillColor(kWhite);
  c5->SetLogx(1);

	
	_NuMuToNuMu3f->GetYaxis()->SetRangeUser(-1.00, 1.00);
	_NuMuToNuMu3f->GetXaxis()->SetTitleSize(0.03);
	_NuMuToNuMu3f->GetXaxis()->SetTitleOffset(2.0);
	_NuMuToNuMu3f->GetYaxis()->SetTitleSize(0.03);
	_NuMuToNuMu3f->GetYaxis()->SetTitleOffset(2.5);
	//_NuMuToNuMu3f->SetTitle->("Prob 3 Flavor NuMu->NuMu");
	_NuMuToNuMu3f->GetXaxis()->SetTitle("Energy [GeV]");
	_NuMuToNuMu3f->GetYaxis()->SetTitle("Cosine Zenith Angle");
	_NuMuToNuMu3f->Draw("cont4z");


  TCanvas *c6 = new TCanvas("c6","",5,5,800,800);
  c6->SetLeftMargin(0.153226);
  c6->SetRightMargin(0.153226);
  c6->SetBottomMargin(0.153226);
  c6->SetTopMargin(0.153226);
  c6->SetRightMargin(0.153226);
  c6->SetFillColor(kWhite);
  c6->SetLogx(1);


	_NuMuToNuTau3f->GetYaxis()->SetRangeUser(-1.00, 1.00);
	_NuMuToNuTau3f->GetXaxis()->SetTitleSize(0.03);
	_NuMuToNuTau3f->GetXaxis()->SetTitleOffset(2.0);
	_NuMuToNuTau3f->GetYaxis()->SetTitleSize(0.03);
	_NuMuToNuTau3f->GetYaxis()->SetTitleOffset(2.5);
	//_NuMuToNuMu3f->SetTitle->("Prob 3 Flavor NuMu->NuTau");
	_NuMuToNuTau3f->GetXaxis()->SetTitle("Energy [GeV]");
	_NuMuToNuTau3f->GetYaxis()->SetTitle("Cosine Zenith Angle");
	_NuMuToNuTau3f->Draw("cont4z");


	
//------------------- End of Raw Probability Plots	




}


















































