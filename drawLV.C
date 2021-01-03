{

   gStyle->SetOptStat(0);

   Double_t navy[3]    = {0, 0, 0.5};
   Double_t blue[3]    = {0, 0.2, 1.0};
   Double_t cyan[3]    = {0, 1, 1};
   Double_t yellow[3]  = {1, 1, 0};
   Double_t orange[3]  = {1, 0.5, 0};
   Double_t magenta[3] = {1, 0, 1};
   Double_t red[3]     = {1, 0, 0};
   Double_t white[3]   = {1, 1, 1};
   Double_t black[3]   = {0, 0, 0};


   Double_t Red[5]    = {black[0], red[0], orange[0], yellow[0], white[0]};
   Double_t Green[5]  = {black[1], red[1], orange[1], yellow[1], white[1]};
   Double_t Blue[5]   = {black[2], red[2], orange[2], yellow[2], white[2]};
   Double_t Length[5] = {   0.00,    0.25,       0.5,      0.75,     1.00};

   TColor::CreateGradientColorTable(5,Length,Red,Green,Blue,50);

   TCanvas * cmain2 = new TCanvas("cmain2", "" , 800, 600 );
   cmain2->SetFillColor(0);

   cmain2->SetLogy(1);
   cmain2->SetLogx(1);

   //TH2D * f = (TH2D*) _file0->Get("NuMuToNuE3f");
   TH2D * fmumu = (TH2D*) _file0->Get("NuMuToNuMu3f");

   fmumu->GetZaxis()->SetRangeUser(0.0,1.0);
   
   fmumu->GetXaxis()->SetTitle("Energy [GeV]");
   fmumu->GetYaxis()->SetTitle("Pathlength [km]");

   fmumu->Draw("col4z");
   cmain2->SaveAs("numu_to_numu_LV.png");


   TH2D * fmue = (TH2D*) _file0->Get("NuMuToNuE3f");
   fmue->GetZaxis()->SetRangeUser(0.0,1.0);
   
   fmue->GetXaxis()->SetTitle("Energy [GeV]");
   fmue->GetYaxis()->SetTitle("Pathlength [km]");

   fmue->Draw("col4z");
   cmain2->SaveAs("numu_to_nue_LV.png");

}
