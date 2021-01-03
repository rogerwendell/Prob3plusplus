#include <iostream>

#include "BargerPropagator.h"
#include "FullSMEPropagator.h"

#include "TFile.h"
#include "TH2D.h"

int main(int argc, char * argv[] )
{

   double dcp_in =  0.;
   double h_in   = 1.0;
   int    v_in   = 1.0;

   if( argc >= 2 ) dcp_in = (double) atof( argv[1] );
   if( argc >= 3 ) h_in   = (double) atof( argv[2] );
   if( argc >= 4 ) v_in   = (int)    atoi( argv[3] );

   h_in = ( h_in > 0 ? 1.0 : -1.0 );

   double PathLength, energy;
   double e_start, e_end, e_step, path_start, path_end, path_step;
   int i, j ;

//// Types of Averaging	
   int NBinsEnergy = 1000;
   int PathLengthNbin = 2000; 


// Pathlength Range [km]
   double PathLengthEdge[PathLengthNbin+1];
   path_start = 25.1 ;
   path_end = 12767.0 ;
   path_step = ( path_end - path_start)/double(PathLengthNbin);


// Energy Range [GeV]
   double EnergyBins[NBinsEnergy+1];
   e_start = 0.110000001;
   e_end = 3000.;
   e_step = log10(e_end/e_start)/double(NBinsEnergy);

   
/// Oscillation Parameters
   bool kSquared = true ;   // are we using sin^2(x) variables?

   int    kNuBar  =  1 * v_in;
   double DM2     =  h_in * 2.5e-3;
   double Theta23 =  0.50   ;
   double Theta13 =  0.0219 ;
   double dm2     =  7.6e-5 ;
   double Theta12 =  0.302  ;
   double delta   =  dcp_in * (3.1415926/180.0);

   std::cout << "Using          " << std::endl
             << "      DM2      " <<  DM2      << std::endl
             << "      Theta23  " <<  Theta23  << std::endl
             << "      Theta13  " <<  Theta13  << std::endl
             << "      dm2      " <<  dm2      << std::endl
             << "      Theta12  " <<  Theta12  << std::endl
             << "      delta    " <<  delta    << std::endl;
   std::cout << "      knubar   " <<  kNuBar   << std::endl;

   std::cout << "From "
	     << " [ " << e_start << " - " << e_end << " ] GeV " << endl;
   
// Methods to Compute Probability
   FullSMEPropagator * SMEnu = new FullSMEPropagator( );

   double Entry = e_start;
   for(i=0; i<NBinsEnergy; i++ )
        {
                Entry = e_start*pow( 10.0 , double(i)*e_step );
                EnergyBins[i] = Entry;
   }
   EnergyBins[NBinsEnergy] = EnergyBins[NBinsEnergy-1]*1.001;

   
   PathLengthEdge[0]= path_start*0.9999;
   for ( i=1; i<PathLengthNbin ; i++ )
   	PathLengthEdge[i] = PathLengthEdge[0] + double(i)*path_step;

   PathLengthEdge[PathLengthNbin] = PathLengthEdge[PathLengthNbin-1]*1.001;

   
   TH2D *NuEToNuE3f 	=  new TH2D("NuEToNuE3f","3 Flavor P_{#nu_{e}#rightarrow#nu_{e}}",
   					NBinsEnergy  -1 , EnergyBins, PathLengthNbin -1, PathLengthEdge);
   TH2D *NuEToNuMu3f 	=  new TH2D("NuEToNuMu3f","3 Flavor P_{#nu_{e}#rightarrow#nu_{#mu}}",
   					NBinsEnergy  -1 , EnergyBins, PathLengthNbin -1, PathLengthEdge);
   TH2D *NuEToNuTau3f 	=  new TH2D("NuEToNuTau3f","3 Flavor P_{#nu_{e}#rightarrow#nu_{#tau}}",
   					NBinsEnergy  -1 , EnergyBins, PathLengthNbin -1, PathLengthEdge);
   TH2D *NuEToNuX3f 	=  new TH2D("NuEToNuX3f","3 Flavor P_{#nu_{e}#rightarrow#nu_{x}}",
   					NBinsEnergy  -1 , EnergyBins, PathLengthNbin -1, PathLengthEdge);
   

   
   TH2D *NuMuToNuE3f 	=  new TH2D("NuMuToNuE3f","3 Flavor P_{#nu_{#mu}#rightarrow#nu_{e}}",
   					NBinsEnergy  -1 , EnergyBins, PathLengthNbin -1, PathLengthEdge);
   TH2D *NuMuToNuMu3f 	=  new TH2D("NuMuToNuMu3f","3 Flavor P_{#nu_{#mu}#rightarrow#nu_{#mu}}",
   					NBinsEnergy  -1 , EnergyBins, PathLengthNbin -1, PathLengthEdge);
   TH2D *NuMuToNuTau3f 	=  new TH2D("NuMuToNuTau3f","3 Flavor P_{#nu_{#mu}#rightarrow#nu_{#tau}}",
   					NBinsEnergy  -1 , EnergyBins, PathLengthNbin -1, PathLengthEdge);
   TH2D *NuMuToNuX3f 	=  new TH2D("NuMuToNuX3f","3 Flavor P_{#nu_{#mu}#rightarrow#nu_{x}}",
   					NBinsEnergy  -1 , EnergyBins, PathLengthNbin -1, PathLengthEdge);

   
   TH2D *NuTauToNuE3f 	=  new TH2D("NuTauToNuE3f","3 Flavor P_{#nu_{#tau}#rightarrow#nu_{e}}",
   					NBinsEnergy  -1 , EnergyBins, PathLengthNbin -1, PathLengthEdge);
   TH2D *NuTauToNuMu3f 	=  new TH2D("NuTauToNuMu3f","3 Flavor P_{#nu_{#tau}#rightarrow#nu_{#mu}}",
   					NBinsEnergy  -1 , EnergyBins, PathLengthNbin -1, PathLengthEdge);
   TH2D *NuTauToNuTau3f 	=  new TH2D("NuTauToNuTau3f","3 Flavor P_{#nu_{#tau}#rightarrow#nu_{#tau}}",
   					NBinsEnergy  -1 , EnergyBins, PathLengthNbin -1, PathLengthEdge);
   TH2D *NuTauToNuX3f 	=  new TH2D("NuTauToNuX3f","3 Flavor P_{#nu_{#tau}#rightarrow#nu_{x}}",
   					NBinsEnergy  -1 , EnergyBins, PathLengthNbin -1, PathLengthEdge);


   
   TH2D *NuMuToNuTau2f 	=  new TH2D("NuMuToNuTau2f","2 Flavor P_{#nu_{#mu}#rightarrow#nu_{#tau}}",
   					NBinsEnergy  -1 , EnergyBins, PathLengthNbin -1, PathLengthEdge);

  
   for ( i = 0 ; i < 3 ; i ++ )
     for ( j = 0 ; j < 3 ; j ++ )
     {
       SMEnu->SetLVMatrixEntry( "A" , i, j, 0.0, 0.0 );
       SMEnu->SetLVMatrixEntry( "C" , i, j, 0.0, 0.0 );
     }

// N.B. Matrix must be hermitian

// SMEnu->SetLVMatrixEntry( "C" , 0, 1, 7.5e-23, 0.0    );
// SMEnu->SetLVMatrixEntry( "C" , 1, 0, 7.5e-23, 0.0    );

// SMEnu->SetLVMatrixEntry( "C" , 0, 2, 7.5e-23, 0.0    );
// SMEnu->SetLVMatrixEntry( "C" , 2, 0, 7.5e-23, 0.0    );

   SMEnu->SetLVMatrixEntry( "C" , 1, 2, 7.5e-23, 0.0    );
   SMEnu->SetLVMatrixEntry( "C" , 2, 1, 7.5e-23, 0.0    );


   // matix must be hermitian
// SMEnu->SetLVMatrixEntry( "A" , 0, 1, 1.0e-22, 0.0    );
// SMEnu->SetLVMatrixEntry( "A" , 1, 0, 1.0e-22, 0.0    );

// SMEnu->SetLVMatrixEntry( "A" , 0, 2, 1.0e-22, 0.0    );
// SMEnu->SetLVMatrixEntry( "A" , 2, 0, 1.0e-22, 0.0    );

// SMEnu->SetLVMatrixEntry( "A" , 1, 2, 1.0e-22, 0.0    );
// SMEnu->SetLVMatrixEntry( "A" , 2, 1, 1.0e-22, 0.0    );


   SMEnu->PrintLVMatrices();
   double total = 0.0;
   // fill the probabilities
   for ( i = 0 ; i <= NBinsEnergy ; i ++ ) 
   {

     energy = e_start*pow(10.0, double(i)*e_step);
     for ( j = 0 ; j <= PathLengthNbin ; j++ )
     { 

        PathLength = path_start + double(j)*path_step;
  
        SMEnu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar ); 
   
        // For propagation in the Earth
        SMEnu->DefinePathFromLength( PathLength, 25.00  );

        // For propagation in constant density matter, here 22 g/cm^3
        //SMEnu->DefineLinearPath( PathLength, 22.0 );

        SMEnu->SetHamiltonian( energy );

        // routine is the same for propagation in the Earth and in 
        // constant density matter
        SMEnu->propagate( kNuBar );

        total = 0.0;
	total += SMEnu->GetProb(1,1);
	total += SMEnu->GetProb(1,2);
	total += SMEnu->GetProb(1,3);

	if( fabs( 1.00 - total ) > 1.0e-7 ) abort();
	
	
        NuEToNuE3f 	->Fill( energy, PathLength, SMEnu->GetProb(1,1) ) ;
        NuEToNuMu3f 	->Fill( energy, PathLength, SMEnu->GetProb(1,2) ) ;
        NuEToNuTau3f 	->Fill( energy, PathLength, SMEnu->GetProb(1,3) ) ;
        NuEToNuX3f 	->Fill( energy, PathLength, 1.0 - SMEnu->GetProb(1,1) ) ;

        NuMuToNuE3f 	->Fill( energy, PathLength, SMEnu->GetProb(2,1) ) ;
        NuMuToNuMu3f 	->Fill( energy, PathLength, SMEnu->GetProb(2,2) ) ;
        NuMuToNuTau3f 	->Fill( energy, PathLength, SMEnu->GetProb(2,3) ) ;
        NuMuToNuX3f 	->Fill( energy, PathLength, 1.0 - SMEnu->GetProb(2,2) ) ;

        NuTauToNuE3f 	->Fill( energy, PathLength, SMEnu->GetProb(3,1) ) ;
        NuTauToNuMu3f 	->Fill( energy, PathLength, SMEnu->GetProb(3,2) ) ;
        NuTauToNuTau3f 	->Fill( energy, PathLength, SMEnu->GetProb(3,3) ) ;
        NuTauToNuX3f 	->Fill( energy, PathLength, 1.0 - SMEnu->GetProb(3,3) ) ;
     }// End of Cosine Z Looping //
   
   } // End Energy Loop //


   TFile *tmp = new TFile("fullSME_Prob.root", "recreate");
   tmp->cd();

   NuEToNuE3f    ->Write();
   NuEToNuMu3f   ->Write();
   NuEToNuTau3f  ->Write();
   NuEToNuX3f    ->Write();
   
   NuMuToNuE3f   ->Write();
   NuMuToNuMu3f  ->Write();
   NuMuToNuTau3f ->Write();
   NuMuToNuX3f   ->Write();

   NuTauToNuE3f  ->Write();
   NuTauToNuMu3f ->Write();
   NuTauToNuTau3f->Write();
   NuTauToNuX3f  ->Write();

   NuMuToNuTau2f ->Write();

   tmp->Close();
   	
   std::cout << std::endl<<"Done LV Cowboy!" << std::endl;
   return 0;
}

