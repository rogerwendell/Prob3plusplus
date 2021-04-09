#include <math.h>
#include <iostream>
#include <fstream>

#include "BargerPropagator.h"

#include "TFile.h"
#include "TH1D.h"

#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>

using namespace std;
int main(int argc, char * argv[] )
{

  double dcp_in =  0.;
  double h_in   = 1.0;
  int    v_in   = 1  ;

  if( argc >= 2 ) dcp_in = (double) atof( argv[1] );
  if( argc >= 3 ) h_in   = (double) atof( argv[2] );
  if( argc >= 4 ) v_in   = (int)    atoi( argv[3] );

  h_in = ( h_in > 0 ? 1.0 : -1.0 );

  double total_prob=0.0;
  double path, energy;
  double e_start, e_end, e_step, path_start, path_end, path_step;
     double d_start, d_end, d_step;
  int i, j ;


//// Binning     
  int NBinsEnergy = 10000;
  int NBinsPath   = 10000; 
  int NBinsDensity = 1000;
  int NBinsLoE = 1000;

//// LoE
  double LoEEdge[ NBinsLoE +1 ];
  double loe_start = 10.0;
  double loe_end = 4.0e3;
  double loe_step = log10( loe_end / loe_start )/double(NBinsLoE);

//// Density
  double DensityEdge[NBinsDensity+1];
  double Density = 2.6;
  d_start = 0.0;
  d_end = 200.0;
  d_step = ( d_end - d_start )/double(NBinsDensity);


// Path Length
  double PathLengthEdge[NBinsPath+1];
  double BasePath = 295.0;
  //double BasePath = 1.5e5;
  path_start = 0.1;
  path_end = 2.0e5;
  path_step = log10( path_end/path_start)/double(NBinsPath);


// Energy Range
  double EnergyBins[NBinsEnergy+1];
  double BaseEnergy = 0.010;
  //e_start = 1.0e-3; //GeV 
  e_start = 1.0; //GeV 
  e_end  =  10.0  ;
  e_step = log10(e_end/e_start)/double(NBinsEnergy);

     
/// Oscillation Parameters
  bool kSquared  = true;   // using sin^2(x) variables?

  int    kNuBar  =  1 * v_in;
  double DM2     =  h_in * 2.4e-3;
  double Theta23 =  0.5    ;
  double Theta13 =  0.0256 ;
  double dm2     =  7.6e-5;
  double Theta12 =  0.312;
  double delta   =  dcp_in * (3.1415926/180.0);

  std::cout << "Using          " << std::endl
            << "      DM2      " <<  DM2      << std::endl
            << "      Theta23  " <<  Theta23  << std::endl
            << "      Theta13  " <<  Theta13  << std::endl
            << "      dm2      " <<  dm2      << std::endl
            << "      Theta12  " <<  Theta12  << std::endl
            << "      dcp      " <<  delta    << std::endl
            << "      nu/nubar " <<  kNuBar   << std::endl;

  BargerPropagator   * bNu; 

  bNu = new BargerPropagator( );
  bNu->UseMassEigenstates( false );


  double Entry = e_start;
  for(i=0; i<NBinsEnergy; i++ )
  {
                Entry = e_start*pow( 10.0 , double(i)*e_step );
                EnergyBins[i] = Entry;
  }
  EnergyBins[NBinsEnergy] = EnergyBins[NBinsEnergy-1]*1.001;

     
  PathLengthEdge[0]= path_start;
  for ( i=1; i<NBinsPath ; i++ ){
  Entry = path_start * pow( 10.0, double(i)*path_step );
  PathLengthEdge[i] = Entry;
  }
  PathLengthEdge[NBinsPath] = PathLengthEdge[NBinsPath-1]*1.001;

  DensityEdge[0]= d_start;
  for ( i=1; i<NBinsDensity ; i++ ){
  Entry = d_start + double(i)*d_step ;
  DensityEdge[i] = Entry;
  }
  DensityEdge[NBinsDensity] = DensityEdge[NBinsDensity-1]*1.001;

  LoEEdge[0]= loe_start;
  for ( i=1; i<NBinsLoE ; i++ ){
  Entry = loe_start*pow(10.0, double(i)*loe_step) ;
  LoEEdge[i] = Entry;
  }
  LoEEdge[NBinsLoE] = LoEEdge[NBinsLoE-1]*1.001;



  stringstream ssE, ssL;
  TH1D * histos[4][2];

  /// mu to E 
  ssE.str(""); ssE <<  "P(#nu_{e} #rightarrow #nu_{e})" << " L = " << BasePath ; 
  ssL.str(""); ssL <<  "P(#nu_{e} #rightarrow #nu_{e})" << " E = " << BaseEnergy; 
  TH1D * le2eE = new TH1D("le2eE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );
  TH1D * le2eL = new TH1D("le2eL", ssL.str().c_str() , NBinsPath  -1 , PathLengthEdge );
  
  ///////////////////////////
  /// mu to E 
  ssE.str(""); ssE <<  "P(#nu_{#mu} #rightarrow #nu_{e})" << " L = " << BasePath ; 
  ssL.str(""); ssL <<  "P(#nu_{#mu} #rightarrow #nu_{e})" << " E = " << BaseEnergy; 
  TH1D * lmu2eE = new TH1D("lmu2eE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );
  TH1D * lmu2eL = new TH1D("lmu2eL", ssL.str().c_str() , NBinsPath  -1 , PathLengthEdge );

  ///////////////////////////
  /// mu to mu 
  ssE.str(""); ssE <<  "P(#nu_{#mu} #rightarrow #nu_{#mu})" << " L = " << BasePath ; 
  ssL.str(""); ssL <<  "P(#nu_{#mu} #rightarrow #nu_{#mu})" << " E = " << BaseEnergy; 
  TH1D * lmu2muE = new TH1D("lmu2muE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );
  TH1D * lmu2muL = new TH1D("lmu2muL", ssL.str().c_str() , NBinsPath  -1 , PathLengthEdge );

  ///////////////////////////
  /// mu to tau 
  ssE.str(""); ssE <<  "P(#nu_{#mu} #rightarrow #nu_{#tau})" << " L = " << BasePath ; 
  ssL.str(""); ssL <<  "P(#nu_{#mu} #rightarrow #nu_{#tau})" << " E = " << BaseEnergy; 
  TH1D * lmu2tauE = new TH1D("lmu2tauE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );
  TH1D * lmu2tauL = new TH1D("lmu2tauL", ssL.str().c_str() , NBinsPath  -1 , PathLengthEdge );

  histos[0][0] = lmu2eE; 
  histos[0][1] = lmu2eL; 

  histos[1][0] = lmu2muE; 
  histos[1][1] = lmu2muL; 

  histos[2][0] = lmu2tauE; 
  histos[2][1] = lmu2tauL; 

  histos[3][0] = le2eE; 
  histos[3][1] = le2eL; 

  // Density
  ssE.str(""); ssE <<  "P(#nu_{e} #rightarrow #nu_{e})" << " L = " << BasePath ; 
  ssE <<  " E = " << BaseEnergy; 
  TH1D * density = new TH1D("Density", ssE.str().c_str() , NBinsDensity  -1 , DensityEdge );
  TH1D * LoE = new TH1D("loe", "", NBinsLoE - 1, LoEEdge );


  for ( i = 0 ; i <= NBinsEnergy ; i ++ ) 
  {
     energy = e_start*pow(10.0, double(i)*e_step);
     bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar ); 
     bNu->propagateLinear( 1*kNuBar, BasePath, Density );

     total_prob = 0.0;
     for(int m=1; m<=3; m++)
       total_prob += bNu->GetProb(2, m); // Normalize the Probabilities //

     if ( total_prob >1.00001 || total_prob<0.99998 )
     { cerr << "ERROR Prob:" << "Energy: "<< energy << " " << endl;   abort();   }


     total_prob = 0.0;
     for(int m=1; m<=3; m++)
       total_prob += bNu->GetProb(1, m); // Normalize the Probabilities //

     if ( total_prob >1.00001 || total_prob<0.99998 )
     { cerr << "ERROR Prob:" << "Energy: "<< energy << " " << endl;   abort();   }


     for( j = 0 ; j < 3 ; j++ )
        histos[j][0]->Fill( energy, bNu->GetProb(2,j+1) );     

     histos[3][0]->Fill( energy, bNu->GetProb(1,1) );     

  } // End Energy Loop //


  for ( i = 0 ; i <= NBinsPath ; i ++ ) 
  {
     path = path_start*pow(10.0, double(i)*path_step ); 
     bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , BaseEnergy, kSquared, kNuBar ); 
     bNu->propagateLinear( 1*kNuBar, path, Density );

     for( j = 0 ; j < 3 ; j++ )
        histos[j][1]->Fill( path , bNu->GetProb(2,j+1) );     

        histos[3][1]->Fill( path , bNu->GetProb(1,1) );     

  } // End Path Loop //

  for ( i = 0 ; i <= NBinsDensity ; i ++ ) 
  {
     Density = d_start +  double(i)*d_step ; 
     bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , BaseEnergy, kSquared, kNuBar ); 
     bNu->propagateLinear( 1*kNuBar, BasePath , Density );

     density->Fill( Density , bNu->GetProb(2,1) );     
  }

  double loe;
  double lSin23 = 1.0;
  for ( i = 0 ; i <= NBinsLoE ; i ++ ) 
  {
     loe = loe_start* pow(10.0, double(i)*loe_step ); 
     LoE->Fill( loe , 1.0 - lSin23*sin( loe * 1.2667 * DM2 )*sin( loe * 1.2667 * DM2 ) );     
  }

  /////
  // Write the output
  TFile *tmp = new TFile("LinearProb.root", "recreate");
  tmp->cd();

  for( j = 0 ; j < 4 ; j++ ){
     histos[j][0]->Write();     
     histos[j][1]->Write();     
  }
  density->Write();
  LoE->Write();

  tmp->Close();

  cout << endl<<"Done Cowboy!" << endl;
  return 0;
}

