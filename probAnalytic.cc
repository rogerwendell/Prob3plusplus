///
// Analytic mu2e probabilities in vacuum
//

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

double GetAnalyticMu2E( double x12, double x13, double x23, double m21, 
                        double m23, double Delta, double Energy_ , bool kSquared, int, double, double density = 0. );

using namespace std;
int main(int argc, char * argv[] )
{

  double dcp_in = 0.;
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

// Path Length
  double PathLengthEdge[NBinsPath+1];
  double BasePath = 295.0;
  path_start = 0.1;
  path_end = 1.0e2;
  path_step = log10( path_end/path_start)/double(NBinsPath);


// Energy Range
  double EnergyBins[NBinsEnergy+1];
  double BaseEnergy = 0.004;
  e_start = 0.010;
  e_end  =  10.0;
  e_step = log10(e_end/e_start)/double(NBinsEnergy);

     
/// Oscillation Parameters
  bool kSquared  = true;   // using sin^2(x) variables?

  int    kNuBar  =  1 * v_in;
  double DM2     =  h_in * 2.4e-3;
  double Theta23 =  0.5   ;
  double Theta13 =  0.025  ;
  double dm2     =  7.6e-5;
  double Theta12 =  0.312;
  double delta   =  dcp_in * (3.1415926/180.0);
  double density =  0.0  ; // g/cm3

  std::cout << "Using          " << std::endl
            << "      DM2      " <<  DM2      << std::endl
            << "      Theta23  " <<  Theta23  << std::endl
            << "      Theta13  " <<  Theta13  << std::endl
            << "      dm2      " <<  dm2      << std::endl
            << "      Theta12  " <<  Theta12  << std::endl;



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

  stringstream ssE, ssL;
  TH1D * histos[3][2];

  ///////////////////////////
  /// mu to E 
  ssE.str(""); ssE <<  "P(#nu_{#mu} #rightarrow #nu_{e})" << " L = " << BasePath ; 
  ssL.str(""); ssL <<  "P(#nu_{#mu} #rightarrow #nu_{e})" << " E = " << BaseEnergy; 
  TH1D * lmu2eE = new TH1D("lmu2eE", ssE.str().c_str() , NBinsEnergy  -1 , EnergyBins );
  TH1D * lmu2eL = new TH1D("lmu2eL", ssL.str().c_str() , NBinsPath  -1 , PathLengthEdge );

  histos[0][0] = lmu2eE; 
  histos[0][1] = lmu2eL; 

  double prob;
  for ( i = 0 ; i <= NBinsEnergy ; i ++ ) 
  {
     energy = e_start*pow(10.0, double(i)*e_step);
     prob = GetAnalyticMu2E( Theta12,  Theta13, Theta23, dm2, DM2, delta , 
                             energy, kSquared, kNuBar, BasePath, density ); 

     for( j = 0 ; j < 1 ; j++ )
        histos[j][0]->Fill( energy, prob );     

  } // End Energy Loop //

  for ( i = 0 ; i <= NBinsPath ; i ++ ) 
  {
    path = path_start*pow(10.0, double(i)*path_step ); 
    prob = GetAnalyticMu2E( Theta12,  Theta13, Theta23, dm2, DM2, delta , 
                             energy, kSquared, kNuBar, path, density ); 

    for( j = 0 ; j < 1 ; j++ )
        histos[j][1]->Fill( path , prob );     

  } // End Path Loop //



  TFile *tmp = new TFile("Analytic.root", "recreate");
  tmp->cd();

  for( j = 0 ; j < 1 ; j++ ){
     histos[j][0]->Write();     
     histos[j][1]->Write();     
  }

  tmp->Close();
         
  cout << endl<<"Done Cowboy!" << endl;
  return 0;
}


// This formula is exact for vacuum oscillations
// And only approximate otherwise 
// Based on B.Richter hep-ph/0008222
double GetAnalyticMu2E( double x12, double x13, double x23, 
                        double m21, double m23, double Delta, 
                        double Energy_ , bool kSquared, int kNuBar, double L, double density )
{
  double Energy = Energy_;

  // 2*sqrt(2)*Gfermi in (eV^2-cm^3)/(mole-GeV) - for e<->[mu,tau] //
  double tworttwoGf = 1.52588e-4;
  double A          = tworttwoGf * density * 0.5 * Energy  ; 
  double Vmat       = A * L / (4.0 * Energy );

  if( kNuBar < 0 ) 
  {  
     Delta *= -1.0 ;
     Vmat  *= -1.0 ;
     A     *= -1.0 ;
  }

  double cd = cos( Delta );
  double sd = sin( Delta );

  double s12;
  double s13;
  double s23;

  double c12;
  double c13;
  double c23;
 
  double s213;
  double s212;
  double s223;

  //if xAB = sin( xAB )^2
  if ( kSquared ){

    s12 = sqrt( x12 );
    s13 = sqrt( x13 );
    s23 = sqrt( x23 );

     s213 = sin( 2.0 * asin( s13 ));
     s212 = sin( 2.0 * asin( s12 ));
     s223 = sin( 2.0 * asin( s23 ));

     c12 = sqrt( 1.0 - x12 );
     c13 = sqrt( 1.0 - x13 );
     c23 = sqrt( 1.0 - x23 );
   }
   else
   {
      //if xAB = sin( 2 xAB )^2
      s12 = sqrt( 0.5*(1 - sqrt(1 - x12 ))  );
      s13 = sqrt( 0.5*(1 - sqrt(1 - x13 ))  );
      s23 = sqrt( 0.5*(1 - sqrt(1 - x23 ))  );

      c12 = sqrt( 1.0 - s12*s12 ); 
      c13 = sqrt( 1.0 - s13*s13 ); 
      c23 = sqrt( 1.0 - s23*s23 ); 

      s213 = sqrt( x13 );
      s212 = sqrt( x12 );
      s223 = sqrt( x23 );

   }

   double d21 = (1.2667 * m21       * L / Energy ); 
   double d32 = (1.2667 * m23       * L / Energy ); 
   double d31 = (1.2667 * (m23+m21) * L / Energy ); 

   double dom =  4.0 * c13 * c13 * s13 * s13 * s23 * s23 * sin( d31) *sin( d31 ); 

   double cpc =  8.0 * c13 * c13 * s12 * s13 * s23 * 
                      ( c12 *c23 * cd - s12*s13*s23 ) * cos( d32 )* sin( d31 ) * sin( d21 ); 
   
   double cpv = -8.0 * c13*c13 * c12 * c23 * s12 * s13 * s23 * sd * sin( d32 ) * sin( d31 ) * sin( d21 );

   double sol =  4.0 * s12 * s12 * c13 * c13 * (    c12 * c12 * c23 * c23   
                                                 + s12 * s12 * s23 * s23 * s13 * s13  
                                         - 2.0   * c12 * c23 * s12 * s23 * s13 * cd ) * sin(d21) * sin( d21 ); 


// This form is only valid when A / m31 is small
   double mat  =   dom * 2.0 * A * ( 1 - 2.0 * s13 * s13 ) / (m23+m21) ; 
          mat += - 8.0 * c13 * c13 * s13 * s13 * s23 * s23 * ( 1 - 2.0 * s13 * s13 ) * 
                            Vmat  * cos( d32 ) * sin( d31 );


   double prob = dom + cpc + cpv + sol + mat ;

   return prob;

}



