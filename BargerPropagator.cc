#include "BargerPropagator.h"

BargerPropagator::BargerPropagator()
{
//Earth = NULL;
//std::cout << "BargerPropagator:: " << Earth << " address " << &Earth <<std::endl;; 
   Earth = new EarthDensity( );	
//std::cout << "BargerPropagator:: " << Earth << " address " << &Earth <<std::endl;; 
   init();
}


BargerPropagator::BargerPropagator( bool k )
{
   Earth = new EarthDensity( );	
   init();
}


BargerPropagator::~BargerPropagator( )
{
   delete Earth;
}

BargerPropagator::BargerPropagator( const char * f )
{
   Earth = new EarthDensity( f );	
   init();
}

void BargerPropagator::init()
{
   kUseMassEigenstates = false;

   //rad earth in [cm] /
   ProductionHeight = 0.0;
   PathLength = 0.0;

   // default is neutral matter
   density_convert = 0.5;

   kAntiMNSMatrix     = false ;
   kSuppressWarnings  = false ;

   kOneDominantMass   = true  ;

   // Default is to choose the first octant when converting from 
   // sin^2 (2x) variables 
   kSx12Octant = 1 ;
   kSx13Octant = 1 ;
   kSx23Octant = 1 ;
}



void BargerPropagator::propagate( int NuFlavor ){

   int    i,j;
   int    Layers;
   double TransitionMatrix[3][3][2];
   //double TransitionProduct[3][3][2]; // Use global one
   double TransitionTemp[3][3][2];	
   double RawInputPsi[3][2];
   double OutputPsi[3][2];


   if( ! kSuppressWarnings )
   if(   
       ( kAntiMNSMatrix && NuFlavor > 0) ||  
       (!kAntiMNSMatrix && NuFlavor < 0)  
     )
   {
      std::cout << " Warning BargerPropagator::propagate - " << std::endl;       
      std::cout << "     Propagating neutrino flavor and MNS matrix definition differ :" << std::endl;       
      std::cout << "     MNS Matrix was defined for : " << ( kAntiMNSMatrix ? " Nubar " : "Nu" )<< std::endl; 
      std::cout << "     Propagation is for         : " << ( NuFlavor < 0   ? " Nubar " : "Nu" )<< std::endl; 
      std::cout << "     Please check your call to BargerPropagator::SetMNS() " << std::endl; 
      std::cout << "     This message can be suppressed with a call to BargerPropagator::SuppressWarnings() " << std::endl;
   
      exit(-1);
   }  

   clear_complex_matrix( TransitionMatrix );
   clear_complex_matrix( TransitionProduct );
   clear_complex_matrix( TransitionTemp );
	
   ClearProbabilities();

	
   Earth->SetDensityProfile( CosineZenith, PathLength, ProductionHeight );
   Layers = Earth->get_LayersTraversed( );

	
   for ( i = 0; i < Layers ; i++ )
   {
      get_transition_matrix(  NuFlavor  , 
                              Energy	,		   // in GeV
                              Earth->get_DensityInLayer(i) * density_convert, 
                              Earth->get_DistanceAcrossLayer(i)/1.0e5,     // in km
                              TransitionMatrix,			           // Output transition matrix
                              0.0  					   // phase offset 
                              );			
      		
      if ( i == 0 )
         copy_complex_matrix( TransitionMatrix , TransitionProduct );

      if ( i >0 ){
         clear_complex_matrix( TransitionTemp );	
         multiply_complex_matrix( TransitionMatrix, TransitionProduct, TransitionTemp ); 
         copy_complex_matrix( TransitionTemp, TransitionProduct );
      }//for other layers
    }// end of layer loop

	
   // loop on neutrino types
   for ( i = 0 ; i < 3 ; i++ )
   {
      for ( j = 0 ; j < 3 ; j++ )
      {  RawInputPsi[j][0] = 0.0; RawInputPsi[j][1] = 0.0;   }
      
      if( kUseMassEigenstates )	      
         convert_from_mass_eigenstate( i+1, NuFlavor,  RawInputPsi );
      else
         RawInputPsi[i][0] = 1.0;

 
      multiply_complex_matvec( TransitionProduct, RawInputPsi, OutputPsi );
      Probability[i][0] += OutputPsi[0][0] * OutputPsi[0][0] + OutputPsi[0][1]*OutputPsi[0][1];
      Probability[i][1] += OutputPsi[1][0] * OutputPsi[1][0] + OutputPsi[1][1]*OutputPsi[1][1];
      Probability[i][2] += OutputPsi[2][0] * OutputPsi[2][0] + OutputPsi[2][1]*OutputPsi[2][1];

   }//end of neutrino loop
	
	
}




void BargerPropagator::ClearProbabilities()
{
   for ( int i = 0 ; i < 3; i++ )
      for ( int j = 0 ; j < 3 ; j++ )
         Probability[i][j] = 0.0;

}

// Reload the same mixing matrix elements , but allow for adjustments 
// based on the neutrino eneryg and flavor
void BargerPropagator::ResetMNS( double energy, int nutype )
{

   SetMNS( fx12, fx13, fx23, fm21, fmAtm, fdelta, energy , fSquared , nutype );
}


void BargerPropagator::SetMNS( double x12, double x13, double x23, 
                               double m21, double mAtm, double delta, 
                               double Energy_ , bool kSquared, int kNuType )
{
   Energy = Energy_;

   fx12      = x12   ;
   fx13      = x13   ;
   fx23      = x23   ;
   fm21      = m21   ;
   fmAtm     = mAtm  ; 
   fdelta    = delta ;
   fSquared  = kSquared ;
 



   double sin12;
   double sin13;
   double sin23;

   double lm32 = mAtm ;
   // Dominant Mixing mode assumes the user 
   // simply changes the sign of the input atmospheric 
   // mixing to invert the hierarchy 
   //  so the input for  NH corresponds to m32  
   // and the input for  IH corresponds to m31
   if( kOneDominantMass ) 
   {
      // For the inverted Hierarchy, adjust the input
      // by the solar mixing (should be positive) 
      // to feed the core libraries the correct value of m32
      if( mAtm < 0.0 )  
          lm32 = mAtm - m21 ;
   }
   else 
   {
       if( !kSuppressWarnings )
       {
         std::cout << " BargerPropagator::SetMNS - " << std::endl;       
         std::cout << "     You have opted to specify the value of m23 by yourself. " << std::endl;       
         std::cout << "     This means you must correct the value of m23 when switching " << std::endl; 
         std::cout << "     between the mass hierarchy options. " << std::endl; 
         std::cout << "     This message can be suppressed with BargerPropagator::SuppressWarnings()"<< std::endl;
      }
   } 



   //if xAB = sin( xAB )^2
   if ( kSquared )
   {
      sin12 = sqrt( x12 );
      sin13 = sqrt( x13 );
      sin23 = sqrt( x23 );
   }
   else
   {
      //if xAB = sin( 2 xAB )^2
      // Default is to specify sin(x) in the first octant 
      // but this may be changed by the user to the second octant 
      // (mostly only an issue for atmospheric mixing) 
      if( kSx12Octant == 1 ) sin12 = sqrt( 0.5*(1 - sqrt(1 - x12 ))  );
      else                   sin12 = sqrt( 0.5*(1 + sqrt(1 - x12 ))  );

      if( kSx13Octant == 1 ) sin13 = sqrt( 0.5*(1 - sqrt(1 - x13 ))  );
      else                   sin13 = sqrt( 0.5*(1 + sqrt(1 - x13 ))  );
   
      if( kSx23Octant == 1 ) sin23 = sqrt( 0.5*(1 - sqrt(1 - x23 ))  );
      else                   sin23 = sqrt( 0.5*(1 + sqrt(1 - x23 ))  );
   
   }

   if ( kNuType < 0 )
   {
     delta *= -1.0 ;
     kAntiMNSMatrix = true ;
   }
   else 
   {
     kAntiMNSMatrix = false ;
   }

   init_mixing_matrix( m21, lm32, sin12, sin23, sin13, delta );

}

void BargerPropagator::DefinePath(double cz, double ProdHeight, bool kSetProfile )
{

   ProductionHeight = ProdHeight*1e5;
   REarth = Earth->GetEarthRadiuskm() * 1.0e5;
   PathLength = sqrt( (REarth + ProductionHeight )*(REarth + ProductionHeight) 
                     - (REarth*REarth)*( 1 - cz*cz)) - REarth*cz;
   CosineZenith = cz;
   if( kSetProfile )
      Earth->SetDensityProfile( CosineZenith, PathLength, ProductionHeight );
	
}


void BargerPropagator::SetMatterPathLength()
{

   int Layers = Earth->get_LayersTraversed( );

   MatterPathLength = 0.0;
   AirPathLength = 0.0;
   for( int i = 1 ; i < Layers ; i++ )
      MatterPathLength +=  Earth->get_DistanceAcrossLayer(i);

   AirPathLength +=  Earth->get_DistanceAcrossLayer(0);

}
 

void BargerPropagator::SetAirPathLength(double x)
{
// argument is [km], convert to [cm]
   AirPathLength = x*1.0e5 - MatterPathLength;
}



double BargerPropagator::GetVacuumProb( int Alpha, int Beta , double Energy, double Path )
{
   // alpha -> 1:e 2:mu 3:tau
   // Energy[GeV]
   // Path[km]
   /// simple referes to the fact that in the 3 flavor analysis 
   //  the solar mass term is zero
   double Probs[3][3];


   get_vacuum_probability( Alpha, Energy, Path, Probs );

   Alpha = abs(Alpha);
   Beta = abs(Beta);

   if ( Alpha > 0 )
      return Probs[Alpha-1][Beta-1]; 

   if ( Alpha < 0 ) // assuming CPT!!!
      return Probs[Beta-1][Alpha-1]; 

   std::cerr << " BargerPropagator::GetVacuumProb neutrino must be non-zero: " << std::endl;
   return -1.0; 
   
}



void BargerPropagator::propagateLinear( int NuFlavor, double pathlength, double Density )
{

   int    i,j;

   double TransitionMatrix[3][3][2];
   //double TransitionProduct[3][3][2]; // use global one
   double TransitionTemp[3][3][2];	
   double RawInputPsi[3][2];
   double OutputPsi[3][2];

   if( ! kSuppressWarnings )
   if(   
       ( kAntiMNSMatrix && NuFlavor > 0) ||  
       (!kAntiMNSMatrix && NuFlavor < 0)     
     )
   {
      std::cout << " Warning BargerPropagator::propagateLinear - " << std::endl;       
      std::cout << "     Propagating neutrino flavor and MNS matrix definition differ :" << std::endl;       
      std::cout << "     MNS Matrix was defined for : " << ( kAntiMNSMatrix ? " Nubar " : "Nu" )<< std::endl; 
      std::cout << "     Propagation is for         : " << ( NuFlavor < 0   ? " Nubar " : "Nu" )<< std::endl; 
      std::cout << "     Please check your call to BargerPropagator::SetMNS() " << std::endl; 
      std::cout << "     This message can be suppressed with a call to BargerPropagator::SuppressWarnings() " << std::endl;
   
      exit(-1);
   }  

   clear_complex_matrix( TransitionMatrix );
   clear_complex_matrix( TransitionProduct );
   clear_complex_matrix( TransitionTemp );
	
   ClearProbabilities();

   
   get_transition_matrix( NuFlavor, 
                  Energy	,		// in GeV
                  Density * density_convert, 
                  pathlength ,      	// in km
                  TransitionMatrix,	// Output transition matrix
                  0.0  			
            );			
		
   copy_complex_matrix( TransitionMatrix , TransitionProduct );

   for ( i = 0 ; i < 3 ; i++ )
   {
      for ( j = 0 ; j < 3 ; j++ )
      {       RawInputPsi[j][0] = 0.0; RawInputPsi[j][1] = 0.0;   }

      if( kUseMassEigenstates )	      
         convert_from_mass_eigenstate( i+1, NuFlavor,  RawInputPsi );
      else
         RawInputPsi[i][0] = 1.0;
   	      
      multiply_complex_matvec( TransitionProduct, RawInputPsi, OutputPsi );

      Probability[i][0] += OutputPsi[0][0] * OutputPsi[0][0] + OutputPsi[0][1]*OutputPsi[0][1];
      Probability[i][1] += OutputPsi[1][0] * OutputPsi[1][0] + OutputPsi[1][1]*OutputPsi[1][1];
      Probability[i][2] += OutputPsi[2][0] * OutputPsi[2][0] + OutputPsi[2][1]*OutputPsi[2][1];

   }// end of loop on neutrino types
	
}


double BargerPropagator::GetPathAveragedDensity( )
{
   unsigned Layers = (unsigned) Earth->get_LayersTraversed( );
	
   double density_sum = 0.;
   double length_sum = 0.;

   for ( unsigned i = 0; i < Layers ; i++ )
   {
      density_sum += Earth->get_DistanceAcrossLayer(i) * Earth->get_DensityInLayer(i) * density_convert ;
      length_sum  += Earth->get_DistanceAcrossLayer(i);
   }

   return density_sum / length_sum ;

}
	

void BargerPropagator::SetDefaultOctant( int var   , int octant )
{

  // don't accept bad octants
  if( octant != 1 && octant != 2 ) return;

  if( var == 12 ) kSx12Octant = octant ;
  if( var == 13 ) kSx13Octant = octant ;
  if( var == 23 ) kSx23Octant = octant ;

  // otherwise do nothing 

}


