#include "FullSMEPropagator.h"

const double FullSMEPropagator::LtoeV_Coeff = 5.0677e9;


FullSMEPropagator::FullSMEPropagator():BargerPropagator()
{

  for(int i=0; i<3; i++)
    for(int j=0; j<i; j++) {
        SetLVMatrixEntry( "A" , i, j, 0.0 , 0.0 );
        SetLVMatrixEntry( "C" , i, j, 0.0 , 0.0 );
   }

}

FullSMEPropagator::FullSMEPropagator(const char * txt):BargerPropagator( txt )
{
  for(int i=0; i<3; i++)
    for(int j=0; j<i; j++) {
        SetLVMatrixEntry( "A" , i, j, 0.0 , 0.0 );
        SetLVMatrixEntry( "C" , i, j, 0.0 , 0.0 );
   }

}

void FullSMEPropagator::DefinePathFromLength(double LinKM, double ProdHeight)
{
   REarth = Earth->GetEarthRadiuskm();
   double minl = ProdHeight;
   double maxl = ProdHeight+2*REarth;

   if (LinKM < minl) {
      std::cout << "Warning: In DefinePathFromLength, " << LinKM << " km is too short, setting to " << minl << std::endl;
      LinKM = minl;
   }
   if (LinKM > maxl) {
      std::cout << "Warning: In DefinePathFromLength, " << LinKM << " km is too long, setting to " << maxl << std::endl; 
      LinKM = maxl;
   }

// ProductionHeight = ProdHeight*1e5;
// PathLength       = LinKM *1e5;

   ProductionHeight = ProdHeight;
   PathLength       = LinKM     ;
   
   CosineZenith = (ProductionHeight*ProductionHeight - PathLength*PathLength + 2*ProductionHeight*REarth)
       / (2*PathLength*REarth);

   fLocalPath = PathLength ;

   // also sets internal density parameter based on CosineZenith
   GetDensityAvg();

}


double FullSMEPropagator::GetDensityAvg()
{
  int    i;
  int    Layers;
  
  Earth->SetDensityProfile( CosineZenith, PathLength*1e5 , ProductionHeight );
  Layers = Earth->get_LayersTraversed( );

  double distance = 0.;
  fAveDensity = 0.;
  
  for ( i = 0; i < Layers ; i++ )
    {
      double thisL = Earth->get_DistanceAcrossLayer(i) / 1.e5; // in km
      distance    += thisL;
      fAveDensity += Earth->get_DensityInLayer(i) * thisL;
    }
    fAveDensity /= distance;

    return fAveDensity;
}

#include <math.h>
#include <TMath.h>

#include <iomanip>


void FullSMEPropagator::SetMNS( double x12, double x13, double x23,
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

   fm31 = mAtm + m21 ;
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
          fm31 = mAtm - m21 ;
   }
   else 
   {
       if( !kSuppressWarnings )
       {
         std::cout << " FullSMEPropagator::SetMNS - " << std::endl;       
         std::cout << "     You have opted to specify the value of m23 by yourself. " << std::endl;       
         std::cout << "     This means you must correct the value of m23 when switching " << std::endl; 
         std::cout << "     between the mass hierarchy options. " << std::endl; 
         std::cout << "     This message can be suppressed with FullSMEPropagator::::SuppressWarnings()"<< std::endl;
      }
   } 

   fDM2[0] = 0    ; 
   fDM2[1] = m21  ;
   fDM2[2] = fm31 ; 


   //if xAB = sin( xAB )^2
   if ( kSquared )
   {
      sin12 =  x12 ;
      sin13 =  x13 ;
      sin23 =  x23 ;
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

   // LIV Routines expect sin^2 
   //  
      sin12 *= sin12 ;
      sin13 *= sin13 ;
      sin23 *= sin23 ;
   
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

   SetLVUPMNSmatrix( sin12, sin23, sin13, delta );

}

void FullSMEPropagator::DefineLinearPath( double l, double rho )
{

  fLocalPath  = l ;
  fAveDensity = rho ;


}


void FullSMEPropagator::SetHamiltonian( double Energy ) 
{
    
   //Matter Effects Hamiltonian
   if ( kAntiMNSMatrix )
   {
     fHLV  = BuildLVHamiltonian( Energy , true );
     fHmat = Calc_Hmat( -1*fAveDensity );
   }
   else
   {
     fHLV  = BuildLVHamiltonian( Energy , false );
     fHmat = Calc_Hmat( fAveDensity );
   }

   //double L = fLocalPath * LtoeV_Coeff / 1.e5; // cm to km
   double L = fLocalPath * LtoeV_Coeff; 
   fHosc = Calc_H0( fPMNS, L , Energy, fDM2);

}

void FullSMEPropagator::DefinePath(double cz, double ProdHeight, bool kSetProfile  )
{

   BargerPropagator::DefinePath(cz, ProdHeight, kSetProfile );
   fLocalPath = PathLength ; 

}



vector< complex<double> > FullSMEPropagator::BuildLVHamiltonian( double E , bool kAnti)
{

  int dimA = 3 ;
  int dimC = 4 ;

  double GeVToeV = 1.0e9;

  complex<double> tmpA;
  complex<double> tmpC;

  double reA, reC;
  double imA, imC;

  //Set now the perturbative (though size may be large) hamiltonian
  vector< complex<double> > dh;
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++) {

      reA  = fALV[i][j].real() * pow(GeVToeV ,4-dimA) * pow(E*GeVToeV ,dimA-3) ;
      imA  = fALV[i][j].imag() * pow(GeVToeV ,4-dimA) * pow(E*GeVToeV ,dimA-3) ;

      reC  = -1. * fCLV[i][j].real() * pow(GeVToeV ,4-dimC) * pow(E*GeVToeV ,dimC-3);
      imC  = -1. * fCLV[i][j].imag() * pow(GeVToeV ,4-dimC) * pow(E*GeVToeV ,dimC-3);

      if( kAnti ) 
      {
         reA *= -1 ;
         imC *= -1 ;
      }
      tmpA = complex<double> ( reA, imA ); 
      tmpC = complex<double> ( reC, imC ); 
     
      dh.push_back( tmpA + (4.0/3.0) * tmpC  );
  }

  return dh;

}

void FullSMEPropagator::propagateLinear( int NuFlavor, double pathlength, double Density )
{

  std::cout << "FullSMEPropagator::propagateLinear -- unused. Instead call DefineLinearPath( pathlength , Density ); " << std::endl;
  std::cout << "                                      then propagate( NuFlavor);" << std::endl;

  exit(-1);
  
}


void FullSMEPropagator::propagate( int knu )
{

   complex<double>** ABCNDeltaE = Calc_ABCNDeltaE( fHosc, fHLV, fHmat);

   // length needs to be in eV
   double L = fLocalPath * LtoeV_Coeff; 

   Probability [0][0] = Calc_Prob_ee     (ABCNDeltaE, L);
   Probability [0][1] = Calc_Prob_emu    (ABCNDeltaE, L);
   Probability [0][2] = Calc_Prob_etau   (ABCNDeltaE, L);

   Probability [1][0] = Calc_Prob_mue    (ABCNDeltaE, L);
   Probability [1][1] = Calc_Prob_mumu   (ABCNDeltaE, L);
   Probability [1][2] = Calc_Prob_mutau  (ABCNDeltaE, L);

   Probability [2][0] = Calc_Prob_taue   (ABCNDeltaE, L);
   Probability [2][1] = Calc_Prob_taumu  (ABCNDeltaE, L);
   Probability [2][2] = Calc_Prob_tautau (ABCNDeltaE, L);
    
   for(int i=0; i<5; i++)
      delete[] ABCNDeltaE[i];
   delete[] ABCNDeltaE;

}
    


////////////////////////////////////////////
//    Added part for Full SME analysis    //
////////////////////////////////////////////

vector< complex<double> > FullSMEPropagator::Calc_H0(complex<double> PMNS[3][3], double L, double E, double DM2_tmp[3]) 
{


  double DM2_12 = DM2_tmp[1]/(2*E*1e9);
  double DM2_23 = DM2_tmp[2]/(2*E*1e9);

  double DM2[3]={ 0.,
		  DM2_12,
		  DM2_23 };
  
  vector< complex<double> > H0;
  complex<double> tmp;

  //Invert i and j so that we fill in term of lines
  for(int j=0; j<3; j++)
    for(int i=0; i<3; i++) {
      tmp = complex<double> (0,0);
      for(int k=0; k<3; k++)
	tmp += PMNS[k][i]*conj(PMNS[k][j])*DM2[k];
      
      H0.push_back(tmp);
    }
  
 return H0;
}


vector< complex<double> > FullSMEPropagator::AddCplxMatrices(vector< complex<double> > M1, vector< complex<double> > M2) {

  vector< complex<double> > H;
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      H.push_back( M1.at(j+3*i) + M2.at(j+3*i) );

  return H;
}


vector< complex<double> > FullSMEPropagator::Calc_ProdMatrices(vector< complex<double> > M1, vector< complex<double> > M2) {

  vector< complex<double> > Prod;
  complex<double> tmp;

  for(int k=0; k<3; k++) {
    for(int j=0; j<3; j++) {
      tmp = complex<double> (0,0);
      for(int i=0; i<3; i++)
        tmp += M1.at(i+3*k) * M2.at(3*i+j);
      
      Prod.push_back(tmp);
    }
  }
  
  
  return Prod;
}


complex<double> FullSMEPropagator::Calc_Trace(vector< complex<double> > H) {
  
  complex<double> tr (0,0);
  for(int i=0; i<3; i++)
    tr += H.at(i+3*i);
  
  return tr;
}


complex<double> FullSMEPropagator::Calc_Determinant(vector< complex<double> > H) {
  
  complex<double> det (0,0);
  det = H.at(0)*(H.at(4)*H.at(8)-H.at(5)*H.at(7)) - H.at(1)*(H.at(3)*H.at(8)-H.at(5)*H.at(6)) + H.at(2)*(H.at(3)*H.at(7)-H.at(4)*H.at(6));
  
  return det;
}



double FullSMEPropagator::Calc_Prob_ee(complex<double>** ABCNDeltaE, double L) {

  complex<double> A[3], B[3], C[3], N[3]; 
  double DeltaE[3];
  for(int i=0; i<3; i++) {
    A[i] = ABCNDeltaE[0][i];
    B[i] = ABCNDeltaE[1][i];
    C[i] = ABCNDeltaE[2][i];
    N[i] = ABCNDeltaE[3][i];
    DeltaE[i] = ABCNDeltaE[4][i].real();
  }

  double Prob = 1.
    - 4.*norm((B[1]*B[0]*C[1]*C[0])/(N[1]*N[0]))*pow(sin(DeltaE[0]*L/2.),2)
    - 4.*norm((B[2]*B[0]*C[2]*C[0])/(N[2]*N[0]))*pow(sin(DeltaE[2]*L/2.),2)
    - 4.*norm((B[2]*B[1]*C[2]*C[1])/(N[2]*N[1]))*pow(sin(DeltaE[1]*L/2.),2) ;

  return Prob;
}


double FullSMEPropagator::Calc_Prob_emu(complex<double>** ABCNDeltaE, double L) {

  complex<double> A[3], B[3], C[3], N[3]; double DeltaE[3];
  for(int i=0; i<3; i++) {
    A[i] = ABCNDeltaE[0][i];
    B[i] = ABCNDeltaE[1][i];
    C[i] = ABCNDeltaE[2][i];
    N[i] = ABCNDeltaE[3][i];
    DeltaE[i] = ABCNDeltaE[4][i].real();
  }
  
  double Prob = -4*(A[1]*conj(A[0])*B[1]*conj(B[0])*norm(C[1]*C[0]/(N[1]*N[0]))).real()*pow(sin(DeltaE[0]*L/2.),2)
    +2*(A[1]*conj(A[0])*B[1]*conj(B[0])*norm(C[1]*C[0]/(N[1]*N[0]))).imag()*sin(DeltaE[0]*L)
    -4*(A[2]*conj(A[0])*B[2]*conj(B[0])*norm(C[2]*C[0]/(N[2]*N[0]))).real()*pow(sin(DeltaE[2]*L/2.),2)
    +2*(A[2]*conj(A[0])*B[2]*conj(B[0])*norm(C[2]*C[0]/(N[2]*N[0]))).imag()*sin(DeltaE[2]*L)
    -4*(A[2]*conj(A[1])*B[2]*conj(B[1])*norm(C[2]*C[1]/(N[2]*N[1]))).real()*pow(sin(DeltaE[1]*L/2.),2)
    +2*(A[2]*conj(A[1])*B[2]*conj(B[1])*norm(C[2]*C[1]/(N[2]*N[1]))).imag()*sin(DeltaE[1]*L);

  return Prob;
}


double FullSMEPropagator::Calc_Prob_etau(complex<double>** ABCNDeltaE, double L) 
{

  complex<double> A[3], B[3], C[3], N[3]; double DeltaE[3];
  for(int i=0; i<3; i++) {
    A[i] = ABCNDeltaE[0][i];
    B[i] = ABCNDeltaE[1][i];
    C[i] = ABCNDeltaE[2][i];
    N[i] = ABCNDeltaE[3][i];
    DeltaE[i] = ABCNDeltaE[4][i].real();
  }

  double Prob = -4*(A[1]*conj(A[0])*conj(C[1])*C[0]*pow((B[1]*conj(B[0])/(N[1]*N[0])),2)).real()*pow(sin(DeltaE[0]*L/2.),2)
    +2*(A[1]*conj(A[0])*conj(C[1])*C[0]*pow((B[1]*conj(B[0])/(N[1]*N[0])),2)).imag()*sin(DeltaE[0]*L)
    -4*(A[2]*conj(A[0])*conj(C[2])*C[0]*pow((B[2]*conj(B[0])/(N[2]*N[0])),2)).real()*pow(sin(DeltaE[2]*L/2.),2)
    +2*(A[2]*conj(A[0])*conj(C[2])*C[0]*pow((B[2]*conj(B[0])/(N[2]*N[0])),2)).imag()*sin(DeltaE[2]*L)
    -4*(A[2]*conj(A[1])*conj(C[2])*C[1]*pow((B[2]*conj(B[1])/(N[2]*N[1])),2)).real()*pow(sin(DeltaE[1]*L/2.),2)
    +2*(A[2]*conj(A[1])*conj(C[2])*C[1]*pow((B[2]*conj(B[1])/(N[2]*N[1])),2)).imag()*sin(DeltaE[1]*L);

  return Prob;
}


double FullSMEPropagator::Calc_Prob_mue(complex<double>** ABCNDeltaE, double L) 
{

  complex<double> A[3], B[3], C[3], N[3]; double DeltaE[3];
  for(int i=0; i<3; i++) {
    A[i] = ABCNDeltaE[0][i];
    B[i] = ABCNDeltaE[1][i];
    C[i] = ABCNDeltaE[2][i];
    N[i] = ABCNDeltaE[3][i];
    DeltaE[i] = ABCNDeltaE[4][i].real();
  }

  double Prob = -4*(conj(A[1])*A[0]*conj(B[1])*B[0]*norm(C[1]*C[0]/(N[1]*N[0]))).real()*pow(sin(DeltaE[0]*L/2.),2)
    +2*(conj(A[1])*A[0]*conj(B[1])*B[0]*norm(C[1]*C[0]/(N[1]*N[0]))).imag()*sin(DeltaE[0]*L)
    -4*(conj(A[2])*A[0]*conj(B[2])*B[0]*norm(C[2]*C[0]/(N[2]*N[0]))).real()*pow(sin(DeltaE[2]*L/2.),2)
    +2*(conj(A[2])*A[0]*conj(B[2])*B[0]*norm(C[2]*C[0]/(N[2]*N[0]))).imag()*sin(DeltaE[2]*L)
    -4*(conj(A[2])*A[1]*conj(B[2])*B[1]*norm(C[2]*C[1]/(N[2]*N[1]))).real()*pow(sin(DeltaE[1]*L/2.),2)
    +2*(conj(A[2])*A[1]*conj(B[2])*B[1]*norm(C[2]*C[1]/(N[2]*N[1]))).imag()*sin(DeltaE[1]*L);

  return Prob;
}


double FullSMEPropagator::Calc_Prob_mumu(complex<double>** ABCNDeltaE, double L) 
{

    
  complex<double> A[3], B[3], C[3], N[3]; double DeltaE[3];
  for(int i=0; i<3; i++) {
    A[i] = ABCNDeltaE[0][i];
    B[i] = ABCNDeltaE[1][i];
    C[i] = ABCNDeltaE[2][i];
    N[i] = ABCNDeltaE[3][i];
    DeltaE[i] = ABCNDeltaE[4][i].real();
  }
  double Prob = 1.
      - 4.*norm((A[1]*A[0]*C[1]*C[0])/(N[1]*N[0]))*pow(sin(DeltaE[0]*L/2.),2)
      - 4.*norm((A[2]*A[0]*C[2]*C[0])/(N[2]*N[0]))*pow(sin(DeltaE[2]*L/2.),2)
      - 4.*norm((A[2]*A[1]*C[2]*C[1])/(N[2]*N[1]))*pow(sin(DeltaE[1]*L/2.),2) ;

  return Prob;
}


double FullSMEPropagator::Calc_Prob_mutau(complex<double>** ABCNDeltaE, double L) 
{

  complex<double> A[3], B[3], C[3], N[3]; double DeltaE[3];
  for(int i=0; i<3; i++) {
    A[i] = ABCNDeltaE[0][i];
    B[i] = ABCNDeltaE[1][i];
    C[i] = ABCNDeltaE[2][i];
    N[i] = ABCNDeltaE[3][i];
    DeltaE[i] = ABCNDeltaE[4][i].real();
  }

  double Prob = -4*(B[1]*conj(B[0])*conj(C[1])*C[0]*norm(A[1]*A[0]/(N[1]*N[0]))).real()*pow(sin(DeltaE[0]*L/2.),2)
    +2*(B[1]*conj(B[0])*conj(C[1])*C[0]*norm(A[1]*A[0]/(N[1]*N[0]))).imag()*sin(DeltaE[0]*L)
    -4*(B[2]*conj(B[0])*conj(C[2])*C[0]*norm(A[2]*A[0]/(N[2]*N[0]))).real()*pow(sin(DeltaE[2]*L/2.),2)
    +2*(B[2]*conj(B[0])*conj(C[2])*C[0]*norm(A[2]*A[0]/(N[2]*N[0]))).imag()*sin(DeltaE[2]*L)
    -4*(B[2]*conj(B[1])*conj(C[2])*C[1]*norm(A[2]*A[1]/(N[2]*N[1]))).real()*pow(sin(DeltaE[1]*L/2.),2)
    +2*(B[2]*conj(B[1])*conj(C[2])*C[1]*norm(A[2]*A[1]/(N[2]*N[1]))).imag()*sin(DeltaE[1]*L);

  return Prob;
}


double FullSMEPropagator::Calc_Prob_taue(complex<double>** ABCNDeltaE, double L) 
{

  complex<double> A[3], B[3], C[3], N[3]; double DeltaE[3];
  for(int i=0; i<3; i++) {
    A[i] = ABCNDeltaE[0][i];
    B[i] = ABCNDeltaE[1][i];
    C[i] = ABCNDeltaE[2][i];
    N[i] = ABCNDeltaE[3][i];
    DeltaE[i] = ABCNDeltaE[4][i].real();
  }

  double Prob = -4*(conj(A[1])*A[0]*C[1]*conj(C[0])*pow(((conj(B[1])*B[0])/(N[1]*N[0])),2)).real()*pow(sin(DeltaE[0]*L/2.),2)
    +2*(conj(A[1])*A[0]*C[1]*conj(C[0])*pow(((conj(B[1])*B[0])/(N[1]*N[0])),2)).imag()*sin(DeltaE[0]*L)
    -4*(conj(A[2])*A[0]*C[2]*conj(C[0])*pow(((conj(B[2])*B[0])/(N[2]*N[0])),2)).real()*pow(sin(DeltaE[2]*L/2.),2)
    +2*(conj(A[2])*A[0]*C[2]*conj(C[0])*pow(((conj(B[2])*B[0])/(N[2]*N[0])),2)).imag()*sin(DeltaE[2]*L)
    -4*(conj(A[2])*A[1]*C[2]*conj(C[1])*pow(((conj(B[2])*B[1])/(N[2]*N[1])),2)).real()*pow(sin(DeltaE[1]*L/2.),2)
    +2*(conj(A[2])*A[1]*C[2]*conj(C[1])*pow(((conj(B[2])*B[1])/(N[2]*N[1])),2)).imag()*sin(DeltaE[1]*L);

  return Prob;
}


double FullSMEPropagator::Calc_Prob_taumu(complex<double>** ABCNDeltaE, double L) 
{

  complex<double> A[3], B[3], C[3], N[3]; double DeltaE[3];
  for(int i=0; i<3; i++) {
    A[i] = ABCNDeltaE[0][i];
    B[i] = ABCNDeltaE[1][i];
    C[i] = ABCNDeltaE[2][i];
    N[i] = ABCNDeltaE[3][i];
    DeltaE[i] = ABCNDeltaE[4][i].real();
  }

  double Prob = -4*(conj(B[1])*B[0]*C[1]*conj(C[0])*norm(A[1]*A[0]/(N[1]*N[0]))).real()*pow(sin(DeltaE[0]*L/2.),2)
    +2*(conj(B[1])*B[0]*C[1]*conj(C[0])*norm(A[1]*A[0]/(N[1]*N[0]))).imag()*sin(DeltaE[0]*L)
    -4*(conj(B[2])*B[0]*C[2]*conj(C[0])*norm(A[2]*A[0]/(N[2]*N[0]))).real()*pow(sin(DeltaE[2]*L/2.),2)
    +2*(conj(B[2])*B[0]*C[2]*conj(C[0])*norm(A[2]*A[0]/(N[2]*N[0]))).imag()*sin(DeltaE[2]*L)
    -4*(conj(B[2])*B[1]*C[2]*conj(C[1])*norm(A[2]*A[1]/(N[2]*N[1]))).real()*pow(sin(DeltaE[1]*L/2.),2)
    +2*(conj(B[2])*B[1]*C[2]*conj(C[1])*norm(A[2]*A[1]/(N[2]*N[1]))).imag()*sin(DeltaE[1]*L);

  return Prob;
}


double FullSMEPropagator::Calc_Prob_tautau(complex<double>** ABCNDeltaE, double L) 
{

  complex<double> A[3], B[3], C[3], N[3]; double DeltaE[3];
  for(int i=0; i<3; i++) {
    A[i] = ABCNDeltaE[0][i];
    B[i] = ABCNDeltaE[1][i];
    C[i] = ABCNDeltaE[2][i];
    N[i] = ABCNDeltaE[3][i];
    DeltaE[i] = ABCNDeltaE[4][i].real();
  }

  double Prob = 1.
    - 4.*norm((A[1]*A[0]*B[1]*B[0])/(N[1]*N[0]))*pow(sin(DeltaE[0]*L/2.),2)
    - 4.*norm((A[2]*A[0]*B[2]*B[0])/(N[2]*N[0]))*pow(sin(DeltaE[2]*L/2.),2)
    - 4.*norm((A[2]*A[1]*B[2]*B[1])/(N[2]*N[1]))*pow(sin(DeltaE[1]*L/2.),2) ;

  return Prob;
}

void FullSMEPropagator::SetCMatrix( double reC[3][3] , double imC[3][3] )
{
     
  for(int i=0; i<3; i++)
    for(int j=i; j<3; j++) {
      SetLVMatrixEntry( "C" , i , j , reC[i][j] , imC[i][j] );
    }

  // check hermicity
  for(int i=0; i<3; i++)
    for(int j=0; j<i; j++) {
      if( fCLV[i][j] != conj( fCLV[j][i] ) )
        std::cout << " FullSMEPropagator::SetCMatrix : C Matrix is not hermitian at entry (" << i << "," << j << ")." << endl; 
    }
}

void FullSMEPropagator::SetAMatrix( double reA[3][3] , double imA[3][3] )
{
     
  for(int i=0; i<3; i++)
    for(int j=i; j<3; j++) {
      SetLVMatrixEntry( "A" , i , j , reA[i][j] , imA[i][j] );
    }

  // check hermicity
  for(int i=0; i<3; i++)
    for(int j=0; j<i; j++) {
      if( fALV[i][j] != conj( fALV[j][i] ) )
        std::cout << " FullSMEPropagator::SetCMatrix : A Matrix is not hermitian at entry (" << i << "," << j << ")." << endl; 
    }

}


void FullSMEPropagator::SetLVMatrixEntry( const char * matrix , int i, int j, double re , double im )
{

   bool proc = true ;
   if ( i < 0 || i > 3 )
      proc = false ;

   if ( i < 0 || i > 3 )
      proc = false ;
 
   if ( ! proc ) 
   {
      std::cout << "FullSMEPropagator::SetLVMatrixEntry : " << i << "," << j << " indices are out of bounds."  << std::endl;
      std::cout << "  Both should be in [0,2] . returning " << std::endl;
      return;
   }

   if      ( strcmp( matrix , "C") == 0  ) fCLV[i][j] = complex<double>( re, im );
   else if ( strcmp( matrix , "A") == 0  ) fALV[i][j] = complex<double>( re, im );
   else 
   {
      std::cout << "FullSMEPropagator::SetMatrixEntry : " << matrix << " is not a valid option " << std::endl;
      std::cout << "  please choose either C or A . " << std::endl;
   }

   
}



vector< complex<double> > FullSMEPropagator::Calc_Hmat(double PathAvgDensity) 
{

  double Gf = 1.166379e-23; // in eV^-2 and Ne is in #e- eV^3
  double Ne = PathAvgDensity * 0.5 * 6.0221413e23 * pow(1.9732697e-5,3); //0.5 for e- being hald the nb of nucleons?

  vector< complex<double> > Hmat;
  Hmat.push_back(sqrt(2)*Gf*Ne); for(int i=0; i<8; i++) Hmat.push_back(0);

  return Hmat;
}

void FullSMEPropagator::PrintLVMatrices()
{
   std::cout << std::endl;
   std::cout << "--- LV Matrices " << std::endl;
   std::cout << "A matrix: " << std::endl;
   PrintMatrix( fALV );

   std::cout << std::endl;

   std::cout << "C matrix: " << std::endl;
   PrintMatrix( fCLV );

}


void FullSMEPropagator::PrintMatrix(  complex<double> H[][3] ) 
{
    
    for (int r = 0; r < 3; r++){
      for (int c = 0; c < 3; c++) {
          cout << "  " << setw(4) << H[r][c] << " ";
      }
      cout << endl;
    }

}


complex<double>** FullSMEPropagator::Calc_ABCNDeltaE(vector< complex<double> > H0, vector< complex<double> > dh, 
							 vector< complex<double> > MatterEffects) 
{

  complex<double>** array2D = new complex<double>*[5];
  for(int i=0; i<5; i++)
    array2D[i] = new complex<double>[3];


  vector < complex<double> > Htmp = AddCplxMatrices(H0,dh);
  vector < complex<double> > H = AddCplxMatrices(Htmp,MatterEffects);

  complex<double> HTrace = Calc_Trace(H);
  complex<double> HDet   = Calc_Determinant(H);
  vector< complex<double> > H2 = Calc_ProdMatrices(H,H);
  complex<double> H2Trace = Calc_Trace(H2);


  //Define the variables a, b, c
  double a = -HTrace.real();
  double b = (1/2.*(HTrace*HTrace-H2Trace)).real();
  double c = (-HDet).real();

  //Define the variables Q, R, Theta
  double Q     = 1/9.*(a*a - 3.*b);
  double R     = 1/54.*(2.*pow(a,3) - 9.*a*b + 27.*c);
  double tmp = R/pow(Q,3/2.);
  if (tmp > 1.) tmp = 1.;
  if (tmp < -1.) tmp = -1.;
  double Theta = acos(tmp);
  //cout<<"Q= "<<Q<<"  R= "<<R<<"  Theta= "<<Theta<<endl;
  
  //Define New eigenvalues
  double E1 = -2.*sqrt(Q)*cos(Theta/3.) - a/3.;
  double E2 = -2.*sqrt(Q)*cos((Theta + 2.*TMath::Pi())/3.) - a/3.;
  double E3 = -2.*sqrt(Q)*cos((Theta - 2.*TMath::Pi())/3.) - a/3.;
  double Evect[3] = { E1, E2, E3 };
  double DeltaE[3] = { E2-E1, E3-E2, E3-E1 };

  //Define the variables Ai, Bi, Ci, Ni
  for(int i=0; i<3; i++) {
    array2D[0][i] = H.at(5)*(H.at(0) - Evect[i]) - H.at(3)*H.at(2); //Ai
    array2D[1][i] = H.at(6)*(H.at(4) - Evect[i]) - H.at(7)*H.at(3); //Bi
    array2D[2][i] = H.at(3)*(H.at(8) - Evect[i]) - H.at(5)*H.at(6); //Ci
    array2D[3][i] = sqrt( norm(array2D[0][i]*array2D[1][i])         //Ni
			+ norm(array2D[0][i]*array2D[2][i])
			+ norm(array2D[1][i]*array2D[2][i]) );
    array2D[4][i] = DeltaE[i];
  }

  return array2D;
}




void FullSMEPropagator::SetLVUPMNSmatrix( double ssq12, double ssq23, double ssq13, double dCP) {

  double th12=asin(sqrt(ssq12)); double c12=cos(th12); double s12=sin(th12);
  double th23=asin(sqrt(ssq23)); double c23=cos(th23); double s23=sin(th23);
  double th13=asin(sqrt(ssq13)); double c13=cos(th13); double s13=sin(th13);
  complex<double> delta= complex<double> (0.,dCP);

  //LV uses transposed matrix!! Do it by hand here - it does?
  fPMNS[0][0]= c12*c13;
  fPMNS[0][1]=-s12*c23 - c12*s23*s13*std::exp(-delta);
  fPMNS[0][2]= s12*s23 - c12*c23*s13*std::exp(-delta);
  fPMNS[1][0]= s12*c13;
  fPMNS[1][1]= c12*c23 - s12*s23*s13*std::exp(-delta);
  fPMNS[1][2]=-c12*s23 - s12*c23*s13*std::exp(-delta);
  fPMNS[2][0]= s13*exp(delta);
  fPMNS[2][1]= s23*c13;
  fPMNS[2][2]= c23*c13;
}


///////////////////////////////////////
// End of LV function implementation //
///////////////////////////////////////


