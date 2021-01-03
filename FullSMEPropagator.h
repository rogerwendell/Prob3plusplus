#ifndef _FullSMEPropagator_
#define _FullSMEPropagator_

/// 
//  Oscillation probabilities for Lorentz invariance violation under 
//  the Standard Model Extension 
//
//  Module computes full three flavor oscillations for an averaged matter 
//  density for neutrinos traveling the Earth ( ::propagate) or 
//  a specified density and pathlength (::propagateLinear) 
//
//  Follows the calculation used in the SK Lorentz Violation Search with 
//  Atmospheric neutrinos:
//     K. Abe et. al Phys.Rev.D91 052003 (2013)
//
//  Oscillation code based on work of  
//  T. Akiri and A. Himmel with help from H. Diaz
//
//


#include "BargerPropagator.h"

#include <iostream>
#include <cstdlib>
#include <complex>


class FullSMEPropagator : public BargerPropagator
{
  public:

    FullSMEPropagator( bool );
    FullSMEPropagator();
    FullSMEPropagator( const char * );
    ~FullSMEPropagator(){;}

    // Main propagation routine. 
    // driving routine for oscillations through the 
    //  -) earth at definite  zenith angle and production height.
    //     Call DefinePathFromLength() before propagate
    //  or 
    //  -) matter of constant densitiy 
    //     Call DefineLinearPath() before propagate()
    //
    // Must be called after SetMNS(), SetLVMatrix, and Define*Path* 
    // specify neutrino type:   +int : neutrino   -int: anti-neutrino 
    // specify Path length in the matter  [km]
    // specify density of the matter [g/cm^3]
    void propagate( int knu );

    // Unused in FullSME calculation 
    // Instead call DefineLinePath and run "propagte" 
    void propagateLinear( int , double, double );

    // Set neutrino propagation through the earth assuming a fixed production height 
    // and pathlength. Internally computes cosine zenith to estimate average 
    // matter density along trajectory. Lengths are in [km]
    void DefinePathFromLength(double LinKM, double ProdHeight);

    // Set pathlength and matter density for single path through matter. Good for 
    // beam experiments. l in [km] , rho in [g/cm^3]
    void DefineLinearPath( double l, double rho );

    // Same as BargerPropagator but also sets fLocalPath
    void DefinePath(double cz, double ProdHeight, bool kSetProfile = true );

    virtual double GetDensityAvg(); 

    void SetMNS( double x12, double x13, double x23,
                 double m21, double mAtm, double delta,
                 double Energy_ , bool kSquared, int kNuType );

    // Build the complete neutrino propagation hamiltonian. 
    // Must be called after SetMNS, DeFinePath*, and LV matrices are set 
    void SetHamiltonian( double Energy ) ;


    // Define LV matrices in bulk. Must be called before SetHamiltonian()
    void SetCMatrix( double reC[3][3] , double imC[3][3] );
    void SetAMatrix( double reA[3][3] , double imA[3][3] );

    // Define LV matrices entry by entry. Matrix is specified as either "A" or "C"
    // Must be called before SetHamiltonian
    void SetLVMatrixEntry( const char * matrix , int i, int j, double re , double im );
   
    void PrintLVMatrices();
    void PrintMatrix( complex<double> H[][3] ) ;

    static const double LtoeV_Coeff;

  private:

    //Standard 3f oscillations
    void SetLVUPMNSmatrix             ( double, double, double, double);

    vector< complex<double> > Calc_H0(complex<double>[][3], double, double, double[3]);
    vector< complex<double> > Calc_Hmat(double PathAvgDensity);
    complex<double>**         Calc_ABCNDeltaE(vector< complex<double> >, vector< complex<double> >,
				          vector< complex<double> >);

    vector< complex<double> > BuildLVHamiltonian( double E , bool kAnti);

    // Complex matrix manipulations 
    vector< complex<double> > AddCplxMatrices(vector< complex<double> >, vector< complex<double> >);
    vector< complex<double> > Calc_ProdMatrices(vector< complex<double> >, vector< complex<double> >);
    complex<double>           Calc_Trace(vector< complex<double> >);
    complex<double>           Calc_Determinant(vector< complex<double> >);

    // Oscillation calculations assuming the matter and LV hamiltonians, 
    // as well as MNS matrix have been set. Should be accessed as the 
    // last step of a calculation
    double Calc_Prob_ee     (complex<double>**, double);
    double Calc_Prob_emu    (complex<double>**, double);
    double Calc_Prob_etau   (complex<double>**, double);
    double Calc_Prob_mue    (complex<double>**, double);
    double Calc_Prob_mumu   (complex<double>**, double);
    double Calc_Prob_mutau  (complex<double>**, double);
    double Calc_Prob_taue   (complex<double>**, double);
    double Calc_Prob_taumu  (complex<double>**, double);
    double Calc_Prob_tautau (complex<double>**, double);


    // Average matter density 
    double fAveDensity;
    double fm31;
    double fDM2[3];

    // Internal variable for neutrino flight length. Used in ::propagate
    double fLocalPath;

    // Internal LV matrices 
    // "A" and "C"
    complex<double> fALV[3][3];
    complex<double> fCLV[3][3];
 
    // Internal PMNS matrix
    complex<double> fPMNS[3][3];

    // Internal Hamiltonians 
    vector< complex<double> > fHosc;
    vector< complex<double> > fHLV;
    vector< complex<double> > fHmat;


};

#endif
