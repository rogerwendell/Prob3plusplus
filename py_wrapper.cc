#include "BargerPropagator.h"

// 
// Simple ctypes-inspired wrappers for 
// BargerPropagator (probably the only thing necessary
// for typical package usage)
//


extern "C"
{

  BargerPropagator * propagator_new() 
  { 
     return new BargerPropagator(); 
  }
 
  void  set_mns( BargerPropagator * p ,
                               double x12, double x13, double x23,
                               double m21, double mAtm, double delta,
                               double Energy_ , bool kSquared, int kNuType )
  {
     p->SetMNS( x12, x13, x23, m21, mAtm, delta, Energy_, kSquared, kNuType );
  }

  void  define_path ( BargerPropagator * p ,  
                      double cz, double ProdHeight, bool kSetProfile ) 
  {
    p->DefinePath( cz, ProdHeight, kSetProfile );
  }

  void propagate    ( BargerPropagator * p ,  int NuFlavor )
  {
    p->propagate( NuFlavor );
  }

  void propagate_linear ( BargerPropagator * p ,  int NuFlavor, double pathlength, double Density )
  {
    p->propagateLinear( NuFlavor, pathlength, Density );
  }

  
  double get_prob( BargerPropagator * p , int nuIn, int nuOut )
  {
    return p->GetProb( nuIn , nuOut );
  }
  
  //////
  // Others
 
  void use_mass_eigenstates( BargerPropagator * p , bool x )
  {
    p->UseMassEigenstates(x);
  }
  
  
  void set_warning_suppression( BargerPropagator * p , bool x )
  {
    p->SetWarningSuppression(x);
  }
  
  void set_one_mass_scale_mode( BargerPropagator * p , bool x )
  {
    p->SetOneMassScaleMode(x);
  }

  
  void set_density_conversion( BargerPropagator * p , double x )
  {
    p->SetDensityConversion( x );
  }
  
  double get_path_length( BargerPropagator * p )
  {
    return p->GetPathLength();
  }

  
  void set_path_length( BargerPropagator * p, double x )
  {
    p->SetPathLength(x);
  }
  
  void set_energy( BargerPropagator * p, double x )
  {
    p->SetEnergy(x);
  }
  
  
  void set_air_path_length( BargerPropagator * p, double x )
  {
    p->SetAirPathLength(x);
  }


  double get_vacuum_prob( BargerPropagator * p , int nuIn, int nuOut, double energy, double path )
  {
    return p->GetVacuumProb( nuIn , nuOut , energy , path );
  }


  void set_default_octant( BargerPropagator * p , int var , int octant)
  {
    return p->SetDefaultOctant( var , octant ) ;
  }
} 
