#include "EarthDensity.h"
#include <iostream>
#include <cstdlib>
#include <sstream>

EarthDensity::EarthDensity( bool x ) 
{
   //cout << "EarthDensity::EarthDensity Using Default density profile  " << endl;

   // radius: [ km ]  density  [ g/cm^3 ]
   _density[ 0 ]       =  13.0 ;
   _density[ 1220.0 ]  =  13.0 ;
   _density[ 3480.0 ]  =  11.3 ;
   _density[ 5701.0 ]  =  5.0 ;
   _density[ 6371.0 ]  =  3.3 ;

   // chemical composition
   _YpMap[ 0 ]       =  0.468 ;
   _YpMap[ 1220.0 ]  =  0.468 ;
   _YpMap[ 3480.0 ]  =  0.497 ;
   _YpMap[ 5701.0 ]  =  0.497 ;
   _YpMap[ 6371.0 ]  =  0.497 ;

   REarth  = 6371.0 ; 

   _TraverseDistance  = NULL;
   _TraverseRhos      = NULL;
   _Yp                = NULL;

   // radius: [ km ]  DensityCoefficients { a, b, c }
   // from http://www.typnet.net/Essays/EarthGrav.htm
   // for conversion, use:
   // [kg/m3] = [(1000g)/(100cm)^3]            = 1e-3 [ g/cm3         ]
   // [kg/m4] = 1e-3 [g/cm3] * [1/(1e-3 km)]   = 1e0  [(g/cm3) (1/km) ]
   // [kg/m5] = 1e-3 [g/cm3] * [1/(1e-3 km)^2] = 1e3  [(g/cm3) (1/km2)]
   _densityCoefficients[ 0 ]       =  { -2.1773e-7,  1.9110e-8, 1.3088e1 } ; // Inner core
   _densityCoefficients[ 1220.0 ]  =  { -2.1773e-7,  1.9110e-8, 1.3088e1 } ; // Inner core
   _densityCoefficients[ 3480.0 ]  =  { -2.4123e-7,  1.3976e-4, 1.2346e1 } ; // Outer core
   _densityCoefficients[ 5701.0 ]  =  { -3.0922e-8, -2.4441e-4, 6.7823e0 } ; // Lower mantle
   _densityCoefficients[ 6371.0 ]  =  {  0.       , -1.2603e-3, 1.1249e1 } ; // Inner transition zone 2


   kUseAverageDensity = x ; 
   
   init();
       
}


EarthDensity::EarthDensity( const char * file )
{
    _TraverseDistance  = NULL;
    _TraverseRhos      = NULL;
    _Yp                = NULL;

    LoadDensityProfile( file );
}

void EarthDensity::LoadDensityProfile( const char * file )
{
        ifstream PREM_dat;
        double r_dist = 0.0;    // radial distance -- map key //
        double rho;             // density at that distance -- map value //
        double lyp;
        double REarth  = 6371.0 ; 

	DensityFileName = file;

        PREM_dat.open(DensityFileName.c_str());
        if(! PREM_dat)
        {
                cerr<<"EarthDensity::Load ERROR OPENING " << DensityFileName << endl;
                exit(1);
        }
        else
           cout << "Loading Density profile from: " << DensityFileName << endl;

        // count the number of columns to separate IO between reading 
        // files like 
        //   radius rho     and  
        //   radius rho  Yp and  
        std::string line;
        getline( PREM_dat , line);         
        std::stringstream s;
        s << line;                   
        int count = 0;
        double dummy; 
        while ( s >> dummy ) count++;

        // rewind the file 
        PREM_dat.clear();
        PREM_dat.seekg(0);
         
        while( !PREM_dat.eof( ) )
        {
                if ( r_dist > REarth ) REarth = r_dist;
            
                if( count == 2 ) 
                { 
                   PREM_dat >> r_dist >> rho ;
                   _density[r_dist] = rho ;
                   _YpMap  [r_dist] = 0.5 ;
                }
                else
                { 
                   PREM_dat >> r_dist >> rho >> lyp ;
                   _density[r_dist] = rho ;
                   _YpMap  [r_dist] = lyp ;
                }

        }
        PREM_dat.close();

        // must be re-initialized after a profile is loaded
        init();        
}


void EarthDensity::init()
{

	Load();
	ComputeMinLengthToLayers();
}


///// Really need to clean this bit up, slow and bulky!
void EarthDensity::SetDensityProfile( double CosineZ, double PathLength , double ProductionHeight)
{
   int i;
   int MaxLayer;
   double km2cm = 1.0e5;
   double TotalEarthLength =  -2.0*CosineZ*REarth*km2cm; // in [cm]
   double CrossThis, CrossNext;

   map<double, double>::iterator _i;


   // path through air
   _TraverseRhos[0] = 0.0 ;
   _Yp[0]           = 0.0 ;
   _TraverseDistance[0] =  PathLength - TotalEarthLength ;

// std::cout << " Earth ... PathLenght: " << PathLength << " totallength " << TotalEarthLength << std::endl;

   if( CosineZ >= 0 )
   {  
       _TraverseDistance[0] =  PathLength;
       Layers = 1; 
       return; 
   }

   	
   Layers = 0;
   for ( _i = _CosLimit.begin(); _i != _CosLimit.end() ; _i++ )
      if( CosineZ < _i->second )
          Layers++;	


   MaxLayer = Layers;

   // the zeroth layer is the air!
   for ( i = 0 ; i< MaxLayer ; i++ ) 
   {
 
      _TraverseRhos[i+1]      = _Rhos[i];
      _Yp          [i+1]      = _Yps [i];

//   CrossThis = 2.0*sqrt( _Radii[i]   * _Radii[i]    - REarth*REarth*( 1 -CosineZ*CosineZ ) );
//   if( i < MaxLayer-1 )
//   {
//    CrossNext = 2.0*sqrt( _Radii[i+1] * _Radii[i+1]  - REarth*REarth*( 1 -CosineZ*CosineZ ) );
//    _TraverseDistance[i+1]  =  0.5*( CrossThis-CrossNext )*km2cm;
//   }
//   else
//    _TraverseDistance[i+1]  =  CrossThis*km2cm;

     double Rmax2 = _Radii[i]*_Radii[i];
     double Rmin2 = REarth*REarth*(1.-CosineZ*CosineZ);
     CrossThis = sqrt( Rmax2 - Rmin2 );

     if ( i < MaxLayer-1 ) {
        CrossNext = sqrt( _Radii[i+1]*_Radii[i+1] - Rmin2 );
     }
     else {
          CrossNext = 0.;
     }

     // calculate average density
     if ( kUseAverageDensity && _Rhos_a.size() > 0)
     {
        double t1 = CrossThis-CrossNext;
        double t0 = 0.;
        double R_t1 = _Radii[i+1];
        double R_t0 = _Radii[i  ];
        if (i == MaxLayer-1) {
            R_t1 = sqrt(Rmin2);
        }
        // R^2 = Rmin^2 + (sqrt(Rmax^2 - Rmin^2) - (t-t0))^2
        //     = A (t-t0)^2 + B (t-t0) + C
        double RA = 1.;
        double RB = -2.*CrossThis;
        double RC = Rmax2;
        double ExpR2 = (RA*pow(t1,3)/3. + RB*pow(t1,2)/2. + RC*t1) / t1; // here we used t0=0
        double ExpR1ss1 = 2.*sqrt(RA)*R_t1;
        double ExpR1ss0 = 2.*sqrt(RA)*R_t0;
        double ExpR1bb1 = 2.*RA*t1 + RB;
        double ExpR1bb0 = 2.*RA*t0 + RB;
        double ExpR1det = RB*RB - 4.*RA*RC;
        double ExpR1 = 1./(8.*sqrt(RA)*(t1-t0)) *
                       ((ExpR1ss1*ExpR1bb1 - ExpR1det*log(ExpR1ss1+ExpR1bb1))
                       -(ExpR1ss0*ExpR1bb0 - ExpR1det*log(ExpR1ss0+ExpR1bb0)));
        double ExpR0 = 1.;
        _TraverseRhos[i+1] = _Rhos_a[i]*ExpR2 + _Rhos_b[i]*ExpR1 + _Rhos_c[i]*ExpR0;
     }

     _TraverseDistance[i+1]  =  (CrossThis-CrossNext)*km2cm;
     if (i == MaxLayer-1) {
         // for center we have two segments
         _TraverseDistance[i+1] *= 2.;
     }
     // end computations for average density






    
     // assumes azimuthal symmetry    
     if( i < MaxLayer ){
        _TraverseRhos    [ 2*MaxLayer - i ] = _TraverseRhos[i];
        _TraverseDistance[ 2*MaxLayer - i ] = _TraverseDistance[i];
        _Yp              [ 2*MaxLayer - i ] = _Yp[i]; 
     }
   }

   Layers = 2*MaxLayer; 

}


// now using Zenith angle to compute minimum conditions...20050620 rvw
void EarthDensity::ComputeMinLengthToLayers()
{
	double x;

	_CosLimit.clear();

	// first element of _Radii is largest radius!
	for(int i=0; i < (int) _Radii.size() ; i++ )
	{
            // Using a cosine threshold instead! //
            x = -1* sqrt( 1 - (_Radii[i] * _Radii[i] / ( REarth*REarth)) );
            if ( i  == 0 ) x = 0;
                _CosLimit[ _Radii[i] ] = x;
	}

}



void EarthDensity::Load( )
{

      int MaxDepth = 0;

      //map<double, double>::iterator _i;
      map<double, double>::reverse_iterator _i;
      map<double, DensityCoefficients>::reverse_iterator _ic;


      if( _TraverseRhos     != NULL ) delete [] _TraverseRhos;
      if( _TraverseDistance != NULL ) delete [] _TraverseDistance;
      if( _Yp               != NULL ) delete [] _Yp; 

      //_i = _density.end();
      //_i--;

      _i = _density.rbegin();
      REarth = _i->first;  // Earth is 6371.0 [km]

      // to get the densities in order of decreasing radii
      for( _i = _density.rbegin() ; _i != _density.rend() ; ++_i )
      {
	 _Rhos .push_back( _i->second );
   	 _Radii.push_back( _i->first  );
	 MaxDepth++;
      } 

      for( _i = _YpMap.rbegin() ; _i != _YpMap.rend() ; ++_i )
   	 _Yps  .push_back( _i->second  );


      for (_ic = _densityCoefficients.rbegin(); _ic != _densityCoefficients.rend(); ++_ic)
      {
         _Rhos_a.push_back(_ic->second.a);
         _Rhos_b.push_back(_ic->second.b);
         _Rhos_c.push_back(_ic->second.c);
      }

      _TraverseRhos      = new double [ 2*MaxDepth + 1 ];
      _TraverseDistance  = new double [ 2*MaxDepth + 1 ];
      _Yp                = new double [ 2*MaxDepth + 1 ];

      return;
                                                                                                                                                             
}


EarthDensity::~EarthDensity( )
{
   if( _TraverseRhos     != NULL ) delete [] _TraverseRhos;
   if( _TraverseDistance != NULL ) delete [] _TraverseDistance;
   if( _Yp               != NULL ) delete [] _Yp; 

}


