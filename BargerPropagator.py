#
# This is a simple interface to the BargerPropagator class 
# There is no fanciness, and blatant disregard for the
# use of default arguments on the C++ side
#


import ctypes as c 
lib = c.cdll.LoadLibrary('libThreeProb.so')
#lib = c.cdll.LoadLibrary('/disk/usr3/raw/Development/atmpd-trunk/src/analysis/Prob3++/libThreeProb.so')

# constructor (This returns a pointer. This is necessary for 64 bit machines.
# The default return type is a 32 bin integer.)
lib.propagator_new.restype = c.c_void_p

# set_mns
lib.set_mns.argtypes = [c.c_void_p,
                        c.c_double, c.c_double , c.c_double , c.c_double , 
                        c.c_double, c.c_double , c.c_double ,             
                        c.c_bool  , c.c_int   
                       ]

# propagate
lib.propagate.argtypes = [c.c_void_p, c.c_int   ]

# propagate_linear
lib.propagate_linear.argtypes = [c.c_void_p, c.c_int   , c.c_double , c.c_double ]

# get_prob
lib.get_prob.restype  =  c.c_double 
lib.get_prob.argtypes = [c.c_void_p , c.c_int , c.c_int ]

# define_path
lib.define_path.argtypes = [c.c_void_p , c.c_double , c.c_double , c.c_bool ]

# use_mass_eigenstates 
lib.use_mass_eigenstates.argtypes = [c.c_void_p , c.c_bool ]

# set_warning_suppression 
lib.set_warning_suppression.argtypes = [c.c_void_p , c.c_bool ]


# set_one_mass_scale_mode 
lib.set_one_mass_scale_mode.argtypes = [c.c_void_p , c.c_bool ]


# set_density_conversion 
lib.set_density_conversion.argtypes = [c.c_void_p , c.c_double ]

# get_path_length 
lib.get_path_length.restype  =  c.c_double 
lib.get_path_length.argtypes = [c.c_void_p ]

# set_path_length 
lib.set_path_length.argtypes = [c.c_void_p, c.c_double ]

# set_energy 
lib.set_energy.argtypes = [c.c_void_p, c.c_double ]

# set_air_path_length 
lib.set_air_path_length.argtypes = [c.c_void_p, c.c_double ]

# get_vacuum_prob
lib.get_vacuum_prob.restype  =  c.c_double 
lib.get_vacuum_prob.argtypes = [c.c_void_p , c.c_int , c.c_int , c.c_double, c.c_double ]


# set_default_octant
lib.set_default_octant.argtypes = [c.c_void_p, c.c_int , c.c_int  ]

class BargerPropagator(object):

    def __init__(self):
        self.obj = lib.propagator_new()

    def SetMNS(self, x12, x13, x23, m21, mAtm, delta, energy, ks, nutype):
        lib.set_mns( self.obj, x12, x13, x23, m21, mAtm, delta, energy, ks, nutype )

    def DefinePath ( self , cz, ProdHeight, kSetProfile ): 
         lib.define_path ( self.obj , cz, ProdHeight, kSetProfile ) 

    def propagate( self, nu ):
        lib.propagate( self.obj , nu ) 

    def propagateLinear( self, nu, path, density ):
        lib.propagate_linear( self.obj , nu , path, density ) 
 
    def GetProb( self , nuIn, nuOut ):
        return lib.get_prob( self.obj , nuIn, nuOut )

    def UseMassEigenstates( self, x ):
        lib.use_mass_eigenstates(self.obj , x )
 
    def SetWarningSuppression( self, x ):
        lib.set_warning_suppression(self.obj , x )

    def SetOneMassScaleMode( self, x ):
        lib.set_one_mass_scale_mode(self.obj, x )

    def SetDensityConversion( self, x):
        lib.set_density_conversion( self.obj , x )

    # The following functions are for expert use
    # please be careful and check BargerPropagator.h 
    # for more information
    def GetPathLength( self  ):
        return lib.get_path_length( self.obj )

    def SetPathLength( self, x ):
        lib.set_path_length( self.obj , x )

    def SetEnergy( self, x ):
        lib.set_energy( self.obj , x )

    def SetAirPathLength( self, x ):
        lib.set_air_path_length( self.obj , x )

    def GetVacuumProb( self, nuin, nuout, energy , path ):
        return lib.get_vacuum_prob( self.obj, nuin, nuout, energy, path )

    def SetDefaultOctant( self, x , y ):
        lib.set_default_octant( self.obj , x , y)

