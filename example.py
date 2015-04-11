#!/usr/bin/env python

import  scipy
from scipy import interpolate

import  qnd

class SmokePoint ( qnd.Point ) :
    
    # Computing the kinetic energy requires the progenitor white dwarf
    # pre-explosion internal energy and gravitational binding energy.  
    # These figures come from Cody Raskin (personal communication), but
    # omit radii he provided (2014-01-14).  The last value listed comes
    # from WKBS 2007 which also provided the equation
    #
    #   E_K = 1.56 M_Ni + 1.74 M_Fe + 1.24 M_IME - E_g + E_int
    #
    # giving kinetic energy in terms of ejecta products and white dwarf
    # energies.  White dwarf energies are determined by interpolating 
    # Cody's table.

    wd_data = scipy.array( [
        #   mass (msun)     E_int (erg)     E_g (erg)
        [   0.39,           1.80e+49,       3.52e+49    ],
        [   0.44,           2.45e+49,       4.74e+49    ],
        [   0.49,           3.26e+49,       6.21e+49    ],
        [   0.54,           4.25e+49,       7.97e+49    ],
        [   0.59,           5.44e+49,       1.00e+50    ],
        [   0.64,           6.87e+49,       1.25e+50    ],
        [   0.69,           8.59e+49,       1.54e+50    ],
        [   0.74,           1.07e+50,       1.88e+50    ],
        [   0.79,           1.31e+50,       2.27e+50    ],
        [   0.84,           1.61e+50,       2.74e+50    ],
        [   0.89,           1.98e+50,       3.29e+50    ],
        [   0.94,           2.42e+50,       3.94e+50    ],
        [   0.99,           2.97e+50,       4.73e+50    ],
        [   1.04,           3.66e+50,       5.68e+50    ],
        [   1.08,           4.55e+50,       6.87e+50    ],
        [   1.13,           5.72e+50,       8.38e+50    ],
        [   1.18,           7.32e+50,       1.04e+51    ],
        [   1.23,           9.66e+50,       1.32e+51    ],
        [   1.28,           1.34e+51,       1.76e+51    ],
        [   1.33,           2.09e+51,       2.58e+51    ],
        [   1.38,           2.89e+51,       3.35e+51    ] ] )

    wd_internal_energy_func  = interpolate.interp1d( wd_data[ :, 0 ], wd_data[ :, 1 ] / 1.0e51 )
    wd_binding_energy_func   = interpolate.interp1d( wd_data[ :, 0 ], wd_data[ :, 2 ] / 1.0e51 )

    @property
    def wd_internal_energy( self ) :
        """Progenitor white dwarf's internal energy (Bethes)."""
        try :
            return self._wd_internal_energy
        except AttributeError :
            self._wd_internal_energy = self.wd_internal_energy_func( self.ejecta_mass )
            return self._wd_internal_energy

    @property
    def wd_binding_energy( self ) :
        """Progenitor white dwarf's binding energy (Bethes)."""
        try :
            return self._wd_binding_energy
        except AttributeError :
            self._wd_binding_energy  = self.wd_binding_energy_func( self.ejecta_mass )
            return self._wd_binding_energy

    def __init__( self, nickel_mass, ime_mass, unburned_mass, nickel_radius, opacity ) :
        self.nickel_mass    = nickel_mass
        self.ime_mass       = ime_mass
        self.unburned_mass  = unburned_mass
        self.nickel_radius  = nickel_radius
        self.opacity        = opacity

    @property
    def ejecta_mass( self ) :
        """Total ejecta mass (solar masses)."""
        try :
            return self._ejecta_mass
        except AttributeError :
            self._ejecta_mass = self.nickel_mass + self.ime_mass + self.unburned_mass
            return self._ejecta_mass
        
    @property
    def kinetic_energy( self ) :
        """Kinetic energy (Bethes)."""
        try :
            return self._kinetic_energy
        except AttributeError :
            self._kinetic_energy  = 1.56 * self.nickel_mass + 1.24 * self.ime_mass
            self._kinetic_energy -= self.wd_binding_energy
            self._kinetic_energy += self.wd_internal_energy
            return self._kinetic_energy

    @property
    def velocity( self ) :
        """Characteristic velocity (10^3 km/s)."""
        try :
            return self._velocity
        except AttributeError :
            self._velocity = scipy.sqrt( 2.0e51 * self.kinetic_energy / self.ejecta_mass / 1.9891e33 ) / 1.0e8
            return self._velocity


    def _is_compliant( self ) :
        """Is this point okay, does it violate our custom constraints?"""

        if self.ejecta_mass < 0.6 :
            return False
        if self.ejecta_mass > 1.38 :
            return False
        if self.nickel_radius < self.nickel_mass:
            return False
        if self.nickel_radius > self.ejecta_mass:
            return False
        if self.velocity < 4.0:
            return False
            
        return True

class SedonaPoint(SmokePoint):
    
    def __init__(self, nickel_mass, ime_mass, unburned_mass, nickel_radius):
        self.nickel_mass    = nickel_mass
        self.ime_mass       = ime_mass
        self.unburned_mass  = unburned_mass
        self.nickel_radius  = nickel_radius


if __name__ == "__main__" :

    lower = scipy.array( [ 0.20, 0.20, 0.00, 0.00 ] )
    upper = scipy.array( [ 1.38, 1.00, 0.10, 1.38 ] )

    design = qnd.random_design( SedonaPoint, lower, upper )
    print design.unscaled
