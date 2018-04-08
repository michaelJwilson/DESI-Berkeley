import  matplotlib;       matplotlib.use('PDF')

from    nbodykit.lab      import *
from    nbodykit          import setup_logging, style

import  os
import  ssl
import  matplotlib.pyplot as plt
import  pylab             as pl
import  requests
import  warnings


plt.style.use(style.notebook)

setup_logging()  ## turn on logging to screen


warnings.simplefilter(action='ignore', category=FutureWarning)

if __name__ == "__main__":
    data_path     = './DR12/galaxy_DR12v5_LOWZ_South.fits'
    randoms_path  = './DR12/random0_DR12v5_LOWZ_South.fits'

    data          = FITSCatalog(data_path)
    randoms       = FITSCatalog(randoms_path)

    print('data columns = ',    data.columns)
    print('randoms columns = ', randoms.columns)

    ZMIN = 0.15
    ZMAX = 0.43

    ## slice the randoms
    valid   = (randoms['Z'] > ZMIN) & (randoms['Z'] < ZMAX)
    randoms =  randoms[valid]

    ## slice the data
    valid   = (data['Z'] > ZMIN) & (data['Z'] < ZMAX)
    data    =  data[valid]

    ## Alam et al. cosmology. 
    cosmo   = cosmology.Cosmology(h=0.676).match(Omega0_m=0.31)

    data['Position']    = transform.SkyToCartesian(data['RA'], data['DEC'], data['Z'], cosmo=cosmo)
    randoms['Position'] = transform.SkyToCartesian(randoms['RA'], randoms['DEC'], randoms['Z'], cosmo=cosmo)

    randoms['WEIGHT']   = 1.0
    data['WEIGHT']      = data['WEIGHT_SYSTOT'] * (data['WEIGHT_NOZ'] + data['WEIGHT_CP'] - 1.0)
    
    ## Combine the data and randoms into a single catalog
    fkp  = FKPCatalog(data, randoms)

    ## Extent of the randoms sets the box size; There is also a BoxSize keyword. 
    mesh = fkp.to_mesh(Nmesh=256, nbar='NZ', fkp_weight='WEIGHT_FKP', comp_weight='WEIGHT', window='tsc')

    r    = ConvolvedFFTPower(mesh, poles=[0,2,4], dk=0.005, kmin=0.)

    for key in r.attrs:
        print("%s = %s" % (key, str(r.attrs[key])))

    poles = r.poles

    for ell in [0, 2, 4]:
        label = r'$\ell=%d$' % (ell)
        P     = poles['power_%d' %ell].real
        
        if ell == 0: 
            P = P - r.attrs['shotnoise']

        plt.plot(poles['k'], poles['k']*P, label=label)
        
        plt.legend(loc=0)
        plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
        plt.ylabel(r"$k \ P_\ell$ [$h^{-2} \ \mathrm{Mpc}^2$]")
        plt.xlim(0.01, 0.25)

    ## Using interlacing (by setting interlaced=True) can also reduce the effects of aliasing on the measured results.
    pl.savefig('DR12_multipoles.pdf')

    print("\n\nDone.")
