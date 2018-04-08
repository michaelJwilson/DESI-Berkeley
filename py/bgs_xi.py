import  matplotlib;                         matplotlib.use('PDF')

from    nbodykit.lab                import  *
from    nbodykit                    import  style, setup_logging
from    astropy.cosmology           import  FlatLambdaCDM
from    astropy.cosmology           import  Planck15              as pycosmo
from    astropy.io                  import  fits
from    sys                         import  argv
from    logmocks                    import  assign_randzs, getopts
from    math                        import  pi
from    scipy.constants             import  speed_of_light        as lightspeed
from    Corrfunc.mocks              import  DDrppi_mocks

import  matplotlib.pyplot           as      plt
import  pylab                       as      pl
import  numpy                       as      np
import  healpy                      as      hp
import  dask.array                  as      da
import  dask
import  os

setup_logging()


args        = getopts(argv)

try:
    seed    = np.int(args['-i'])

except:
    seed    = 1                                     ## Default argument

redshift    = 0.3

home_dir    = os.environ['HOME']
scratch_dir = os.environ['SCRATCH']

mock_output = 'Scratch'
output_dirs = {'Home': home_dir, 'Scratch': scratch_dir}

midchi     = pycosmo.comoving_distance([redshift])  ##  Mpc.                                                                                         
midchi    *= pycosmo.h                              ##  Mpc / h.                                                                                           
                                     
midchi     = midchi.value

cosmo      = cosmology.Planck15
Plin       = cosmology.LinearPower(cosmo, redshift, transfer='EisensteinHu')


b1         = 1.0e+0

res        = 'hi'
chunks     =  100
interlaced =  False

lores      = {'nbar': 1.0e-4, 'boxsize': 1.0e+3, 'fftsize': 128}
hires      = {'nbar': 1.0e-4, 'boxsize': 1.5e+3, 'fftsize': 128}

res_opts   = {'lo': lores, 'hi': hires}

nbar       = res_opts[res]['nbar']
boxsize    = res_opts[res]['boxsize']
fftsize    = res_opts[res]['fftsize']

if __name__ == "__main__":
    rmin             =  1.0
    rmax             = 20.0
    nbins            =   20

    r                = np.linspace(rmin, rmax, nbins + 1)
    logr             = np.logspace(np.log10(rmin), np.log10(rmax), nbins + 1)

    '''
    data, randoms    = assign_randzs(ztype = "LZEE", num = seed)
    data['Weight']   = 1.0

    data['Position'] = transform.SkyToCartesian(data['RA'], data['DEC'], data['RZEE'], cosmo, degrees=True)
    
    print('data columns = ',       data.columns)
    print('randoms columns = ', randoms.columns)
    
    print data['RA']   ##.compute()
    print data['DEC']  ##.compute()
    print data['LZEE'] ##.compute()
    
    ra      = data['RA'].compute().astype(np.float64) 
    dec     = data['DEC'].compute().astype(np.float64)

    czee    = data['LZEE'].compute().astype(np.float64)
    czee   *= lightspeed

    weight  = data['Weight'].compute().astype(np.float64)
    '''
    '''
    ngal    = 5000
    
    ra      = 360. *  np.random.uniform(low=0.0, high=1.0, size=ngal)
    dec     = 180. * (np.random.uniform(low=0.0, high=1.0, size=ngal) - 0.5)  

    zee     =         np.random.uniform(low=0.0, high=1.0, size=ngal)
    czee    =  zee *  lightspeed  

    weight  =         np.random.uniform(low=0.0, high=1.0, size=ngal)

    print  ra
    print  dec
    print  czee
    print  weight

    print  logr

    ## Valid cosmology (integer) are 1 -> LasDamas and 2 -> Planck cosmology.
    results          =   DDrppi_mocks(autocorr=1,                 cosmology=1,                       nthreads=1,\
                                      binfile=logr,               RA1=ra,                            DEC1=dec,\
                                      CZ1=czee,                   weights1=weight,                   is_comoving_dist=False,\
                                      output_rpavg=False,         fast_divide=False,                 pimax=20,\
                                      xbin_refine_factor=2,       ybin_refine_factor=2,              zbin_refine_factor=1,\
                                      max_cells_per_dim=100,      verbose=True,\
                                      c_api_timer=False,          isa=u'fastest',                    weight_type="pair_product")

    
    print("\n\nDone.\n\n")


