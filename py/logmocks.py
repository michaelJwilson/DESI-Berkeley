import  matplotlib;                matplotlib.use('PDF')
 
from    nbodykit.lab       import  *
from    nbodykit           import  style, setup_logging
from    astropy.cosmology  import  FlatLambdaCDM
from    astropy.cosmology  import  Planck15              as pycosmo
from    astropy.io         import  fits
from    sys                import  argv
from    nbodykit           import  algorithms

import  matplotlib.pyplot  as      plt
import  pylab              as      pl
import  numpy              as      np
import  healpy             as      hp
import  dask.array         as      da
import  dask 
import  os

setup_logging()


def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.

    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.

        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.

    return opts

def null_power(k):
    return  1e-10

args       = getopts(argv)

try:
    seed   = np.int(args['-i'])    
except:
    seed   = 42                    ## Default argument  

## How to: collect complete array from distribution amongst allocated cores.
## data    = numpy.concatenate(catalog.comm.allgather(catalog['Position'].compute()), axis=0)

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

print("\n\nRedshift: %.3lf, b1: %.3lf, seed: %d" % (redshift, b1, seed))
print("Resolution params: %s-res, nbar: %.3le, boxsize: %.3lf, fft size: %d" % (res, nbar, boxsize, fftsize))


def create_catalog(plot_cubemultipoles = False, num=seed):
  ## Generate lognormal cube with observer at the centre.  Add real, global and local pp redshift space and a flag for in / out of BGS.
  
  cat      = LogNormalCatalog(Plin=Plin,       nbar=nbar, BoxSize=boxsize, Nmesh=fftsize, bias=b1,  seed = num)
  ## rand  = LogNormalCatalog(Plin=null_power, nbar=nbar, BoxSize=boxsize, Nmesh=fftsize, bias=1.0, seed = 142, cosmo=cosmo, redshift=redshift)

  print("Created lognormal cube.")
  
  ## Possibility of a uniform cat. 
  ## uniform          = UniformCatalog(nbar=150, BoxSize=1.0, seed=42)

  ## Place (0,0,0) at the centre.  Note: Dask array calls are placed on a list of commands, all evaluated on .compute()                                                                           
  cat['Position']    -= boxsize / 2.

  ## Place in redshift-space according to the global plane-parallel approx.                                                                               
  RSD_prefactor       = (1. + redshift) / (100. * cosmo.efunc(redshift))
    
  cat['norm']         = da.sqrt(da.sum(cat['Position']**2., axis=1))

  GLOS                = [0, 0, 1]            ## Assumes a Kaiser pp approx. along the z-axis                                                             
  
  ## Unit vector in the direction of the galaxy by broadcasting.                                                                                                             
  cat['LLOS']         = cat['Position'] / cat['norm'][:, None]

  cat['GLOS_pos']     = cat['Position'] + RSD_prefactor * cat['Velocity'] * GLOS
  cat['LLOS_pos']     = cat['Position'] + RSD_prefactor * cat['Velocity'] * cat['LLOS']

  ## RA and DEC will be returned in degrees, with RA in the range [0,360] and DEC in the range [-90, 90]
  ## Could just add dz explicitly.
  (cat['ra'], cat['dec'], cat['RZEE']) = transform.CartesianToSky(cat['Position'], cosmo, velocity=None, observer = [0, 0, 0], zmax = 100.0)

  ## Now apply a zee cut according to local line-of-sight.                                                                                                  
  ZMIN             = 0.05
  ZMAX             = 0.20

  valid            = (cat["RZEE"] > ZMIN) & (cat["RZEE"] < ZMAX)
  cat              =  cat[valid]

  print("Applied %.3lf < (real-space) z < %.3lf cut to catalogue." % (ZMIN, ZMAX))

  (cat['ra'], cat['dec'], cat['GZEE']) = transform.CartesianToSky(cat['GLOS_pos'], cosmo, velocity=None, observer = [0, 0, 0], zmax = 100.0)
  (cat['ra'], cat['dec'], cat['LZEE']) = transform.CartesianToSky(cat['LLOS_pos'], cosmo, velocity=None, observer = [0, 0, 0], zmax = 100.0)

  print("Created redshifts for real, global and local pp redshift space.")
  
  if plot_cubemultipoles is True:
    pl.clf()

    ## Convert the catalog to the mesh, with CIC interpolation                                                                                               
    mesh            = cat.to_mesh(compensated=True, window='cic', position='Position', interlaced=interlaced)                                                     
    r               = FFTPower(mesh, mode='2d', dk=0.005, kmin=0.01, Nmu=20, los = GLOS, poles=[0, 2, 4])                                            
    poles           = r.poles                                                                                                                               

    for ell in [0, 2]:                                                                                                                                    
        label  = r'$\ell=%d$'     % ell                                                                                                                    
        P      = poles['power_%d' % ell].real                                                                                                                

        if ell == 0:                                                                                                                                        
            P  = P - poles.attrs['shotnoise']                                                                                                               

    pl.loglog(poles['k'], P, label=label)                                                                                                                   
    
    ## Plot real-space power.                                                                                                                           
    k = np.logspace(-2, 0, 512)                                                                                                                         
    plt.loglog(k, b1**2 * Plin(k), c='k', label=r'$b_1^2 P_\mathrm{lin}$')                                                                                  

    plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")                                                                                                          
    plt.ylabel(r"$P_\ell(k)$ [$(h^{-1} \ \mathrm{Mpc})^3$]")                                                                                                

    plt.xlim(0.01, 0.30)                                                                                                                                 
    plt.ylim( 1e3,  1e5)                                                                                                                                    

    pl.savefig('lognormalcube_multipoles.pdf')
  
    print("Plotted lognormal cube multipoles.")

  ## DESI BGS Healpix nside defined footprint; hardcoded nside.  
  nside            = 128
  imap             = np.loadtxt("desibgs_imap_nside_%d.txt" % nside)

  cat['theta']     = 0.5 * np.pi - da.deg2rad(cat['dec'])
  cat['phi']       =               da.deg2rad(cat['ra'])

  ipix             = hp.pixelfunc.ang2pix(nside, cat['theta'].compute(), cat['phi'].compute(), nest=False, lonlat=False)
  cat['in_bgs']    = imap[ipix]

  ngal             = cat['in_bgs'].shape[0]
  bgs_ngal         = da.sum(cat['in_bgs']).compute()

  ## Accepted by BGS mask.
  print("Number of galaxies created: %d, Number in BGS: %d, Percentage accepted: %.3lf" % (ngal, bgs_ngal, 100. * bgs_ngal / ngal))

  ## Cut by in BGS.
  valid            = cat['in_bgs'] > 0.0
  cat              = cat[valid]

  print("Applied BGS footprint cut to catalogue.")
  
  print("Writing lognormal BGS mock to fits.")
  
  ## Create .fits
  ra               = fits.Column(name='RA',     format='D', array= cat['ra'].compute())
  dec              = fits.Column(name='DEC',    format='D', array=cat['dec'].compute())

  rzee             = fits.Column(name='RZEE',   format='D', array=cat['RZEE'].compute())
  gzee             = fits.Column(name='GZEE',   format='D', array=cat['GZEE'].compute())
  lzee             = fits.Column(name='LZEE',   format='D', array=cat['LZEE'].compute())

  cols             = fits.ColDefs([ra, dec, rzee, gzee, lzee])
  hdr              = fits.Header()

  hdr['Creator']   = 'M. J. Wilson'
  hdr['COMMENT']   = "Lognormal BGS mocks in real, global and local pp redshift space."
    
  hdu              = fits.BinTableHDU.from_columns(cols, header=hdr)

  print(num)

  hdu.writeto(output_dirs[mock_output] + '/desi/logmocks/lognormal_bgs_seed-%03d.fits' % num, overwrite=True)
  

def assign_randzs(ztype = "LZEE", num = seed):
    ## Add named columns to Martin's randoms.fits, with redshifts drawn from the data.  
    dfname          =  output_dirs[mock_output] + '/desi/logmocks/lognormal_bgs_seed-%03d.fits' % num
    rfname          =  "/global/homes/m/mjwilson/desi/randoms/randoms.fits"
    
    print("Loading:  ", dfname, rfname)
    
    data            =  FITSCatalog(dfname)
    rand            =  FITSCatalog(rfname)  ## Martin's DESI / BGS randoms.
    
    ngal            =  len(data)

    ## 20 x randoms as galaxies.                                                                                                                           
    rand            =  rand.gslice(0, 20 * ngal, redistribute=False)
    nrand           =  len(rand)

    ncopy           =  np.int(np.floor(1.0 * nrand / ngal))
    
    ## Damp removal of intrinsic radial structure.                                                                                               
    data['blur']    =  0.05 * (da.random.uniform(low=0.0, high=1.0, size=data['GZEE'].shape, chunks=chunks) - 0.5) 

    print("Calculating randoms redshifts for z type: %s" % ztype)
    
    shuf            =  np.arange(ngal)
    
    np.random.shuffle(shuf)                           ## Make sure there's no clustering in the redshift assignment
                                                      ## Would be a problem if randoms are ordered on the sky. 
    
    ##  Check that blurred redshifts are positive.  
    array           =  da.tile(data[ztype][shuf] + data['blur'], ncopy)  

    rand[ztype]     =  da.from_array(array, chunks = chunks)  ## da.random.choice(data[ztype] + data['blur'], size = rand['RA'].shape, chunks=chunks) 

    return data, rand


def multipoles(DR12=False, ztype="LZEE", num = seed):
    if DR12 is True:
        print("Loading BOSS DR12.")

        ## BOSS DR12 FITS for comparison.                                                                                                         
        data    = FITSCatalog("/global/homes/m/mjwilson/cov/DR12/data/galaxy_DR12v5_LOWZ_South.fits")                                             
        randoms = FITSCatalog("/global/homes/m/mjwilson/cov/DR12/data/random0_DR12v5_LOWZ_South.fits")                                                

        ZMIN    = 0.15
        ZMAX    = 0.43

        # slice the randoms
        valid   = (randoms['Z'] > ZMIN)&(randoms['Z'] < ZMAX)
        randoms = randoms[valid]

        # slice the data
        valid   = (data['Z'] > ZMIN)&(data['Z'] < ZMAX)
        data    = data[valid]

        data['WEIGHT']        = data['WEIGHT_SYSTOT'] * (data['WEIGHT_NOZ'] + data['WEIGHT_CP'] - 1.0)
        randoms['WEIGHT']     = 1.0
        
        randoms['WEIGHT_FKP'] = 1.0
        data['WEIGHT_FKP']    = 1.0

        ztype                 = "Z"

    else:
        print("Loading BGS lognormal mock with redshift type:  %s." % ztype)
                                                                                                                                 
        data, randoms  = assign_randzs(ztype = ztype, num = num)

        data['NZ']     = nbar
        randoms['NZ']  = nbar 

        ## Add weights;                                                                                                                                    
        ## Next, we specify the completeness weights. By construction, there are no systematic variations in the number density of , so the completenesss   
        ## weights are set to unity for all objects. For  catalog, the completeness weights are computed as defined in eq. 48 of Reid et al. 2016. These 
        ## weights account for systematic issues, redshift failures, and missing objects due to close pair collisions on the fiber plate.

        randoms['WEIGHT']     = 1.0
        data['WEIGHT']        = 1.0

        randoms['WEIGHT_FKP'] = 1.0
        data['WEIGHT_FKP']    = 1.0

    print('data columns = ',       data.columns)
    print('randoms columns = ', randoms.columns)

    data['Position']      = transform.SkyToCartesian(   data['RA'],    data['DEC'],    data[ztype], cosmo=cosmo)
    randoms['Position']   = transform.SkyToCartesian(randoms['RA'], randoms['DEC'], randoms[ztype], cosmo=cosmo)

    fkp                   = FKPCatalog(data, randoms)
    
    print("Writing BGS data / mock to mesh.")

    mesh                  = fkp.to_mesh(Nmesh=fftsize, nbar='NZ', fkp_weight='WEIGHT_FKP', comp_weight='WEIGHT', window='tsc', BoxSize=boxsize, \
                                        interlaced=interlaced)

    print("Calculating multipoles.")

    rpoles                 = ConvolvedFFTPower(mesh, poles=[0, 2], dk=0.005, kmin=0.005, use_fkp_weights=False)

    for key in rpoles.attrs:
        print("%s = %s" % (key, str(rpoles.attrs[key])))

    return rpoles

def plot_multipoles(rpoles, ztype = "LZEE"):
    poles            = rpoles.poles

    for ell in [0]:
        label        =      r'$\ell=%d$' % ell
        P            = poles['power_%d'  % ell].real

        if ell == 0: 
            P        = P - rpoles.attrs['shotnoise']

        plt.loglog(poles['k'], np.abs(P), label=label)

    ## Plot real-space power.                                                                                           
    k = np.logspace(-2, 0, 512)                         

    plt.loglog(k, b1**2 * Plin(k), c='k', label=r'$b_1^2 P_\mathrm{lin}$') 

    plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
    plt.ylabel(r"$P_\ell(k)$ [$(h^{-1} \ \mathrm{Mpc})^3$]")

    plt.xlim(0.01, 0.30)
    plt.ylim( 1e3,  1e5)

    pl.savefig('bgs_mockmultipoles_%s.pdf' % ztype, bbox_inches="tight")
    

## Correlation functions.                                                                                                                            
## SimulationBoxPairCount and SurveyDataPairCount.                                                                                                          
 
## Users can compute the correlation function of a Catalog using SimulationBox2PCF and SurveyData2PCF.                                                  
## Available modes: '1d', '2d', 'projected' or 'angular'.                                                                                                   

## xi = nbodykit.algorithms.paircount_tpcf.tpcf.SurveyData2PCF(mode, data, randoms, edges, cosmo=cosmo, Nmu=20, pimax=None, data2=None, \                  
##                                                             randoms2=None, ra='RA', dec='DEC', redshift='Redshift', weight='Weight', show_progress=False)


if __name__ == "__main__":
    nmocks      =  5

    ztypes      =  ["RZEE", "GZEE", "LZEE"]

    create_mock =   False
    calc_poles  =   False

    list_poles  = {x: [] for x in ztypes}

    for num in np.arange(0, nmocks, 1):
        if create_mock is True:
          create_catalog(num = num)
    
        for ztype in ztypes:
            fname  = "./poles/poles_%s_%d.json" % (ztype, num)

            if calc_poles is True:
                ## assign_randzs(ztype = ztype)

                rpoles = multipoles(DR12=False, ztype = ztype, num = num)

                rpoles.save(fname)

            else:
                rpoles = ConvolvedFFTPower.load(fname)
                
            list_poles[ztype].append(rpoles)
            

    for ztype in ztypes:
        pl.clf()

        for rpoles in list_poles[ztype]: 
            plot_multipoles(rpoles, ztype = ztype)

    print("\n\nDone.\n\n")
