import  matplotlib;             matplotlib.use('PDF')
 
import  math
import  numpy           as      np
import  pylab           as      pl

from    nbodykit.lab    import  *
from    astropy.io      import  fits


cosmo              = cosmology.Cosmology(h=0.7).match(Omega0_m=0.31)

ftiles             = fits.open('./desi-tiles.fits')
ftiles             = ftiles[1]


tiles              = {}
tiles['META']      = {}

## tile radius [degree].                                                                                            
tiles['META']['DRADIUS'] = 1.605
tiles['META']['RRADIUS'] = np.radians(tiles['META']['DRADIUS'])

for tile in ftiles.data[::1]:
  if tile[4] == 1: ## In DESI;  Tiled whole 4PI steradians incase of footprint changes.    
    id                = tile[0]

    tiles[id]         = {}

    ## RA, DEC in DEG
    tiles[id]['RA']   = tile[1] 
    tiles[id]['DEC']  = tile[2]

    ## print id, tiles[id]['RA'], tiles[id]['DEC']
    
    tiles[id]['THETA'] = np.pi/2. - np.radians(tiles[id]['DEC'])

    tiles[id]['SINT']  = np.sin(tiles[id]['THETA'])
    tiles[id]['COST']  = np.cos(tiles[id]['THETA'])

    pl.plot(tiles[id]['RA'], tiles[id]['DEC'], 'k.')

'''
## RAND GENERATION
NRANDPERTILE       = 1
 
rands              = {}

for id in tiles:
    if id is not 'META':
        print id, tiles[id]['COST']
    
        radraws              =  np.random.uniform(low =  0.0, high = 1.0, size = NRANDPERTILE)
        cosdraws             =  np.random.uniform(low = -0.5, high = 0.5, size = NRANDPERTILE)

        rands[id]            = {}
        
        rands[id]['RAS']     =  2. * np.pi * radraws
        rands[id]['RAS']     =  np.degrees(rands[id]['RAS'])

        rands[id]['COSTS']   =  tiles[id]['COST'] ## + tiles[id]['SINT'] * tiles['METAx']['RRADIUS'] * cosdraws

        rands[id]['THETAS']  =  np.arccos(rands[id]['COSTS'])

        rands[id]['DECS']    =  np.pi/2. - np.arccos(rands[id]['THETAS'])
        rands[id]['DECS']    =  np.degrees(rands[id]['DECS'])

        pl.plot(rands[id]['RAS'], rands[id]['DECS'], 'k.')
'''   
edges              = np.logspace(0.0, 3.0, num=100)

tiles['Position']  = transform.SkyToCartesian(tiles[10]['RA'], tilse[10]['DEC'], np.ones_like(tiles[10]['RA']), cosmo=cosmo)

## http://nbodykit.readthedocs.io/en/latest/api/_autosummary/nbodykit.algorithms.paircount_tpcf.tpcf.html
## nbodykit.algorithms.paircount_tpcf.tpcf.SurveyData2PCF: 
##
## Default:  'RA', 'DEC', 'Redshift'   
## print SurveyData2PCF('2d', data, randoms, edges, cosmo=cosmo, Nmu=20)
     
pl.xlim(360.0,  0.0)
pl.ylim(-20.0, 90.0)

pl.savefig('footprint.pdf')

