import  matplotlib;                              matplotlib.use('PDF')

import  itertools
import  numpy                                 as np
import  healpy                                as hp
import  astropy.io.fits                       as fits
import  pylab                                 as pl 

from    astropy              import units     as u
from    astropy.coordinates  import Galactic
from    astropy.coordinates  import SkyCoord


## -- Generate fractional areas that each tile contributes to each pixel, via a high resolution hp.map -- 
nside         = 1024
npix          = hp.nside2npix(nside)
imap          = np.zeros(npix)
pix_area      = hp.nside2pixarea(nside, degrees = False) ## [radians].

ftiles        = fits.open('./desi-tiles.fits')
ftiles        = ftiles[1]

## tile radius [degree]; radius [radians]                                                                                                                                                                                               
dradius       = 1.605
radius        = np.radians(dradius)

vecs          = []

for tile in ftiles.data[::1]:
  if tile[4]  == 1: ## In DESI;  Tiled whole 4PI steradians incase of footprint changes.                                                                                                                                                 
    id         = tile[0]

    ra         = tile[1]  ## RA, DEC in DEG 
    dec        = tile[2]

    colat      = np.radians(90. - dec)
    longitude  = np.radians(ra)

    vecs.append(hp.pixelfunc.ang2vec(colat, longitude, lonlat=False))

## -- QUERY DISC --                                                                                                                                                                                                                       
## Returns (nested) pixels whose centers lie within the disk defined by vec and radius (in radians) if inclusive is False,                                                                                                                 
## or which overlap with this disk if inclusive is True.                                                                                                                                                                                   
##                                                                                                                                                                                                                                         
## -- Args --                                                                                                                                                                                                                              
## radius:  The radius [radians] of the disk                                                                                                                                                                                               
## fact:    When inclusive is True, overlapping test will be done at resolution fact*nside; For NESTED,                                                                                                                                    
##          fact must be positive integer satisfying 2^n where n <= 30; Default = 4.                                                                                                                                                       
## vec:     The coordinates of unit vector defining the disk center.

## Inclusive =  True implies that a pixel can overlap multiple tiles (pointings).
## query_disc:  returns list of pixels whose centre falls in given tile, sorted by increasing order. 
lnested_ipix = [hp.query_disc(nside, vec, radius, inclusive=True, fact=2, nest=False) for vec in vecs]

## Convert list of unequal length lists to rows in a np.array;  Python 3:  itertools.zip_longest, Python 2:  itertools.izip_longest.
nested_ipix  =  np.array(list(itertools.izip_longest(*lnested_ipix, fillvalue = -99)), dtype=np.int32).T

##  Find unique pixels (Not necessarily imap unless tiles tesselate the full sky).
unique       =  np.unique(nested_ipix)
unique       =  unique[1:]              ##  First in unique is -99

imap[unique] =  1.0

np.savetxt("desi_footprint_imap.dat", imap)

hp.mollview(imap, title="DESI", coord=['G','E'])

pl.savefig("footprint.pdf")

'''
for u in unique[1:]:
    ## Find the 'first' occurence of a given pixel in the array, set all the rest with that value to nan.
    args                  = np.argwhere(nested_ipix == u)
    
    ## Keep first, set all others to -99.
    nested_ipix[args[1:]] = -99

## Convert back to list, removing -99s (repeated pixels).
result       = [x[x >= 0] for x in nested_ipix] 

for x in result:
    print x 
    print
'''


