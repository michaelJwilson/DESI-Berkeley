from    astropy.io   import  fits

import  numpy        as      np
import  collections
import  glob
import  os


DMOCKS     = os.environ['DMOCKS']
CHALLENGE  = os.environ['CHALLENGE']
DR12       = os.environ['DR12']

infiles    = open('../infiles.txt', 'r')

allfits    = collections.OrderedDict()
roots      = {'DMOCKS': DMOCKS, 'CHALLENGE': CHALLENGE, 'DR12': DR12}


def get_filenames(printit = True):
  for line in infiles:
    if not line.strip().startswith("#"):
      [name, root, path]   = line.split() 
         
      allfits[name]        = roots[root] + path
         
      if printit:
        print allfits[name]

  return allfits

def print_header(allfits, name, printhead = False, tile = '00098'):  
  fname                    = allfits[name]
  
  try:
    dat                    =  fits.open(fname)

  except:
    if name == "FIBERS":
      fname                = "/project/projectdirs/desi/datachallenge/dc17b/fiberassign/tile_%s.fits" % tile
      dat                  = fits.open(fname)

    elif name == "TILES":
      print("\n\nLoading tiles.\n\n")

    else:
      exit(1)

  hdrs                     = [x.header for x in dat]
  hdrkeys                  = list(dat[1].header.keys()) 
  
  print("\nGetting:  %s\n" % fname)

  if printhead:
    print(repr(dat[1].header))

  for i, x in enumerate(hdrkeys):
    if "TTYPE" in x:
      hdr = dat[1].header[i]

      if "COEFF" not in hdr:
        print "".join(np.str(x).ljust(25) for x in [i, hdr, dat[1].data[0][hdr], dat[1].data[1][hdr], dat[1].data[2][hdr], dat[1].data[3][hdr]])

'''
Available keys:  
Input:      ELG, LRG, QSO, Lya
Imaging:    TARGETS, TARGETSTRUTH, FIBERS, CONDITIONS, IN2OUT
DR12:       DR12SDATA, DR12SRAND
DarkSky:    DSELG, DSLRG, DSQSO, DSRAN
Pointings:  TILES
'''

if __name__ == "__main__":
  print("\n\nWelcome to seefits.\n\n")

  allfits          = get_filenames(printit = False)
  allfits['TILES'] = '../desi-tiles.fits'

  print("Available keys:  ")

  for x in allfits.keys():
    print x

  name   = None

  while name is not "Exit":
    name = input("\n\nEnter the key for your chosen file:\n")
    name = str(name)

    try: 
      print_header(allfits, name)

    except:
      raiseValueError("\n\nChosen key is not available.")

  print("\n\nDone.\n\n")
