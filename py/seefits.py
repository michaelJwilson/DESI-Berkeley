from    astropy.io   import  fits

import  numpy        as      np
import  collections
import  sys
import  os


DMOCKS     = os.environ['DMOCKS']
CHALLENGE  = os.environ['CHALLENGE']
DR12       = os.environ['DR12']

infiles    = open('../infiles.txt', 'r')

allfits    = collections.OrderedDict()
roots      = {'DMOCKS': DMOCKS, 'CHALLENGE': CHALLENGE, 'DR12': DR12}


def get_filenames(printit = True):
  ## From ../infiles.txt, read the relative and absolute path for 
  ## available keys.

  for line in infiles:
    if not line.strip().startswith("#"):
      [name, root, path]   = line.split() 
         
      allfits[name]        = roots[root] + path
         
      if printit:
        print allfits[name]

  return allfits

def print_header(allfits, name, printhead = False, tile = '00098'):  
  ## Given a dictionary allfits, with each key corresponding to a .fits file
  ## print out the available headers and the a number of instances as columns. 

  fname                    = allfits[name]

  ## Currently hard coded fibers file;                                                                                                                      
  if name == "FIBERS":
    fname                  = "/project/projectdirs/desi/datachallenge/dc17b/fiberassign/tile_%s.fits" % tile
  
  print("\nGetting:  %s\n" % fname)
  
  try:
    dat                    =  fits.open(fname)

  except:
      raiseValueError("Provided key was not available in provided dictionary of .fits files.")
  
  hdrs                     = [x.header for x in dat]
  hdrkeys                  = list(dat[1].header.keys()) 
  
  if printhead:
    print(repr(dat[1].header))
  
  for i, x in enumerate(hdrkeys):
    if "TTYPE" in x:
      hdr = dat[1].header[i]

      if "COEFF" not in hdr:
        print "".join(np.str(x).ljust(25) for x in [i, hdr, dat[1].data[0][hdr], dat[1].data[1][hdr], dat[1].data[2][hdr], dat[1].data[3][hdr]])
  

if __name__ == "__main__":
  print("\n\nWelcome to seefits.\n\n")

  allfits          = get_filenames(printit = False)
  allfits['TILES'] = '../desi-tiles.fits'

  print("Available keys:  ")

  '''                                                                                                                                                       
  Input:      ELG, LRG, QSO, Lya                                                                                                                            
  Imaging:    TARGETS, TARGETSTRUTH, FIBERS, CONDITIONS, IN2OUT                                                                                             
  DR12:       DR12SDATA, DR12SRAND                                                                                                                             DarkSky:    DSELG, DSLRG, DSQSO, DSRAN                                                                                                                       Pointings:  TILES                                                                                                                                          
  '''

  for x in allfits.keys():
    print x

  name   = ""

  while name is not "Exit":
    name = raw_input("\n\nEnter the key for your chosen file or Exit to finish:\n")

    if name == "Exit":
      print("\n\nDone.\n\n")

      sys.exit(0)

    if name.endswith('.fits'):
      ## If user input ends in .fits, assume user is providing absolute path to a fits file. 
      ## Write headers of this file to screen.
      allfits['USER'] = name

      name            = 'USER'

    print_header(allfits, name)
