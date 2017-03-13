import os
import sys
import numpy as np
import glob
from subprocess import call
import pickle

# Parse input, if possible
if len(sys.argv) != 4:
  print 'ERROR format for calling batch.py is'
  print '  python batch.py [folder] [N] [tmin]'
  print '  where'
  print '    [folder]: (string) path to dump files'
  print '    [N]:      (float)  superphoton resolution parameter'
  print '    [tmin]:   (float)  time to begin averaging'
  sys.exit()
folder = os.path.join(sys.argv[1], '')
N      = float(sys.argv[2])
tmin   = float(sys.argv[3])

# Get dumps
dumps = np.sort(glob.glob(folder+'dump_fluid*'))
Nd = len(dumps)

print '\n  **************************************************'                   
print '  ******************** BATCH.PY ********************'                     
print '  **************************************************'                     
print ''                                                                         
print '  ****************** CONFIGURATION *****************'                     
print '    Folder: %s' % folder                                                  
print '    N:      %g' % N                                                  
print '    tmin:   %g' % tmin                                               
print '    dumps:  %i' % Nd                                                      
print ''

# Get dimensions and set arrays
EBINS = 200
THBINS = 6
VAR = 7
DTd = 5.

nu = np.zeros(EBINS)
nuLnu = np.zeros([THBINS, EBINS])
navg = 0

ME  = 9.1093897e-28
CL  = 2.99792458e10
HPL = 6.6260755e-27

# Call grmonty for each dump
for n in xrange(Nd):
  f = open(dumps[n], 'r')
  t = float(f.readline().split()[0])
  print t
  f.close()
  if t >= tmin:

  #if n*DTd >= tmin:
    navg += 1
    print dumps[n]
    print ['./igrmonty2d', '%g' % N, dumps[n]]
    call(['./igrmonty2d', '%g' % N, dumps[n]])

    # Read in grmonty output
    spec = np.loadtxt('grmonty.spec', skiprows=1)
    for k in xrange(THBINS):
      for n in xrange(EBINS):
        nuLnu[k][n] += spec[n][k*VAR+1]
    for n in xrange(EBINS):
      nu[n] = 10.**(spec[n][0])*ME*CL*CL/HPL

for k in xrange(THBINS):
  for n in xrange(EBINS):
    nuLnu[k][n] /= navg

avgs = {}
avgs['THBINS'] = THBINS
avgs['NUBINS'] = EBINS
avgs['nu'] = nu
avgs['nuLnu'] = nuLnu

pickle.dump(avgs, open('batch_out.p','wb'))

