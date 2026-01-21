#!/usr/bin/python3
import numpy as np
from os import listdir
from sys import argv

def GetParamFromString(string, keyword): return string.split(keyword+'_')[1].split('_')[0]
def StripField(string, keyword): 
  parts    = string.split(keyword+'_')
  parts[1] = '_'.join(parts[1].split('_')[1:])
  return parts[0]+parts[1]
def ExtractVal(keyword, string):
  return string.split(keyword+'=')[1].split()[0]

# returns the log of a sum of exponentials, whose exponents
# are given in array "exponents"
# (needed when exponents are large enough that the direct sum
# is problematic)
def add_exp(exponents):
    emax = np.max(exponents)
    return emax + np.log(np.sum(np.exp(exponents-emax)))

# read in a free energy profile
def readhf(filename, rcutoff):
    infile = open(filename,"r")
    lines = infile.readlines()
    infile.close()

    hvals = []
    Fvals = []
    for line in lines:
        numbers = line.split()
        if numbers[0] != '#' and float(numbers[0]) >= rcutoff:            
          hvals.append( float(numbers[0]) )
          Fvals.append( float(numbers[2]) )
    
    return [np.array(hvals), np.array(Fvals), len(hvals)], lines[0]

def ExtractFeatures(rawDat):
  hvals, Fvals, ndat = rawDat
  imax  = np.argmax(Fvals)
  Fmax  = Fvals[imax]
  imin  = np.argmin(Fvals)
  Fmin  = Fvals[imin]
  if hvals[imax] < 30:
    Fddag = Fmax
  else:
    Fddag = np.max(Fvals[:imin])
  #Fddag = Fmax if hvals[imax] < 30. else np.amax(Fvals[:imax])
  #FvalsLeft = Fvals[:imax]
  #Fmax =np.max(FvalsLeft)
  return [Fmax, Fmin, Fddag]

def ComputeEnhancement( rawDat, Features, Params):
  hvals, Fvals, ndat = rawDat
  Fmax, Fmin, Fddag  = Features
  alpha, R           = Params 
  print('FDag is ' + str(Fddag) +"\n" + "FMax is "+ str(Fmax) + "\n" + "Fmin is " + str(Fmin))
  #print('Intsum is ' + str(R))
  dh     = hvals[1:] - hvals[:-1]
  ExpdF  = np.exp(Fvals[:-1] - Fmax)
  intsum = sum(dh*ExpdF)
  print("Intsum is " + str(intsum))
  A      = Fmax + np.log(intsum/R)
  
  intsum2 = 0 
  for i in range(hvals.shape[0] - 1):
    dh2     = hvals[1:i+1] - hvals[:i]
    ExpdF   = np.exp(-Fvals[i] + Fmin + Fvals[:i] - Fddag)
    intsum2 = intsum2 + dh[i]*np.sum(dh2*ExpdF)
  #print("Intsum2 is " + str(intsum2))
  B = Fddag - Fmin + np.log(alpha/R**2) + np.log(intsum2)
  #B = Fddag - Fmin + np.log(alpha) + np.log(intsum2)
  lograte = -add_exp(np.array([0,A,B]))
  #print(R)
  print(np.exp(lograte))
  print('\n')
  return [lograte, A, B, intsum, intsum2, Fddag, Fmin, Fmax] 

def ObtainRates(R, alpha, InDir, fields, rcutoff):
  AllFiles = [InDir+f for f in filter(lambda x: 'ThetaAveFE_' in x, listdir(InDir))]
  results  = np.zeros((len(fields) + 8,))
  for i, fname in enumerate(AllFiles):
    #print(fname)
    rawDat, firstline = readhf(fname, rcutoff)
    res  = []
    Features = ExtractFeatures(rawDat)
    for field in fields:
      res.append(ExtractVal(field, firstline))
    res += ComputeEnhancement(rawDat, Features, [alpha, R])
    results  = np.vstack((results, res))
 
  return results[1:, :]

def StoreResults(results, alpha, OutDir, fields):
  ind = np.lexsort([results[:, -i].astype(float) for i in range(results.shape[1] - 1) ])
  results = results[ind]

  fp = open(OutDir+'kineticEnhancement_alpha_{}_.dat'.format(alpha), 'w')
  line ='#'
  for i, p in enumerate(fields):
    line = line + p + '('+str(i+1)+') '
  line += 'enh({}) '.format(i + 2)
  line += 'A({}) '.format(i + 3)
  line += 'B({}) '.format(i + 4)
  line += 'intsum1({}) '.format(i + 5)
  line += 'intsum2({}) '.format(i + 6)
  line += 'Fddag({}) '.format(i + 7)
  line += 'Fmin({}) '.format(i + 8)
  line += 'Fmax({}) '.format(i + 9)
  fp.write( line + '\n')
  
  for res in results:
    for r in res:
      fp.write(str(r)+' ')
    fp.write('\n')
  fp.close()

InDir  = 'Free_Energies/'
OutDir = 'Kinetic_Rates/'

if len(argv)>1:
  InDir = argv[1]
  if InDir[-1] != '/':InDir += '/'

# geometric parameters for rate calculation
#R = 1e3

R=1000
alphalist = [1e-11]
#alphalist = [1e-20]
fields    = ['LSBT', 'rhoSBT', 'LSBB', 'rhoSBB', 'LBYT', 'rhoBYT', 'LBYB', 'rhoBYB', 'LLBT', 'rhoLBT', 'LLBB', 'rhoLBB']
rcutoff   = 15.

for alpha in alphalist:
  results = ObtainRates(R, alpha, InDir, fields, rcutoff)
  StoreResults(results, alpha, OutDir, fields)
  #MakePlots(resultset, keys, alpha)



