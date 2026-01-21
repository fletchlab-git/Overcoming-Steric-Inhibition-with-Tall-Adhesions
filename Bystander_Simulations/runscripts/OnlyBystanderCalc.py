import os
import numpy as np
import subprocess as sp

#os.chdir("C:\\Users\\jared\\Downloads\\bystander_simulations\\bystander_simulations\\runscripts\\")
def PrintCMDS(cmds):
  for c in cmds:
    print(c.split(';')[0], c.split('-w')[1].split('--')[0])

def ConstructCMD(Params):
  CMD = '''LatticeMembrane --Tag='JaredExperimentHarmonicNoBind_rhoBYT_{P[rhoBYT]}_rhoLBT_{P[rhoLBT]}_ULB_{P[ULB]}_percent_{P[perc]}'  '''.format(P = Params)
  CMD = CMD + "--MCSweeps={P[MCSweeps]} --BiasStr={P[Bst]} --Equil='No' -v 'Minimum' ".format(P = Params)
  CMD = CMD + "--MeshDisp={P[MDisp]} --z0={P[zc]} -w {P[zc]} --BC='Frame' ".format(P = Params)
  CMD = CMD + "--Lx={P[Lx]} --Ly={P[Lx]} --kc={P[kc]} --Nx {P[Nx]} --Ny {P[Nx]} ".format(P = Params)
  CMD = CMD + "--NumPType=3 --BindPType='{P[USB]} 0.0 {P[ULB]}' ".format(P = Params)
  CMD = CMD + "--LPType='{P[LSBT]} {P[LBYT]} {P[LLBT]} {P[LSBB]} {P[LBYB]} {P[LLBB]}' ".format(P = Params)
  CMD = CMD + "--RanPType='4.0 0.0 4.0' --Densities='{P[rhoSBT]} {P[rhoBYT]} {P[rhoLBT]} {P[rhoSBB]} {P[rhoBYB]} {P[rhoLBB]}' ".format(P = Params)
  CMD = CMD + "--OutDir='{P[OutDir]}' --framenums={P[fnums]} --datnums={P[dnums]} ".format(P = Params)
  CMD = CMD + " --SEED=2"

  return CMD

def ComputeZcValsParams(Params):
  Params['wnums'] = int( (Params['zmax'] - Params['zmin'])/Params['res'] ) + 1
  #Params['zcvals'] = np.linspace(Params['zmin'], Params['zmax'], Params['RepNums']*Params['wnums'], endpoint=True)
  Params['zcvals'] = [15.1,15.5,15.9]
  Params['BiasStrs'] = [50. for z in Params['zcvals']] 
  return Params

rhoLBlist   = [0.0]
rhoBYlist   = [100,200,300,500, 1000.,2000.,4000.,6000., 8000.,10000]
#rhoBYlist   = [2000.,4000.]
#rhoBYlist   = [100,200,300]

Params = {}
Params['cores']    = 12
Params['res']      = 2.0
Params['RepNums']  = 1

Params['MCSweeps'] = 1000000

Params['MDisp']    = 1.5

Params['LSBT']     = 3.5
Params['LSBB']     = 7.5
Params['rhoSBT']   = Params['rhoSBB'] = 200.0

Params['LBYT']     = 17.5
Params['LBYB']     = 17.5

Params['LLBT']     = 30.0
Params['LLBB']     = 30.0
Params['rhoLBB']   = Params['rhoLBT'] = 0.0

Params['USB']      = 10.5
Params['ULB']      = 4.75

#Params['z0']       = (Params['LBYT'] + Params['LBYB']) + 10.0

Params['Nx']       = 70
Params['Lx']       = 12*Params['Nx']
Params['kc']       = 20.0

Params['BiasStrFloor'] = 20.0
Params['BiasStrMax']   = 200.0

Params['OutDir']     = '../Simulation_Results/'
Params['fnums']      = 0
Params['dnums']      = 1000000

#Params['MCSweeps']   = 200
#percentlist          = percentlist[:2]
#lenlist              = [lenlist[0]]

cmds = []
for i in range(len(rhoBYlist)):
  Params['rhoBYT'] = Params['rhoBYB'] = rhoBYlist[i]
  Params['perc']    = rhoBYlist[i]
        
  Params['res']      = 2.0
  Params['zmin']     = 5.0
  Params['zmax']     = (Params['LBYT'] + Params['LBYB']) + 2.0 + 10.0 
  Params             = ComputeZcValsParams(Params)
  for j in range(len(Params['zcvals'])):
    Params['zc']  = Params['zcvals'][j]
    Params['z0']  = Params['zc']
    Params['Bst'] = Params['BiasStrs'][j]
    CMD = ConstructCMD(Params)
    cmds.append(CMD)        

PrintCMDS(cmds)
with open('status.dat', 'w') as fp: print()
print(os.getcwd())

i=0;processes = [];
while cmds != []:
  while len(processes)<Params['cores'] and cmds!=[]:
      processes.append(sp.Popen(cmds.pop(), stdout = sp.PIPE, shell = True))
      i = i + 1
  exitval = [p.wait() for p in processes]
  processes = []
  with open('status.dat', 'a') as fp: fp.write('sims remaining: {} {}\n'.format(len(cmds), i))



