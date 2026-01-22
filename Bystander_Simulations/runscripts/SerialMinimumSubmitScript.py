import os
import numpy as np
import subprocess as sp

def PrintCMDS(cmds):
  for c in cmds:
    print(c.split(';')[0], c.split('-w')[1].split('--')  [0])

def ConstructCMD(Params):
  CMD = '''LatticeMembrane --Tag='_Z_{P[zc]}'  '''.format(P = Params)
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
  Params['zcvals'] = np.linspace(Params['zmin'], Params['zmax'], Params['RepNums']*Params['wnums'], endpoint=True)
  Params['zcvals'] = [15.1,15.5,15.9]
  zcutoff  = (Params['LLBT'] + Params['LLBB']) - 4.0
  compBval = lambda z: max(Params['BiasStrFloor'], Params['BiasStrMax']*(Params['perc']/100.)) if (z > zcutoff) else Params['BiasStrFloor']
  #Params['BiasStrs'] = [compBval(z) for z in Params['zcvals']]
  Params['BiasStrs'] = [30.0]
  return Params

def CompleteParams(Params):
  Params['rhoBYB'] = (1. - Params['perc']/100)*Params['rhototB']
  #Params['rhoBYB'] = 7600.0
  Params['rhoLBB'] = Params['SwitchLBB']*Params['rhototB']*Params['perc']/100
  #Params['rhoLBB'] = 400.0
  return Params

percentlist = [6.0]
#percentlist = [2.0, 4.0, 6.0, 8.0]
#percentlist = [5.0]
#percentlist = [11.0,12.0]

#lenlist     = [20., 30., 40.]
#lenlist     = [10.0,15.0,20.0,25.0,30.0,35.0,40.0]
lenlist      = [20.0]
rhoLBlist   = [500.0]

#rhoBYlist   = [5000.0]
rhoBYlist = [4000.0]

Params = {}
Params['cores']    = 8
#Params['res']      = 2.0
Params['res']      = 2.0
Params['RepNums']  = 1

Params['MCSweeps'] = 100000
#Params['MCSweeps'] = 500000

Params['MDisp']    = 1.5

Params['LSBT']     = 3.5
Params['LSBB']     = 7.5
Params['rhoSBT']   = Params['rhoSBB'] = 200.0

Params['LBYT']     = 17.5

Params['LBYB']     = 17.5
Params['LLBB']     = 30.0

Params['rhototB']  = 8000.0

Params['SwitchLBB'] = 1

Params['USB']      = 10.5
#Params['ULB']      = 4.75
LBlist = [4.75]
Params['z0']       = (lenlist[-1] + Params['LLBB']) + 10.0
#Params['z0']       = [15.1,15.5,15.9]

Params['Nx']       = 70
Params['Lx']       = 12*Params['Nx']
Params['kc']       = 20.0

Params['BiasStrFloor'] = 20.0
Params['BiasStrMax']   = 200.0

Params['OutDir']     = '../Simulation_Results/'
Params['fnums']      = 0
Params['dnums']      = 50000

#Params['MCSweeps']   = 200
#percentlist          = percentlist[:2]
#lenlist              = [lenlist[0]]

cmds = []
for rhoBY in rhoBYlist:
    Params['rhoBYT'] = rhoBY
    for rhoLB in rhoLBlist:
        Params['rhoLBT'] = rhoLB
        for LB in LBlist:
          Params['ULB']      = LB
          for l in lenlist[:]:
            Params['LLBT'] = l
            for p in percentlist[:]:
                Params['perc']     = p
                Params             = CompleteParams(Params)
        
                Params['res']      = 2.0
                Params['zmin']     = 15.0
                Params['zmax']     = (Params['LLBT'] + Params['LLBB']) + 2.0 + 10.0 
                Params             = ComputeZcValsParams(Params)

                for i in range(len(Params['zcvals'])):
                    Params['zc']  = Params['zcvals'][i]
                    #Params['Bst'] = Params['BiasStrs'][i]
                    Params['Bst'] = 20.0
                    CMD = ConstructCMD(Params)
                    cmds.append(CMD)        

                #Params['res']      = 0.25
                #Params['zmin']     = (Params['LLBB'] + Params['LLBT']) - 2.0
                #Params['zmax']     = (Params['LLBB'] + Params['LLBT']) + 4.0  
                #Params             = ComputeZcValsParams(Params)
                        
                #for i in range(len(Params['zcvals'])):
                #    Params['zc']  = Params['zcvals'][i]
                #    Params['Bst'] = Params['BiasStrs'][i]
                #    CMD = ConstructCMD(Params)
                #    cmds.append(CMD)

PrintCMDS(cmds)
with open('status.dat', 'w') as fp: print()

i=0;processes = [];
while cmds != []:
  while len(processes)<Params['cores'] and cmds!=[]:
      processes.append(sp.Popen(cmds.pop(), stdout = sp.PIPE, shell = True))
      i = i + 1
  exitval = [p.wait() for p in processes]
  processes = []
  with open('status.dat', 'a') as fp: fp.write('sims remaining: {} {}\n'.format(len(cmds), i))



