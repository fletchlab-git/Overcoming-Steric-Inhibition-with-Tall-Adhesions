import os
import numpy as np
import subprocess as sp

def PrintCMDS(cmds):
  for c in cmds:
    print(c.split(';')[0], c.split('-w')[1].split('--')[0])

ParamTemplate = {
    'cores'    : 1, 
    'res'      : 2.0,
    'RepNums'  : 1,
    'MCSweeps' : 200000,
    #'MCSweeps' : 500000,
    
    'MDisp'    : 1.5,

    'LSBT'     : 3.5,
    'LSBB'     : 7.5,
    'rhoSBT'   : 200.0,
    'rhoSBB'   : 200.0,

    'LBYT'     : 17.5,

    'LBYB'     : 17.5,
    'LLBB'     : 30.0,

    'SwitchLBB' : 1,

    'USB'      : 10.5,
    'ULB'      : 4.75,

    'z0'       : 45.0,

    'Nx'       : 70,
    'Lx'       : 12*70,
    'kc'       : 20.0,

    'BiasStrFloor' : 50.0,
    'BiasStrMax'   : 200.0,

    'OutDir'     : '../DataTest/',
    'fnums'      : 0,
    'dnums'      : 10000,
}

def ConstructCMD(Params):
  CMD = '''LatticeMembrane --Tag='AriExperimentMin{P[TagAdd]}_rhoBYT_{P[rhoBYT]}_rhoLBT_{P[rhoLBT]}_Ll_{P[LLBT]}_ULB_{P[ULB]}_percent_{P[perc]}'  '''.format(P = Params)
  CMD = CMD + "--MCSweeps={P[MCSweeps]} --BiasStr={P[BiasStrs]} --Equil='No' -v 'Minimum' ".format(P = Params)
  CMD = CMD + "--MeshDisp={P[MDisp]} --z0={P[z0]} -w {P[zc]} --BC='Frame' ".format(P = Params)
  CMD = CMD + "--Lx={P[Lx]} --Ly={P[Lx]} --kc={P[kc]} --Nx {P[Nx]} --Ny {P[Nx]} ".format(P = Params)
  CMD = CMD + "--NumPType=3 --BindPType='{P[USB]} 0.0 {P[ULB]}' ".format(P = Params)
  CMD = CMD + "--LPType='{P[LSBT]} {P[LBYT]} {P[LLBT]} {P[LSBB]} {P[LBYB]} {P[LLBB]}' ".format(P = Params)
  CMD = CMD + "--RanPType='4.0 0.0 4.0' --Densities='{P[rhoSBT]} {P[rhoBYT]} {P[rhoLBT]} {P[rhoSBB]} {P[rhoBYB]} {P[rhoLBB]}' ".format(P = Params)
  CMD = CMD + "--OutDir='{P[OutDir]}' --framenums={P[fnums]} --datnums={P[dnums]} ".format(P = Params)
  CMD = CMD + " --SEED=${{RANDOM}} ".format(P = Params)

  return CMD

def CompleteParams(Params):
  Params['rhoBYB'] = (1. - Params['perc']/100)*Params['rhototB']
  Params['rhoLBB'] = Params['SwitchLBB']*Params['rhototB']*Params['perc']/100
  return Params

def CheckParamSets(ParamSets):
    length = len(list(ParamSets.values())[0])
    
    for Plist in list(ParamSets.values())[1:]:
        if len(Plist) != length:
            print("Error parameter lengths do not match")
            for key, Plist in ParamSets.items():
                print(key, length, len(Plist), Plist)
            exit(1)
        else:
            length = len(Plist)
    return length


def MakeCMDList(Params, ParamSets):
    SimNums = CheckParamSets(ParamSets)
    cmds = []
    for s in range(SimNums):
        for key in ParamSets.keys():
            Params[key] = ParamSets[key][s]
        Params             = CompleteParams(Params)
        Params['z0']       = max(Params['z0'], Params['LLBT'] + Params['LLBB'] + 10.)
        print("setting z0 to {}, u cool with that?".format(Params['z0']))
        
        zcutoff            = (Params['LLBT'] + Params['LLBB']) - 4.0
        compBval           = lambda z: max(Params['BiasStrFloor'], Params['BiasStrMax']*(Params['perc']/100.)) if (z > zcutoff) else Params['BiasStrFloor']
        Params['BiasStrs'] = compBval(Params['zc'])  
        
        ###############
        #Params['rhoBYB'] = Params['rhoBYT']
        #Params['rhoLBB'] = 0.0
        ###############
        
        CMD                = ConstructCMD(Params)
        cmds.append(CMD)

    return cmds

ParamSets = {}
#ParamSets['zc']      = [15.5]*4          
#ParamSets['rhoBYT']  = [4000.]*4             
#ParamSets['rhoLBT']  = [200.]*4              
#ParamSets['LLBT']    = [20.]*3 + [40.]       
#ParamSets['rhototB'] = [8000.0]*4            
#ParamSets['perc']    = [5.0, 0.0, 10.0, 5.0] 
#ParamSets['TagAdd']  = ['']*4          


#ParamSets['zc']      = [15.5]*2          
#ParamSets['rhoBYT']  = [4000.]*2             
#ParamSets['rhoLBT']  = [200.]*2              
#ParamSets['LLBT']    = [20., 40.]
#ParamSets['rhototB'] = [8000.0]*2            
#ParamSets['perc']    = [2.0, 5.0] 
#ParamSets['TagAdd']  = ['']*2  

ParamSets['zc']      = [15.5]          
ParamSets['rhoBYT']  = [4000.]             
ParamSets['rhoLBT']  = [200.]              
ParamSets['LLBT']    = [40.]
ParamSets['rhototB'] = [8000.0]            
ParamSets['perc']    = [4.0] 
ParamSets['TagAdd']  = ['']


#ParamSets['zc']      = [13.0]*3          
#ParamSets['rhoBYT']  = [4000.0, 6000.0, 8000.0]
#ParamSets['rhoLBT']  = [0.0]*3
#ParamSets['LLBT']    = [30.0]*3
#ParamSets['rhototB'] = [4000.0, 6000.0, 8000.0]
#ParamSets['perc']    = [4000.0, 6000.0, 8000.0]
#ParamSets['TagAdd']  = ['NoBind']*3

Params = ParamTemplate.copy()
cmds   = MakeCMDList(Params, ParamSets)

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



