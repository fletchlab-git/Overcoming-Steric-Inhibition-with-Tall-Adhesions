import os
import numpy as np

def Preamb(Params):
  return '''#!/bin/bash -l
#SBATCH --qos={P[Queue]}     #Use the shared QOS
#SBATCH --array=1-{P[wnums]}
#SBATCH --ntasks=1
#SBATCH -t {P[rtime]}
#SBATCH -J Harmonic
#SBATCH -A m1415
#SBATCH --constraint=haswell

zcs={P[zcvals]}
BiasStrs={P[BiasStrs]}

'''.format(P = Params)

def MakeCommand(Params, i):
  return '''##Run No. {index}
Index=$((${{SLURM_ARRAY_TASK_ID}}-1+{shift:d}))
Tag="AriExperimentMin_rhoBYT_{P[rhoBYT]}_rhoLBT_{P[rhoLBT]}_Ll_{P[LLBT]}_ULB_{P[ULB]}_percent_{P[perc]}"
zc=${{zcs[${{Index}}]}}
BiasStr=${{BiasStrs[${{Index}}]}}

args="--MCSweeps={P[MCSweeps]} --Tag=${{Tag}} --BiasStr=${{BiasStr}} --Equil='No' -v 'Minimum' "
args=${{args}}"--MeshDisp={P[MDisp]} --z0={P[z0]} -w ${{zc}} --BC='Frame' "
args=${{args}}"--Lx={P[Lx]} --Ly={P[Lx]} --kc={P[kc]} --Nx {P[Nx]} --Ny {P[Nx]} "
args=${{args}}"--NumPType=3 --BindPType='{P[USB]} 0.0 {P[ULB]}' "
args=${{args}}"--LPType='{P[LSBT]} {P[LBYT]} {P[LLBT]} {P[LSBB]} {P[LBYB]} {P[LLBB]}' "
args=${{args}}"--RanPType='4.0 0.0 4.0' --Densities='{P[rhoSBT]} {P[rhoBYT]} {P[rhoLBT]} {P[rhoSBB]} {P[rhoBYB]} {P[rhoLBB]}' "
args=${{args}}"--OutDir='{P[OutDir]}' --framenums={P[fnums]} --datnums={P[dnums]}"

echo ../../LatticeMembrane ${{args}}
time bash -c "../../LatticeMembrane ${{args}} --SEED=${{RANDOM}}"

'''.format(P = Params, shift=i*Params['wnums'], index=i)

def ConstructScript(Params):
  cmd = Preamb(Params)
  for i in range(Params['RepNums']):cmd=cmd+MakeCommand(Params, i)
  return cmd

def AssignIndex(z, zmin, dz):return int(round( (z-zmin)/dz, 0 ))

def ComputeZcValsParams(Params):
  if Params['Queue'] == 'debug': Params['wnums'] = 2 
  else:                          Params['wnums'] = int( (Params['zmax'] - Params['zmin'])/Params['res'] ) + 1
    
  zcvals           = np.linspace(Params['zmin'], Params['zmax'], Params['RepNums']*Params['wnums'], endpoint=True)
  Params['zcvals'] = '('+' '.join(map(str, zcvals))+')'

  zcutoff  = (Params['LLBT'] + Params['LLBB']) - 4.0
  compBval = lambda z: max(Params['BiasStrFloor'], Params['BiasStrMax']*(Params['perc']/100.)) if (z > zcutoff) else Params['BiasStrFloor']
  bvals              = [compBval(z) for z in zcvals] 
  Params['BiasStrs'] = '('+' '.join(map(str, bvals))+')'

  #zcvals             = [z for z in zcvals if z>zcutoff]
  #Params['zcvals']   = '('+' '.join(map(str, zcvals))+')'
  #bvals              = [compBval(z) for z in zcvals if z>zcutoff]
  #Params['BiasStrs'] = '('+' '.join(map(str, bvals))+')'
  #Params['wnums']    = len(zcvals)

  return Params

def CompleteParams(Params):
  Params['rhoBYB'] = (1. - Params['perc']/100)*Params['rhototB']
  Params['rhoLBB'] = Params['SwitchLBB']*Params['rhototB']*Params['perc']/100
  return Params

#percentlist = [0.0, 1.0, 10.0, 25.0, 50.0, 75.0, 100.0]
percentlist = [0.0, 1.0, 5.0, 10.0, 13.0, 15.0, 18.0, 20.0, 25.0, 50.0, 75.0, 100.]
lenlist     = [20., 40.]
rhoLBlist   = [1000., 2000.]
rhoBYlist   = [4000.]

#lenlist     = [20., 30., 40.]
#rhoLBlist   = [200., 1000., 2000.]
#rhoBYlist   = [4000., 10000.]
#lenlist     = [40.]

Params = {}

Params['Queue']    = 'shared'
Params['res']      = 2.0
Params['RepNums']  = 1

Params['MCSweeps'] = 1000000
Params['rtime']    = '4:00:00'
Params['MDisp']    = 1.5

Params['LSBT']     = 3.5
Params['LSBB']     = 7.5
Params['rhoSBT']   = Params['rhoSBB'] = 200.0

Params['LBYT']     = 17.5

Params['LBYB']     = 17.5
Params['LLBB']     = 30.0
Params['rhototB']  = 8000.0

Params['SwitchLBB'] = 1
#Params['rhoLBT']   = 0.0 #Perhaps needs to be increased b/c the longbinders on the macrophage is unknown
#Params['SwitchLBB'] = 0

Params['USB']      = 10.5
Params['ULB']      = 4.75

Params['z0']       = (lenlist[-1] + Params['LLBB']) + 10.0

Params['Nx']       = 70
Params['Lx']       = 12*Params['Nx']
Params['kc']       = 20.0

Params['BiasStrFloor'] = 20.0
Params['BiasStrMax']   = 200.0

Params['OutDir']     = '/global/homes/j/jhasnain/Data/FBMembrane/'
Params['fnums']      = 0
Params['dnums']      = 10000

#Params['MCSweeps']   = 200
#Params['rtime']      = '00:02:00'
#Params['Queue']      = 'debug'
#percentlist          = percentlist[:1]
#lenlist              = [lenlist[0]]
for rhoBY in rhoBYlist:
    Params['rhoBYT'] = rhoBY
    for rhoLB in rhoLBlist:
        Params['rhoLBT'] = rhoLB
        for l in lenlist[:]:
            Params['LLBT'] = l
            for p in percentlist[:]:
                Params['perc']     = p
                Params             = CompleteParams(Params)
        
                Params['res']      = 2.0
                Params['zmin']     = 5.0
                Params['zmax']     = (Params['LLBT'] + Params['LLBB']) + 2.0 + 10.0 
                Params             = ComputeZcValsParams(Params)
                open('MinSim.sbatch', 'w').write(ConstructScript(Params))
                os.system('sbatch MinSim.sbatch; rm MinSim.sbatch;') 
        
                Params['res']      = 0.25
                Params['zmin']     = (Params['LLBB'] + Params['LLBT']) - 2.0
                Params['zmax']     = (Params['LLBB'] + Params['LLBT']) + 4.0  
                Params             = ComputeZcValsParams(Params)
                open('MinSim.sbatch', 'w').write(ConstructScript(Params))
                os.system('sbatch MinSim.sbatch; rm MinSim.sbatch;') 

