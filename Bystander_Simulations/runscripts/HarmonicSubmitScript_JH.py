import os
import numpy as np

def Preamb(Params):
  return '''#!/bin/bash -l
#SBATCH --qos={P[Queue]}     #Use the shared QOS
#SBATCH --array=1-{P[wnums]}
#SBATCH --ntasks=1
#SBATCH -t {P[rtime]}
#SBATCH -J Harmonic
#SBATCH -A m3012
#SBATCH --constraint=haswell

zcs={P[zcvals]}

'''.format(P = Params)

def MakeCommand(Params, i):
  return '''##Run No. {index}
Index=$((${{SLURM_ARRAY_TASK_ID}}-1+{shift:d}))
Tag="AriExperiment_percent_{P[perc]}_zc_${{zcs[${{Index}}]}}"
zc=${{zcs[${{Index}}]}}

args="--MCSweeps={P[MCSweeps]} --Tag=${{Tag}} --BiasStr={P[BiasStr]} --Equil='Yes' -v 'Harmonic' "
args=${{args}}"--MeshDisp={P[MDisp]} --z0=${{zc}} -w ${{zc}} "
args=${{args}}"--Lx={P[Lx]} --Ly={P[Lx]} --kc={P[kc]} --Nx {P[Nx]} --Ny {P[Nx]} "
args=${{args}}"--NumPType=3 --BindPType='{P[USB]} 0.0 {P[ULB]}' "
args=${{args}}"--LPType='{P[LSBinder]} {P[LNonBinder]} {P[LLBinder]} {P[LSBinder]} {P[LNonBinder]} {P[LLBinder]}' "
args=${{args}}"--RanPType='4.0 0.0 4.0' --Densities='{P[rhoSBinder]} {P[rhononBinder]} {P[rhoLBinder]} {P[rhoSBinder]} {P[rhononBinder]} {P[rhoLBinder]}' "
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
  zcvals           = np.linspace(Params['zmin'], Params['zmax'], Params['RepNums']*Params['wnums'], endpoint=True)
  Params['zcvals'] = '('+' '.join(map(str, zcvals))+')'
  return Params

def CompleteParams(Params):
  Params['rhononBinder'] = (1. - Params['perc']/100)*Params['rhotot']
  Params['rhoLBinder']   = Params['rhotot']*Params['perc']/100
  return Params

percentlist = [0., 5., 10., 25., 30., 40.0, 100.]

Params = {}

Params['Queue']      = 'shared'
Params['wnums']      = 15
Params['RepNums']    = 1

Params['MCSweeps']   = 2200000
Params['rtime']      = '24:00:00'
Params['MDisp']      = 0.4

Params['zmin']       = 10.0
Params['zmax']       = 25.0

Params['LSBinder']   = 5.5
Params['LNonBinder'] = 17.5
Params['LLBinder']   = 8.0
Params['USB']        = 10.5
Params['ULB']        = 4.0
Params['rhoSBinder'] = 200.0
Params['rhotot']     = 8000.0

Params['Nx']         = 50
Params['Lx']         = 4*Params['Nx']
Params['kc']         = 5.0
Params['BiasStr']    = 16.0

Params['OutDir']     = '/global/homes/j/jhasnain/Data/FBMembrane/'
Params['fnums']      = 20
Params['dnums']      = 10000

#Params['wnums']      = 3
#Params['MCSweeps']   = 200
#Params['rtime']      = '00:02:00'
#Params['Queue']      = 'debug'


Params = ComputeZcValsParams(Params)

for p in percentlist[:]:
  Params['perc'] = p
  Params         = CompleteParams(Params)
  open('HarmonicSim.sbatch', 'w').write(ConstructScript(Params))
  os.system('sbatch HarmonicSim.sbatch;rm HarmonicSim.sbatch;') 

