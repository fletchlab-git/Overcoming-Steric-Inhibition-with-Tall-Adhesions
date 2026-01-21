import os

def ConstructScript(Params):
  return '''#!/bin/bash -l
#SBATCH --qos={P[Queue]}     #Use the shared QOS
#SBATCH --array=1-{P[wnums]}
#SBATCH --ntasks=1
#SBATCH -t {P[rtime]}
#SBATCH -J Entropic
#SBATCH -A m3012
#SBATCH --constraint=haswell

zmins={P[zminvals]}
zmaxs={P[zmaxvals]}
zcs={P[zcenters]}
bins={P[binnums]}

Index=$((${{SLURM_ARRAY_TASK_ID}}-1))
Tag="AriExperiment_WL_percent_{P[perc]}_zc_${{zcs[${{Index}}]}}"
zmin=${{zmins[${{Index}}]}}
zmax=${{zmaxs[${{Index}}]}}
zc=${{zcs[${{Index}}]}}
binnums=${{bins[${{Index}}]}}

args="--MCSweeps={P[MCSweeps]} --Tag=${{Tag}} --BiasStr={P[BiasStr]} --Equil='Yes' -v 'WangLandau' "
args=${{args}}"--ABSFmin=${{zmin}} --ABSFmax=${{zmax}} --laps=${{binnums}} "
args=${{args}}"--MeshDisp={P[MDisp]} --z0=${{zc}} -w ${{zc}} "
args=${{args}}"--Lx={P[Lx]} --Ly={P[Lx]} --kc={P[kc]} --Nx {P[Nx]} --Ny {P[Nx]} "
args=${{args}}"--NumPType=3 --BindPType='{P[USB]} 0.0 {P[ULB]}' "
args=${{args}}"--LPType='{P[LSBinder]} {P[LNonBinder]} {P[LLBinder]} {P[LSBinder]} {P[LNonBinder]} {P[LLBinder]}' "
args=${{args}}"--RanPType='4.0 0.0 4.0' --Densities='{P[rhoSBinder]} {P[rhononBinder]} {P[rhoLBinder]} {P[rhoSBinder]} {P[rhononBinder]} {P[rhoLBinder]}' "
args=${{args}}"--OutDir='{P[OutDir]}' --framenums={P[fnums]} --datnums={P[dnums]}"

echo ../../LatticeMembrane ${{args}}
time bash -c "../../LatticeMembrane ${{args}} --SEED=${{RANDOM}}"
'''.format(P = Params)

def AssignIndex(z, zmin, dz):return int(round( (z-zmin)/dz, 0 ))

def ComputeWLParams(Params):
  Params['wnums'] = Params['wnums'] + (1 - Params['wnums']%2) 
  totbnums        = int((Params['gzmax'] - Params['gzmin'])/Params['res'] + 1)
  globalbinvals   = [ Params['gzmin'] + i*Params['res'] for i in range(totbnums) ]

  wnumon2  = (Params['wnums'] + 1)/2
  dw       = (Params['gzmax'] - Params['gzmin'])/wnumon2
  dwon2    = dw/2

  zminvals = [ Params['gzmin'] + i*dwon2 for i in range(Params['wnums']) ]  
  zmaxvals = [ Params['gzmin'] + dw + i*dwon2 for i in range(Params['wnums']) ]  

  zminvals = [ globalbinvals[AssignIndex(z, Params['gzmin'], Params['res'])] for z in zminvals ]
  zmaxvals = [ globalbinvals[AssignIndex(z, Params['gzmin'], Params['res'])] for z in zmaxvals ]
  zcvals   = [ (zmaxvals[i] + zminvals[i])*0.5 for i in range(Params['wnums']) ]
  zcvals   = [ globalbinvals[AssignIndex(z, Params['gzmin'], Params['res'])] for z in zcvals ]
  bnums    = [ int(round((zmaxvals[i] - zminvals[i])/Params['res'] + 1, 0)) for i in range(Params['wnums']) ]
 
  Params['zminvals'] = '('+' '.join(map(str, zminvals))+')'
  Params['zcenters'] = '('+' '.join(map(str, zcvals))+')'
  Params['zmaxvals'] = '('+' '.join(map(str, zmaxvals))+')'
  Params['binnums']  = '('+' '.join(map(str, bnums))+')'

  return Params

def CompleteParams(Params):
  Params['rhononBinder'] = (1. - Params['perc']/100)*Params['rhotot']
  Params['rhoLBinder']   = Params['rhotot']*Params['perc']/100
  return Params

percentlist = [0., 5., 25.0, 40.0, 50., 80.0, 100.]

Params = {}
Params['Queue']      = 'shared'
Params['wnums']      = 21
Params['res']        = 0.1
Params['gzmin']      = 10.5
Params['gzmax']      = 120
Params['MCSweeps']   = 4000000
Params['rtime']      = '20:00:00'
Params['MDisp']      = 0.4

Params['LSBinder']   = 5.5
Params['LNonBinder'] = 17.5
Params['LLBinder']   = 40.0
Params['USB']        = 10.5
Params['ULB']        = 4.0
Params['rhoSBinder'] = 200.0
Params['rhotot']     = 8000.0

Params['Nx']         = 80
Params['Lx']         = 4*Params['Nx']
Params['kc']         = 20.0
Params['BiasStr']    = 4.0

Params['OutDir']     = '/global/homes/j/jhasnain/Data/FBMembrane/'
Params['fnums']      = 20
Params['dnums']      = 10000

#Params['wnums']      = 3
#Params['MCSweeps']   = 200
#Params['rtime']      = '00:02:00'
#Params['Queue']      = 'debug'


Params = ComputeWLParams(Params) 

for p in percentlist[2:3]:
  Params['perc'] = p
  Params = CompleteParams(Params)
  open('WLSim.sbatch', 'w').write(ConstructScript(Params))
  os.system('sbatch WLSim.sbatch;rm WLSim.sbatch;') 

