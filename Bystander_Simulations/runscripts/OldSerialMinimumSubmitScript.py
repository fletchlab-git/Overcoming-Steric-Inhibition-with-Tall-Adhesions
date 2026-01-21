import os
import numpy as np


def ConstructSim(Params, i):
  args = "time ../LatticeMembrane --MCSweeps={P[MCSweeps]} --Tag="+Params['Tag']+" --BiasStr={P[BiasStr]} --Equil='No' -v 'Minimum' "
  args = args + "--MeshDisp={P[MDisp]} --z0={P[zmax]} -w {zc} --BC='Frame' "
  args = args + "--Lx={P[Lx]} --Ly={P[Lx]} --kc={P[kc]} --Nx {P[Nx]} --Ny {P[Nx]} "
  args = args + "--NumPType=3 --BindPType='{P[USB]} 0.0 {P[ULB]}' "
  args = args + "--LPType='{P[LSBinder]} {P[LNonBinder]} {P[LLBinder]} {P[LSBinder]} {P[LNonBinder]} {P[LLBinder]}' "
  args = args + "--RanPType='4.0 0.0 4.0' --Densities='{P[rhoSBinder]} {P[rhononBinder]} {P[rhoLBinder]} {P[rhoSBinder]} {P[rhononBinder]} {P[rhoLBinder]}' "
  args = args + "--OutDir='{P[OutDir]}' --framenums={P[fnums]} --datnums={P[dnums]} --SEED=${{RANDOM}}"
  return args.format(P = Params, zc = Params['zcvals'][i])

def ComputeZcValsParams(Params):
  zcvals           = np.linspace(Params['zmin'], Params['zmax'], Params['wnums'], endpoint=True)
  Params['zcvals'] = map(str, zcvals)
  return Params

def CompleteParams(Params):
  Params['rhononBinder'] = (1. - Params['perc']/100)*Params['rhotot']
  Params['rhoLBinder']   = Params['rhotot']*Params['perc']/100
  return Params

percentlist = [0., 10., 25., 40.0, 100.]

Params = {}

Params['wnums']      = 40

Params['MCSweeps']   = 50000
Params['MDisp']      = 0.4

Params['LSBinder']   = 5.5
Params['LNonBinder'] = 17.5
Params['LLBinder']   = 10.0
Params['USB']        = 10.5
Params['ULB']        = 4.0
Params['rhoSBinder'] = 200.0
Params['rhotot']     = 8000.0

Params['zmin']       = 0.0
Params['zmax']       = 2*Params['LLBinder'] + 2.0 + 20.0

Params['Nx']         = 50
Params['Lx']         = 4*Params['Nx']
Params['kc']         = 20.0
Params['BiasStr']    = 3.0

Params['OutDir']     = '../Data/'
Params['fnums']      = 20
Params['dnums']      = 10000

Params['Tag']        = "AriExperimentMin_Ll_{P[LLBinder]}_percent_{P[perc]}"

#Params['wnums']      = 3
#Params['MCSweeps']   = 200
#percentlist          = percentlist[:1]

Paramc = ComputeZcValsParams(Params)

for p in percentlist[1:3]:
  Params['perc'] = p
  Params         = CompleteParams(Params)
  for i in range(Params['wnums'])[:-1]: 
    print ConstructSim(Params, i)
    os.system(ConstructSim(Params, i))
   

