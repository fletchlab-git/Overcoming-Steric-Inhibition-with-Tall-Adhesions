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

Restarts={P[Restarts]}

Index=$((${{SLURM_ARRAY_TASK_ID}}-1))
Restart="${{Restarts[${{Index}}]}}"

echo ../../LatticeMembrane --Input ${{Restart}}
time bash -c "../../LatticeMembrane --Input ${{Restart}} "
'''.format(P = Params)

def GetValFromStr(string, keyword): return string.split(keyword+'_')[1].split('_')[0]

def ListUniquePercentages(InDir):
  Rlist       = filter(lambda x: 'Restart' in x, os.listdir(InDir+'/'))
  preamb      = Rlist[0].split('percent')[0]
  percentlist = list(set( [GetValFromStr( Rfile, 'percent' ) for Rfile in Rlist] ))
  percentlist = sorted( percentlist, key=lambda x: float(x) )
  return [ preamb+'percent_'+f for f in percentlist]

def MakeRestarts(InDir, p):
  Rs = [ x for x in os.listdir(InDir+'/') if p in x]
  Rs = sorted( [InDir+'/'+r for r in Rs], key = lambda x: float(GetValFromStr(x, 'zc')))
  return '('+' '.join(Rs)+')'

InDir = 'Restarts'

Params = {}
Params['Queue']      = 'shared'
Params['wnums']      = 21
Params['rtime']      = '28:00:00'

#Params['wnums']      = 3
#Params['MCSweeps']   = 200
#Params['rtime']      = '00:02:00'
#Params['Queue']      = 'debug'

percentlist = ListUniquePercentages(InDir)

for p in percentlist[:]:
  Params['Restarts'] = MakeRestarts(InDir, p)
  open('WLSim.sbatch', 'w').write(ConstructScript(Params))
  os.system('sbatch WLSim.sbatch') 

