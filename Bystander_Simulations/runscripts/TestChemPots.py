import os
import glob
import numpy as np

def ConstructSimulationCommand(c, BoxLength, NodeNums, PRate, GCRate, MCSweepsRun):
  return "../MembraneFusion --NumPType 1 --ChemPots '{} {}' --Lx {} --Ly {} --Nx {} --Ny {} --MeshRate 0.0 --PRate {} --GCRate {} --BiasStr 0 --MCSweeps {} --MeshRate 0.0 ".format(
    c[0], c[1],
    BoxLength, BoxLength, NodeNums, NodeNums, PRate, GCRate, MCSweepsRun
    )

MCSweepsRun=100000
GCRate=0.5;
PRate=1.0-GCRate;

minprob=0.01;
maxprob=10.0;
windownums=11;

BoxLength=100.0;
NodeNums=20;

ChemPots=np.array([ np.log(np.linspace(minprob, 1.0, windownums)), np.log(np.linspace(1.0, maxprob, windownums)) ])
ChemPots=np.swapaxes(ChemPots, 0, 1);

for c in ChemPots:
  command = ConstructSimulationCommand(c, BoxLength, NodeNums, GCRate, PRate, MCSweepsRun)
  print command
  os.system(command);
  
  
  
  
  
  
  
  
#Code to restart simulation  
#list_of_files = glob.glob('/Data/FletcherMembrane/Restarts/*.rest') 
#latest_file = max(list_of_files, key=os.path.getctime)
#os.system("sed -i 's/^\( *MC Sweeps: *\) [0-9]\+/\\1 {}/' {}".format(MCSweepsRun, latest_file))
#os.system("mv {} ./MyRestart.rest".format(latest_file))
#os.system("../MembraneFusion --Input MyRestart.rest")
