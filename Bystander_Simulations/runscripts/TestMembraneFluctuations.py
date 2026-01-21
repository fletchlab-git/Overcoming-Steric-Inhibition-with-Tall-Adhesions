import os
import glob


def ConstructSimulationCommand(BoxLength, NodeNums, MeshDisp, MCSweeps, kc):
  return "../MembraneFusion --NumPType 0 --kc {} --Lx {} --Ly {} --Nx {} --Ny {} --MeshRate 1.0 --PRate 0.0 --GCRate 0.0 --BiasStr 0 --MCSweeps {} --z0 50.0 --BiasCenter 0.0 --MeshDisp {} --framenums 1000".format(kc, BoxLength, BoxLength, NodeNums, NodeNums, MCSweeps, MeshDisp)

BoxLength=10.0;
NodeNums=20;
kcs = [1.0];
MeshDisp=0.18

MCSweepsEqui = 10000;MCSweepsRun=2000000

for kc in kcs:
  command = ConstructSimulationCommand(BoxLength, NodeNums, MeshDisp, MCSweepsEqui, kc)
  os.system(command);
  list_of_files = glob.glob('/Data/FletcherMembrane/Restarts/*.rest') 
  latest_file = max(list_of_files, key=os.path.getctime)
  os.system("sed -i 's/^\( *MC Sweeps: *\) [0-9]\+/\\1 {}/' {}".format(MCSweepsRun, latest_file))
  os.system("mv {} ./MyRestart.rest".format(latest_file))
  os.system("../MembraneFusion --Input MyRestart.rest")
