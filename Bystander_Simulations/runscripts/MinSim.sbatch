#!/bin/bash -l
#SBATCH --qos=shared     #Use the shared QOS
#SBATCH --array=1-42
#SBATCH --ntasks=1
#SBATCH -t 12:00:00
#SBATCH -J Harmonic
#SBATCH -A m3012
#SBATCH --constraint=haswell

zcs=(10.0 12.0 14.0 16.0 18.0 20.0 22.0 24.0 26.0 28.0 30.0 32.0 34.0 36.0 38.0 40.0 42.0 44.0 46.0 48.0 50.0 52.0 54.0 56.0 58.0 60.0 62.0 64.0 66.0 68.0 70.0 72.0 74.0 76.0 78.0 80.0 82.0 84.0 86.0 88.0 90.0 92.0)

##Run No. 0
Index=$((${SLURM_ARRAY_TASK_ID}-1))
Tag="MyRun"
zc=${zcs[${Index}]}

args="--MCSweeps=4500000 --Tag=${Tag} --BiasStr=20.0 --Equil='No' -v 'Minimum' "
args=${args}"--MeshDisp=1.2 --z0=82.0 -w ${zc} --BC='Frame' "
args=${args}"--Lx=840 --Ly=840 --kc=20.0 --Nx 70 --Ny 70 "
args=${args}"--NumPType=3 --BindPType='10.5 0.0 4.0' "
args=${args}"--LPType='5.5 17.5 40.0 5.5 17.5 40.0' "
args=${args}"--RanPType='4.0 0.0 4.0' --Densities='200.0 6400.0 1600.0 200.0 6400.0 1600.0' "
args=${args}"--OutDir='./MyOutDir' --framenums=5 --datnums=10000"

echo ./myprogram ${args}
#./myprogram ${args}
