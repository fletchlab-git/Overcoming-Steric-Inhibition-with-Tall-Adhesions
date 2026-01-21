#!/bin/bash -l
#SBATCH --qos=shared     #Use the shared QOS
#SBATCH --array=1-20
#SBATCH --ntasks=1
#SBATCH -t 2:00:00
#SBATCH -J Entropic
#SBATCH -A m3012
#SBATCH --constraint=haswell

zcs=(13.0 15.5641025641 18.1282051282 20.6923076923 23.2564102564 25.8205128205 28.3846153846 30.9487179487 33.5128205128 36.0769230769 38.641025641 41.2051282051 43.7692307692 46.3333333333 48.8974358974 51.4615384615 54.0256410256 56.5897435897 59.1538461538 61.7179487179 64.2820512821 66.8461538462 69.4102564103 71.9743589744 74.5384615385 77.1025641026 79.6666666667 82.2307692308 84.7948717949 87.358974359 89.9230769231 92.4871794872 95.0512820513 97.6153846154 100.179487179 102.743589744 105.307692308 107.871794872 110.435897436 113.0)

Index=$((${SLURM_ARRAY_TASK_ID}-1+0))
Tag="AriExperiment_percent_40.0_zc_${zcs[${Index}]}"
zc=${zcs[${Index}]}

args="--MCSweeps=4000000 --Tag=${Tag} --BiasStr=8.0 --Equil='Yes' "
args=${args}"--MeshDisp=0.4 --z0=${zc} -w ${zc} "
args=${args}"--Lx=120 --Ly=120 --kc=20.0 --Nx 30 --Ny 30 "
args=${args}"--NumPType=3 --BindPType='10.5 0.0 4.0' "
args=${args}"--LPType='5.5 17.5 40.0 5.5 17.5 40.0' "
args=${args}"--RanPType='4.0 0.0 4.0' --Densities='200.0 4800.0 3200.0 200.0 4800.0 3200.0' "
args=${args}"--OutDir='/global/homes/j/jhasnain/Data/FBMembrane/' --framenums=20 --datnums=10000"

echo ../../LatticeMembrane ${args}
time bash -c "../../LatticeMembrane ${args} --SEED=${RANDOM}"
Index=$((${SLURM_ARRAY_TASK_ID}-1+20))
Tag="AriExperiment_percent_40.0_zc_${zcs[${Index}]}"
zc=${zcs[${Index}]}

args="--MCSweeps=4000000 --Tag=${Tag} --BiasStr=8.0 --Equil='Yes' "
args=${args}"--MeshDisp=0.4 --z0=${zc} -w ${zc} "
args=${args}"--Lx=120 --Ly=120 --kc=20.0 --Nx 30 --Ny 30 "
args=${args}"--NumPType=3 --BindPType='10.5 0.0 4.0' "
args=${args}"--LPType='5.5 17.5 40.0 5.5 17.5 40.0' "
args=${args}"--RanPType='4.0 0.0 4.0' --Densities='200.0 4800.0 3200.0 200.0 4800.0 3200.0' "
args=${args}"--OutDir='/global/homes/j/jhasnain/Data/FBMembrane/' --framenums=20 --datnums=10000"

echo ../../LatticeMembrane ${args}
time bash -c "../../LatticeMembrane ${args} --SEED=${RANDOM}"
