#!/bin/bash
#SBATCH --job-name="Syn"               # Job name
#SBATCH -n 32                    # Number of cores
#SBATCH -N 1                     # Number of nodes
#SBATCH -o slurmlog_%j.out       # File to which STDOUT will be written
#SBATCH -e slurmlog_%j.err       # File to which STDERR will be written
#SBATCH --gres=gpu:4             # require 4 gpus
#SBATCH -p faster

##SBATCH --nodelist=node[056]


## create hostfile for mpirun
#scontrol show hostnames $SLURM_JOB_NODELIST|cat | awk ' { print $1, str} ' str="slots=10" > hostfile_$SLURM_JOB_ID
#env > fullenv
#gpuid=`env | grep CUDA_VISIBLE_DEVICES`
#echo "gpuid is $gpuid"


module unload openmpi/gcc/64/1.10.1
module unload cuda80/toolkit/8.0.61

export PATH=/home/zgjia/Software/openmpi/306/bin/:$PATH:
export  LD_LIBRARY_PATH=/home/zgjia/Software/openmpi/306/lib:$LD_LIBRARY_PATH:
export PATH=$PATH:/home/zgjia/Software/gromacs/plumed262/bin/
export  LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/zgjia/Software/gromacs/plumed262/lib/

module load cuda/10.1.243

 export PATH=/home/zgjia/Software/gcc/7.5/bin/:$PATH
 export LD_LIBRARY_PATH=/home/zgjia/Software/gcc/7.5/lib/:/home/zgjia/Software/gcc/7.5/lib64/:/home/zgjia/Software/gcc/7.5/libexec/:$LD_LIBRARY_PATH

 export PATH=/home/zgjia/Software/gromacs/gromacs2020_plumed262/bin/:$PATH
 export LD_LIBRARY_PATH=/home/zgjia/Software/gromacs/gromacs2020_plumed262/lib64/:$LD_LIBRARY_PATH


###  job  ** after each step, regenrate the mdp and  continue

  gmx_mpi grompp -f Protein_0-1000ns_1.0_steer.mdp -o  Protein_0-1000ns_2020.tpr -c initial.gro  -r  initial.gro -n md.ndx -p Protein.top

  gmx_mpi mdrun   -s Protein_0-1000ns_2020.tpr  -deffnm  Protein_0-1000ns -ntomp 16 -nb gpu  -pme gpu -bonded gpu  -gpu_id 0 -pin on  -pinoffset 0 # -cpi Protein_0-1000ns.cpt -append

  gmx_mpi grompp -f Protein_0-1000ns_1.0.mdp -o  Protein_0-1000ns_2020.tpr -c Protein_0-1000ns.gro  -r initial.gro -n md.ndx -p Protein.top

  gmx_mpi mdrun   -s Protein_0-1000ns_2020.tpr  -deffnm  Protein_0-1000ns -ntomp 16 -nb gpu  -pme gpu -bonded gpu  -gpu_id 0 -pin on  -pinoffset 0 -cpi Protein_0-1000ns.cpt -append -px Protein_0-1000ns_pullx.xvg -pf Protein_0-1000ns_pullf.xvg


















