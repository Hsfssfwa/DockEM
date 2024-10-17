#!/bin/bash
#SBATCH -J dock
#SBATCH -o dock.log
#SBATCH -e dock.e
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -t 10-00:00:00

#输入要执行的命令，do your jobs
../../DockRMSD ./af2_edock/test_0_F.mol2 lig_charge_trans.mol2