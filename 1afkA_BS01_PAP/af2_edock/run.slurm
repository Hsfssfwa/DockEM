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
../../../edock_ligand_flex_2023new/edock_flexible af2_protein.mol2 super_lig_af2.mol2 bindingsite_pre.dat test1225 rec_keypoint_pre rec_keypoint_pre.sph 1afkA_BS01_PAP_map.mrc res 8.21 sam 0