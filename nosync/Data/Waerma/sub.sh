#! /bin/bash
#
#SBATCH --job-name=xj
#SBATCH --ntasks=1
#SBATCH --output=slurm_%j.out
echo cp ~/SHUD/shud .
cp ~/SHUD/shud .
chmod +x shud
./shud waerma

