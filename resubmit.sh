#!/bin/bash
#SBATCH -J gizmo
#SBATCH -t 24:00:00
#SBATCH -N 6 -n 144

export OMP_NUM_THREADS=2
mpiexec -np 144 ../gizmo/GIZMO ./params.txt 1 1>gizmo.out 2>gizmo.err

#aprun -S 4 -cc numa_node -N 16 -d 2 -n 16 ./GIZMO ./params.txt 1>gizmo.out 2>gizmo.err

