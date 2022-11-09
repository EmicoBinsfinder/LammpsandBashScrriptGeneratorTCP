"""
Author: Egheosa Ogbomo

Making QSUB Files (Young's) for all of the different experiments

"""

SlidingSpeeds = ['20ms', '30ms', '40ms', '50ms']
Temperatures = ['700K']
Pressures = ['1GPa', '2GPa', '3GPa', '4GPa', '5GPa']

for SlidingSpeed in SlidingSpeeds:
    for Temperature in Temperatures:
        for Pressure in Pressures:
            text = f"""
#!/bin/bash -l

# Batch script to run an MPI parallel job under SGE with Intel MPI.

# Request ten minutes of wallclock time (formathours:minutes:seconds).
#$ -l h_rt=96:00:0

# Request 1 gigabyte of RAM per process (must be an integer followed by M, G, or T)and budgets.
#$ -P GoldLong
#$ -A Imperial_MEng
#$ -l mem=1G

# Set the name of the job.
#$ -N TCP{SlidingSpeed}{Temperature}{Pressure}

# Select the MPI parallel environment.
#$ -pe mpi 120

# Set the working directory to somewhere in your scratch space.
#$ -wd /home/mmm1058/Scratch/TCP_Templates/DifferentSlidingSpeeds/{SlidingSpeed}/{Temperature}/{Pressure} 
mkdir results
cp * results

# Run our MPI job.  GERun is a wrapper that launchesMPI jobs on our clusters.

gerun /lustre/home/mmm1058/LAMMPS/lammps-install/bin/lmp -l log.lammps -in TCP.lammps
"""

            file = open(f"D:/PhD/TCPDecompositionExperiments/Templates/DifferentSlidingSpeeds/QSUBFILES/TCP{SlidingSpeed}{Temperature}{Pressure}.txt", "w")
            file.write(text)
            file.close()