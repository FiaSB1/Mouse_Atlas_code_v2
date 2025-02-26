#PBS -N mouseatlas_tf
#PBS -l ncpus=30
#PBS -l mem=400gb
#PBS -l walltime=120:00:00

# Send an email when this job aborts, begins or ends.
#PBS -m abe
#PBS -M fia.boedijono@uts.edu.au

# Create a unique /scratch directory.
SCRATCH="/scratch/${USER}_${PBS_JOBID%.*}"
mkdir ${SCRATCH}
mkdir PISCES
cd PISCES
mkdir ARACNe_input
mkdir ARACNe_output
mkdir Regs

# Change to the PBS working directory where qsub was started from.
cd ${PBS_O_WORKDIR}

# Copy your input data to this scratch directory.
cp /shared/ci/kyle/macrophage.rds ${SCRATCH}
cp ~/work/R/data/scmouseCOPD/PISCES/MouseTF.txt ${SCRATCH}


# Change directory to the scratch directory and run your program.
# my_program uses input.dat and creates an output file called output.dat
cd ${SCRATCH}
Rscript ${PBS_O_WORKDIR}/TFanalysis.R

# Copy results back to your working directory.
mv ${SCRATCH} ${PBS_O_WORKDIR}/

# Clean up
cd ${PBS_O_WORKDIR}
rm -r ${SCRATCH}
