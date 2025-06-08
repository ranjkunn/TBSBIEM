#!/bin/bash
#SBATCH --job-name=GPU_Test
#SBATCH --output=./dats/gpuTest_%j.out
#SBATCH --error=./dats/gpuTest_%j.err
#SBATCH --ntasks-per-node=16 # Total CPUs
#SBATCH --time=20-23:59:59 # format HH:MM:SS
#SBATCH --mem=32000 # RAM Size
#SBATCH --partition=gpu_scholar
#SBATCH --gres=gpu:a100:1

echo "Total CPUs assigned=>"
nproc
echo "Total GPUs assigned with ID=>"
echo $CUDA_VISIBLE_DEVICES
 

# Move to the HybPiper directory
cd ./

echo "Running..."
rm ./dats/*.dat
rm ./main
# Your code below this line
nvfortran -fopenmp -fast -g -Minfo=accel -acc -gpu=managed -O3 -o main ./src/*.cuf -lcufft

if [ $? -eq 0 ]; then
    echo "Compilation successful. Running program..."
# Run the executable	
else
    echo "Compilation failed."
fi
echo  "Execution Started"


# Remove intermediate object and module files generated during the compilation
rm *.o 
rm *.mod 

# Running the executable
./main

# for i in {0001..0115}; do    gnuplot -e "set terminal jpeg; set hidden3d; set xlabel 'x1 (km)'; set ylabel 'x3 (km)'; set zlabel 'Slip (m)'; set xrange [-15:15]; set yrange [-7.5:7.5];   set zrange [0.0:0.5]; set cbrange [0.0:0.5]; set view map; splot 'Out$i.dat' u 2:3:4 ps 0.1 palette" > Slip_Top$i.jpeg; done

# ffmpeg -r 10 -i Slip_Top%04d.jpeg  -vf "fps=10" Slip_Top.mp4


echo "SLURM JOB Finished!"


