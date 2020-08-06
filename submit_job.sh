#!/bin/bash

#SBATCH -J ISM_Simulation
#SBATCH -n 2
#SBATCH --mem=32G
#SBATCH -t 1:00:00

module load ffmpeg/4.0.1
module load matlab/R2017b
matlab-threaded -nodisplay -r "geom_torus; exit"
ffmpeg -r 10 -i sim_%03d.png -vcodec libx264 -y -an sim_movie.mp4 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"
