# Example script: sample_script.sh
#!/bin/bash
# Set current working directory
#$ -cwd
#Use current environment variables/ modules
#$ -V
#Request one hour of runtime
#$ -l h_rt=48:00:00
#Email at the beginning and end of the job
#$ -m be
#Run the executable 'myprogram' from the current working directory

mkdir -p ../output/coupled-2/
./main-coupled-2 fem.prefix:../output/coupled-2/ > ../output/coupled-2/main.out
