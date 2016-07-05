#!/bin/bash

mkdir -p ../output/surface-3/
screen -dm ./main-surface-3 fem.prefix:../output/surface-3/ > ../output/surface-3/main.out

mkdir -p ../output/coupled-2/
screen -dm ./main-coupled-2 fem.prefix:../output/coupled-2/ > ../output/coupled-2/main.out
