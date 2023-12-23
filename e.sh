#!/bin/sh

cd model
CMD='mpirun -np 3 python3 KINSHIP_DYNAMICS_AND_SELECTION.py'
start=$(date '+%Y-%m-%d %H:%M:%S')
echo "\nSTART $start"
eval $CMD
end=$(date '+%Y-%m-%d  %H:%M:%S')
echo "DONE $end ($(($(date -d "$end" +%s) - $(date -d "$start" +%s))) sec)\n"
rm -rf "__pycache__"

cd ../plot
rwd="$(cd "$(dirname $0)" && pwd)"
Rscript "./Fig1_Fig2_Fig3.R" ${rwd}
Rscript "./FigA1_FigA2.R" ${rwd}