#!sh
#$ -q krt,bio,pub64
#$ -t 1-250
#$ -o /dev/null
#$ -e /dev/null
module load boost/1.53.0
module load R
cd /bio/krthornt/revisions_K/Nc100.N500/s0.05/K5
if [ ! -s out.$SGE_TASK_ID.1000.scores.txt.gz ] || [ ! -s out.$SGE_TASK_ID.100.scores.txt.gz ] || [ ! -s out.$SGE_TASK_ID.500.scores.txt.gz ]
then
SEED=`echo "$SGE_TASK_ID*$RANDOM"|bc -l`
echo "$SGE_TASK_ID $SEED" > seedfile.$SGE_TASK_ID
~/src/forward_src/expevol_region_breakpoints_Ksel /kevin/krthornt/expevol_project/coalescent/macfile.$SGE_TASK_ID.gz 100 500 0.05 0.5 5 0.025 100 500 1000 25 25 summary.$SGE_TASK_ID.gz out.$SGE_TASK_ID $SEED
for NGEN in 100 500 1000
do
FN=out.$SGE_TASK_ID.$NGEN.gz
OFN=out.$SGE_TASK_ID.$NGEN.scores.txt
R --no-save --slave <<EOF
fn="$FN"
ofn="$OFN"
Nc=100
source("/kevin/krthornt/expevol_project/scores/Rcode/calc_scores.R")
EOF
if [ -e $OFN.gz ]
then
rm -f $OFN.gz
fi
gzip $OFN
rm -f $FN
done
fi
