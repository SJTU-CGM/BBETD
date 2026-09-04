#!/bin/bash
#SBATCH -J BBETD_Job
#SBATCH -p cpu
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --output=%J.out
#SBATCH --error=%J.err

module load miniconda3
eval "$(conda shell.bash hook)"
conda activate bbetd


if [ ! -f ../test_data/simulated_genome.fa ];then
    wget -c https://cgm.sjtu.edu.cn/test/hzxue/simulated_genome.fa.gz
    mv simulated_genome.fa.gz ../test_data/
    zcat ../test_data/simulated_genome.fa.gz > ../test_data/simulated_genome.fa
fi
if [ ! -f ../test_data/simulated_annotation.gff3 ];then
    wget -c https://cgm.sjtu.edu.cn/test/hzxue/simulated_annotation.gff3.gz
    mv simulated_annotation.gff3.gz ../test_data/
    zcat ../test_data/simulated_annotation.gff3.gz > ../test_data/simulated_annotation.gff3
fi
if [ ! -f ../test_data/simulated_genome.fa.fai ];then
    wget -c https://cgm.sjtu.edu.cn/test/hzxue/simulated_genome.fa.fai
    mv simulated_genome.fa.fai ../test_data/
fi


bbetd \
  --genome ../test_data/simulated_genome.fa \
  --gff ../test_data/simulated_annotation.gff3 \
  --outdir ./results \
  --prefix simulated_bbetd \
  --threads 8 \
  --max-gap 5