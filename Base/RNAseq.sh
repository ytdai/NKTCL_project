####################################################################################
######## RNA-seq data analysis [run on HPC]
####################################################################################


#!/bin/bash
#SBATCH -p small # partition (queue)
#SBATCH --job-name=rnaseq
#SBATCH -n 8
#SBATCH -t 7-00:00 # time (D-HH:MM)
#SBATCH -o _log/stat.%N.%A_%a.out # STDOUT
#SBATCH -e _log/stat.%N.%A_%a.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=XXX # send-to address

fq_path=${base_dir}
out_path=${output_dir}

id=${sample_id}

fq1=${fq_path}/${id}_R1.fq.gz
fq2=${fq_path}/${id}_R2.fq.gz
gtf_file=${salmon_index}/gencode.v37.annotation.gtf

salmon quant -p 10 -l IU -i ${salmon_index} -o ${out_path}/${id} -1 ${fq1} -2 ${fq2} -g ${gtf_file} --gcBias --validateMappings





#!/bin/bash
#SBATCH -p small # partition (queue)
#SBATCH --job-name=rnaseq
#SBATCH -n 40
#SBATCH -t 7-00:00 # time (D-HH:MM)
#SBATCH -o _log/stat.%N.%A_%a.out # STDOUT
#SBATCH -e _log/stat.%N.%A_%a.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=XXX # send-to address

fq_path=${base_dir}
out_path=${output_dir}

id=${sample_id}

fq1=${fq_path}/${id}_R1.fq.gz
fq2=${fq_path}/${id}_R2.fq.gz
star_index=${human_ebv_star_index}
gtf_file=${human_ebv_anno_gtf}

echo ${id}
mkdir ${out_path}/${id}
STAR --runThreadN 40 --genomeDir ${star_index} \
  --readFilesIn ${fq1} ${fq2} \
  --readFilesCommand zcat \
  --outFileNamePrefix ${out_path}/${id}/${id}. \
  --sjdbGTFfile ${gtf_file} 

samtools view -@ 40 -bF -h ${out_path}/${id}/${id}.Aligned.out.sam -o ${out_path}/${id}/${id}.Aligned.out.bam

samtools sort -@ 40 -n ${out_path}/${id}/${id}.Aligned.out.bam -o ${out_path}/${id}/${id}.sort_by_name.bam

samtools sort -@ 40 ${out_path}/${id}/${id}.Aligned.out.bam -o ${out_path}/${id}/${id}.sort.bam

samtools index ${out_path}/${id}/${id}.sort.bam

rm ${out_path}/${id}/${id}.Aligned.out.sam ${out_path}/${id}/${id}.Aligned.out.bam





