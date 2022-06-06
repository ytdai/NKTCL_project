####################################################################################
######## Single-cell RNA sequencing data analysis [run on HPC]
####################################################################################

#!/bin/bash
#SBATCH -p small # partition (queue)
#SBATCH --job-name=scRNA
#SBATCH -n 40
#SBATCH -t 7-00:00 # time (D-HH:MM)
#SBATCH -o _log/stat.%N.%A_%a.out # STDOUT
#SBATCH -e _log/stat.%N.%A_%a.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=XXX # send-to address

id=${sample_id}

fq_path=${base_dir}
index_path=${cellranger_index}/refdata-gex-GRCh38-2020-A
index_path_vdj=${cellranger_index}/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0

cellranger count --id=${id} \
  --fastqs=${fq_path}/${id} \
  --sample=${id} \
  --transcriptome=${index_path}

cellranger vdj --id=${id} \
  --fastqs=${fq_path}/${id}-B \
  --sample=${id}-B \
  --reference=${index_path_vdj} 

cellranger vdj --id=${id} \
  --fastqs=${fq_path}/${id}-T \
  --sample=${id}-T \
  --reference=${index_path_vdj} 





######## virus analysis
samtools view -@ 40 -h -bf 4 possorted_genome_bam.bam -o possorted_genome_bam.f4.bam
samtools index possorted_genome_bam.f4.bam
samtools view -@ 40  possorted_genome_bam.f4.bam | awk '{print ">"$1"_"$19"\n"$10}' | gzip - > possorted_genome_bam.f4.fasta.gz 


star_index=${human_ebv_index}
gtf_file=${human_ebv_anno_gtf_file}

STAR --runThreadN 40 --genomeDir ${star_index} \
  --readFilesIn possorted_genome_bam.f4.fasta.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix hs.V01555.2. \
  --sjdbGTFfile ${gtf_file} 
awk '/V01555/' hs.V01555.2.Aligned.out.sam | awk '$1 !~ /@/' > hs.V01555.2.extract.sam









