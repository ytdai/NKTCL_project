####################################################################################
######## Chip-seq data analysis [run on HPC]
####################################################################################

#!/bin/bash
#SBATCH -p SVC # partition (queue)
#SBATCH --job-name=chipseq
#SBATCH -n 40
#SBATCH -t 7-00:00 # time (D-HH:MM)
#SBATCH -o _log/stat.%N.%A_%a.out # STDOUT
#SBATCH -e _log/stat.%N.%A_%a.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=XXX # send-to address

fq_path=${base_dir}
out_path=${output_dir}

id=${sample_id}

bowtie2 -x ${bowtie2_index} -1 ${fq_path}/${id}_1.fq.gz -2 ${fq_path}/${id}_2.fq.gz -X 2000 -p 50 -3 75 -S ${out_path}/${id}.bowtie2.sam
samtools sort ${out_path}/${id}.bowtie2.sam -o ${out_path}/${id}.sort.bam -@ 40
samtools index ${out_path}/${id}.sort.bam -@ 40

macs2 callpeak -t ${out_path}/${id}.sort.bam -n ${id} -q 0.2 -g hs --keep-dup all --nomodel --extsize 200 --outdir .

igvtools count -z 5 -w 20 -e 200 ${out_path}/${id}.sort.bam ${out_path}/${id}.igv.tdf ${bowtie2_index}/hg19.chrom.sizes


######### call peaks, run on compute node
macs2 callpeak -t ${id}.sort.bam -c ${id}_input.sort.bam -n ${id} -f BAM -g hs -q 1 -g hs --keep-dup all --nomodel --extsize 200 --outdir macs2_q1 &


######### generate bigwig file
bamCoverage -b ${out_path}/${id}.sort.bam -o ${out_path}/${id}.bw -bs 20 --smoothLength 0 --ignoreDuplicates --normalizeUsing RPKM -p 100


######### visualization
computeMatrix reference-point \
  -R GRCh37_genes_coding.bed \
  -S ${id}.bw ${id}_input.bw \
  -a 3000 -b 3000 \
  --binSize 20 \
  --skipZeros \
  -o computeMatrix/GRCh37_genes_coding.gz 
plotHeatmap --heatmapHeight 18 --colorMap Reds -m computeMatrix/GRCh37_genes_coding.gz -out computeMatrix/GRCh37_genes_coding.pdf


######### Find motif
findMotifsGenome.pl peak.txt hg19 motif.${id} -size 200 -len 8,10,12 &












