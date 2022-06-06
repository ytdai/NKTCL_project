####################################################################################
######## WGS/WES data analysis [run on HPC]
####################################################################################

#!/bin/bash
#SBATCH -p small # partition (queue)
#SBATCH --job-name=WXSmapping
#SBATCH -n 40
#SBATCH -t 7-00:00 # time (D-HH:MM)
#SBATCH -o _log/stat.%N.%A_%a.out # STDOUT
#SBATCH -e _log/stat.%N.%A_%a.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=XXX # send-to address

fq_path=${base_dir}
out_path=${output_dir}

id=${sample_id}

fq1=${fq_path}/${id}_1.fq.gz
fq2=${fq_path}/${id}_2.fq.gz
index_path=${bwa_index}
fa=${human_reference}

echo ${id}
mkdir ${out_path}/${id}
bwa mem -t 40 -M ${fa} ${fq1} ${fq2} > ${out_path}/${id}/${id}.sam

samtools view -@ 40 -bF -h ${out_path}/${id}/${id}.sam -o ${out_path}/${id}/${id}.bam

samtools sort -@ 40 ${out_path}/${id}/${id}.bam -o ${out_path}/${id}/${id}.sort.bam

samtools index ${out_path}/${id}/${id}.sort.bam

rm ${out_path}/${id}/${id}.bam ${out_path}/${id}/${id}.sam


dbsnp=${db_dir}/dbsnp_138.hg19.vcf
dbindel=${db_dir}/1000G_phase1.indels.hg19.sites.vcf
dbsnp_1000G=${db_dir}/1000G_phase1.snps.high_confidence.hg19.sites.vcf

# AddGroup
java -jar $PICARD AddOrReplaceReadGroups I=${out_path}/${id}/${id}.sort.bam O=${out_path}/${id}/${id}.AddGroup.bam RGID=1 RGLB=NKTCL RGPL=ILLUMINA RGPU=novaseq RGSM=${id}

# MarkDuplicates
java -jar $PICARD MarkDuplicates I=${out_path}/${id}/${id}.AddGroup.bam O=${out_path}/${id}/${id}.MarkDup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${out_path}/${id}/${id}.MarkDup.bam.metrics

# BaseRecalibrator
gatk --java-options "-Xmx128G" BaseRecalibrator -I ${out_path}/${id}/${id}.MarkDup.bam -R ${fa} -O ${out_path}/${id}/${id}.recal_data.grp -known-sites ${dbindel} -known-sites ${dbsnp_1000G} -known-sites ${dbsnp}

# ApplyBQSR
gatk --java-options "-Xmx128G" ApplyBQSR --add-output-sam-program-record -R ${fa} -I ${out_path}/${id}/${id}.MarkDup.bam --use-original-qualities -O ${out_path}/${id}/${id}.BQSR.bam --bqsr-recal-file ${out_path}/${id}/${id}.recal_data.grp 

# HaplotypeCaller
gatk --java-options "-Xmx128G" HaplotypeCaller -R ${fa} -I ${out_path}/${id}/${id}.BQSR.bam -O ${out_path}/${id}/${id}.vcf.gz --dbsnp ${dbsnp} --dont-use-soft-clipped-bases







#!/bin/bash
#SBATCH -p small # partition (queue)
#SBATCH --job-name=EBVanalysis
#SBATCH -n 10
#SBATCH -t 7-00:00 # time (D-HH:MM)
#SBATCH -o _log/stat.%N.%A_%a.out # STDOUT
#SBATCH -e _log/stat.%N.%A_%a.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=XXX # send-to address

source activate R

fq_path=${base_dir}
out_path=${output_dir}

id=${sample_id}

ebv_index=${ebv_reference}
hsebv_index=${human_ebv_reference}

samtools view -@ 8 -h -bf 4 ${bam_path}/${id}/${id}.MarkDup.bam -o ${out_path}/${id}.f4.bam
samtools index ${out_path}/${id}.f4.bam

samtools view -@ 8 ${out_path}/${id}.f4.bam | awk '{print $1}' > ${out_path}/${id}.f4.name.txt
sort -n ${out_path}/${id}.f4.name.txt | uniq > ${out_path}/${id}.f4.name.uniq.txt

java -jar $PICARD FilterSamReads FILTER=includeReadList I=${bam_path}/${id}/${id}.MarkDup.bam RLF=${out_path}/${id}.f4.name.uniq.txt WRITE_READS_FILES=False O=${out_path}/${id}.unmapped.bam

java -jar $PICARD SamToFastq I=${out_path}/${id}.unmapped.bam F=${out_path}/${id}.unmapped_1.fq F2=${out_path}/${id}.unmapped_2.fq

gzip ${out_path}/${id}.unmapped_1.fq
gzip ${out_path}/${id}.unmapped_2.fq


######### map virus
bwa mem -t 5 ${ebv_index} ${out_path}/${id}.unmapped_1.fq.gz ${out_path}/${id}.unmapped_2.fq.gz > ${out_path}/${id}.ebv.sam
samtools view -@ 5 -b -F 4 -h ${out_path}/${id}.ebv.sam -o ${out_path}/${id}.unsort.ebv.bam
samtools sort -@ 5 ${out_path}/${id}.unsort.ebv.bam -o ${out_path}/${id}.ebv.bam
samtools index ${out_path}/${id}.ebv.bam
rm ${out_path}/${id}.ebv.sam ${out_path}/${id}.unsort.ebv.bam
samtools depth {out_path}/${id}.ebv.bam > {out_path}/${id}.ebv.depth.txt


######### map human and virus
bwa mem -t 5 ${hsebv_index} ${out_path}/${id}.unmapped_1.fq.gz ${out_path}/${id}.unmapped_2.fq.gz > ${out_path}/${id}.hsebv.sam
samtools view -h -Sb ${out_path}/${id}.hsebv.sam -F 12 -o ${out_path}/${id}.hsebv.unsort.bam
samtools sort -@ 5 ${out_path}/${id}.hsebv.unsort.bam -o ${out_path}/${id}.hsebv.bam
samtools index ${out_path}/${id}.hsebv.bam
rm ${out_path}/${id}.hsebv.unsort.bam ${out_path}/${id}.hsebv.sam
samtools view -h ${out_path}/${id}.hsebv.bam > ${out_path}/${id}.hsebv.sam


awk '$7 !~ /=/' ${out_path}/${id}.hsebv.sam > ${out_path}/${id}.hsebv.RP.sam
awk '/V01555/' ${out_path}/${id}.hsebv.RP.sam | awk '$1 !~ /@/' > ${out_path}/${id}.hsebv.RP.extract.sam




#!/bin/bash
#SBATCH -p small # partition (queue)
#SBATCH --job-name=CNV
#SBATCH -n 10
#SBATCH -t 7-00:00 # time (D-HH:MM)
#SBATCH -o _log/stat.%N.%A_%a.out # STDOUT
#SBATCH -e _log/stat.%N.%A_%a.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=XXX # send-to address

source activate gistic
tumor_bam=`cat nktcl_tumor.txt`
normal_bam=`cat nktcl_normal.txt`
cnvkit.py batch ${tumor_bam} --normal ${normal_bam} -p 40 \
    --method wgs \
    --annotate ${db_dir}/hg19.refFlat.txt \
    --fasta ${db_dir}/ucsc.hg19.fasta \
    --access ${db_dir}/access-5k-mappable.hg19.bed \
    --output-reference nktcl.wgs.cnvkit.cnn --output-dir nktcl.wgs.cnvkit/ \
    --diagram --scatter

tumor_bam=`cat nktcl_all_tumor.txt`
cnvkit.py batch ${tumor_bam} -p 40 \
    --method wgs \
    -r nktcl.wgs.cnvkit.cnn \
    --output-dir nktcl.wgs.all.tumor \
    --diagram --scatter



##################  GISTIC
echo --- running GISTIC ---
## input file definitions
basedir=nktcl.wes.all.tumor
segfile=gistic.input.NKTCL.txt
refgenefile=hg19.mat
## call script that sets MCR environment and calls GISTIC executable 
gistic2 -b $basedir -seg $segfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.98 -conf 0.9 -armpeel 1 -savegene 1 -gcm extreme -ta 0.1 -td 0.1 -cap 1.5 -maxseg 4000 -js 4
















