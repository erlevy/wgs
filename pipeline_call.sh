Pipeline: /storage/home.hburkhardt/kd_pipeline.py
GATK: /opt/bin/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar
Picard: /opt/bin/picard-tools-1.91/
snpEff: /opt/bin/snpEff/
Ref: /storage/refDB/hg19_reordered.fa
Fastq 1: /data/vol1/raw_fastq/PG0000???.1.fastq
Fastq 2: /data/vol1/raw_fastq/PG0000???.2.fastq
Known sites: /storage/home/hburkhardt/reorderedHg19/common_all_prefix_re-reordered.vcf

python kd_pipeline_v3.py /opt/bin/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar /opt/bin/picard-tools-1.91/ /opt/bin/snpEff/ /storage/refDB/hg19_reordered.fa  /data/vol1/raw_fastq/PG0000536.1.fastq  /data/vol1/raw_fastq/PG0000536.2.fastq /storage/refDB/dbsnp/common_all_prefix.vcf
python kd_pipeline_v3.py /opt/bin/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar /opt/bin/picard-tools-1.91/ /opt/bin/snpEff/ /storage/refDB/hg19_reordered.fa  /data/vol1/raw_fastq/PG0000536.1.fastq  /data/vol1/raw_fastq/PG0000536.2.fastq /storage/refDB/dbsnp/common_all_prefix.vcf redo_with_pipeline_v3
python kd_pipeline_v3_1.py /opt/bin/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar /opt/bin/picard-tools-1.96/ /opt/bin/snpEff/ /storage/refDB/hg19_reordered.fa  /data/vol1/raw_fastq/PG0000538.1.fastq  /data/vol1/raw_fastq/PG0000538.2.fastq /storage/refDB/dbsnp/common_all_prefix.vcf rerun-3-1