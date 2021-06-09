export FQDIR="/home/ec2-user/sctl/Hoke.scRNAseq/FASTQ"
export REF="/home/ec2-user/sctl/references/cellranger/refdata-gex-mm10-2020-A"

cellranger count \
--id=082120DRG \
--transcriptome=${REF} \
--fastqs=${FQDIR} \
--sample=082120DRG 

cellranger count \
--id=101920DRG \
--transcriptome=${REF} \
--fastqs=${FQDIR} \
--sample=101920DRG 

cellranger count \
--id=103020DRG \
--transcriptome=${REF} \
--fastqs=${FQDIR} \
--sample=103020DRG 

cellranger count \
--id=110220DRG \
--transcriptome=${REF} \
--fastqs=${FQDIR} \
--sample=110220DRG 

cellranger count \
--id=111220DRG \
--transcriptome=${REF} \
--fastqs=${FQDIR} \
--sample=111220DRG 
