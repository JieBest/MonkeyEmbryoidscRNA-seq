#Code script for cellranger pipeline which is executed on Linux System. 
#Pipeline can be built accodrding to the tutorials provided by 10xgenomics at https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials.
#Example is taken on the basis of embryoid on D17.
#/data2/software/cellranger-6.1.2/bin/cellranger: cellranger path
#/data2/ref/10x/mfas5ncbi: reference file location
#/data2/mfasEmbryoid/raw/D17: raw fastq files location

#Aligns sequencing reads in FASTQ files to generate expression matrix files for further analysis

/data2/software/cellranger-6.1.2/bin/cellranger count --id=D17 \
                   --transcriptome=/data2/ref/10x/mfas5ncbi \
                   --fastqs=/data2/mfasEmbryoid/raw/D17 \
                   --sample=D17 \
                   --expect-cells=10000 \
                   --localcores=8 \
                   --nosecondary
