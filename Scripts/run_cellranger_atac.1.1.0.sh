#!/bin/bash
#$ -S /bin/bash
# Version, 1.1.0
cellranger-atac count --id=NA104 --reference=refdata-cellranger-atac-mm10-1.1.0 --fastqs=NA104 --sample=NA104 --localcores=16 --localmem=15
