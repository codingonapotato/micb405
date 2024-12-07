#!/bin/bash

STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles GCF_013372085.1_ASM1337208v1_genomic.fna --sjdbGTFfile GCF_013372085.1_ASM1337208v1_genomic.gtf --runThreadN 8 --genomeSAindexNbases 10 --sjdbGTFfeatureExon CDS
