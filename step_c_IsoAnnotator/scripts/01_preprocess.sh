#!/bin/bash

## This script extracts (experimentally) verified from GENCODE v43 GTF (hg38).
## Specifically, we extract entries annotated as level 1 (verified) and gene type != TEC (To be Experimentally Confirmed)
## Contact: Louis SHE (snakesch@connect.hku.hk)

GTF_PATH="../resources/gencode.v43.annotation.gtf"

grep -wv CDS "${GTF_PATH}" | grep -w 'level 1' | grep -vw TEC > "../resources/gencode.v43.annotation.no_cds.level_1.gtf"

cut -f1,4,5,3,7,9- "../resources/gencode.v43.annotation.no_cds.level_1.gtf" > "../resources/cleaned_annotations.gtf"