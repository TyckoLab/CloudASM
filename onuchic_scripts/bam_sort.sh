#!/bin/bash


PICARD="/genomics-packages/picard-tools-1.97"

TMP_DIR="/mnt/data/tmp/"
mkdir -p ${TMP_DIR}
             
java \
    -Djava.io.tmpdir=${TMP_DIR} \
    -Xmx52g \
    -jar $PICARD/SortSam.jar \
    I=${BAM} \
    O=${BAM_SORTED} \
    SORT_ORDER=coordinate \
    MAX_RECORDS_IN_RAM=9000000