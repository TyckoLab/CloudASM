########################## Simple test for dsub ################################

# This should not lead to any error
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --regions $REGION_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --output OUT=gs://$OUTPUT_B/output/out.txt \
  --command 'echo "Hello World" > "${OUT}"' \
  --wait

dsub \
  --provider local \
  --logging $WD/logging/ \
  --output $WD/out.txt \
  --command 'echo "Hello World" > "${OUT}"' \
  --wait

# Delete the output
gsutil rm -r gs://$OUTPUT_B/output



########################## Make an extract on ENCODE to test the pipeline ################################

# Type this in extract.tsv

--input ZIPPED	--output EXTRACT
gs://encode-wgbs/A549/ENCFF327QCK.fastq.gz	gs://encode-wgbs/A549-extract/ENCFF327QCK.fastq.gz
gs://encode-wgbs/A549/ENCFF986UWM.fastq.gz	gs://encode-wgbs/A549-extract/ENCFF986UWM.fastq.gz
gs://encode-wgbs/gm12878/ENCFF798RSS.fastq.gz	gs://encode-wgbs/gm12878-extract/ENCFF798RSS.fastq.gz
gs://encode-wgbs/gm12878/ENCFF113KRQ.fastq.gz	gs://encode-wgbs/gm12878-extract/ENCFF113KRQ.fastq.gz

dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --logging gs://$OUTPUT_B/logging/ \
  --disk-size 800 \
  --preemptible \
  --image $DOCKER_IMAGE \
  --command 'gunzip ${ZIPPED} && \
             head -12000000 ${ZIPPED%.gz} > ${EXTRACT%.gz} && \
             gzip ${EXTRACT%.gz}' \
  --tasks extract.tsv \
  --wait