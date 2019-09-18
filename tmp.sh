# Upload the database of SNPs to Big Query. (the first 57 rows are part of the header)
bq --location=US load \
               --replace=true \
               --source_format=CSV \
               --skip_leading_rows=57 \
               --field_delimiter "\t" \
               ${DATASET_ID}."grch38_db151" \
               gs://${REF_DATA_B}/grch38/variants/All_20180418.vcf \
               chr:STRING,pos:INTEGER,snp_id:STRING,ref:STRING,alt:STRING,qual:STRING,filter:STRING,info:STRING


# Upload the database of SNPs to Big Query. (the first 57 rows are part of the header)
bq --location=US load \
               --replace=true \
               --source_format=CSV \
               --skip_leading_rows=57 \
               --field_delimiter "\t" \
               ${DATASET_ID}."grch37_db150" \
               gs://${REF_DATA_B}/dbSNP150_grc37_GATK/no_chr_dbSNP150_GRCh37.vcf \
               chr:STRING,pos:INTEGER,snp_id:STRING,ref:STRING,alt:STRING,qual:STRING,filter:STRING,info:STRING
