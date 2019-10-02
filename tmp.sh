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




###################################################################

## Catherine: compress 20k txt files.

dos2unix Non_cpg_list.txt

echo -e "--input TXT\t--env FILE\t--output ZIPPED" > noncpg.tsv


# Test
#echo -e "gs://em-scratch/Non_cpg_list.txt\tgs://em-scratch/Non_cpg_list.txt\tgs://em-scratch/Non_cpg_list.txt.gz" >> noncpg.tsv


awk 'BEGIN { FS=OFS="\t" } { if (NR!=1) print $1, $1, $1".gz" }'  Non_cpg_list.txt > noncpg.tsv

split -l 1000 --numeric-suffixes --suffix-length=2 --additional-suffix=.tsv noncpg.tsv split_noncpg

for FILE in split_noncpg*; do
    sed -i '1i --input TXT\t--env FILE\t--output ZIPPED' $FILE
done 

# Final test
head -2 split_noncpg00.tsv > test.tsv

# Working
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_GCP \
  --logging gs://$OUTPUT_B/logging/ \
  --preemptible \
  --machine-type n1-standard-4 \
  --disk-size 175 \
  --command 'gzip ${TXT} && \
             mv ${TXT}.gz ${ZIPPED} && \
             gsutil rm ${FILE}' \
  --tasks split_noncpg13.tsv \
  --wait



		
		

split_noncpg14.tsv		
split_noncpg15.tsv		



# DONE
split_noncpg00_update.tsv
  dstat --provider google-v2 --project hackensack-tyco --jobs 'gzip--emmanuel--190919-192918-47' --users 'emmanuel' --status '*'
split_noncpg00_update_list.txt

# ALL THESE FILES WERE EMPTY TO BEGIN WITH ON GOOGLE CLOUD
gzip           975  worker was terminated  2019-09-19 19:42:18
gzip           878  worker was terminated  2019-09-19 19:41:00
gzip           672  worker was terminated  2019-09-19 19:38:31
gzip           546  worker was terminated  2019-09-19 19:33:22
gzip           514  worker was terminated  2019-09-19 19:37:07
gzip           465  worker was terminated  2019-09-19 19:42:18
gzip           305  worker was terminated  2019-09-19 19:32:44
gzip           288  worker was terminated  2019-09-19 19:34:55
gzip           286  worker was terminated  2019-09-19 19:35:13
gzip            78  worker was terminated  2019-09-19 19:37:28
gzip            35  worker was terminated  2019-09-19 19:41:02
gzip            32  worker was terminated  2019-09-19 19:39:51

########
Curent : split_noncpg01.tsv
dstat --provider google-v2 --project hackensack-tyco --jobs 'gzip--emmanuel--190919-215608-09' --users 'emmanuel' --status '*'
6 empty files
########
Curent : split_noncpg02.tsv
dstat --provider google-v2 --project hackensack-tyco --jobs 'gzip--emmanuel--190919-221916-69' --users 'emmanuel' --status '*' > split_noncpg02_list.txt
cat split_noncpg02_list.txt | grep -v Success
sed -n '933p' < split_noncpg02.tsv
########
Curent : split_noncpg03.tsv
########
Curent : split_noncpg04.tsv
########
Curent : split_noncpg05.tsv
The canceled jobs were canceled by mistake. SHould check the 04 file
re-run
sed -n '1p;676p;669p;652p;648p;541p;502p;472p;455p;412p' < split_noncpg05.tsv > rerun_05.tsv
Ran successfully the re-run
########
Curent : split_noncpg06.tsv
1 error because of an empty file.
########
Curent : split_noncpg07.tsv
2 errors: empty file
########
Curent : split_noncpg08.tsv
3 errors. all empty files.
########
Curent : split_noncpg09.tsv
5 errors. all empty
########
Curent : split_noncpg10.tsv
6 failed. One empty (442)
sed -n '1p;491p;480p;422p;414p;395p' < split_noncpg10.tsv > rerun_10.tsv
rerun10: one failed.
dstat --provider google-v2 --project hackensack-tyco --jobs 'gzip--emmanuel--190920-143852-25' --users 'emmanuel' --status '*'
all completed.

#########
Current: split_noncpg18.tsv
  dstat --provider google-v2 --project hackensack-tyco --jobs 'gzip--emmanuel--190920-145528-54' --users 'emmanuel' --status '*'
Need to check the results

########
Current: split_noncpg19.tsv
  dstat --provider google-v2 --project hackensack-tyco --jobs 'gzip--emmanuel--190920-145528-54' --users 'emmanuel' --status '*'

split_noncpg20.tsv
  dstat --provider google-v2 --project hackensack-tyco --jobs 'gzip--emmanuel--190921-093905-25' --users 'emmanuel' --status '*'

split_noncpg17.tsv
  dstat --provider google-v2 --project hackensack-tyco --jobs 'gzip--emmanuel--190921-105623-17' --users 'emmanuel' --status '*'

split_noncpg16.tsv
  dstat --provider google-v2 --project hackensack-tyco --jobs 'gzip--emmanuel--190921-124928-61' --users 'emmanuel' --status '*'

split_noncpg21.tsv
  dstat --provider google-v2 --project hackensack-tyco --jobs 'gzip--emmanuel--190921-133450-90' --users 'emmanuel' --status '*'

split_noncpg05.tsv
  dstat --provider google-v2 --project hackensack-tyco --jobs 'gzip--emmanuel--190921-145720-04' --users 'emmanuel' --status '*'

split_noncpg11.tsv
  dstat --provider google-v2 --project hackensack-tyco --jobs 'gzip--emmanuel--190921-191648-93' --users 'emmanuel' --status '*'

split_noncpg12.tsv
  dstat --provider google-v2 --project hackensack-tyco --jobs 'gzip--emmanuel--190921-193302-04' --users 'emmanuel' --status '*'


# Limits
200 T @ 175G: 1,142 at the same time
10k CPU: 4


disk-size: 175

Batch of 1,000 files. 175G + 4 CPUS + preemptible