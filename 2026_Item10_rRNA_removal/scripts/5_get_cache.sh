# 检查 nf-core 缓存目录
ls -la /home/gao/.nextflow/assets/nf-core/rnaseq/assets/

# 或者检查工作目录中的 staged 文件
find /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work -name "*.fasta" | grep -i rrna


# ls -la /home/gao/.nextflow/assets/nf-core/rnaseq/assets/
# total 112
# drwxrwxr-x  2 gao gao  4096 Mar 22 15:06 .
# drwxrwxr-x 13 gao gao  4096 Mar 22 15:06 ..
# -rw-rw-r--  1 gao gao  2698 Mar 18 23:25 adaptivecard.json
# -rw-rw-r--  1 gao gao  3643 Mar 22 15:06 email_template.html
# -rw-rw-r--  1 gao gao  1459 Mar 22 15:06 email_template.txt
# -rw-rw-r--  1 gao gao 74201 Mar 22 15:06 nf-core-rnaseq_logo_light.png
# -rw-rw-r--  1 gao gao   735 Mar 18 23:25 samplesheet.csv
# -rw-rw-r--  1 gao gao  1841 Mar 18 23:25 schema_input.json
# -rw-rw-r--  1 gao gao  1155 Mar 18 23:25 sendmail_template.txt
# -rw-rw-r--  1 gao gao  1788 Mar 18 23:25 slackreport.json
# (regular_bioinfo) (base) [21:49:25] [/home/gao/projects/2026_Item10_rRNA_removal/scripts]:
# gao@us1 $ find /Wofind /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work -name "*.fasta" | grep -i rrna
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/stage-ed1842a3-b949-471e-867f-a4d8392915b6/b6/1e6751bb4aebb8f8ed8ef319412ef3/rfam-5.8s-database-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/stage-ed1842a3-b949-471e-867f-a4d8392915b6/5f/b0d594299dbfd460290c80927bb220/silva-bac-16s-id90.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/stage-ed1842a3-b949-471e-867f-a4d8392915b6/ac/249ad5603f0be3ed7f7cdcb50cc633/silva-euk-18s-id95.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/stage-ed1842a3-b949-471e-867f-a4d8392915b6/78/bcc7d42cf597d9d023782c4f4b50da/silva-euk-18s-id95.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/stage-ed1842a3-b949-471e-867f-a4d8392915b6/16/ee14d36cbaf863257e9149e82e8120/rfam-5.8s-database-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/stage-ed1842a3-b949-471e-867f-a4d8392915b6/f5/ae962fa513641b0cfd1afd416e49cf/silva-euk-28s-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/stage-ed1842a3-b949-471e-867f-a4d8392915b6/20/63ef30e0038a4b117d18cecfa3235a/rfam-5s-database-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/stage-ed1842a3-b949-471e-867f-a4d8392915b6/d2/5d905d9b2bbd21dba34c9474d6d6bf/silva-arc-23s-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/stage-ed1842a3-b949-471e-867f-a4d8392915b6/0d/c2eb20c49dd481e386e53dd843fe9e/silva-bac-23s-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/stage-ed1842a3-b949-471e-867f-a4d8392915b6/dd/190ea89673c3b239e1b5ca0131f13b/rfam-5s-database-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/stage-ed1842a3-b949-471e-867f-a4d8392915b6/75/35d901009867f2434834f1dc9bdc43/silva-arc-16s-id95.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/d2/e3ee4afac4bb296a8c7bd62eef07c4/silva-euk-28s-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/d2/e3ee4afac4bb296a8c7bd62eef07c4/silva-bac-23s-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/d2/e3ee4afac4bb296a8c7bd62eef07c4/silva-arc-23s-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/d2/e3ee4afac4bb296a8c7bd62eef07c4/rfam-5.8s-database-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/d2/e3ee4afac4bb296a8c7bd62eef07c4/silva-arc-16s-id95.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/d2/e3ee4afac4bb296a8c7bd62eef07c4/rfam-5s-database-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/d2/e3ee4afac4bb296a8c7bd62eef07c4/silva-bac-16s-id90.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/d2/e3ee4afac4bb296a8c7bd62eef07c4/silva-euk-18s-id95.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/c1/e545f7ba3f75d589a8e2edbea7f618/silva-euk-28s-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/c1/e545f7ba3f75d589a8e2edbea7f618/silva-bac-23s-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/c1/e545f7ba3f75d589a8e2edbea7f618/silva-arc-23s-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/c1/e545f7ba3f75d589a8e2edbea7f618/rfam-5.8s-database-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/c1/e545f7ba3f75d589a8e2edbea7f618/silva-arc-16s-id95.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/c1/e545f7ba3f75d589a8e2edbea7f618/rfam-5s-database-id98.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/c1/e545f7ba3f75d589a8e2edbea7f618/silva-bac-16s-id90.fasta
# /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/work/c1/e545f7ba3f75d589a8e2edbea7f618/silva-euk-18s-id95.fasta