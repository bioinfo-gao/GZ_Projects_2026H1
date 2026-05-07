# conda config --add channels conda-forge
# conda config --add channels bioconda # Nextflow is available via the Bioconda channel.
# mamba install nextflow 


# Check if rRNA databases exist locally
RIBO_MANIFEST="/Work_bio/gao/projects/2026_Item12_rRNAal/scripts/rRNA_databases/sortmerna_database_manifest.txt"
if [ ! -f "$RIBO_MANIFEST" ]; then
    echo "ERROR: rRNA database manifest not found at $RIBO_MANIFEST!"
    echo "Please run download_rRNA_databases.sh on a machine with internet access,"
    echo "then transfer the rRNA_databases folder to this directory."
    exit 1
fi
echo "All rRNA database files verified. Starting pipeline..."

nextflow -version
tmux new -s RNA5  
# tmux a

cd /Work_bio/gao/projects/2026_Item12_rRNAal/scripts

# 限制 Nextflow 自身的内存开销，确保它不被 Killed
export NXF_OPTS="-Xms512m -Xmx2g"
# 设置 Singularity 缓存目录以避免重复下载容器镜像
export NXF_SINGULARITY_CACHEDIR="/home/gao/.singularity/nf-core"

# Check if rRNA databases exist locally
RIBO_MANIFEST="/Work_bio/gao/projects/2026_Item12_rRNAal/scripts/rRNA_databases/sortmerna_database_manifest.txt"

# 2026-05-05 添加  改回同时跑两个样本，MEM 90G

nextflow run nf-core/rnaseq \
    -r 3.15.1 \
    -profile singularity \
    -c local_optimized.config \
    -c avoid_download.config \
    --input /home/gao/projects/2026_Item12_rRNAal/scripts/Samples.csv \
    --outdir ../output_results \
    --fasta /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/GRCh38.primary_assembly.genome.fa \
    --gtf /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/gencode.v45.annotation.gtf \
    --star_index /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/star_index \
    --gencode \
    --aligner star_salmon \
    --remove_ribo_rna \
    --ribo_database_manifest "$RIBO_MANIFEST" \
    --save_non_ribo_reads \
    --max_cpus 28 \
    --max_memory '90.GB' \
    -resume  #参数         # 暂时移除 -resume # 参数, 还有上一行 末尾的 \
    