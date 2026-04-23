conda config --add channels conda-forge
conda config --add channels bioconda # Nextflow is available via the Bioconda channel.

mamba install nextflow 
nextflow -version
# export NXF_SINGULARITY_CACHEDIR=/Work_bio/gao/projects/singularity_cache
# mkdir -p $NXF_SINGULARITY_CACHEDIR
# # 如果这个目录为空或不存在，说明镜像尚未完全下载。不用担心，下一次运行时会自动下载。
# # 首先确定你的工作目录正确不正确， 且其中有无 
# ls -lh $NXF_SINGULARITY_CACHEDIR

tmux new -s RNA  
#tmux a -t rnaseq                # 如果改环境已存在 == tmux attach -t rnaseq
#  SortMeRNA 的开发者修改或删除了 GitHub 上的库文件（v4.3.4 标签下的数据库文件），导致 nf-core/rnaseq 流程中硬编码的下载链接失效了。
#  必须升级版本，本处使用 3.15.1 版本
#  下一次似乎可以继续提高CPU 和  Memory
# SortMeRNA 报错是因为之前的版本硬编码了一个指向 GitHub master 分支的链接，而该文件最近被移动或重命名了。
# 在 3.14.1 版本中，开发者将链接指向了固定的提交（commit）地址，解决了下载失败的问题。
# 在 3.15.1 版本中，不仅修复了此问题，还包含了一些其他的性能改进和 Bug 修复。


# 重要提示：当你提供了 star_index 时，Nextflow 会自动跳过生成索引的任务，极大节省时间和内存。
# 移除了 --save_reference: 已经有了现成的索引，通常不需要再把索引文件拷贝到输出目录（节省空间）。如果你确实想在 output 里存一份，可以加回去。

# 人类基因组的 STAR Index 读入内存大约需要 30GB - 32GB。
# 如果你同时跑 2 个样本，STAR 会启动两个进程，总共消耗约 60GB - 64GB。
# 加上 Nextflow 其他中间步骤（如 FastQC, Salmon, Sortmerna）的开销，80GB 的限制刚好可以支撑 2 个样本同时进行比对。
# 建议：在这种情况下，设置 80.GB 是非常合理的。


cd /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts

# 限制 Nextflow 自身的内存开销，确保它不被 Killed
export NXF_OPTS="-Xms512m -Xmx2g"

# WARN: Singularity cache directory has not been defined -- Remote image will be stored in the path: /Work_bio/gao/projects/2026_Item9_gc/scripts/work/singularity
# -- Use the environment variable NXF_SINGULARITY_CACHEDIR to specify a different location

# 这样做的好处：

# 避免重复下载：nf-core 容器镜像只会下载一次并缓存在指定位置
# 节省磁盘空间：不会在每次运行时都在 work/singularity 目录中存储重复的镜像
# 提高性能：后续运行可以直接使用缓存的镜像，加快启动速度
# 消除警告信息：不会再看到关于 Singularity 缓存目录未定义的警告
# 修改后的脚本现在会将 Singularity 容器镜像缓存在 /home/gao/.singularity/nf-core 目录中，这是一个标准且合理的缓存位置。

# 设置 Singularity 缓存目录以避免重复下载容器镜像

export NXF_SINGULARITY_CACHEDIR="/home/gao/.singularity/nf-core"


nextflow run nf-core/rnaseq \
    -r 3.15.1 \
    -profile singularity \
    -c local_optimized.config \
    --input /home/gao/projects/2026_Item10_rRNA_removal/scripts/Sample_Sheet2.csv \
    --outdir ../output_results \
    --fasta      /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/GRCh38.primary_assembly.genome.fa \
    --gtf        /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/gencode.v45.annotation.gtf \
    --star_index /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/star_index \
    --gencode \
    --aligner star_salmon \
    --remove_ribo_rna \
    --save_non_ribo_reads \
    --max_cpus 28 \
    --max_memory '90.GB' \
#    -resume # 如果这个程序运行过，下次运行时，会自动跳过已经完成的步骤，但如果使用了错误参数，会继续使用错误参数，所有必须谨慎，确认一切正确之后，再加上这个resume 参数


