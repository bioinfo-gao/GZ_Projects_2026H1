conda config --add channels conda-forge
conda config --add channels bioconda # Nextflow is available via the Bioconda channel.

mamber install nextflow 
nextflow -version
# export NXF_SINGULARITY_CACHEDIR=/Work_bio/gao/projects/singularity_cache
# mkdir -p $NXF_SINGULARITY_CACHEDIR
# # 如果这个目录为空或不存在，说明镜像尚未完全下载。不用担心，下一次运行时会自动下载。
# # 首先确定你的工作目录正确不正确， 且其中有无 
# ls -lh $NXF_SINGULARITY_CACHEDIR

tmux new -s rnaseq  
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

nextflow run nf-core/rnaseq \
    -r 3.15.1 \
    -profile singularity \
    -resume \
    --input /home/gao/projects/2026_Item10_rRNA_removal/scripts/Sample_Sheet2.csv \
    --outdir /home/gao/projects/2026_Item10_rRNA_removal/output \
    --star_index /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/star_index \
    --fasta /Work_bio/references/Homo_sapiens/GRCh38/GENCODE/human_gencode_v45/GRCh38.primary_assembly.genome.fa \
    --gtf /Work_bio/references/Homo_sapiens/GRCh38/GENCODE/human_gencode_v45/gencode.v45.annotation.gtf \
    --gencode \
    --aligner star_salmon \
    --remove_ribo_rna \
    --save_non_ribo_reads \
    --max_cpus 32 \
    --max_memory '80.GB'

# 修改进程资源分配（高级用法）：
# 如果发现任务因为申请不到内存而排队，你可以在当前目录下创建一个名为 local.config 的文件，内容如下：

# process {
#     withName: 'STAR_ALIGN' {
#         memory = '35 GB'
#         cpus = 12
#     }
# }