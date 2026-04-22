# 运行方法
# cd /home/gao/projects/2026_Item7_LJZ/scripts
# python 1_produce_nf-core_Samplesheet.py

import pandas as pd
import glob
import os

# 配置路径
fastq_base_dir = "/home/gao/Dropbox/P2026_04_20/Guangcan/"
# fastq_base_dir = "/home/gao/Dropbox/LZJ/01.RawData/"

original_csv = "/home/gao/projects/2026_Item7_LZJ/scripts/Analysis_LZJ.csv"
output_samplesheet = "/home/gao/projects/2026_Item7_LZJ/scripts/nf_core_samplesheet.csv"

# 读取原始表格
df = pd.read_csv(original_csv)
samples = df['Name in File'].unique() # 获取所有样本名称, ATH 使用这个固定格式
samples 

data = []
for sample in samples:
    # 构建样本目录路径
    sample_dir = os.path.join(fastq_base_dir, sample)
    
    if not os.path.exists(sample_dir):
        print(f"Warning: 找不到样本 {sample} 的目录 {sample_dir}")
        continue
    
    # 查找 R1 和 R2 文件 (支持两种命名格式)
    # 格式1: {sample}_R1.fq.gz / {sample}_R2.fq.gz
    r1_pattern1 = os.path.join(sample_dir, f"{sample}_R1.fq.gz")
    r2_pattern1 = os.path.join(sample_dir, f"{sample}_R2.fq.gz")
    
    # 格式2: {sample}_CKDL..._1.fq.gz / {sample}_CKDL..._2.fq.gz
    r1_files = glob.glob(os.path.join(sample_dir, f"{sample}_*_1.fq.gz"))
    r2_files = glob.glob(os.path.join(sample_dir, f"{sample}_*_2.fq.gz"))
    
    r1 = None
    r2 = None
    
    # 优先检查格式1
    if os.path.exists(r1_pattern1) and os.path.exists(r2_pattern1):
        r1 = r1_pattern1
        r2 = r2_pattern1
    # 否则使用格式2
    elif r1_files and r2_files:
        r1 = r1_files[0]
        r2 = r2_files[0]
    
    # 检查文件是否存在
    if r1 and r2:
        # strandedness 设为 reverse (如果是标准库)
        # strandedness 设为 auto 让流程自动检测，或者设为 reverse (如果是标准库)
        data.append([sample, r1, r2, 'reverse'])
    else:
        print(f"Warning: 找不到样本 {sample} 的 fastq 文件于 {sample_dir}")

# 生成新表格
final_df = pd.DataFrame(data, columns=['sample', 'fastq_1', 'fastq_2', 'strandedness'])
final_df.to_csv(output_samplesheet, index=False)
print(f"成功生成 nf-core 专用 Samplesheet: {output_samplesheet}")
print(f"共处理 {len(data)} 个样本")