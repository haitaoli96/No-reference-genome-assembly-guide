# No-reference-genome-assembly-guide
# Illumina Genome Assembly & Decontamination Pipeline

[![Conda](https://img.shields.io/badge/Conda-env-blue.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)]()
[![Platform](https://img.shields.io/badge/Platform-Linux_x86_64-lightgrey.svg)]()
[![Status](https://img.shields.io/badge/Status-Stable-success.svg)]()

本项目提供一条龙的 **Illumina 二代测序基因组无参组装流程**，从原始 reads 自动完成：

**质控 → 去宿主 → 去细菌污染 → k-mer survey → 组装 → 抛光 → 去冗余 → QUAST/BUSCO 评估**

适用于：
- 寄生虫（Cryptosporidium / Giardia / Myxozoa）
- 鱼类寄生虫 / 腔肠动物
- 所有 5–500 Mb 的中小基因组

---

# 📦 Repository Structure

流程脚本执行后将自动生成如下目录：

```
project/
├── genome-illumina.yml                  # Conda 环境文件
├── run_spades_bowtie_kraken_cdhit_busco.sh
├── KrakenTools/                         # Git clone KrakenTools（用于去污染）
├── 00_qc/                               # fastp, fastqc, multiqc
├── 01_hostsub/                          # Bowtie2 去宿主
├── 02_kraken/                           # Kraken2 分类
├── 03_clean/                            # 按 taxid 清洗后的 clean reads
├── 04_kmer/                             # Jellyfish k-mer 分析
├── 05_asm/                              # SPAdes 组装
├── 06_polish/                           # Pilon 抛光（可选）
├── 06_cdhit/                            # CD-HIT 去冗余
├── 07_eval/                             # QUAST / BUSCO 评估
└── logs/                                # 流程日志
```

---

# 🔄 Pipeline Flowchart (Mermaid)

```mermaid
flowchart TD

A[Raw Reads
R1/R2] --> B[Step0: fastp + fastqc + multiqc]
B --> C[Step1: Bowtie2 Host Subtraction]
C --> D[Step2: Kraken2 Classification]
D --> E[Step2b: extract_kraken_reads.py
Exclude TaxIDs (Bacteria=2)]
E --> F[Step3: Jellyfish k-mer Survey]
F --> G[Step4: SPAdes Assembly
(k=21,33,55,77,99,127)]
G --> H{Pilon Polishing?}
H -->|Yes| I[Step5: Pilon Polishing
Fix snps,indels]
H -->|No| J[Skip]
I --> K[≥500bp Contig Filtering]
J --> K
K --> L[Step6: CD-HIT-EST
Redundancy Removal]
L --> M[Step7a: QUAST Assembly Evaluation]
L --> N[Step7b: BUSCO Genome Completeness]
```

---

# 🧬 1. Conda Environment (`genome-illumina.yml`)

```yaml
name: genome-illumina
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python=3.10
  - fastp
  - fastqc
  - multiqc
  - kraken2
  - bracken
  - jellyfish
  - bwa
  - samtools
  - spades
  - megahit
  - pilon
  - openjdk=17
  - quast
  - busco
  - blast
  - kmc
  - bbmap
  - pigz
  - parallel
  - git
  - wget
  - pandas
  - cd-hit
```

## Create Environment

```bash
mamba env create -f genome-illumina.yml
conda activate genome-illumina
```

---

# 🧬 2. Install KrakenTools （必需）

```bash
git clone https://github.com/jenniferlu717/KrakenTools.git
export PATH="$PWD/KrakenTools:$PATH"

echo 'export PATH="/path/to/KrakenTools:$PATH"' >> ~/.bashrc
which extract_kraken_reads.py
```

---

# 🚀 3. Main Pipeline Script

脚本：  
```
run_spades_bowtie_kraken_cdhit_busco.sh
```

功能包括：

- QC（fastp/fastqc/multiqc）
- 去宿主（Bowtie2）
- 去细菌污染（Kraken2 + extract_kraken_reads.py）
- Jellyfish k-mer 调查
- SPAdes 组装 + ≥500bp 过滤
- Pilon（可选）抛光（snps/indels）
- CD-HIT 去冗余
- QUAST / BUSCO 评估
- 支持断点续跑

---

# ▶️ 4. Usage

```bash
conda activate genome-illumina
export PATH="/path/to/KrakenTools:$PATH"

chmod +x run_spades_bowtie_kraken_cdhit_busco.sh
./run_spades_bowtie_kraken_cdhit_busco.sh
```

最终推荐使用：

```
06_cdhit/assembly.cdhit0.98.fa
```

---

# 📊 5. Output Summary

| Step | Description | Output Directory |
|------|-------------|------------------|
| 00 | QC (fastp/fastqc/multiqc) | `00_qc/` |
| 01 | Bowtie2 去宿主 | `01_hostsub/` |
| 02 | Kraken2 分类 | `02_kraken/` |
| 02b | Taxid-based filtering | `03_clean/` |
| 03 | Jellyfish k-mer survey | `04_kmer/` |
| 04 | SPAdes 组装（≥500bp） | `05_asm/` |
| 05 | Pilon 抛光（可选） | `06_polish/` |
| 06 | CD-HIT 去冗余 | `06_cdhit/` |
| 07 | QUAST / BUSCO | `07_eval/` |

---

# ⚠️ Notes

- SPAdes、Pilon、CD-HIT 对内存要求较高，建议 **≥64GB RAM**。  
- 若 Pilon 内存不足，可设置：  
  ```bash
  DO_POLISH=false
  ```
- CD-HIT 阈值：  
  - 更严格：`-c 0.99`  
  - 更宽松：`-c 0.97`
- 若样本杂合度高，可进一步使用：  
  - GenomeScope  
  - Smudgeplot

---

# 📚 Citation

请引用本流程依赖的工具：

fastp • Bowtie2 • Kraken2 • Jellyfish • SPAdes • Pilon • CD-HIT • QUAST • BUSCO

---

# 📄 License

MIT License

---

# 🎉 完成！
