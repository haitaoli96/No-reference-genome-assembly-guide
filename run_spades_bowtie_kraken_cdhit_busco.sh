#!/usr/bin/env bash
set -euo pipefail

########################################
# 0. 用户可修改参数
########################################

# ===== 输入 reads（原始或已初步clean都可以）=====
R1="QDLSDPC-DNA_FDSW190647919-1a_1.clean.fq"
R2="QDLSDPC-DNA_FDSW190647919-1a_2.clean.fq"

# ===== 宿主参考基因组（Bowtie2 去宿主）=====
HOST_FASTA="/data/ab1/genome_asemmble/M_cephalus/M_cephalus.fna"

# ===== Kraken2 数据库路径 =====
KRAKEN_DB="/data/ab1/Kraken2_db/standard_db"

# ===== 要排除的 taxid（细菌=2）=====
EXCLUDE_TAXIDS="2"

# ===== 线程 / 内存等 =====
THREADS=48
SPADES_MEM_GB=250          # SPAdes 内存（GB）
JAVA_MEM="128g"            # Pilon -Xmx 内存
JELLY_SIZE="1G"            # Jellyfish hash size
K_LIST="21,33,55,77,99,127"

# ===== 是否执行可选步骤 =====
DO_POLISH=true             # 是否执行 Pilon 抛光
RUN_BUSCO=true             # 是否执行 BUSCO
BUSCO_LINEAGE="/data/ab1/busco_db/eukaryota_odb10"

########################################
# 1. 环境检查 & 目录准备
########################################
need(){ command -v "$1" >/dev/null 2>&1 || { echo "缺少命令：$1"; exit 1; }; }

for c in fastp fastqc multiqc bowtie2-build bowtie2 kraken2 extract_kraken_reads.py \
         jellyfish bwa samtools spades.py quast.py pilon busco; do
  need "$c"
done

# cd-hit 有两个常见名字：cd-hit-est / cdhit_est
CDHIT_BIN=""
if command -v cd-hit-est >/dev/null 2>&1; then
  CDHIT_BIN="cd-hit-est"
elif command -v cdhit_est >/dev/null 2>&1; then
  CDHIT_BIN="cdhit_est"
else
  echo "缺少命令：cd-hit-est 或 cdhit_est"
  exit 1
fi

mkdir -p 00_qc 01_hostsub 02_kraken 03_clean 04_kmer 05_asm 06_polish 06_cdhit 07_eval logs

log(){ echo "[$(date '+%F %T')] $*" | tee -a logs/run.log; }

# 统计 fq / fq.gz reads 数
count_reads(){
  local f="$1"
  local n
  if [[ "$f" == *.gz ]]; then
    n=$(zcat "$f" | wc -l)
  else
    n=$(wc -l < "$f")
  fi
  echo $(( n / 4 ))
}

########################################
# Step0: fastp + fastqc + multiqc
########################################
if [[ ! -s 00_qc/R1.clean.fq.gz || ! -s 00_qc/R2.clean.fq.gz ]]; then
  log "Step0: fastp 质控"

  RAW_R1_READS=$(count_reads "$R1")
  RAW_R2_READS=$(count_reads "$R2")
  {
    echo -e "Raw_R1_reads\t${RAW_R1_READS}"
    echo -e "Raw_R2_reads\t${RAW_R2_READS}"
  } > 00_qc/raw_read_counts.txt

  fastp -i "$R1" -I "$R2" \
        -o 00_qc/R1.clean.fq.gz \
        -O 00_qc/R2.clean.fq.gz \
        -h 00_qc/fastp.html \
        -j 00_qc/fastp.json \
        -w "$THREADS"

  CLEAN_R1_READS=$(count_reads 00_qc/R1.clean.fq.gz)
  CLEAN_R2_READS=$(count_reads 00_qc/R2.clean.fq.gz)
  {
    echo -e "Clean_R1_reads\t${CLEAN_R1_READS}"
    echo -e "Clean_R2_reads\t${CLEAN_R2_READS}"
  } > 00_qc/clean_read_counts.txt

  {
    echo -e "Type\tR1_reads\tR2_reads"
    echo -e "Raw\t${RAW_R1_READS}\t${RAW_R2_READS}"
    echo -e "After_fastp\t${CLEAN_R1_READS}\t${CLEAN_R2_READS}"
  } > 00_qc/read_counts_summary.txt

  fastqc -t "$THREADS" 00_qc/R1.clean.fq.gz 00_qc/R2.clean.fq.gz || true
  multiqc -o 00_qc 00_qc || true

else
  log "Step0: 已存在 00_qc/R1.clean.fq.gz / R2.clean.fq.gz，跳过 fastp"
fi

R1C=00_qc/R1.clean.fq.gz
R2C=00_qc/R2.clean.fq.gz

########################################
# Step1: Bowtie2 去宿主
########################################
R1HS=01_hostsub/hostsub.clean.R1.fq.gz
R2HS=01_hostsub/hostsub.clean.R2.fq.gz

if [[ ! -s "$R1HS" || ! -s "$R2HS" ]]; then
  log "Step1: Bowtie2 去宿主"

  [[ -s "$HOST_FASTA" ]] || { echo "HOST_FASTA 不存在：$HOST_FASTA"; exit 1; }

  if [[ ! -s 01_hostsub/hostidx.1.bt2 ]]; then
    log "  构建宿主 Bowtie2 索引"
    bowtie2-build "$HOST_FASTA" 01_hostsub/hostidx 1>>logs/run.log 2>&1
  fi

  HS_BEFORE_R1=$(count_reads "$R1C")
  HS_BEFORE_R2=$(count_reads "$R2C")

  bowtie2 -x 01_hostsub/hostidx \
          -1 "$R1C" -2 "$R2C" \
          -p "$THREADS" \
          --very-sensitive \
          --no-unal \
          --un-conc-gz 01_hostsub/hostsub.clean.R%.fq.gz \
          -S /dev/null \
          1>>logs/run.log \
          2>01_hostsub/bowtie2_host.log

  HS_AFTER_R1=$(count_reads "$R1HS")
  HS_AFTER_R2=$(count_reads "$R2HS")

  HS_REMOVED=$(( HS_BEFORE_R1 - HS_AFTER_R1 ))
  if [[ "$HS_BEFORE_R1" -gt 0 ]]; then
    HS_PCT=$(awk -v b="$HS_BEFORE_R1" -v r="$HS_REMOVED" 'BEGIN{printf "%.2f", r*100.0/b}')
  else
    HS_PCT="0.00"
  fi

  {
    echo -e "Bowtie2 host subtraction summary"
    echo -e "Before_R1_reads\t${HS_BEFORE_R1}"
    echo -e "After_R1_reads\t${HS_AFTER_R1}"
    echo -e "Removed_R1_reads\t${HS_REMOVED}"
    echo -e "Host_contamination_pct_R1\t${HS_PCT}"
    echo
    echo "原始 Bowtie2 日志：01_hostsub/bowtie2_host.log"
  } > 01_hostsub/hostsub_stats.txt

  log "  去宿主完成：R1 ${HS_BEFORE_R1} → ${HS_AFTER_R1}，移除 ${HS_REMOVED}（约 ${HS_PCT}%）"

else
  log "Step1: 去宿主结果已存在，跳过"
fi

########################################
# Step2: Kraken2 分类（在去宿主后的 reads 上）
########################################
if [[ ! -s 02_kraken/kraken.out || ! -s 02_kraken/kraken.report ]]; then
  log "Step2: Kraken2 分类（快速看细菌污染，给后面去细菌做准备）"

  kraken2 \
    --db "$KRAKEN_DB" \
    --threads "$THREADS" \
    --paired \
    --report 02_kraken/kraken.report \
    --output 02_kraken/kraken.out \
    --classified-out 02_kraken/classified#.fq \
    --unclassified-out 02_kraken/unclassified#.fq \
    "$R1HS" "$R2HS"

  # 生成 kraken 总结：总 reads、unclassified vs classified、Top20
  awk 'NR==1{
          un_pct=$1; un=$2;
          total=int(un*100/un_pct+0.5);
          class=total-un;
          class_pct=100-un_pct;
          printf "Total_fragments\t%d\nUnclassified_pct\t%.2f\nUnclassified_fragments\t%d\nClassified_pct\t%.2f\nClassified_fragments\t%d\n\n",
                 total, un_pct, un, class_pct, class;
          print "Top20 taxa from kraken.report:";
       }
       NR<=21{print $0}' 02_kraken/kraken.report > 02_kraken/kraken_summary.txt

  log "  Kraken2 完成，概要见 02_kraken/kraken_summary.txt"
else
  log "Step2: Kraken2 结果已存在，跳过"
fi

########################################
# Step2b: Kraken2 去细菌（taxid=2）→ clean_reads
########################################
CR1=03_clean/clean_R1.fq
CR2=03_clean/clean_R2.fq

if [[ ! -s "$CR1" || ! -s "$CR2" ]]; then
  log "Step2b: 利用 extract_kraken_reads.py 删除所有细菌 reads（taxid=2）"

  readarray -t TAX_ARR < <(echo "$EXCLUDE_TAXIDS" | tr ',' '\n' | sed 's/^[ \t]*//;s/[ \t]*$//')

  extract_kraken_reads.py \
    -k 02_kraken/kraken.out \
    -r 02_kraken/kraken.report \
    -s "$R1HS" \
    -s2 "$R2HS" \
    -t "${TAX_ARR[@]}" \
    --exclude \
    --include-children \
    --fastq-output \
    -o "$CR1" \
    -o2 "$CR2"

  KR_BEFORE_R1=$(count_reads "$R1HS")
  KR_AFTER_R1=$(count_reads "$CR1")
  KR_REMOVED=$(( KR_BEFORE_R1 - KR_AFTER_R1 ))
  if [[ "$KR_BEFORE_R1" -gt 0 ]]; then
    KR_PCT=$(awk -v b="$KR_BEFORE_R1" -v r="$KR_REMOVED" 'BEGIN{printf "%.2f", r*100.0/b}')
  else
    KR_PCT="0.00"
  fi

  {
    echo -e "Kraken-based bacterial decontamination summary"
    echo -e "Before_R1_reads\t${KR_BEFORE_R1}"
    echo -e "After_R1_reads\t${KR_AFTER_R1}"
    echo -e "Removed_R1_reads\t${KR_REMOVED}"
    echo -e "Kraken_removed_pct_R1\t${KR_PCT}"
    echo -e "Excluded_taxids\t${EXCLUDE_TAXIDS} (Bacteria)"
  } > 03_clean/kraken_decontam_stats.txt

  log "  Kraken 去细菌完成：R1 ${KR_BEFORE_R1} → ${KR_AFTER_R1}，移除 ${KR_REMOVED}（约 ${KR_PCT}%）"
else
  log "Step2b: clean reads 已存在，跳过"
fi

########################################
# Step3: Jellyfish k-mer survey（用 clean_reads）
########################################
if [[ ! -s 04_kmer/clean.histo ]]; then
  log "Step3: Jellyfish k-mer (k=21)"

  mkdir -p 04_kmer
  jellyfish count -C -m 21 -s "$JELLY_SIZE" -t "$THREADS" \
           <(cat "$CR1") <(cat "$CR2") \
           -o 04_kmer/clean.jf

  jellyfish histo -t "$THREADS" 04_kmer/clean.jf > 04_kmer/clean.histo

  awk '{sum+=$1*$2; cnt+=$2} END{if(cnt>0){printf "Mean_kmer_cov\t%.2f\n", sum/cnt}else{print "Mean_kmer_cov\tNA"}}' \
      04_kmer/clean.histo > 04_kmer/kmer_stats.txt

  log "  Jellyfish 完成（clean.histo / kmer_stats.txt）"
else
  log "Step3: Jellyfish 结果已存在，跳过"
fi

########################################
# Step4: SPAdes 组装，用 clean_reads
########################################
GE500=05_asm/spades/scaffolds.ge500.fasta
if [[ ! -s "$GE500" ]]; then
  mkdir -p 05_asm
  if [[ ! -s 05_asm/spades/scaffolds.fasta ]]; then
    log "Step4: SPAdes 组装 (-k $K_LIST --careful, 使用 clean_reads)"
    spades.py \
      --pe1-1 "$CR1" \
      --pe1-2 "$CR2" \
      -k $K_LIST \
      --careful \
      -t "$THREADS" \
      -m "$SPADES_MEM_GB" \
      -o 05_asm/spades
  else
    log "Step4: 发现已有 scaffolds.fasta，跳过组装"
  fi

  log "Step4b: 过滤 ≥500bp"
  awk 'BEGIN{FS=" "}/^>/{if(seqlen>=500){print header"\n"seq} header=$0; seq=""; seqlen=0; next}
       {seq=seq$0; seqlen+=length($0)}
       END{if(seqlen>=500){print header"\n"seq}}' \
       05_asm/spades/scaffolds.fasta > "$GE500"

  {
    echo -e "Metric\tValue"
    echo -n "Num_contigs_ge500bp\t"
    grep -c "^>" "$GE500"
    echo -n "Total_bases_ge500bp\t"
    awk '/^>/{next}{len+=length($0)}END{print len}' "$GE500"
  } > 05_asm/assembly_ge500_stats.txt

  log "  组装完成，统计见 05_asm/assembly_ge500_stats.txt"
else
  log "Step4: ≥500bp 组装已存在，跳过"
fi

CUR="$GE500"

########################################
# Step5: Pilon 抛光（可选，用 clean_reads）
########################################
if [[ "$DO_POLISH" == "true" ]]; then
  if [[ ! -s 06_polish/polished.ge500.fasta ]]; then
    log "Step5: Pilon 抛光"

    mkdir -p 06_polish
    bwa index "$CUR"

    if [[ ! -s 06_polish/aln.bam ]]; then
      bwa mem -t "$THREADS" "$CUR" "$CR1" "$CR2" \
        | samtools sort -@ "$THREADS" -o 06_polish/aln.bam
      samtools index 06_polish/aln.bam
    fi

    if [[ ! -s 06_polish/aln.slim.bam ]]; then
      samtools view -@ "$THREADS" -b -f 2 -F 0x904 -q 20 06_polish/aln.bam \
        | samtools sort -@ "$THREADS" -o 06_polish/aln.slim.bam
      samtools index 06_polish/aln.slim.bam
    fi

    PILON_JAVA_OPTS="-Xmx${JAVA_MEM}" pilon \
      --genome "$CUR" \
      --frags 06_polish/aln.slim.bam \
      --outdir 06_polish \
      --output polished \
      --fix snps,indels \
      --changes

    awk 'BEGIN{FS=" "}/^>/{if(seqlen>=500){print header"\n"seq} header=$0; seq=""; seqlen=0; next}
         {seq=seq$0; seqlen+=length($0)}
         END{if(seqlen>=500){print header"\n"seq}}' \
         06_polish/polished.fasta > 06_polish/polished.ge500.fasta

    {
      echo -e "Metric\tValue"
      echo -n "Num_contigs_ge500bp\t"
      grep -c "^>" 06_polish/polished.ge500.fasta
      echo -n "Total_bases_ge500bp\t"
      awk '/^>/{next}{len+=length($0)}END{print len}' 06_polish/polished.ge500.fasta
    } > 06_polish/polished_ge500_stats.txt

    log "  Pilon 抛光完成，统计见 06_polish/polished_ge500_stats.txt"

  else
    log "Step5: polished.ge500.fasta 已存在，跳过"
  fi
  CUR="06_polish/polished.ge500.fasta"
else
  log "Step5: 未开启 DO_POLISH，跳过"
fi

########################################
# Step6: CD-HIT-EST 去冗余
########################################
CDHIT_OUT=06_cdhit/assembly.cdhit0.98.fa
mkdir -p 06_cdhit

if [[ ! -s "$CDHIT_OUT" ]]; then
  log "Step6: CD-HIT-EST 去冗余 (-c 0.98 -aS 0.9)"

  "$CDHIT_BIN" \
    -i "$CUR" \
    -o "$CDHIT_OUT" \
    -c 0.98 \
    -aS 0.9 \
    -G 0 \
    -g 1 \
    -M 0 \
    -T "$THREADS"

  {
    echo -e "Metric\tValue"
    echo -n "Num_contigs_after_cdhit\t"
    grep -c "^>" "$CDHIT_OUT"
    echo -n "Total_bases_after_cdhit\t"
    awk '/^>/{next}{len+=length($0)}END{print len}' "$CDHIT_OUT"
  } > 06_cdhit/cdhit_stats.txt

  log "  CD-HIT 完成，统计见 06_cdhit/cdhit_stats.txt"

else
  log "Step6: CD-HIT 结果已存在，跳过"
fi

FINAL="$CDHIT_OUT"

########################################
# Step7: QUAST & BUSCO 评估
########################################
if [[ ! -s 07_eval/quast/report.tsv ]]; then
  log "Step7: QUAST 评估"
  mkdir -p 07_eval/quast
  quast.py -o 07_eval/quast -t "$THREADS" "$FINAL"
else
  log "Step7: QUAST 已存在，跳过"
fi

if [[ "$RUN_BUSCO" == "true" ]]; then
  if [[ ! -d 07_eval/busco/run_busco ]]; then
    log "Step7: BUSCO 评估 ($BUSCO_LINEAGE)"
    busco \
      -i "$FINAL" \
      -l "$BUSCO_LINEAGE" \
      -m genome \
      -o 07_eval/busco \
      -c "$THREADS"
  else
    log "Step7: BUSCO 已存在，跳过"
  fi
else
  log "Step7: 未开启 RUN_BUSCO，跳过"
fi

log "全部完成！最终基因组：$FINAL"

