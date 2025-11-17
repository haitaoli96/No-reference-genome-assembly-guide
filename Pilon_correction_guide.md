# Pilon 基因组抛光（Polishing）使用指南
## 1. 基本流程

```bash
bwa index assembly.fasta
bwa mem -t 12 assembly.fasta illumina_1.fastq.gz illumina_2.fastq.gz > assembly_illumina.sam
samtools sort -@ 12 -O bam -o assembly_illumina.sorted.bam assembly_illumina.sam
pilon --genome assembly.fasta --fix all --changes --frags assembly_illumina.sorted.bam --output pilon --outdir pilon_result --vcf
```

## 2. 常见内存错误解决
Pilon version 1.24 Thu Jan 28 13:00:45 2021 -0500
Genome: 05_asm/spades/scaffolds.ge500.fasta
Exception in thread "main" java.lang.OutOfMemoryError: Java heap space
        at org.broadinstitute.pilon.GenomeRegion.<init>(GenomeRegion.scala:57)
        at org.broadinstitute.pilon.GenomeFile.$anonfun$contigRegions$1(GenomeFile.scala:73)
        at org.broadinstitute.pilon.GenomeFile.$anonfun$contigRegions$1$adapted(GenomeFile.scala:73)
        at org.broadinstitute.pilon.GenomeFile$$Lambda$30/0x00007f71ec0ee378.apply(Unknown Source)
        at scala.collection.immutable.Range.map(Range.scala:59)
        at org.broadinstitute.pilon.GenomeFile.contigRegions(GenomeFile.scala:73)
        at org.broadinstitute.pilon.GenomeFile.$anonfun$regions$1(GenomeFile.scala:53)
        at org.broadinstitute.pilon.GenomeFile$$Lambda$29/0x00007f71ec0e7c00.apply(Unknown Source)
        at scala.collection.immutable.List.map(List.scala:250)
        at org.broadinstitute.pilon.GenomeFile.<init>(GenomeFile.scala:53)
        at org.broadinstitute.pilon.Pilon$.main(Pilon.scala:108)
        at org.broadinstitute.pilon.Pilon.main(Pilon.scala)
这是 Pilon 内存不够（Java heap OOM）。不需要重跑前面步骤，直接在抛光这一步做“加内存 / 降负载 / 分批跑”三选一（或叠加）即可。给你可直接复制的修正方案——按顺序尝试，通常 A 就能解决。

编辑 Pilon 启动脚本：
```
which pilon
vim /path/to/pilon
```

增大：

```
from os import access, getenv, X_OK
jar_file = 'pilon.jar'

default_jvm_mem_opts = ['-Xms32g', '-Xmx128g']

# !!! End of parameter section. No user-serviceable code below this line !!!

def real_dirname(path):
    """Return the symlink-resolved, canonicalized directory-portion of path."""
    return os.path.dirname(os.path.realpath(path))


```
