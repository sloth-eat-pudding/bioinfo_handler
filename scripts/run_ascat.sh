#!/bin/bash

# 參數：腫瘤BAM、正常BAM、輸出前綴
TUMOR_BAM=$1
NORMAL_BAM=$2
OUTPUT_PREFIX=$3
gender=$4

# 指定 alleleCounter 路徑與前綴
ALLELE_COUNTER="/bip8_disk/zhenyu112/alleleCount/bin/alleleCounter"
ALLELES_PREFIX="/bip8_disk/zhenyu112/ascat_data/alleles/G1000_alleles_hg38_chr"
LOCI_PREFIX="/bip8_disk/zhenyu112/ascat_data/loci/G1000_loci_hg38_chr"

# 建立一個暫時性的 R script
R_SCRIPT=$(mktemp --suffix=.R)

cat <<EOF > $R_SCRIPT
library(ASCAT)

setwd("${OUTPUT_PREFIX}")

ascat.prepareHTS(
  tumourseqfile = "${TUMOR_BAM}",
  normalseqfile = "${NORMAL_BAM}",
  tumourname = "tumor",
  normalname = "normal",
  allelecounter_exe = "${ALLELE_COUNTER}",
  skip_allele_counting_normal = FALSE,
  skip_allele_counting_tumour = FALSE,
  alleles.prefix = "${ALLELES_PREFIX}",
  loci.prefix = "${LOCI_PREFIX}",
  gender = "${gender}",
  genomeVersion = "hg38",
  nthreads = 50,
  tumourLogR_file = "Tumor_LogR.txt",
  tumourBAF_file = "Tumor_BAF.txt",
  normalLogR_file = "Germline_LogR.txt",
  normalBAF_file = "Germline_BAF.txt",
  loci_binsize = 500,
  min_base_qual= 10,
  additional_allelecounter_flags="-f 0")

ascat.bc = ascat.loadData(
  Tumor_LogR_file = "Tumor_LogR.txt",
  Tumor_BAF_file = "Tumor_BAF.txt",
  Germline_LogR_file = "Germline_LogR.txt",
  Germline_BAF_file = "Germline_BAF.txt",
  gender = "${gender}",
  genomeVersion = "hg38"
)

ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")

ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, write_segments = TRUE)
EOF

# 執行 R script
Rscript $R_SCRIPT

# 刪除暫時的 R script（除非你想保留）
rm $R_SCRIPT
