args <- commandArgs(trailingOnly = TRUE)

# 引数の数を確認し、引数がない場合はエラーメッセージを表示する
if (length(args) == 0) {
  stop("引数が必要です。")
}

# 引数を表示する
library(readxl) 
#BiocManager::install("ChIPpeakAnno") ←インストール済
library("ChIPpeakAnno")
#ファイル名を入れる↓↓↓
file_path = args
# パスとファイル名に分割する関数を定義
split_path <- function(file_path) {
  # 最後のスラッシュを探す
  last_slash <- max(gregexpr("/", file_path)[[1]])
  
  # パスとファイル名に分割
  directory <- substr(file_path, 1, last_slash)
  file_name <- substr(file_path, last_slash + 1, nchar(file_path))
  
  return(list(directory = directory, file_name = file_name))
}
fname <- split_path(file_path)

peak.df = read_excel(paste0(file_path, ".xlsx"))
peak.df$peak <- 1:nrow(peak.df)
peak.rd <- GRanges(seqnames=peak.df$chr, ranges=IRanges(start=peak.df$start, end=peak.df$end, names= peak.df$peak))
#head(peak.rd)
gene_reference.df = read.table(paste0(fname$directory, "Athaliana_gene_ref.txt"), sep="\t", header=FALSE)
#head(gene_reference.df)
colnames(gene_reference.df) <- c("name", "txStart", "txEnd", "chrom", "strand")
gene_reference.rd <- GRanges(seqnames=Rle(gene_reference.df$chrom), ranges=IRanges(start=gene_reference.df$txStart, end=gene_reference.df$txEnd, name=gene_reference.df$name), strand=Rle(gene_reference.df$strand))
#head(gene_reference.rd)
peak.anno <- annotatePeakInBatch(peak.rd, AnnotationData=gene_reference.rd, featureType="TSS", output="both", multiple=TRUE, maxgap=10, PeakLocForDistance="middle", FeatureLocForDistance="TSS", select="first")
peak.add <- peak.df[,-c(1,2,3)]
peak.merge <- merge(peak.anno, peak.add, by.x="peak", by.y="peak", all.x=T)

#v出力ファイル
outfile = paste0(fname$directory, "annotation_", fname$file_name, ".csv")
write.csv(peak.merge, file = outfile)
#同じ階層にcsv出力

