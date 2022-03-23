# my RNA seq pipline

# initial
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ ll
total 539704
drwxr-xr-x 1 xuzihe xuzihe       512 Mar  5 20:59 ./
drwxrwxrwx 1 xuzihe xuzihe       512 Mar  5 15:05 ../
-rwxrwxrwx 1 xuzihe xuzihe 253935557 Jan 24  2014 chr1.fa*
-rwxrwxrwx 1 xuzihe xuzihe  51834845 Jan 24  2014 chr22.fa*
-rwxrwxrwx 1 xuzihe xuzihe 122621321 Nov 11 19:48 test_R1.fq.gz*
-rwxrwxrwx 1 xuzihe xuzihe 124261935 Nov 11 19:48 test_R2.fq.gz*

(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ mkdir 1fastqc_Result


# step I. Fastqc
# fastqc 的使用:t参数代表线程数/核心数，-o表示输出结果放置的位置，后面跟的.gz文件表示需要质检的文件，最后的&表示后台运行
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ fastqc -t 8 -o ./1fastqc_Result -q test_*.gz

# 查看结果
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ cd 1fastqc_Result
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData/1fastqc_Result$ ll
total 2100
drwxr-xr-x 1 xuzihe xuzihe    512 Mar  5 21:05 ./
drwxr-xr-x 1 xuzihe xuzihe    512 Mar  5 21:04 ../
-rw-r--r-- 1 xuzihe xuzihe 642924 Mar  5 21:05 test_R1_fastqc.html
-rw-r--r-- 1 xuzihe xuzihe 377146 Mar  5 21:05 test_R1_fastqc.zip
-rw-r--r-- 1 xuzihe xuzihe 639654 Mar  5 21:05 test_R2_fastqc.html
-rw-r--r-- 1 xuzihe xuzihe 374597 Mar  5 21:05 test_R2_fastqc.zip

# step II. cutadapt
# cutadapt 的使用，j参数表示线程数/核心数，--times表示默认一条read只去一次adapter，-e表示容许的错误率0.1则是容许10%的错误；
# -O表示与参考序列的对齐数其中，3表示3bp与参考序列对齐， --quality-cutoff表示去除3’端的低质量read参数一般选25；
# -m 55表示仅需要保留去除adapter后大于55bp的序列，-a表示Read1的adapt序列 -A表示Read2的adapt序列，illuma公司获取（illuma公司的双端测序）。
# -o表示Read1的输出文件名.fq.gz结尾，-p表示Read2的输出文件名.fq.gz结尾；后续紧跟对应的Read1和Read2对应的原始测序文件
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ mkdir 2cutAd_Result

RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ cutadapt \
> -j 6 --times 1  -e 0.1  -O 3  --quality-cutoff 25  -m 55 \
> -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
> -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
> -o ./2cutAd_Result/test_R1_cutadapt.temp.fq.gz \
> -p ./2cutAd_Result/test_R2_cutadapt.temp.fq.gz \
> test_R1.fq.gz \
> test_R2.fq.gz > ./2cutAd_Result/test_cutadapt.temp.log 2>&1 &

# 查看结果
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ cd 2cutAd_Result
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData/2cutAd_Result$ ll
total 229000
drwxr-xr-x 1 xuzihe xuzihe       512 Mar  5 21:14 ./
drwxr-xr-x 1 xuzihe xuzihe       512 Mar  5 21:11 ../
-rw-r--r-- 1 xuzihe xuzihe 115727423 Mar  5 21:14 test_R1_cutadapt.temp.fq.gz
-rw-r--r-- 1 xuzihe xuzihe 118259336 Mar  5 21:14 test_R2_cutadapt.temp.fq.gz
-rw-r--r-- 1 xuzihe xuzihe      4763 Mar  5 21:14 test_cutadapt.temp.log

# Step III. 参考序列拼接
# 下载染色体序列文件chr1~22 chrX chrY UCSC网站
# 合并
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ cat chr1.fa chr22.fa > ref_hg38.fa &

# step IV. build index与alignment and Analysis
# 方法很多，有bowtie/bowtie2/tophat/tophat2/star/hisat/hisat2
# 主要使用star和hisat2方法

# ###########################################################################################
# STAR -> feature Counts -> Deseq2 -> GO/KEGG -> heatmap/WGCNA
# 需要提前下载hg38_refseq_from_ucsc.rm_XM_XR.fix_name.gtf GTF文件
# ###########################################################################################
# build index
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ STAR --runThreadN 4 --runMode genomeGenerate \
> --genomeDir 3-1star_index \
> --genomeFastaFiles ref_hg38.fa \
> --sjdbGTFfile hg38_NCBI_RefSeq.gtf \
> --sjdbOverhang 150 &
# alignment
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ mkdir 3-2star_test
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ STAR --genomeDir 3-1star_index \
> --runThreadN 6 \
> --readFilesIn ./2cutAd_Result/test_R1_cutadapt.temp.fq.gz ./2cutAd_Result/test_R2_cutadapt.temp.fq.gz \
> --readFilesCommand zcat \
> --outFileNamePrefix 3-2star_test \
> --outSAMtype BAM Unsorted \
> --outSAMstrandField intronMotif \
> --outSAMattributes All \
> --outFilterIntronMotifs RemoveNoncanonical > ./3-2star_test/test_STAR.log 2>&1 &
# samtools排序
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ samtools sort -O BAM 、
> -o 3-2star_testAligned.sort.out.bam \
> -@ 6 -m 2G \
> -T 3-2star_testAligned.sort.out.bam.temp 3-2star_testAligned.out.bam &
# feature Counts
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ featureCounts -t exon -g gene_id \
> -Q 6 --primary -s 0 -p -T 1 \
> -a hg38_NCBI_RefSeq.gtf \
> -o ./3-2star_test/star_featureCounts.featureCounts \
> 3-2star_testAligned.sort.out.bam > ./3-2star_test/star_featureCounts.log 2>&1 &
# 使用feature Counts文件进行DEseq2差异分析-> GO/KEGG -> heatmap/WGCNA (R)

# ###########################################################################################
# hisat2 -> feature Counts -> Deseq2 -> GO/KEGG -> heatmap/WGCNA
# 需要提前下载hg38_refseq_from_ucsc.rm_XM_XR.fix_name.gtf GTF文件
# ###########################################################################################
# 提前下载SNP文件(基因突变文件)、GTF文件(可变剪切文件) UCSC网站
# 以下三条处理是hisat2自带的操作
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ mkdir 4-1hisat2_index
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ mkdir 4-2hisat2_test
# make exon 
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ hisat2_extract_exons.py hg38_NCBI_RefSeq.gtf > ./4-1hisat2_index/hg38_refseq.exon &
# make splice site
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ hisat2_extract_splice_sites.py hg38_NCBI_RefSeq.gtf > ./4-1hisat2_index/hg38_refseq.ss &
# make snp and haplotype
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ hisat2_extract_snps_haplotypes_UCSC.py ref_hg38.fa snp151Common.txt snp151Common &
# build index
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ hisat2-build -p 6 \
> --snp snp151Common.snp \
> --haplotype snp151Common.haplotype \
> --exon ./4-1hisat2_index/hg38_refseq.exon \
> --ss ./4-1hisat2_index/hg38_refseq.ss ref_hg38.fa  ./4-1hisat2_index/ref_hg38.fa.snp_gtf > ./4-1hisat2_index/hisat2_build.log 2>&1 &
# alignmnet
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ hisat2 -p 6 \
> -x ./4-1hisat2_index/ref_hg38.fa.snp_gtf \
> -1 ./2cutAd_Result/test_R1_cutadapt.temp.fq.gz \
> -2 ./2cutAd_Result/test_R2_cutadapt.temp.fq.gz \
> -S ./4-2hisat2_test/hisat2.sam > ./4-2hisat2_test/hisat2.log 2>&1 &
# samtools排序并转换成bam文件
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData/4-2hisat2_test$ samtools sort -O BAM -o hisat2.sort.sam -@ 6 -m 2G -T hisat2.sort.bam.temp hisat2.sam &
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData/4-2hisat2_test$ samtools view -b -S hisat2.sort.sam > hisat2.sort.bam
# feature Counts
(RNAseq) xuzihe@LAPTOP-R8KHR20U:~/RNA-seq/testData$ featureCounts -t exon -g gene_id \
> -Q 6 --primary -s 0 -p -T 1 \
> -a hg38_NCBI_RefSeq.gtf \
> -o ./4-2hisat2_test/hisat2_featureCounts.featureCounts \
> ./4-2hisat2_test/hisat2.sort.bam > ./4-2hisat2_test/hisat2.featureCounts.log  2>&1 &
# 使用feature Counts文件进行DEseq2差异分析-> GO/KEGG -> heatmap/WGCNA (R)
