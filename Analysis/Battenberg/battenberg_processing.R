#!/usr/bin/env Rscript

###############################
library(optparse)
library(Battenberg)
library(doParallel)

# read in options
option_list = list(
  make_option(c("-r", "--ref"), type="character", default=NULL,
              help="Name of the reference sample. Must be supplied.", metavar="character"),
  make_option(c("-t", "--tumour"), type="character", default=NULL,
              help="Name of the tumour sample. Must be supplied.", metavar="character"),
  make_option(c("-c", "--cutoff"), type="numeric", default=10,
              help="Minimal coverage in the normal samples to keep SNP. [default = %default]", metavar="numeric"),
  make_option(c("-p", "--threads"), type="numeric", default=1,
              help="Number of threads for the parallel part of the pipeline. [default = %default]", metavar="numeric"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="Output directory. Only full path is accepted. Must be supplied.", metavar="character"),
  make_option(c("-g", "--gender"), type="character", default="XY",
              help="Patient's gender [default = %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$ref) | is.null(opt$tumour) | is.null(opt$out)){
  print_help(opt_parser)
  stop("All arguments except for cutoff (default = 10 reads), gender (default = XY) and number of threads (default = 1) must be supplied.", call.=FALSE)
}

gender <- opt$gender
# gender <- "XY"

# # Identify all the files
normalname <- opt$ref
tumourname <- opt$tumour
# normalname <- "LP6008268-DNA_B04"
# tumourname <- "SLX-16213.DNAA012"

# prepare output directory
directory<-paste0(opt$out,"/", tumourname, "_vs_", normalname)
# directory<-paste0("/mnt/data/Karol/BE_origin/Newest/battenberg/out/","/", tumourname, "_vs_", normalname)

if(!dir.exists(directory)) {
  dir.create(directory, recursive = TRUE)
}

setwd(directory)
sink(file = paste0(directory,"/log.txt"))
# sink(file = paste0(directory,"/log.txt"), type = "message", append = TRUE)

MIN_NORMAL_DEPTH = opt$cutoff
# MIN_NORMAL_DEPTH = 10

chrom_names.real<-c(1:22, "X")
chrom_names.number<-c(1:23)


ALLELECOUNTER = "/home/karolno/Packages/alleleCounter/bin/alleleCounter"
IMPUTEINFOFILE = "/home/karolno/Packages/battenberg/impute_info.txt"
G1000PREFIX = "/home/karolno/Packages/battenberg/battenberg_1000genomesloci2012_v3/1000genomesAlleles2012_chr"
G1000PREFIX_AC = "/home/karolno/Packages/battenberg/battenberg_1000genomesloci2012_v3/1000genomesloci2012_chr"
GCCORRECTPREFIX = "/home/karolno/data/Karol/mutREAD/data/battenberg_wgs_gc_correction_1000g_v3/1000_genomes_GC_corr_chr_"
REPLICCORRECTPREFIX = "/home/karolno/data/Karol/mutREAD/data/battenberg_wgs_replic_correction_1000g_v3/1000_genomes_replication_timing_chr_"
IMPUTE_EXE = "/home/karolno/Packages/bin/impute2"
PLATFORM_GAMMA = 1
PHASING_GAMMA = 1
SEGMENTATION_GAMMA = 10
SEGMENTATIIN_KMIN = 3
PHASING_KMIN = 1
CLONALITY_DIST_METRIC = 0
ASCAT_DIST_METRIC = 1
MIN_PLOIDY = 1.6
MAX_PLOIDY = 4.8
MIN_RHO = 0.1
MIN_GOODNESS_OF_FIT = 0.63
BALANCED_THRESHOLD = 0.51
MIN_BASE_QUAL = 20
MIN_MAP_QUAL = 35
CALC_SEG_BAF_OPTION = 3
NTHREADS = opt$threads
PROBLEMLOCI = "/home/karolno/Packages/battenberg/probloci_270415.txt.gz"

# tumour sample
all.counts.tumour<-read.delim(paste0("/home/karolno/data/Karol/BE_origin/Newest/barretts_origin_g1k_snps_strelka_snv_counts/",tumourname, "_1000gp_readcounts_processed.txt" ))
colnames(all.counts.tumour)<-c("#CHR", "POS", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth")
all.counts.normal<-read.delim(paste0("/home/karolno/data/Karol/BE_origin/Newest/barretts_origin_g1k_snps_strelka_snv_counts/",normalname, "_1000gp_readcounts_processed.txt" ))
colnames(all.counts.normal)<-c("#CHR", "POS", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth")


# process the samples into counts per individual chromosome. Only print data for SNP that were identified in both samples. 
for (chr in 1:23) {
  # chr <- 1
  tmp.tumour<-all.counts.tumour[all.counts.tumour$`#CHR`==chrom_names.real[chr],]
  tmp.normal<-all.counts.normal[all.counts.normal$`#CHR`==chrom_names.real[chr],]
  write.table(tmp.tumour[tmp.tumour$POS %in% intersect(tmp.tumour$POS, tmp.normal$POS),],file = paste0(tumourname,"_alleleFrequencies_chr",chr, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(tmp.normal[tmp.normal$POS %in% intersect(tmp.tumour$POS, tmp.normal$POS),],file = paste0(normalname,"_alleleFrequencies_chr",chr, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
}


# run the script from battenberg that crates BAF and LogR files.
getBAFsAndLogRs(tumourAlleleCountsFile.prefix=paste(tumourname,"_alleleFrequencies_chr", sep=""),
                normalAlleleCountsFile.prefix=paste(normalname,"_alleleFrequencies_chr", sep=""),
                figuresFile.prefix=paste(tumourname, "_", sep=''),
                BAFnormalFile=paste(tumourname,"_normalBAF.tab", sep=""),
                BAFmutantFile=paste(tumourname,"_mutantBAF.tab", sep=""),
                logRnormalFile=paste(tumourname,"_normalLogR.tab", sep=""),
                logRmutantFile=paste(tumourname,"_mutantLogR.tab", sep=""),
                combinedAlleleCountsFile=paste(tumourname,"_alleleCounts.tab", sep=""),
                chr_names=chrom_names.real,
                g1000file.prefix=G1000PREFIX,
                minCounts=MIN_NORMAL_DEPTH,
                samplename=tumourname)


gc.correct.wgs(Tumour_LogR_file=paste(tumourname,"_mutantLogR.tab", sep=""),
               outfile=paste(tumourname,"_mutantLogR_gcCorrected.tab", sep=""),
               correlations_outfile=paste(tumourname, "_GCwindowCorrelations.txt", sep=""),
               gc_content_file_prefix=GCCORRECTPREFIX,
               replic_timing_file_prefix=REPLICCORRECTPREFIX,
               chrom_names=chrom_names.real)


battenberg(tumourname=tumourname, 
           normalname=normalname, 
           # tumour_data_file=list.files(bamlocation, pattern=tumour.sample, all.files=TRUE, full.names=TRUE)[1], 
           # normal_data_file=list.files(bamlocation, pattern=ref.sample, all.files=TRUE, full.names=TRUE)[1], 
           ismale=(gender == "XY"), 
           imputeinfofile=IMPUTEINFOFILE, 
           g1000prefix=G1000PREFIX, 
           g1000allelesprefix=G1000PREFIX_AC, 
           gccorrectprefix=GCCORRECTPREFIX, 
           repliccorrectprefix=REPLICCORRECTPREFIX, 
           problemloci=PROBLEMLOCI, 
           data_type="wgs",
           impute_exe=IMPUTE_EXE,
           allelecounter_exe=ALLELECOUNTER,
           nthreads=NTHREADS,
           platform_gamma=PLATFORM_GAMMA,
           phasing_gamma=PHASING_GAMMA,
           segmentation_gamma=SEGMENTATION_GAMMA,
           segmentation_kmin=SEGMENTATIIN_KMIN,
           phasing_kmin=PHASING_KMIN,
           clonality_dist_metric=CLONALITY_DIST_METRIC,
           ascat_dist_metric=ASCAT_DIST_METRIC,
           min_ploidy=MIN_PLOIDY,
           max_ploidy=MAX_PLOIDY,
           min_rho=MIN_RHO,
           min_goodness=MIN_GOODNESS_OF_FIT,
           uninformative_BAF_threshold=BALANCED_THRESHOLD,
           min_normal_depth=MIN_NORMAL_DEPTH,
           min_base_qual=MIN_BASE_QUAL,
           min_map_qual=MIN_MAP_QUAL,
           calc_seg_baf_option=CALC_SEG_BAF_OPTION,
           skip_allele_counting=TRUE,
           skip_preprocessing=TRUE,
           skip_phasing=FALSE,
           prior_breakpoints_file=NULL)


sessionInfo()

