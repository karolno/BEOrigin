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
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="Output directory. Only full path is accepted. Must be supplied.", metavar="character"),
  make_option(c("-c", "--rho"), type="numeric", default=NULL,
              help="Predefine rho (cellularity). Must be supplied.", metavar="numeric"),
  make_option(c("-p", "--psi"), type="numeric", default=NULL,
              help="Predefine psi (ploidy). Must be supplied.", metavar="numeric")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$ref) | is.null(opt$tumour) | is.null(opt$out) | is.null(opt$psi) | is.null(opt$rho)){
  print_help(opt_parser)
  stop("All arguments must be supplied.", call.=FALSE)
}

# # Identify all the files
normalname <- opt$ref
tumourname <- opt$tumour
# normalname <- "LP6008264-DNA_B06"
# tumourname <- "LP6008280-DNA_C06"

# prepare output directory
directory.in<-paste0(opt$out,"/", tumourname, "_vs_", normalname)
# directory.in<-paste0("/mnt/data/Karol/BE_origin/Newest/battenberg/out/","/", tumourname, "_vs_", normalname)

directory<-paste0(opt$out,"/", tumourname, "_vs_", normalname, "_final")
# directory<-paste0("/mnt/data/Karol/BE_origin/Newest/battenberg/out/","/", tumourname, "_vs_", normalname, "_final")

if(!dir.exists(directory)) {
  dir.create(directory, recursive = TRUE)
}

setwd(directory)
# sink(file = paste0(directory,"/log.txt"))
# sink(file = paste0(directory,"/log.txt"), type = "message", append = TRUE)

# MIN_NORMAL_DEPTH = opt$cutoff
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

# Fit a clonal copy number profile
fit.copy.number(samplename=tumourname,
                outputfile.prefix=paste(tumourname, "_", sep=""),
                inputfile.baf.segmented=paste(directory.in, "/",tumourname, ".BAFsegmented.txt", sep=""),
                inputfile.baf=paste(directory.in, "/",tumourname,"_mutantBAF.tab", sep=""),
                inputfile.logr=paste(directory.in, "/",tumourname,"_mutantLogR.tab", sep=""),
                dist_choice=CLONALITY_DIST_METRIC,
                ascat_dist_choice=ASCAT_DIST_METRIC,
                min.ploidy=MIN_PLOIDY,
                max.ploidy=MAX_PLOIDY,
                min.rho=MIN_RHO,
                min.goodness=MIN_GOODNESS_OF_FIT,
                uninformative_BAF_threshold=BALANCED_THRESHOLD,
                gamma_param=PLATFORM_GAMMA,
                use_preset_rho_psi=T,
                preset_rho=opt$rho, # cellularity
                preset_psi=opt$psi, # ploidy
                # preset_rho=0.72, # cellularity
                # preset_psi=2, # ploidy
                read_depth=30)

# Go over all segments, determine which segements are a mixture of two states and fit a second CN state
callSubclones(sample.name=tumourname,
              baf.segmented.file=paste(directory.in, "/",tumourname, ".BAFsegmented.txt", sep=""),
              logr.file=paste(directory.in, "/",tumourname,"_mutantLogR.tab", sep=""),
              rho.psi.file=paste(tumourname, "_rho_and_psi.txt",sep=""),
              output.file=paste(tumourname,"_subclones.txt", sep=""),
              output.figures.prefix=paste(tumourname,"_subclones_chr", sep=""),
              output.gw.figures.prefix=paste(tumourname,"_BattenbergProfile", sep=""),
              masking_output_file=paste(tumourname, "_segment_masking_details.txt", sep=""),
              chr_names=chrom_names.real,
              gamma=PLATFORM_GAMMA,
              segmentation.gamma=NA,
              siglevel=0.05,
              maxdist=0.01,
              noperms=1000)

# Prepare some figures to make it compatible with the initial output of the battenberg
make_posthoc_plots(samplename=tumourname, 
                   logr_file=paste(directory.in, "/",tumourname, "_mutantLogR_gcCorrected.tab", sep=""), 
                   subclones_file=paste(tumourname, "_subclones.txt", sep=""), 
                   rho_psi_file=paste(tumourname, "_rho_and_psi.txt", sep=""), 
                   bafsegmented_file=paste(directory.in, "/",tumourname, ".BAFsegmented.txt", sep=""), 
                   logrsegmented_file=paste(tumourname, ".logRsegmented.txt", sep=""), 
                   allelecounts_file=paste(directory.in, "/", tumourname, "_alleleCounts.tab", sep=""))


sessionInfo()

