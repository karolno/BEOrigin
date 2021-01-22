#!/usr/bin/env Rscript

###############################################################################
# load libraries
library(matrixStats)
library(data.table)
library(dpclust3p)
library(DPClust)
library(optparse)

# read in options
option_list = list(
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="Name of the patient to be run. Must be supplied.", metavar="character"),
  make_option(c("-d", "--dirout"), type="character", default=NULL,
              help="Directory to save the files. Must be supplied.", metavar="character"),
  make_option(c("-g", "--gender"), type="character", default=NULL,
              help="Sex of the sample. Must be supplied (one of male or female).", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$sample) | is.null(opt$dirout) | is.null(opt$gender)){
  print_help(opt_parser)
  stop("Sample name must be supplied.", call.=FALSE)
}


directory.out <- opt$dirout
sample <- opt$sample
gender <- opt$gender

if (!gender %in% c("male", "female")) {
  stop("Provide male or female as gender")
}


# get the directory with SNV counts
directory<-"/home/karolno/data/Karol/BE_origin/Newest/barretts_origin_strelka_snvs_case_counts/"

# get the list of samples with metadata
samples <- read.delim("/home/karolno/Dropbox/Postdoc/git/BE_WGS/Samples", stringsAsFactors = FALSE)

# suffix with the name of the files
suffix <- "_snvs_readcounts.txt"



directory.out<-paste0(directory.out, sample, "_final")
if(!dir.exists(directory.out)) {
  dir.create(directory.out, recursive = TRUE)
}

setwd(directory.out)

###############################################################################
# process the data to make them compatible with dpclust
###############################################################################
print(sample)


# read all data for this sample
all.data <- read.delim(paste0(directory,sample, suffix ), stringsAsFactors = FALSE)
# keep only 23 chromosomes
all.data <- all.data[all.data$Chromosome %in% c(1:22,"X"),]
# all.data$tracking <- paste(all.data$Chromosome, all.data$Position, sep = "_")

# Calculate the count of each nucleotide
all.data2<-data.frame(#tracking = all.data$tracking, 
  Chromosome = all.data$Chromosome,
  Position = all.data$Position,
  Reference_base = all.data$Reference_base,
  Alternate_allele = all.data$Alternate_allele,
  Sample = all.data$Sample, 
  Depth = all.data$Depth,
  A = all.data[,8] + all.data[,9], 
  C = all.data[,10] + all.data[,11], 
  G = all.data[,12] + all.data[,13], 
  T = all.data[,14] + all.data[,15]
)

# Get tissue IDs
all.data2$Sample <- samples$Pool[samples$Patient == sample][match(all.data2$Sample, samples$Tissue[samples$Patient == sample])]
tissues<- unique(all.data2$Sample)

# cast data to get information for missing lines (counter was not counting the coverage and nucleotides if coverage was 0 in a given samples)
all.data3<-data.table::dcast(setDT(all.data2), 
                             formula = Chromosome + Position + Reference_base + Alternate_allele ~ Sample, 
                             value.var = c("Depth", "A", "C", "G", "T"), 
                             fill = 0) 
all.data3<-as.data.frame(all.data3)

if(!dir.exists(paste0(directory.out,"/counts/"))) {
  dir.create(paste0(directory.out,"/counts/"), recursive = TRUE)
}

write.table(all.data3[,1:4],file = paste0(directory.out,"/counts/","loci", ".txt"), row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)

# write the count data for each sample in the dpclust compatible format
for (tissue in tissues) {
  # tissue<-tissues[1]
  tmp<-all.data3[,c(1,2, grep(tissue, colnames(all.data3)))]
  tmp<-tmp[,c(1,2,4:7,3)]
  colnames(tmp)<-c("#CHR", "POS", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth")
  write.table(tmp,file = paste0(directory.out,"/counts/",tissue, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
}

# Save some memory
rm(all.data, all.data2, all.data3)

###############################################################################
# run dpclust3d to prepare all files for dpclust functions
###############################################################################

settings.file<-list()

# ID of the D2 sample
D2.sample<-samples$Pool[samples$Patient == sample & samples$Tissue == "D2"]

if(!dir.exists(paste0(directory.out,"/dpclust3d/"))) {
  dir.create(paste0(directory.out,"/dpclust3d/"), recursive = TRUE)
}


# prepare the INfo file for each samples together wiht infomration needed for the settings file
for(tmp in samples$Pool[samples$Patient == sample & samples$Tissue == "GC"]) {
  # tmp<-samples$Pool[samples$Patient == sample & samples$Tissue != "D2"][1]
  runGetDirichletProcessInfo(loci_file=paste0(directory.out,"/counts/","loci", ".txt"),
                             allele_frequencies_file=paste0(directory.out,"/counts/",tmp, ".txt"),
                             cellularity_file=paste0("/home/karolno/data/Karol/BE_origin/Newest/battenberg/out/",tmp, "_vs_", D2.sample, "_final/", "/", tmp, "_rho_and_psi.txt"),
                             subclone_file=paste0("/home/karolno/data/Karol/BE_origin/Newest/battenberg/out/",tmp, "_vs_", D2.sample, "_final/", "/", tmp, "_subclones.txt"),
                             gender=gender,
                             SNP.phase.file="NA",
                             mut.phase.file="NA",
                             output_file=paste0(directory.out,"/dpclust3d/",tmp, "_dpclust.txt"))
  settings.file[[tmp]]<-c("sample" = sample, "subsample" = tmp, "datafile" = paste0(tmp, "_dpclust.txt"), "cellularity" = read.delim(paste0("/home/karolno/data/Karol/BE_origin/Newest/battenberg/out/",tmp, "_vs_", D2.sample, "_final/", "/", tmp, "_rho_and_psi.txt"))["FRAC_GENOME", "rho"])
}

# save setting file
settings.file<-as.data.frame(do.call(rbind, settings.file))
write.table(settings.file,file = paste0(directory.out,"/dpclust3d/","dpclust_settings", ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")



###############################################################################
# 2018-11-01
# DPClust v2.2.6 pipeline implementation.
# sd11 [at] sanger.ac.uk
# This was copied from https://github.com/Wedge-lab/dpclust/blob/master/inst/example/dpclust_pipeline.R
###############################################################################

#####################################################################################
# Command line options i.e. variable parameters
#####################################################################################



if(!dir.exists(paste0(directory.out,"/dpclust"))) {
  dir.create(paste0(directory.out,"/dpclust"), recursive = TRUE)
}


run = 1
datpath = paste0(directory.out,"/dpclust3d")
outdir = paste0(directory.out,"/dpclust")
purity_file = paste0(directory.out,"/dpclust3d/","dpclust_settings", ".txt")
analysis_type = "nd_dp" # Analysis type to run
no.iters = 2000 # Number of iterations to run the MCMC chain
no.iters.burn.in = 1000 # Number of iterations to discard as burnin
mut.assignment.type = 1 # Mutation assignment method
num_muts_sample = 50000 # Number of mutations from which downsampling starts
bin.size = NULL # Binsize to use when constructing multi-dimensional density - only used when number of samples > 1
seed = 123 # Provide a seed
assign_sampled_muts = TRUE # Whether to assign mutations that have been removed during sampling
keep_temp_files = T # Keep intermediate output files
min_muts_cluster = -1 # Minimum number of mutations per cluster required for it to be kept in the final output, set to -1 to disable (default), see also --min_frac_muts_cluster
min_frac_muts_cluster = -1 # Minimum fraction of mutations per cluster required for it to be kept in the final output, set to -1 to disable, see also --min_muts_cluster

if (is.null(outdir)) { outdir = getwd() }

#####################################################################################
# Fixed parameters
#####################################################################################
is.male = (gender == "male")
#not used: min_sampling_factor = 1.1
#not used: is.vcf = F
# assign_sampled_muts = T
sample.snvs.only = T # Perform sampling on just the SNVs and not on CNAs
remove.snvs = F # Clear all SNVs, to perform clustering on CNAs only - This needs a better solution
generate_cluster_ordering = T
species = "human" # mouse also supported, just changes the chromosomes on which mutations are kept, has not effect on functionality

# Cocluster CNA parameters
co_cluster_cna = F
add.conflicts = F # Make the conflicts matrix in a dataset - Flag pertains to both copy number and mut2mut phasing
cna.conflicting.events.only = F # Add only those CNAs that are conflicting
num.clonal.events.to.add = 1 # Add this many clonal CNA events to the clustering
min.cna.size = 100 # Minim size in 10kb for a CNA event to be included

#####################################################################################
# Process input
#####################################################################################
# Parse the input file and obtain the required data for this run
sample2purity = read.table(purity_file, header=T, stringsAsFactors=F)
samplename = unique(sample2purity$sample)[run]
datafiles = sample2purity[sample2purity$sample==samplename,]$datafile
subsamples = sample2purity[sample2purity$sample==samplename,]$subsample
cellularity = sample2purity[sample2purity$sample==samplename,]$cellularity
cndatafiles = NA
# if ("sex" %in% colnames(sample2purity)) {
#   is.male = (sample2purity[sample2purity$sample==samplename,]$sex=="male")[1]
#   cndatafiles = sample2purity[sample2purity$sample==samplename,]$cndatafile
# } else {
#   is.male = T
#   cndatafiles = NA
# }

if ("mutphasing" %in% colnames(sample2purity)) {
  mutphasingfiles = sample2purity[sample2purity$sample==samplename,]$mutphasing
} else {
  mutphasingfiles = NA
}


#####################################################################################
# Create output path
#####################################################################################
outdir = file.path(outdir, paste(samplename, "_DPoutput_", no.iters,"iters_", no.iters.burn.in, "burnin_seed", seed, "/", sep=""))

#####################################################################################
# Setup parameters
#####################################################################################
run_params = make_run_params(no.iters, no.iters.burn.in, mut.assignment.type, num_muts_sample, is.male=is.male, min_muts_cluster=min_muts_cluster, min_frac_muts_cluster=min_frac_muts_cluster, species=species, assign_sampled_muts=assign_sampled_muts, keep_temp_files=keep_temp_files, generate_cluster_ordering=generate_cluster_ordering)
sample_params = make_sample_params(datafiles, cellularity, is.male, samplename, subsamples, mutphasingfiles)
advanced_params = make_advanced_params(seed)

#####################################################################################
# Run clustering
#####################################################################################
RunDP(analysis_type=analysis_type, 
      run_params=run_params,
      sample_params=sample_params, 
      advanced_params=advanced_params,
      outdir=outdir, 
      cna_params=NULL, 
      mutphasingfiles=NULL
)


sessionInfo()