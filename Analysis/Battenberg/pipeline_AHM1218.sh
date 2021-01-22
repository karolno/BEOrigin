#!/usr/bin/bash
strelkadir=/home/karolno/data/Karol/BE_origin/Newest/barretts_origin_g1k_snps_strelka_snv_counts/
dirout=/mnt/data/Karol/BE_origin/Newest/battenberg/out/
sample=AHM1218
gender=male
reference=LP6008338-DNA_G03 # D2
tumour1=SLX-16213.DNAA007 # NE
tumour2=SLX-16213.DNAA013 # GC
tumour3=LP6008337-DNA_G03 # BE


#######################################################
# Split the data from the counting Ginny did into individual chromosomes
#######################################################
# perl -alne '$F[0]=~s/chr//; print $F[0],"\t", $F[1], "\t", $F[6]+$F[7],  "\t", $F[8]+$F[9], "\t", $F[10]+$F[11], "\t", $F[12]+$F[13], "\t", $F[5] ' <$strelkadir/${reference}_1000gp_readcounts.txt >$strelkadir/${reference}_1000gp_readcounts_processed.txt
# perl -alne '$F[0]=~s/chr//; print $F[0],"\t", $F[1], "\t", $F[6]+$F[7],  "\t", $F[8]+$F[9], "\t", $F[10]+$F[11], "\t", $F[12]+$F[13], "\t", $F[5] ' <$strelkadir/${tumour1}_1000gp_readcounts.txt >$strelkadir/${tumour1}_1000gp_readcounts_processed.txt
# perl -alne '$F[0]=~s/chr//; print $F[0],"\t", $F[1], "\t", $F[6]+$F[7],  "\t", $F[8]+$F[9], "\t", $F[10]+$F[11], "\t", $F[12]+$F[13], "\t", $F[5] ' <$strelkadir/${tumour2}_1000gp_readcounts.txt >$strelkadir/${tumour2}_1000gp_readcounts_processed.txt
# perl -alne '$F[0]=~s/chr//; print $F[0],"\t", $F[1], "\t", $F[6]+$F[7],  "\t", $F[8]+$F[9], "\t", $F[10]+$F[11], "\t", $F[12]+$F[13], "\t", $F[5] ' <$strelkadir/${tumour3}_1000gp_readcounts.txt >$strelkadir/${tumour3}_1000gp_readcounts_processed.txt

#######################################################
# Run main buttenberg scripts
#######################################################

# mkdir -p ${dirout}/${tumour1}_vs_${reference}
# # Run battenberg script
# echo "Running sample:" $tumour1 >${dirout}/${tumour1}_vs_${reference}/stdout.log
# echo "Running sample:" $tumour1 >${dirout}/${tumour1}_vs_${reference}/stderr.log
# Rscript ~/Dropbox/Postdoc/git/BEOrigin/Analysis/Battenberg/battenberg_processing.R -r $reference -t $tumour1 -c 10 -o ${dirout} -p 10 -g $gender >>${dirout}/${tumour1}_vs_${reference}/stdout.log 2>>${dirout}/${tumour1}_vs_${reference}/stderr.log 
# 
# mkdir -p ${dirout}/${tumour2}_vs_${reference}
# # Run battenberg script
# echo "Running sample:" $tumour2 >${dirout}/${tumour2}_vs_${reference}/stdout.log
# echo "Running sample:" $tumour2 >${dirout}/${tumour2}_vs_${reference}/stderr.log
# Rscript ~/Dropbox/Postdoc/git/BEOrigin/Analysis/Battenberg/battenberg_processing.R -r $reference -t $tumour2 -c 10 -o ${dirout} -p 10 -g $gender >>${dirout}/${tumour2}_vs_${reference}/stdout.log 2>>${dirout}/${tumour2}_vs_${reference}/stderr.log 
# 
# mkdir -p ${dirout}/${tumour3}_vs_${reference}
# # Run battenberg script
# echo "Running sample:" $tumour3 >${dirout}/${tumour3}_vs_${reference}/stdout.log
# echo "Running sample:" $tumour3 >${dirout}/${tumour3}_vs_${reference}/stderr.log
# Rscript ~/Dropbox/Postdoc/git/BEOrigin/Analysis/Battenberg/battenberg_processing.R -r $reference -t $tumour3 -c 10 -o ${dirout} -p 10 -g $gender >>${dirout}/${tumour3}_vs_${reference}/stdout.log 2>>${dirout}/${tumour3}_vs_${reference}/stderr.log 

#######################################################
# Perform refiting of normal samples if required. 
# The samples should be called with cellularity close to 1 and ploidy 2. 
# If this is not the case, they should be recalled using hardcoded values for this specimen.
#######################################################
#
#######################################################
# Both NE and GC have to be refitted
#######################################################
# 
# # Run battenberg script
# mkdir -p ${dirout}/${tumour1}_vs_${reference}_refit
# # Run battenberg script
# echo "Running sample:" $tumour1 >${dirout}/${tumour1}_vs_${reference}_refit/stdout.log
# echo "Running sample:" $tumour1 >${dirout}/${tumour1}_vs_${reference}_refit/stderr.log
# Rscript ~/Dropbox/Postdoc/git/BEOrigin/Analysis/Battenberg/battenberg_processing_rerun.R -r $reference -t $tumour1 -o ${dirout} -c 0.98 -p 2  >>${dirout}/${tumour1}_vs_${reference}_refit/stdout.log 2>>${dirout}/${tumour1}_vs_${reference}_refit/stderr.log 
# 
# # Run battenberg script
# mkdir -p ${dirout}/${tumour2}_vs_${reference}_refit
# # Run battenberg script
# echo "Running sample:" $tumour2 >${dirout}/${tumour2}_vs_${reference}_refit/stdout.log
# echo "Running sample:" $tumour2 >${dirout}/${tumour2}_vs_${reference}_refit/stderr.log
# Rscript ~/Dropbox/Postdoc/git/BEOrigin/Analysis/Battenberg/battenberg_processing_rerun.R -r $reference -t $tumour2 -o ${dirout} -c 0.98 -p 2  >>${dirout}/${tumour2}_vs_${reference}_refit/stdout.log 2>>${dirout}/${tumour2}_vs_${reference}_refit/stderr.log 

#######################################################
# Run single samples DPclust
#######################################################

DPout="/mnt/data/Karol/BE_origin/Newest/DPClust_BE/"
mkdir -p ${DPout}
mkdir -p ${DPout}/${sample}
# Run battenberg script
echo "Running sample:" $sample >${DPout}/${sample}/stdout.log
echo "Running sample:" $sample >${DPout}/${sample}/stderr.log
Rscript ~/Dropbox/Postdoc/git/BEOrigin/Analysis/Battenberg/DPClust_final_BE.R -s $sample -d $DPout -g $gender >>${DPout}/${sample}/stdout.log 2>>${DPout}/${sample}/stderr.log 


DPout="/mnt/data/Karol/BE_origin/Newest/DPClust_GC/"
mkdir -p ${DPout}
mkdir -p ${DPout}/${sample}_refit
# Run battenberg script
echo "Running sample:" $sample >${DPout}/${sample}_refit/stdout.log
echo "Running sample:" $sample >${DPout}/${sample}_refit/stderr.log
Rscript ~/Dropbox/Postdoc/git/BEOrigin/Analysis/Battenberg/DPClust_final_GC_refit.R -s $sample -d $DPout -g $gender >>${DPout}/${sample}_refit/stdout.log 2>>${DPout}_refit/${sample}/stderr.log 


DPout="/mnt/data/Karol/BE_origin/Newest/DPClust_NE/"
mkdir -p ${DPout}
mkdir -p ${DPout}/${sample}_refit
# Run battenberg script
echo "Running sample:" $sample >${DPout}/${sample}_refit/stdout.log
echo "Running sample:" $sample >${DPout}/${sample}_refit/stderr.log
Rscript ~/Dropbox/Postdoc/git/BEOrigin/Analysis/Battenberg/DPClust_final_NE_refit.R -s $sample -d $DPout -g $gender >>${DPout}/${sample}_refit/stdout.log 2>>${DPout}/${sample}_refit/stderr.log 

#######################################################
# Check the solution of the above runs.
# The major peak of the BE sample should be aruond 1
# Since the cellularity of NE and GC was hardcoded to be close to 1, the major clone will not be around 1 for these samples
# This means these samples will have to rurun to make sure that they are realigned with and with new cellularities.
# First I have to run battenberg again to refit the data around the values obtained in the above scripts.
# For this patient, the cellularities are:
#######################################################
NEcell=0.43
GCcell=0.13

# Run battenberg script
mkdir -p ${dirout}/${tumour1}_vs_${reference}_final
# Run battenberg script
echo "Running sample:" $tumour1 >${dirout}/${tumour1}_vs_${reference}_final/stdout.log
echo "Running sample:" $tumour1 >${dirout}/${tumour1}_vs_${reference}_final/stderr.log
Rscript ~/Dropbox/Postdoc/git/BEOrigin/Analysis/Battenberg/battenberg_processing_final.R -r $reference -t $tumour1 -o ${dirout} -c $NEcell -p 2  >>${dirout}/${tumour1}_vs_${reference}_final/stdout.log 2>>${dirout}/${tumour1}_vs_${reference}_final/stderr.log

# Run battenberg script
mkdir -p ${dirout}/${tumour2}_vs_${reference}_final
# Run battenberg script
echo "Running sample:" $tumour2 >${dirout}/${tumour2}_vs_${reference}_final/stdout.log
echo "Running sample:" $tumour2 >${dirout}/${tumour2}_vs_${reference}_final/stderr.log
Rscript ~/Dropbox/Postdoc/git/BEOrigin/Analysis/Battenberg/battenberg_processing_final.R -r $reference -t $tumour2 -o ${dirout} -c $GCcell -p 2  >>${dirout}/${tumour2}_vs_${reference}_final/stdout.log 2>>${dirout}/${tumour2}_vs_${reference}_final/stderr.log


#######################################################
# Finally, rerun single sample DPclust to see if the values cluster around 1 for NE and GC
#######################################################
DPout="/mnt/data/Karol/BE_origin/Newest/DPClust_NE/"
mkdir -p ${DPout}/${sample}_final
# Run battenberg script
echo "Running sample:" $sample >${DPout}/${sample}_final/stdout.log
echo "Running sample:" $sample >${DPout}/${sample}_final/stderr.log
Rscript ~/Dropbox/Postdoc/git/BEOrigin/Analysis/Battenberg/DPClust_final_NE_final.R -s $sample -d $DPout -g $gender >>${DPout}/${sample}_final/stdout.log 2>>${DPout}/${sample}_final/stderr.log &

DPout="/mnt/data/Karol/BE_origin/Newest/DPClust_GC/"
mkdir -p ${DPout}/${sample}_final
# Run battenberg script
echo "Running sample:" $sample >${DPout}/${sample}_final/stdout.log
echo "Running sample:" $sample >${DPout}/${sample}_final/stderr.log
Rscript ~/Dropbox/Postdoc/git/BEOrigin/Analysis/Battenberg/DPClust_final_GC_final.R -s $sample -d $DPout -g $gender >>${DPout}/${sample}_final/stdout.log 2>>${DPout}/${sample}_final/stderr.log &

