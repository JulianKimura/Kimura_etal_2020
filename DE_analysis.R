#convert salmon output to kallisto format.  Each file listed here is folder straight from salmon output, with a .sf file inside
library(wasabi)
library(dplyr)
setwd("/Users/JulianKimura/Documents/Lab/Broseq/Broseq_2019_analysis/Salmon_output/")

sfdirs <-  c("EP1/", "EP2/","EP3/","Morula1/", "Morula2/", "Morula3/", "Dimple1/", "Dimple2/", "Dimple3/", "X55_751/","X55_752/", "X55_753/", "X75_951/", "X75_952/", "X75_953/", "Pill1/", "Pill2/", "Pill3/", "PH1/", "PH2/", "PH3/", "PPH1/", "PPH2/", "PPH3/")

prepare_fish_for_sleuth(sfdirs)

sample_id <- dir(file.path("/Users/JulianKimura/Documents/Lab/Broseq/Broseq_2019_analysis/Salmon_output/"))

kal_dirs <- file.path("/Users/JulianKimura/Documents/Lab/Broseq/Broseq_2019_analysis/Salmon_output/", sample_id)

kal_dirs

#move s2c.csv file into folder before running this command
s2c<-read.table("/Users/JulianKimura/Documents/Lab/Broseq/Broseq_2019_analysis/Salmon_output/Workbook2.txt",header = TRUE, stringsAsFactors=FALSE)

s2c <- mutate(s2c, path = kal_dirs)

####this will be the master table from which different sleuth objects can be made from. For each run where we construct a new sleuth object, make sure to create a table in s2c format with only TWO different conditions.
write.table(s2c, "/Users/JulianKimura/Documents/Lab/Broseq/Broseq_2019_analysis/s2c_MASTER.txt", sep = "\t")


library(sleuth)

####MAke Sleuth Object
####Read in the s2c table which tells sleuth where the quant files are as well as the metadata for stage and conditions
s2c<-read.table("/Users/JulianKimura/Documents/Lab/Broseq/Broseq_2019_analysis/s2c_Matrices/s2c.txt",header = TRUE, stringsAsFactors=FALSE)

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE,read_bootstrap_tpm=TRUE)

#####Differential Expression
so <- sleuth_fit(so, ~condition, 'full')

so <- sleuth_fit(so, ~1, 'reduced')

so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = TRUE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant)

####Adding up vs downreg
PPH<-s2c$sample[s2c$condition=="Pigmented_Prehatchling"]
EE<-s2c$sample[s2c$condition=="Everything_else"]
hiraw<-subset(so$obs_raw,so$obs_raw$sample %in% PPH)
lowraw<-subset(so$obs_raw,so$obs_raw$sample %in% EE)
highmeans<-aggregate(tpm~target_id,data=hiraw,FUN=function(x) c(mean=mean(x)))
colnames(highmeans)<-c("target_id","PPH")
lowmeans<-aggregate(tpm~target_id,data=lowraw,FUN=function(x) c(mean=mean(x)))
colnames(lowmeans)<-c("target_id","EE")
merged_means<-merge(lowmeans,highmeans,by=c("target_id"))
sleuth_table_wTPM<-left_join(sleuth_table,merged_means)

write.csv(sleuth_table_wTPM, "/Users/JulianKimura/Documents/Lab/Broseq/Broseq_2019_analysis/sleuth_PPH_v_EE_wTMP_FINAL.csv")