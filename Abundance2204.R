library('microclass')
library('readr')
args <- commandArgs(trailingOnly = TRUE)
data<-kraken2_read_report(args[1])

count_abundance_func<-function(data){
  data_filtered <- data[data$rank == "S" & data$clade_count>=50,] #рассматриваем виды  с потерей информации
  new_abundance<-data_filtered$clade_count/sum(data_filtered$clade_count)
  data_filtered <- cbind(data_filtered,new_abundance)
  names(data_filtered)[7]<-"Abundance_exact" #аbundance на уровне видов без G1 и S1
  Abundance_exact100<-(data_filtered$Abundance_exact*100)
  data_filtered<-cbind(data_filtered, Abundance_exact100)
  return(data_filtered)
}

count_ab<-count_abundance_func(data) # rename to abundances
set.seed(1001)
    
#growing_org <- "" 

data_binomial<-data.frame()
nbinom_distr_ab_func<-function(count_ab){
  for (r in 1:nrow(count_ab)){
    if (count_ab[r,]$name=="Streptococcus sp. S5") {
      abundance_with_noise <-rep(0, 7) # distr <- abundance_with_noise
      for (i in 0:6){
        growing_mu<-rnbinom(n=1,size=3,mu=count_ab[r,]$clade_count+i*20000)#growing_mu
        abundance_with_noise[i+1]<-growing_mu
      }
    }
    else{
      abundance_with_noise<-rnbinom(n=7,size=3,mu=count_ab[r,]$clade_count)
    }
    data_binomial<-rbind(data_binomial,abundance_with_noise)
  }
names(data_binomial) <- c("clade_nbin1","clade_nbin2","clade_nbin3","clade_nbin4","clade_nbin5","clade_nbin6","clade_nbin7")
abundance_nbin_sum<- data.frame(matrix(ncol = ncol(data_binomial), nrow = nrow(data_binomial)))
for (l in 1:ncol(data_binomial)){
  abundance_nbin<-data_binomial[[l]]/sum(data_binomial[[l]])
  abundance_nbin_sum<-cbind(abundance_nbin_sum,abundance_nbin)
}
abundance_nbin_sum<-abundance_nbin_sum[,8:ncol(abundance_nbin_sum)]
names(abundance_nbin_sum)<-c("Abundance_nbin1","Abundance_nbin2","Abundance_nbin3","Abundance_nbin4","Abundance_nbin5","Abundance_nbin6","Abundance_nbin7")
result_df<-cbind(count_ab,abundance_nbin_sum)
return(result_df)
}

nbinom_distr<-nbinom_distr_ab_func(count_ab)
write_delim(nbinom_distr,"SRR12183113_filtered_rnbinom_size3.txt")

