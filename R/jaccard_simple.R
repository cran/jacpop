#'Calculate pairwise Jaccard similarity matrix
#'
#'Computes pairwise Jaccard similarity matrix from sequencing data and performs
#'PCA on it. The function is specifically useful to detect population stratification in
#'rare variant sequencing data.
#'
#'In order to account for population structure in sequencing data we propose
#'to calculate the Jaccard similarity instead of the genetic covariance
#'between individuals.
#'\deqn{Jac(A,B)=\frac{|A \cap B|}{|A \cup B|}}{Jac(A,B)=intersection(A,B)/union(A,B)}
#'This similarity index is most suitable for sparse data,
#'which is the case, when we restrict our analysis to variants with low minor
#'allele frequencies. The pairwise Jaccard similarity matrix can be further used
#'in Principal Component Analysis.
#'
#'Although the function does basic filtering (singletons, SNPs with missing entries),
#'we recommend to extract a subset of possibly independent SNPs (500k - 1M should be enough)
#' from your initial dataset for population structure identification. You could
#' either extract a random subset of variants or prune your dataset.
#'@param geno A numeric nxm matrix of genotypes. Rows are individuals and columns are variants.
#'The genotypes should be coded 0, 1 and 2. Missing entries are coded as NA. The natural input would be
#'a matrix produced by PLINK using the option --recodeA and removing the first row and the first
#'6 columns.
#'@param pop.label (Optional) A numeric or character vector of n population labels(if known).
#'It is used for plotting purposes.
#'@param plot_it A logical scalar. Should the first 2 principal components be plotted?
#'@param n.pcs A numeric scalar. Number of principal components to extract from the Jaccard
#'similarity matrix. Set to NULL, if you want just the similarity matrix.
#'
#'@return A list of 2 elements:
#'\itemize{
#'  \item A nxn numeric matrix, where the entries are Jaccard similarity indicies between a pair of individuals. The order of individuals corresponds to the order in the input genotype matrix.
#'  \item A data.frame of principal components, which can be further used in an association analysis. The order of individuals corresponds to the order in the input genotype matrix.
#'}
#'@references Prokopenko, D., Hecker, J., Silverman, E., Pagano, M., Noethen, M. M., Dina, C., Lange, C., Fier, H. L. (2015). Utilizing the Jaccard index to reveal population stratification in sequencing data: A simulation study and an application to the 1000 Genomes Project. \emph{Bioinformatics}, 32, 1366-1372.
#'@importFrom graphics legend plot
#'@importFrom stats rbeta runif
#'@examples
#'#####1)Toy example
#'#Simulate genotypes in 2 populations
#'nsnps=10000
#'fst=0.01
#'nind=20
#'maffilter=0.05
#'p<-runif(nsnps,0,maffilter)
#'freq1<-sapply(1:length(p),function(x) rbeta(1,p[x]*(1-fst)/fst,(1-p[x])*(1-fst)/fst))
#'freq2<-sapply(1:length(p),function(x) rbeta(1,p[x]*(1-fst)/fst,(1-p[x])*(1-fst)/fst))
#'
#'pop1<-sapply(1:nsnps, function(x) sample(c(0,1,2),nind,replace=TRUE,
#'  prob=c(((1-freq1[x])^2),(2*freq1[x]*(1-freq1[x])),(freq1[x]^2))))
#'pop2<-sapply(1:nsnps, function(x) sample(c(0,1,2),nind,replace=TRUE,
#'  prob=c(((1-freq2[x])^2),(2*freq2[x]*(1-freq2[x])),(freq2[x]^2))))
#'all<-as.matrix(rbind(pop1,pop2))
#'
#'#Generate the Jaccard similarity index and plot the first 2 principal components
#'res<-generate_pw_jaccard(geno=all,pop.label=c(rep(1,nind),rep(2,nind)))
#'\dontrun{
#'#####2)PLINK files
#'#If you are working with plink files after filtering the dataset consider
#'#to create a genotype count file by using the option --recodeA.
#'#After that remove the first row and the first 6 columns. Now you can
#'# read it in in the following way:
#'geno<-matrix(scan('sample.geno'),nrow=nind,byrow=T)
#'# nind is the number of individuals(rows)
#'}
#'@export
#TODO: add example how to read in big data
###MAIN Function: generates pairwise Jaccard similarity matrix###
#input: genotype matrix
#Optional: 1)number of principal components to output, default=10
#2)plot the first 2 PCs with labels from pop.label, if included
generate_pw_jaccard<-function(geno,pop.label=NULL,n.pcs=10,plot_it=TRUE){
 #jaccard[i,j]=#matches/(#of ones in i + # of ones in j - #matches)


##################Remove SNPS with missing entries
  bad.missing<-sum(is.na(colSums(geno)))
  if(sum(is.na(colSums(geno)))!=0){
    #stop('SNPs with MAF==0')
    #remove SNPs with 0 MAF

    geno<-geno[,-which(is.na(colSums(geno)))]


  }
  cat(bad.missing,"SNPs with missing entries removed.\n")
#################Remove SNPS with 0 MAF
bad<-sum(colSums(geno)==0)
if(sum(colSums(geno)==0)!=0){
	#stop('SNPs with MAF==0')
#remove SNPs with 0 MAF

	geno<-geno[,-which(colSums(geno)==0)]


}
cat(bad,"SNPs with MAF==0 removed.\n")
################Remove singletons
singletons<-sum(colSums(geno)==1)
if(sum(colSums(geno)==1)!=0){
  geno<-geno[,-which(colSums(geno)==1)]
}
cat(singletons,"singletons removed.\n")

cat(ncol(geno),"SNPs left.\n")
 n<-nrow(geno) #number of ind
  m<-ncol(geno) #number of snps
 geno.mod<-geno
 geno.mod[which(geno.mod==2)]<-1

  #A<-geno%*%t(geno)
  A<-tcrossprod(geno.mod) #faster
  if (sum(is.na(A))>0){
	stop('ERROR: check NA')
  }
  ind<-which(!is.na(A),arr.ind=T)
	#expand.grid
  b<- rowSums(geno.mod)
  res<-matrix(A / (b[ind[,1]] + b[ind[,2]] - A),nrow=n,byrow=TRUE)
#  pc.jac<-prcomp(res)$rotation[,1:n.pcs]
###raw PCA on centered JAC:
  if (!is.null(n.pcs)){
    pc.jac<-svd(scale(res,scale=FALSE,center=TRUE),nu=0)[["v"]][,1:n.pcs]
  }else{
    pc.jac<-NULL
  }

  if (plot_it & !is.null(n.pcs)){
	if (!is.null(pop.label)){
		pop.label<-as.factor(pop.label)
		plot(pc.jac[,c(1,2)],col=pop.label,main='Principal components on Jaccard matrix',xlab="PC1",ylab="PC2")
		legend('bottomleft',legend=unique(pop.label),fill=unique(pop.label))
	}else{
		plot(pc.jac[,c(1,2)],main='Principal components on Jaccard matrix',xlab="PC1",ylab="PC2")
	}
  }
  return(list(Jac=res,pcs=pc.jac))
}

