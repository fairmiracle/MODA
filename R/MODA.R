#' Illustration of partition density
#' 
#' Calculate the average density of all resulting modules from a partition. The 
#' density of each module is defined as the average adjacency of the module 
#' genes.
#' 
#' @param ADJ gene similarity matrix 
#' @param PartitionSet vector indicates the partition label for genes
#' 
#' @return partition density, defined as average density of all modules
#' 
#' @references Langfelder, Peter, and Steve Horvath. "WGCNA: an R package for 
#' weighted correlation network analysis." BMC bioinformatics 9.1 (2008): 1.
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords density
#' 
#' @examples
#' data(synthetic)
#' ADJ1=abs(cor(datExpr1,use="p"))^10
#' dissADJ=1-ADJ1
#' hierADJ=hclust(as.dist(dissADJ), method="average" )
#' groups <- cutree(hierADJ, h = 0.8)
#' pDensity <- PartitionDensity(ADJ1,groups) 
#' @export
#' 
PartitionDensity <- function(ADJ, PartitionSet){
    labels <- unique(PartitionSet)
    numP <- length(labels)
    pdensity <- 0
    for (i in seq_len(numP)){
        idx = which(PartitionSet == labels[i])
        
        if (length(idx) == 1){
            pdensity <- pdensity + 0
        }else{
            Aq = ADJ[idx,idx]
            sum = sum(Aq)
            nrow = dim(Aq)[1]
            pdensity <- pdensity + sum*(sum-nrow)/(nrow*nrow-nrow)
        }
        
    }
    pdensity <- pdensity*2/sum(ADJ)
    pdensity
}

#' Illustration of modularity density
#'
#' Calculate the average modularity of a partition. The modularity of each 
#' module is defined from a natural generalization of unweighted case.
#'
#' @param ADJ gene similarity matrix 
#' @param PartitionSet vector indicates the partition label for genes
#' 
#' @return partition modularity, defined as average modularity of all modules
#' 
#' @references Newman, Mark EJ. "Analysis of weighted networks." Physical 
#' review E 70.5 (2004): 056131.
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords modularity
#' 
#' @examples
#' data(synthetic)
#' ADJ1=abs(cor(datExpr1,use="p"))^10
#' dissADJ=1-ADJ1
#' hierADJ=hclust(as.dist(dissADJ), method="average" )
#' groups <- cutree(hierADJ, h = 0.8)
#' pDensity <- PartitionModularity(ADJ1,groups) 
#' 
#' @export
#'
PartitionModularity <- function(ADJ, PartitionSet){
    labels <- unique(PartitionSet)
    numP <- length(labels)
    m = sum(ADJ)
    K = colSums(ADJ)
    pModularity <- 0
    for (i in seq_len(numP)){
        idx = which(PartitionSet == labels[i])
        
        if (length(idx) == 1){
            pModularity <- pModularity + 0
        }else{
            Aq = ADJ[idx,idx]
            deg = K[idx]
            pModularity = pModularity + sum(Aq)- sum(deg*deg)/m
            
        }
    }
    pModularity <- pModularity/m
    pModularity
}


#' Modules detection by hierarchical clustering
#' 
#' Module detection based on the optimal cutting height of dendrogram, which is 
#' selected to make the average density or modularity of resulting partition
#' maximal. The clustering and visulization function are from WGCNA.
#'
#' @param datExpr gene expression profile, rows are samples and columns genes
#' @param foldername where to store the clusters
#' @param indicatename normally a specific tag of condition
#' @param cutmethod cutting the dendrogram based on maximal average Density 
#' or Modularity
#' @param power the power parameter of WGCNA, W_{ij}=|cor(x_i,x_j)|^power
#'
#' @return The number of clusters
#' 
#' @references Langfelder, Peter, and Steve Horvath. "WGCNA: an R package for 
#' weighted correlation network analysis." BMC bioinformatics 9.1 (2008): 1.
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords cutting dendrogram
#' @seealso \code{\link{PartitionDensity}}
#' @seealso \code{\link{PartitionModularity}}
#' 
#' @import WGCNA
#' @import igraph
#' @importFrom dynamicTreeCut cutreeDynamic
#' @examples
#' data(synthetic)
#' ResultFolder = 'ForSynthetic' # where middle files are stored
#' CuttingCriterion = 'Density' # could be Density or Modularity
#' indicator1 = 'X'     # indicator for data profile 1
#' indicator2 = 'Y'      # indicator for data profile 2
#' specificTheta = 0.1 #threshold to define condition specific modules
#' conservedTheta = 0.1#threshold to define conserved modules
#' intModules1 <- WeightedModulePartitionHierarchical(datExpr1,ResultFolder,
#' indicator1,CuttingCriterion) 
#' mymodule <- getPartition(ResultFolder)
#' randIndex(table(mymodule,truemodule),adjust=F)

#' @export
#' 
WeightedModulePartitionHierarchical <- function(datExpr,foldername,indicatename,
                        cutmethod=c('Density','Modularity'), power=10){
    dir.create(file.path('./', foldername), showWarnings = FALSE)
    
    ADJ1=abs(cor(datExpr,use="p"))^power
    dissADJ=1-ADJ1
    hierADJ=hclust(as.dist(dissADJ), method="average" )
    #NumCutHeights = 27440 for DaphIX, too much, just change
    #NumCutHeights <- length(hierADJ$height)
    cutHeights <- seq(0.05,1,by=0.05)
    NumCutHeights <- length(cutHeights)
    pDensity <- numeric(length <- NumCutHeights)
    maxpDensity <- 0
    maxHeight <- 0
    for (i in 1:NumCutHeights) {
        #groups <- cutree(hierADJ, h = hierADJ$height[i]) # cut tree into 5
        groups <- cutree(hierADJ, h = cutHeights[i])
        if (cutmethod == 'Density'){
            pDensity[i] <- PartitionDensity(ADJ1,groups)
        }else{
            pDensity[i] <- PartitionModularity(ADJ1,groups)
        }
        
        if (pDensity[i] > maxpDensity){
            maxpDensity <- pDensity[i]
            maxHeight <- cutHeights[i]
        }
    }

    # We like large modules, so we set the minimum module size relatively high:
    minModuleSize = 30;
    # Module identification using dynamic tree cut:
    dynamicMods = cutreeDynamic(dendro = hierADJ, distM = dissADJ, cutHeight = 
                            maxHeight,deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
    
    dynamicColors = labels2colors(dynamicMods)
    intModules = table(dynamicColors)
     
    #######Visualization#############
    #pdf(paste(foldername,"/Partitions_",indicatename,".pdf",sep=""),width = 11, 
    #    height = 8)
    png(paste(foldername,"/Partitions_",indicatename,".png",sep=""))
    marAll = c(1, 5, 3, 1)
    layout(matrix(c(1,2,3,0), 2, 2, byrow = TRUE), widths=c(0.8,0.2),
           heights=c(0.8,0.2))
    par(mar = c(0, marAll[2], marAll[3], marAll[4]))
    plot(hierADJ, labels = FALSE, xlab="", sub="", ylim=c(0,1),hang = -1,
         axes=FALSE,main = "Gene hierarchical clustering dendrogram")
    abline(h = maxHeight, col = 'red')
    axis(side = 2, at = seq(1,0,-0.2),labels =seq(1,0,-0.2))
    par(mar = c(0, 1, marAll[3], marAll[4]),mgp = c(0, 1, 0))
    plot(pDensity, cutHeights, type = 'n', ylim=c(0,1), axes=FALSE,ylab = "")
    title(paste('Average',cutmethod), font.main = 1, cex.main = 1,line  = 0)
    axis(side = 1, at	= seq(0,(max(pDensity)+0.1),0.1),
         labels = seq(0,(max(pDensity)+0.1),0.1))
    axis(side = 2, at	= seq(1,0,-0.2),labels =seq(1,0,-0.2))
    lines(pDensity, cutHeights, lwd = 2, col = 'blue')
    abline(h = maxHeight, col = 'red')
    par(mar = c(marAll[1], marAll[2], 0, marAll[4]))
    plotColorUnderTree(hierADJ, dynamicColors, groupLabels = NULL, 
                       rowText = NULL, rowLabels = 'Clusters')
    dev.off()
    
    for(J in 1:length(intModules)){
        idx <- which(dynamicColors == names(intModules)[J])
        DenseGenes = colnames(datExpr)[idx]
        densegenefile <- paste(foldername,"/DenseModuleGene_",
                               indicatename,"_",J,".txt",sep="")
        write.table(DenseGenes,densegenefile,sep = "\n",col.names = FALSE,
                    row.names = FALSE,quote = FALSE)
    }
    
    return (length(intModules))
}

#' Modules detection by spectral clustering
#' 
#' Module detection based on the spectral clustering algorithm, which mainly
#' solve the eigendecomposition on Laplacian matrix
#'
#' @param datExpr gene expression profile, rows are samples and columns genes
#' @param foldername where to store the clusters
#' @param indicatename normally a specific tag of condition
#' @param GeneNames normally the gene official names to replace the colnames of datExpr
#' @param power the power parameter of WGCNA, W_{ij}=|cor(x_i,x_j)|^power
#' @param nn the number of nearest neighbor, used to construct the affinity matrix
#' @param k the number of clusters(modules)
#'
#' @return None
#' 
#' @references Von Luxburg, Ulrike. "A tutorial on spectral clustering." 
#' Statistics and computing 17.4 (2007): 395-416.
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords cutting dendrogram
#' 
#'
#' @import igraph
#' @import cluster
#' @examples
#' data(synthetic)
#' ResultFolder <- 'ForSynthetic' # where middle files are stored
#' indicator <- 'X'     # indicator for data profile 1
#' GeneNames <- colnames(datExpr1)
#' WeightedModulePartitionSpectral(datExpr1,ResultFolder,indicator,
#' GeneNames,k=5)
#' truemodule <- c(rep(1,100),rep(2,100),rep(3,100),rep(4,100),rep(5,100))
#' mymodule <- getPartition(ResultFolder)
#' randIndex(table(mymodule,truemodule),adjust=F)
#' @export
#' 
WeightedModulePartitionSpectral <- function(datExpr, foldername, indicatename, 
                                    GeneNames, power=6, nn=10, k=2){
    dir.create(file.path('./', foldername), showWarnings = FALSE)
    ADJ1=abs(cor(datExpr,use="p"))^power
 
    W = TOMsimilarity(ADJ1)   # similarity matrix
    A = make.affinity(W,nn)
    d <- apply(A, 1, sum)
    L <- diag(d)-A                        # unnormalized version
    L <- diag(d^-0.5)%*%L%*% diag(d^-0.5) # normalized version
    evL <- eigen(L,symmetric=TRUE)  # evL$values is decreasing sorted when symmetric=TRUE
    Z <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
    spc <- pam(Z,k)
    colorSpectralTom <- labels2colors(spc$cluster)
    intModules = table(colorSpectralTom)
    for(J in 1:length(intModules)){
        idx <- which(colorSpectralTom == names(intModules)[J])
        DenseGenes = GeneNames[idx]
        densegenefile <- paste(foldername,"/DenseModuleGene_",
                               indicatename,"_",J,".txt",sep="")
        write.table(DenseGenes,densegenefile,sep = "\n",col.names = FALSE,
                    row.names = FALSE,quote = FALSE)
    }                                    
}

#' Modules detection by Louvain algorithm
#' 
#' Module detection based on the Louvain algorithm, which tries to maximize 
#' overall modularity of resulting partition.
#'
#' @param datExpr gene expression profile, rows are samples and columns genes
#' @param foldername where to store the clusters
#' @param indicatename normally a specific tag of condition
#' @param GeneNames normally the gene official names to replace the colnames of datExpr
#' @param maxsize the maximal nodes allowed in one module
#' @param minsize the minimal nodes allowed in one module
#' @param power the power parameter of WGCNA, W_{ij}=|cor(x_i,x_j)|^power
#' @param tao the threshold to cut the adjacency matrix
#'
#' @return The number of clusters
#' 
#' @references Blondel, Vincent D., et al. "Fast unfolding of communities in 
#' large networks." Journal of statistical mechanics: theory and experiment 
#' 2008.10 (2008): P10008.
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords cutting dendrogram
#' 
#'
#' @import igraph
#' @examples
#' data(synthetic)
#' ResultFolder <- 'ForSynthetic' # where middle files are stored
#' indicator <- 'X'     # indicator for data profile 1
#' GeneNames <- colnames(datExpr1)
#' intModules1 <- WeightedModulePartitionLouvain(datExpr1,ResultFolder,indicator,GeneNames)
#' truemodule <- c(rep(1,100),rep(2,100),rep(3,100),rep(4,100),rep(5,100))
#' mymodule <- getPartition(ResultFolder)
#' randIndex(table(mymodule,truemodule),adjust=F)
#' @export
#' 
WeightedModulePartitionLouvain <- function(datExpr,foldername,indicatename,GeneNames,
                                  maxsize=200, minsize=30, power=6, tao=0.2){
    ADJ <- abs(cor(datExpr,use="p"))^power
    ADJ[ADJ < tao] <- 0
    g <- graph_from_adjacency_matrix(ADJ,mode='undirected',weighted=TRUE)
    V(g)$name=1:length(V(g))
    dir.create(foldername, showWarnings = FALSE)
    tmfile <- paste(foldername,'tmp.txt',sep='')
    recursiveigraph(g,tmfile,'louvain',maxsize, minsize)
    num <- modulesRank(foldername,indicatename,GeneNames)
    num
}

#' Modules detection by AMOUNTAIN algorithm
#' 
#' Module detection based on the AMOUNTAIN algorithm, which tries to find the  
#' optimal module every time and use a modules extraction way
#'
#' @param datExpr gene expression profile, rows are samples and columns genes
#' @param foldername where to store the clusters
#' @param GeneNames normally the gene official names to replace the colnames of datExpr
#' @param Nmodule the number of clusters(modules)
#' @param maxsize the maximal nodes allowed in one module
#' @param minsize the minimal nodes allowed in one module
#' @param power the power parameter of WGCNA, W_{ij}=|cor(x_i,x_j)|^pwr
#' @param tao the threshold to cut the adjacency matrix
#' 
#' @return None
#' 
#' @references Blondel, Vincent D., et al. "Fast unfolding of communities in 
#' large networks." Journal of statistical mechanics: theory and experiment 
#' 2008.10 (2008): P10008.
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords optimization
#' 
#' @import AMOUNTAIN
#' @examples
#' data(synthetic)
#' ResultFolder <- 'ForSynthetic' # where middle files are stored
#' GeneNames <- colnames(datExpr1)
#' intModules1 <- WeightedModulePartitionAmoutain(datExpr1,5,ResultFolder,
#' GeneNames,maxsize=100,minsize=50)
#' truemodule <- c(rep(1,100),rep(2,100),rep(3,100),rep(4,100),rep(5,100))
#' mymodule <- getPartition(ResultFolder)
#' randIndex(table(mymodule,truemodule),adjust=F)
#' @export
#' 
WeightedModulePartitionAmoutain <- function(datExpr,Nmodule,foldername,GeneNames,
                                    maxsize=200, minsize=3, power=6, tao=0.2){
    W <- abs(cor(datExpr,use="p"))^power
    W[W < tao] <- 0
    N = dim(W)[1]
    z <- rep(0,N)
    idx <- 1
    dir.create(foldername, showWarnings = FALSE)
    saveAtomfile <- paste(foldername,'/AtomModule',sep='')
    for (ii in 1:Nmodule) {
        abegin = 0.01
        aend = 0.9
        for (i in 1:20) {
            #x <- CGPFixSS(W,z,rep(1/N,N),a=(abegin+aend)/2,lambda=0,maxiter = 50)
            x <- moduleIdentificationGPFixSS(W,z,rep(1/N,N),a=(abegin+aend)/2,
                                             lambda=0,maxiter = 50)
            predictedid <- which(x[[2]]!=0)
            if(length(predictedid) > maxsize){
                abegin = (abegin+aend)/2
            }else if (length(predictedid) < minsize){
                aend = (abegin+aend)/2
            }else
                break
        }
        
        if(length(predictedid) <= maxsize){
            modulescore = sum(W[predictedid,predictedid])
            write.table(GeneNames[predictedid],file = paste(foldername,'/',
                        floor(modulescore),'-moduleid-',idx,'.txt',sep=''),
                        quote = FALSE, row.names = FALSE, col.names = FALSE)
            idx <- idx+1
        }
        else if(length(predictedid) > maxsize){
            modulescoreW = W[predictedid,predictedid]
            print(paste('Atom! with size',length(predictedid),sep=' '))
            tmpstr = as.numeric(GeneNames[predictedid])-1
            forTotalcompletegraph(tmpstr,modulescoreW,saveAtomfile)
        }
        W = W[-predictedid,-predictedid]
        GeneNames = GeneNames[-predictedid]
        z = z[-predictedid]
        N = length(GeneNames)
        print(paste('Finishing module ',ii,sep=''))
        
        if(N < 3 | sum(W)==0)
            break
    }
}

#' Illustration of two networks comparison
#' 
#' Compare the background network and a condition-specific network. A Jaccard
#' index is used to measure the similarity of two sets, which represents the 
#' similarity of each module pairs from two networks.
#'
#' @param sourcehead prefix of where to store results
#' @param nm1 how many modules in the background network
#' @param nm2 how many modules in the condition-specific network
#' @param ind1 indicator of the background network
#' @param ind2 indicator of the condition-specific network
#'
#' @return A matrix where each entry is the Jaccard index of corresponding 
#' modules from two networks
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords module comparison
#' 
#' @examples
#' data(synthetic)
#' ResultFolder = 'ForSynthetic' # where middle files are stored
#' CuttingCriterion = 'Density' # could be Density or Modularity
#' indicator1 = 'X'     # indicator for data profile 1
#' indicator2 = 'Y'      # indicator for data profile 2
#' intModules1 <- WeightedModulePartitionHierarchical(datExpr1,ResultFolder,
#' indicator1,CuttingCriterion) 
#' intModules2 <- WeightedModulePartitionHierarchical(datExpr2,ResultFolder,
#' indicator2,CuttingCriterion) 
#' JaccardMatrix <- comparemodulestwonets(ResultFolder,intModules1,intModules2,
#' paste('/DenseModuleGene_',indicator1,sep=''),
#' paste('/DenseModuleGene_',indicator2,sep=''))
#' 
#' @export
#' 
comparemodulestwonets <- function(sourcehead,nm1,nm2,ind1,ind2){
    my.array<-array(0,dim=c(nm1,nm2))
    for(i1 in 1:nm1){
        for (i2 in 1:nm2){
            # sourcehead = "Networks/DenseModuleGene_"
            densegenefile1 <- paste(sourcehead,ind1,"_",i1,".txt",sep="")
            densegenefile2 <- paste(sourcehead,ind2,"_",i2,".txt",sep="")
            list1 <- readLines(densegenefile1)
            list2 <- readLines(densegenefile2)
            my.array[i1,i2] = length(intersect(list1,list2))/
                length(union(list1,list2))
        }
    }
    return (my.array)
}

#' Illustration of network comparison
#' 
#' Compare the background network and a set of condition-specific network. 
#' Conserved or condition-specific modules are indicated by the plain files, 
#' based on the statistics 
#'
#' @param ResultFolder where to store results
#' @param intModules how many modules in the background network
#' @param indicator identifier of current profile, served as a tag in name
#' @param intconditionModules a numeric vector, each of them is the number 
#' of modules in each condition-specific network. Or just single number
#' @param conditionNames a character vector, each of them is the name 
#' of condition. Or just single name
#' @param specificTheta the threshold to define min(s)+specificTheta, 
#' less than which is considered as condition specific module. 
#' s is the sums of rows in Jaccard index matrix. See supplementary file. 
#' @param conservedTheta The threshold to define max(s)-conservedTheta, 
#' greater than which is considered as condition conserved module. 
#' s is the sums of rows in Jaccard index matrix. See supplementary file.
#'
#' @return None
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @seealso \code{\link{WeightedModulePartitionHierarchical}},
#' \code{\link{comparemodulestwonets}}
#' @keywords module differential
#' 
#' @examples
#' data(synthetic)
#' ResultFolder = 'ForSynthetic' # where middle files are stored
#' CuttingCriterion = 'Density' # could be Density or Modularity
#' indicator1 = 'X'     # indicator for data profile 1
#' indicator2 = 'Y'      # indicator for data profile 2
#' specificTheta = 0.1 #threshold to define condition specific modules
#' conservedTheta = 0.1#threshold to define conserved modules
#' intModules1 <- WeightedModulePartitionHierarchical(datExpr1,ResultFolder,
#' indicator1,CuttingCriterion) 
#' intModules2 <- WeightedModulePartitionHierarchical(datExpr2,ResultFolder,
#' indicator2,CuttingCriterion) 
#' CompareAllNets(ResultFolder,intModules1,indicator1,intModules2,indicator2,
#' specificTheta,conservedTheta)
#' 
#' @export
#' 
CompareAllNets <- function(ResultFolder,intModules,indicator,
                          intconditionModules,conditionNames,specificTheta,
                          conservedTheta){
    for (i in 1:length(conditionNames)) {
        ArrayGroup1 <- comparemodulestwonets(ResultFolder,intModules,
                    intconditionModules[i],paste('/DenseModuleGene_',
                        indicator,sep=''),paste('/DenseModuleGene_',
                    conditionNames[i],sep=''))
        dir.create(paste(ResultFolder,'/',conditionNames[i],sep=''), showWarnings = FALSE)
        fileprefix <- paste(ResultFolder,'/',conditionNames[i],'/',sep='')

        #pdf(paste(fileprefix,'module_overlap_remove',conditionNames[i],
        #          '.pdf',sep=''),width = 10, height = 8)
        png(paste(fileprefix,'module_overlap_remove',conditionNames[i],
                  '.png',sep=''))
        plot(1:intModules,rowSums(ArrayGroup1),xlim = c(0,(intModules + 1)),
             ylim = c(0,max(rowSums(ArrayGroup1)) + 0.1),
             xlab="Module ID",ylab = 'RowSums of jaccard matrix', 
             type = 'n',main="Overlapped ratio")
        lines(1:intModules, rowSums(ArrayGroup1), col = 'black', type = "b")
        abline(h = min(rowSums(ArrayGroup1)) + 0.1, col='red')
        abline(h = max(rowSums(ArrayGroup1)) - 0.1, col='green')
        text(5, min(rowSums(ArrayGroup1)) - 0.05, "condition response", 
             col = "red",cex = 0.75)
        text(5, max(rowSums(ArrayGroup1)) + 0.05, "conserved response", 
             col = "green",cex = 0.75)
        dev.off()
        specificmoduleid <- which(rowSums(ArrayGroup1) <= 
                                min(rowSums(ArrayGroup1)) + specificTheta)
        conservedmoduleid <- which(rowSums(ArrayGroup1) >= 
                                max(rowSums(ArrayGroup1)) - conservedTheta)
        write.table(specificmoduleid,file = paste(fileprefix,
                    'sepcificModuleid.txt',sep=''),
                    row.names = FALSE, col.names = FALSE)
        write.table(conservedmoduleid,file = paste(fileprefix,
                    'conservedModuleid.txt',sep=''),
                    row.names = FALSE, col.names = FALSE)
    }
}

#' Statistics of all conditions
#' 
#' Statistics of all conditions. To highlight conserved or condition-specific 
#' by counting how frequent each module is lablelled as which, and then 
#' visualize the frequency by bar plot. 
#'
#' @param ResultFolder where to store results
#' @param intModules how many modules in the background network
#' @param indicator identifier of current profile, served as a tag in name
#' @param conditionNames a character vector, each of them is the name 
#' of condition. Or just single name
#' @return None
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @seealso \code{\link{WeightedModulePartitionHierarchical}},
#' \code{\link{WeightedModulePartitionLouvain}},
#' \code{\link{WeightedModulePartitionSpectral}},
#' \code{\link{WeightedModulePartitionAmoutain}},
#' \code{\link{CompareAllNets}}
#' @keywords module differential Statistics
#' 
#' @import RColorBrewer
#' 
#' @export
#' 
ModuleFrequency <- function(ResultFolder,intModules, conditionNames, indicator){
    Ncon <- length(conditionNames)
    wide <- matrix (0,nrow = Ncon, ncol = intModules)
    for (i in 1:Ncon) {
        fileprefix <- paste(ResultFolder,'/',conditionNames[i],'/',sep='')
        tad <- read.table(paste(fileprefix,'sepcificModuleid.txt',sep=''),header = TRUE, sep = "")
        moduleid <-tad[,1]
        wide[i,moduleid] <- 1
    }
    
    # sequential <- brewer.pal(length(conditionNames), "Greens")
    # sequential <- brewer.pal(length(conditionNames), col=c('black', 'white', 'blue', 'red', 'yellow', 'purple', 'green'))
    # sequential <- c('black', 'white', 'blue', 'red', 'yellow', 'purple', 'green')
    sequential <- brewer.pal(Ncon, "Set1")
    pdf(paste(ResultFolder,'/specificmembership.pdf',sep=''),width = 10, height = 5)
    #  png('specificmembership.png',width = 1000, height = 500)
    barplot(wide,
            names.arg = 1:intModules,
            cex.names = 0.7, # makes x-axis labels small enough to show all
            col = sequential, # colors
            xlab = "module id",
            ylab = "if in condition specific network",
            xlim = c(0,intModules+1), # these two lines allow space for the legend
            main = 'Frequency of condition specific module',
            width = 0.75) # these two lines allow space for the legend
    legend("topright", 
           legend = conditionNames, #in order from top to bottom
           fill = sequential, # 6:1 reorders so legend order matches graph
           title = "conditons",cex = 0.75)
    dev.off()
    
    for (ispec in 1:Ncon) {
        idx1 <- intersect(which(colSums(wide)!=0), which(colSums(wide) <= ispec))
        uniqueaffiliation <- vector('list',Ncon)
        names (uniqueaffiliation) <- conditionNames
        for (i in 1:length(idx1)) {
            tmp <- which(wide[,idx1[i]] == 1)
            for (j in 1:length(tmp)) {
                uniqueaffiliation[[tmp[j]]] <- union(uniqueaffiliation[[tmp[j]]],idx1[i])
            }
        }
        write(paste('If we allow at most ', ispec ,' share one condition\n',sep=''),
              file = paste(ResultFolder,'/specificmembership.txt',sep=''),append=TRUE)
        for (i in 1:Ncon) {
            write(conditionNames[i],file = paste(ResultFolder,'/specificmembership.txt',sep=''),append=TRUE)
            if(is.null(uniqueaffiliation[[i]])){
                write('NULL',file = paste(ResultFolder,'/specificmembership.txt',sep=''),append = TRUE)
            }else{
                write(uniqueaffiliation[[i]],file = paste(ResultFolder,'/specificmembership.txt',sep=''),append = TRUE)
            }
            write('\n',file = paste(ResultFolder,'/specificmembership.txt',sep=''),append=TRUE)
        }
    }
    idx1 <- which(colSums(wide)!=0)
    
    wide <- matrix (0,nrow = Ncon, ncol = intModules)
    for (i in 1:length(conditionNames)) {
        fileprefix <- paste(ResultFolder,'/',conditionNames[i],'/',sep='')
        tad <- read.table(paste(fileprefix,'conservedModuleid.txt',sep=''),header = TRUE, sep = "")
        moduleid <- tad[,1]
        wide[i,moduleid] <- 1
    }
    #sequential <- c('black', 'white', 'blue', 'red', 'yellow', 'purple', 'green')
    sequential <- brewer.pal(Ncon, "Set1")
    pdf(paste(ResultFolder,'/conservedmembership.pdf',sep=''),width = 10, height = 5)
    #  png('specificmembership.png',width = 1000, height = 500)
    barplot(wide,
            names.arg = 1:intModules,
            cex.names = 0.7, # makes x-axis labels small enough to show all
            col = sequential, # colors
            xlab = "module id",
            ylab = "if in condition specific network",
            xlim = c(0,intModules+1), # these two lines allow space for the legend
            main = 'Frequency of conserved module',
            width = 0.75) # these two lines allow space for the legend
    legend("topright", 
           legend =conditionNames, #in order from top to bottom
           fill = sequential, # 6:1 reorders so legend order matches graph
           title = "conditons",cex = 0.75)
    dev.off()
    
    idx2 <- which(colSums(wide)!=0)
    #write(idx2,file = paste(ResultFolder,'/conservedmembership.txt',sep=''),append = FALSE)
    write.table(idx2,file = paste(ResultFolder,'/conservedmembership.txt',sep=''),
                row.names = FALSE, col.names = FALSE)
    dir.create(paste(ResultFolder,'/','interestedModules',sep=''))
    for (i in idx1) {
        file.copy(paste(ResultFolder,"/DenseModuleGene_",indicator,"_",i,".txt",sep=""), 
                  paste(ResultFolder,'/','interestedModules',sep=''))
    }
    for (i in idx2) {
        file.copy(paste(ResultFolder,"/DenseModuleGene_",indicator,"_",i,".txt",sep=""), 
                  paste(ResultFolder,'/','interestedModules',sep=''))
    }
}

#' datExpr1
#' 
#' Synthetic gene expression profile with 20 samples and 500 genes.
#'
#' @name datExpr1
#' @docType data
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords data
#' @format A matrix with 20 rows and 500 columns.
#' @examples
#' data(synthetic)
#' ## plot the heatmap of the correlation matrix ...
#' \dontrun{heatmap(cor(as.matrix(datExpr1)))}
NULL

#' datExpr2
#' 
#' Synthetic gene expression profile with 25 samples and 500 genes.
#'
#' @name datExpr2
#' @docType data
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords data
#' @format A matrix with 25 rows and 500 columns.
#' @examples
#' data(synthetic)
#' ## plot the heatmap of the correlation matrix ...
#' \dontrun{heatmap(cor(as.matrix(datExpr2)))}
NULL