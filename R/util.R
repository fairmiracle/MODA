##  utility

#' Illustration of network comparison by NMI
#' 
#' Compare the background network and a set of condition-specific network. 
#' returning a pair-wise matrix to show the normalized mutual information between
#' each pair of networks in terms of partitioning
#'
#' @param ResultFolder where to store results
#' @param intModules how many modules in the background network
#' @param indicator identifier of current profile, served as a tag in name
#' @param intconditionModules a numeric vector, each of them is the number 
#' of modules in each condition-specific network. Or just single number
#' @param conditionNames a character vector, each of them is the name 
#' of condition. Or just single name
#' @param Nsize The number of genes in total
#' @param legendNames a character vector, each of them is the condition name 
#' showing up in the similarity matrix plot if applicable
#' @param plt a boolean value to indicate whether plot the similarity matrix
#'
#' @return NMI matrix indicating the similarity between each two networks
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @seealso \code{\link{CompareAllNets}}
#' @keywords module differential NMI
#' 
#' 
NMImatrix <- function(ResultFolder,intModules,indicator,intconditionModules,
                      conditionNames, Nsize, legendNames=NULL, plt=FALSE){
    Ncon <- length(conditionNames)
    matnmi <- matrix(0, nrow = Ncon+1,ncol = Ncon+1)
    
    sourcehead <- paste(ResultFolder,'/DenseModuleGeneID_',sep='')
    comm <- numeric(length=Nsize)
    
    for(i1 in 1:intModules){
        densegenefile1 <- paste(sourcehead,indicator,"_",i1,".txt",sep="")
        list1 <- readLines(densegenefile1)
        comm[as.numeric(list1)] <- i1
    }
    
    for(i in 1:Ncon){
        tmpcomm <- numeric(length=Nsize)
        for(i1 in 1:intconditionModules[i]){
            densegenefile1 <- paste(sourcehead,conditionNames[i],"_",i1,".txt",sep="")
            list1 <- readLines(densegenefile1)
            tmpcomm[as.numeric(list1)] <- i1
        }
        matnmi[1,i+1] <- compare(comm, tmpcomm,"nmi")
        matnmi[i+1,1] <- matnmi[1,i+1]
        assign(paste('comm', i, sep=''), tmpcomm)
    }
    
    for(i in 1:(Ncon-1)){
        for (j in (i+1):Ncon) {
            matnmi[i+1,j+1] <- compare(get(paste('comm', i, sep='')), 
                                       get(paste('comm', j, sep='')),"nmi")
            matnmi[j+1,i+1] <- matnmi[i+1,j+1]
        }
    }
    # range01 <- function(x){(x-min(x[x > 0]))/(max(x)-min(x[x > 0]))}
    diag(matnmi) <- 1
    if(plt){
        textMatrix =  paste(signif(matnmi))
        pdf(file = paste(ResultFolder,"/NMIsimilarity.pdf",sep=''), width = 10, height = 7);
        par(mar = c(6, 8.8, 3, 2.2));
        labeledHeatmap(Matrix = matnmi,
                       xLabels = c(indicator,legendNames),
                       yLabels = c(indicator,legendNames),
                       ySymbols = c(indicator,legendNames),
                       colorLabels = FALSE,
                       colors = blueWhiteRed(50),
                       textMatrix = textMatrix,
                       setStdMargins = FALSE,
                       cex.text = 0.5,
                       zlim = c(-1,1),
                       cex.lab.x = 0.75,
                       cex.lab.y = 0.75,
                       main = "Network similarities by NMI")
        dev.off()
    }
    return(matnmi)
}

#' Modules detection by each condition
#' 
#' Module detection on each condition-specific network, which is constructed from
#' all samples but samples belonging to that condition
#'
#' @param datExpr gene expression profile, rows are samples and columns genes, 
#' rowname should contain condition specifier
#' @param conditionNames character vector, each as the condition name
#' @param ResultFolder where to store the clusters
#' @param GeneNames normally the gene official names to replace the colnames of datExpr
#' @param maxsize the maximal nodes allowed in one module
#' @param minsize the minimal nodes allowed in one module
#' 
#' @return a numeric vector, each entry is the number of modules in condition-specific network
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords multiplecondition
#' 
#'
MIcondition <- function(datExpr,conditionNames,ResultFolder,GeneNames,maxsize=100, minsize=30){
    intconditionModules <- numeric(length = length(conditionNames))
    for (i in 1:length(conditionNames)) {
        # user can specify rows to remove each time here according to the data itself
        removeid <- c()
        for (j in 1:length(rownames(datExpr))){
            if(grepl(conditionNames[i],rownames(datExpr)[j]))
                removeid <- union(removeid,j)
        }
        
        # datExprsConditionRemoved is the rest data without a specific stage
        datExprsConditionRemoved <- datExpr[-removeid,]
        ids <- which(duplicated(t(datExprsConditionRemoved))==TRUE)
        if(length(ids) > 0)
            datExprsConditionRemoved <- datExprsConditionRemoved[,-ids]
        
        #clean genes that have identical values across all the samples, 
        #which makes standard deviation zero, leading NAs in correlation matrix
        colSD <- apply(datExprsConditionRemoved, 2, sd)
        ids <- which(colSD==0)
        if(length(ids) > 0)
            datExprsConditionRemoved <- datExprsConditionRemoved[,-ids]
        
        #intconditionModules[i] = WeightedModulePartitionHierarchical(datExprsConditionRemoved,ResultFolder,conditionNames[i],CuttingCriterion)
        #intconditionModules[i] = WeightedModulePartitionDensity(datExprsConditionRemoved,ResultFolder,conditionNames[i],CuttingCriterion)
        intconditionModules[i] <- WeightedModulePartitionLouvain(datExprsConditionRemoved,ResultFolder,conditionNames[i],
                                                                 GeneNames,maxsize=maxsize, minsize=minsize,power=pwd,tao=0.2)
    }
    intconditionModules
}

#' Modules identification by recursive community detection
#' 
#' Modules detection using igraph's community detection algorithms, when the
#' resulted module is larger than expected, it is further devided by the same program
#'
#' @param g igraph object, the network to be partitioned
#' @param savefile plain text, used to store module, each line as a module
#' @param method specify the community detection algorithm
#' @param maxsize maximal module size
#' @param minsize minimal module size
#' @return None
#' @references Blondel, Vincent D., et al. "Fast unfolding of communities in 
#' large networks." Journal of statistical mechanics: theory and experiment 
#' 2008.10 (2008): P10008.
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords community detection
#' 
#' @import igraph

recursiveigraph <- function(g, savefile, method = c('fastgreedy','louvain'),
                            maxsize=200, minsize=30){
    
    if (method == "fastgreedy")
        fc <- cluster_fast_greedy(g)
    else if (method == "louvain")
        fc <- cluster_louvain(g)
    
    memfc <- membership(fc)
    msize <- sizes(fc)
    
    if(length(msize) > 1){
        
        for (i in 1:length(msize)) {
            if(msize[i] <= maxsize & msize[i] >= minsize){
                mgeneids <- V(g)$name[which(memfc==i)]
                #mgeneids <- as.numeric(mgeneids) - 1
                cp = c()
                for (k in 1:length(mgeneids)){
                    cp = paste(cp,mgeneids[k],sep='\t')
                }
                write(cp,file = savefile,append = TRUE)
            } else if (msize[i] > maxsize) {
                #large modules, in recursive way
                ids = which(memfc==i)
                g2 <- induced.subgraph(graph=g,vids=ids)
                recursiveigraph(g2,savefile,method,maxsize,minsize)
            } else {
                next
            }
        }
    }
    else{
        print(paste('Atom! with size',length(V(g)),sep=' '))
        edges <- get.edgelist(g)
        write.table(edges,paste(savefile,'Atomsize',length(V(g)),sep='_'),
                    sep = "\t",row.names = FALSE,col.names = FALSE,
                    quote = FALSE,append = TRUE)      
    }
}

#' Modules rank from recursive communities detection
#' 
#' Assign the module scores by weights, and rank them from highest to lowest
#'
#' @param foldername folder used to save modules
#' @param indicator normally a specific tag of condition
#' @param GeneNames Gene symbols, sometimes we need them instead of probe ids
#' @return The numeber of modules
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @seealso \code{\link{recursiveigraph}}
#' 
#' @import igraph

modulesRank <- function(foldername,indicator,GeneNames){
    
    tmfile <- paste(foldername,'tmp.txt',sep='')
    rlines <- readLines(tmfile)
    file.remove(tmfile)
    for (i in 1:length(rlines)) {
        ap <- strsplit(rlines[i],'\t')[[1]]
        ap <- as.numeric(ap[2:length(ap)])
        #mscore <- sum(W[ap,ap])
        #write.table(GeneNames[match(ap,rownames(W))],file = paste(foldername,'/',floor(mscore),'-moduleid-',i,'.txt',sep=''),
        #            quote = FALSE, row.names = FALSE, col.names = FALSE)
        write.table(GeneNames[ap],file = paste(foldername,'/DenseModuleGene_',indicator,'_',i,'.txt',sep=''),
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
        write.table(ap,file = paste(foldername,'/DenseModuleGeneID_',indicator,'_',i,'.txt',sep=''),
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    length(rlines)
}

# Purity
ClusterPurity <- function(clusters, classes) {
    sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

# Gaussian distance
GaussianDis <- function(x1, x2, alpha=1) {
    exp(- alpha * norm(as.matrix(x1-x2), type="F"))
}

# k-nearest neighorhood to filter affinity
make.affinity <- function(S, n.neighboors=3) {
    N <- length(S[,1])
    if (n.neighboors >= N) {  # fully connected
        A <- S
    } else {
        A <- matrix(rep(0,N^2), ncol=N)
        for(i in 1:N) { # for each line
            # only connect to those points with larger similarity 
            best.similarities <- sort(S[i,], decreasing=TRUE,index.return = TRUE)
            for (j in best.similarities$ix[1:n.neighboors]) {
                A[i,j] <- S[i,j]
                A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
            }
        }
    }
    A  
}

#' Get numeric partition from folder
#' 
#' Get identified partitionAssignment, only for synthetic data where gene names are numbers
#'
#' @param ResultFolder folder used to save modules
#' @return Number of partitions
#' 
getPartition <- function(ResultFolder){
    mymodule <- rep(0,500)
    ResultFiles <- list.files(ResultFolder)
    ResultFiles <- ResultFiles[grepl('.txt',ResultFiles)]
    for (i in 1:length(ResultFiles)){
        ap <- as.numeric(readLines(paste(ResultFolder,'/',ResultFiles[i],sep='')))
        mymodule[ap] <- i
    }
    mymodule
}

# a large clique, even AMOUNTAIN cannot devided it
forTotalcompletegraph <- function(predictedid,modulescoreW,savefile1){
    fixedlen = 60
    for (i in 1:(floor(length(predictedid)/fixedlen))) {
        cp = c()
        for (k in (fixedlen*(i-1)+1):(fixedlen*i)){
            cp = paste(cp,predictedid[k],sep='\t')
        }
        mscore = sum(modulescoreW[(fixedlen*(i-1)+1):(fixedlen*i),(fixedlen*(i-1)+1):(fixedlen*i)])
        write(paste(mscore,cp,sep=''),file = savefile1,append = TRUE)
    }
    
    if( (length(predictedid) - fixedlen*i)>=3){
        cp = c()
        for (k in (fixedlen*i+1):length(predictedid)){
            cp = paste(cp,predictedid[k],sep='\t')
        }
        mscore = sum(modulescoreW[(fixedlen*i+1):length(predictedid),(fixedlen*i+1):length(predictedid)])
        write(paste(mscore,cp,sep=''),file = savefile1,append = TRUE)
    }
}

# in case the program stopped, count the number of modules
countintconditionModules <- function(conditionNames,ResultFolder){
    intconditionModules <- numeric(length = length(conditionNames))
    fileNames <- list.files(ResultFolder)
    for(i in 1:length(conditionNames)){
        intconditionModules[i] <- length(which(grepl(paste('DenseModuleGene_',conditionNames[i],'_',sep=''),
                                                     fileNames)==TRUE))
    }
    intconditionModules
}

# see which module contains DE
HitGenes <- function(degenelist,interesFolder){
    allgenes <-c()
    fileNames <- list.files(interesFolder)
    hitfrequency <- numeric(length=length(fileNames)/3) #DE in module
    uniquenumincharactervector <- unique(gsub("\\D", "", fileNames))
    names(hitfrequency) <- as.character(sort(as.numeric(uniquenumincharactervector)))
    for (fs in fileNames) {
        if (!grepl('Symbol',fs) & !grepl('Entrezid',fs)& !grepl('blast2go',fs)){
            ids <- readLines(paste(interesFolder,'/',fs,sep=''))
            #if(goalGene %in% ids)
            #print(fs)
            hitfrequency[gsub("\\D", "", fs)] <- length(intersect(ids,degenelist))
            allgenes <- c(allgenes,ids)
        }
    }
    return (hitfrequency)
}

# when you have 50k genes, requring package Matrix
# from http://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r