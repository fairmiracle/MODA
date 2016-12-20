##  utility

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
#'
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
                recursiveigraph(g2,savefile,method)
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

# get identified partitionAssignment, only for synthetic data where gene names are numbers
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