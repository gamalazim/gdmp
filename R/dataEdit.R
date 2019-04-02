#' @import methods
#' @import utils
#' @useDynLib gdmp


#' @export
toArray <- function(finalRep) {

finalRep[,1] <- as.character(finalRep[,1]); finalRep[,2] <- as.character(finalRep[,2]);
snpNo <- nlevels(as.factor(finalRep[,1]))
animalNo <- nlevels(as.factor(finalRep[,2]))
k <- ifelse(animalNo > 2, sample(2:animalNo, 1), 2)

if(animalNo > 1) {
  if(!all(finalRep[1:snpNo, 1] == finalRep[(k-1)*snpNo + (1:snpNo), 1]))
    stop("Sequence of SNPs is not unique across animals")
}

if(sum(diff(table(finalRep[,2])))) stop("Variable number of SNPs per animal")
snps <- finalRep[1:snpNo,1]; animals <- unique(finalRep[,2])

t( array(paste(finalRep[,3],finalRep[,4],sep=""), dim=c(snpNo, animalNo),
      dimnames = list(snps, animals)) )

} # End toArray

## Append two arrays -- add new data
## #################################
#' @export
arrayAppend = function(topA, bottomA, missingVal = 5) {
#
# Old Data = topA
# New Data = bottomA
# missingVal = value to use for missing SNPs in new data, default = 5
# For differing SNPs, go by the SNP list of the old data
#

  if( ncol(topA) == ncol(bottomA) ) {
    ##if(all.equal(colnames(topA), colnames(bottomA))==TRUE)
    if(setequal(colnames(topA), colnames(bottomA)))
      topA <- rbind(topA, bottomA)
    else {
      ## 1. Files that follow may contain 4275 SNPs without the 'ARS-' prefix
      ## The prefix exists in Illumina files such as allele report and final
      ## SNP list. Therefore, the prefix will be added.
      u <- which((!is.element(colnames(bottomA), colnames(topA))))
      colnames(bottomA)[u] <- paste("ARS-", colnames(bottomA)[u], sep="")
      if(sum(is.element(colnames(bottomA), colnames(topA))) != dim(topA)[2])
        stop("Something other than ARS-, chech and interfere manually!")

      ## 2. Order of SNPs is not the same
      topA = rbind(topA, bottomA[,colnames(topA), drop=FALSE])
    }
  }
  else {
    ## Exclude SNPs in bottomA (new data) that are not in topA (old data)
    u <- which((is.element(colnames(bottomA), colnames(topA))))
    bottomA <- bottomA[,u, drop=FALSE]

    ## if bottomA is for a single animal
    if(dim(as.matrix(bottomA))[2] == 1) bottomA <- t(as.matrix(bottomA))

    ## Add SNPs in topA (old data) that are not in bottomA (new data)
    u <- which((!is.element(colnames(topA), colnames(bottomA))))
    tA <- array(missingVal, c(nrow(bottomA), length(u)), dimnames=list(rownames(bottomA), colnames(topA)[u]))
    bottomA <- cbind( bottomA, tA )
    topA = rbind(topA, bottomA[,colnames(topA), drop=FALSE])
  }

  topA
}

## Get 2 unique bases present at a SNP
## ###################################
getBases <- function(snpG, all.char = FALSE) {
## Get 2 unique bases across all animals present at a SNP
## snpG is a column vector in genotype array 'ga' with "AA", "AG", "GA", "-A", "--", whaterver ...
## Option 'all.char = TRUE' prints out all characters not only elements of (A,G,C,T)
  x <- paste(unique(snpG), collapse ='')
  nx <- nchar(x)
  y <- character(nx)
  for(i in 1:nx) y[i] <- substr(x, i, i)
  if(!all.char) sort(unique(y[is.element(y, c("A","G","C","T"))]))
  else sort(unique(y))
} # End getBases

##
##
## WARNING:
## VERY SLOW VERSION
## Not in NAMESPACE
## Not documented
## Recode base genotypes as 0, 1, 2
## #################################
snpRecode.nodesig <- function(snpG) {
# nodesig: no designations needed
## Recode snp genotypes by counting number of allele A
## snpG is a column vector in genotype array 'ga' with "AA", "AG", "GA", "-A", "--", whaterver ...
## Unknown genotypes are those with non A/G/C/T bases, code those as 5
## Takes about 2.1 seconds per animal for the 50K Chip, about 1/2 hr is needed to recode 800 animals!

  b <- getBases(snpG)
  n <- length(snpG)
  x <- integer(n)

  if(length(b) == 0) { ## missing all snp genotypes
    x <- rep.int(5, n)
  }

  if(length(b) == 1) { ## X-linked, or non polymorphic snps, code as 2
    for(i in 1:n) {
      if(substr(snpG[i],1,1) == b && substr(snpG[i],2,2) == b) x[i] <- 2
      else x[i] <- 5  ## unknown genotype
    }
  }

  if(length(b) == 2) { ## regular autosomal snp
    for(i in 1:n) {
      alleleA <- substr(snpG[i],1,1)
      alleleB <- substr(snpG[i],2,2)
      if(!is.element(alleleA, b) || !is.element(alleleB, b)) x[i] <- 5
      else { if(alleleA == b[1]) x[i] <- 1; if(alleleB == b[1]) x[i] <- x[i]+1; }
    }
  }
  x
} # Slow snpRecode


## Recode SNP genotypes as 0, 1, 2
## #################################
#' @export
snpRecode <- function(snpG, designat) {

  b <- getBases(snpG)
  n <- length(snpG)
  x <- rep.int(5, n)

  if(length(b) == 1) { ## X-linked, or non polymorphic snps, code as 0 or 2
    if(b == designat[2]) { BB <- paste(b, b, sep=""); x[snpG == BB] <- 0 }
    if(b == designat[1]) { AA <- paste(b, b, sep=""); x[snpG == AA] <- 2 }
  }

  if(length(b) == 2) { ## regular autosomal snp
    if(b[1] != designat[1]) {
      tempb = b; b[1] = tempb[2]; b[2] = tempb[1]
      ## if Designated and actual bases mismatch use unkown '--'
      if(b[1] != designat[1]) b[1] <- b[2] <- '-'
    }
    AA <- paste(b[1], b[1], sep="")
    AB <- paste(b[1], b[2], sep="")
    BA <- paste(b[2], b[1], sep="")
    BB <- paste(b[2], b[2], sep="")

    x[snpG == AA] <- 2
    x[snpG == AB] <- 1
    x[snpG == BA] <- 1
    x[snpG == BB] <- 0
  }
  x
} # End snpRecode

#' @export
is.hwEq <- function(snpG, diff) {

  g <- table(snpG)
  g <- g[c("0","1","2")]
  g <- g[!is.na(g)] ## just in case of a missing 0, 1, or 2
  if(length(g) < 3) hwEq <- FALSE
  else {
    g <- g/sum(g)
    if(abs(sqrt(4*g[1]*g[3]) - g[2]) > diff) hwEq <- FALSE
    else hwEq <- TRUE
  }
  hwEq
} # End is.hwEq

#' @export
snpSelect <- function(ga.r, select.method=c("call.rate", "heterozygosity", "HW.Eq", "MAF"),
  call.rate.min = .85, hz.min = .01, hz.diff = .15, MAF.min = .005) {

  select.method <- match.arg(select.method)

  if(select.method == "call.rate") {
    poorThreshold = round((1 - call.rate.min)*(dim(ga.r)[2]))
    poor <- apply(ga.r, 1, function(x) table(x)["5"])
    ga.r[poor < poorThreshold,]
  }
  else if(select.method == "heterozygosity") {
    hzThreshold = round(hz.min*(dim(ga.r)[1]))
    hzNumber <- apply(ga.r, 2, function(x) table(x)["1"])
    hzNumber[is.na(hzNumber)] <- 0
    ga.r[,(hzNumber > hzThreshold)]
  }
  else if(select.method == "HW.Eq") {
    keep <- apply(ga.r, 2, function(x) is.hwEq(x, hz.diff))
    ga.r[,keep]
  }
  else {
    maf <-  apply(ga.r, 2, getMAF)
    ga.r[,(maf > MAF.min)]
  }
} # End snpSelect

#' @export
is.identical <- function(x, y, allow = .005) {

  n <- length(x)
  k <- round(allow*n)
  result <- 0
  if(.C("cIdentical", as.integer(x), as.integer(y), as.integer(k), as.integer(n),
  result = as.integer(result), PACKAGE="gdmp")$result) TRUE
  else FALSE
}

#' @export
GetHCS <- function(ga.r, Exclude=1:ncol(ga.r), allow = .005) {

  ga.r <- ga.r[,Exclude]
  n <- dim(ga.r)[2]
  l <- dim(ga.r)[1]
  k <- round(allow*l)
  iex <- .C("cGetHCS", as.integer(ga.r),iex=as.integer(rep.int(1,n)),as.integer(l),as.integer(n),as.integer(k), PACKAGE="gdmp")$iex
  Exclude[iex == 0]
}

#' @export
getMAF <- function(snpG) {
## Calculate Minor Allele Frequency for a single SNP
  maf <- 0
  g <- table(snpG)[c("0","1","2")]
  g <- g[!is.na(g)] ## just in case of a missing 0, 1, or 2
  if(length(g) == 3) maf <- (2*g[1] + g[2])
  else if(length(g) == 2) {
    if(!is.na(g["1"])) maf <- g["1"]
    if(!is.na(g["0"])) maf <- maf + (2*g["0"])
    else if(!is.na(g["2"])) maf <- maf + (2*g["2"])
  }

  maf <- maf/(2*sum(g))
  ifelse(maf <= .5, maf, 1-maf)
} # End getMAF



##
## The following functions are not documented and not included in NAMESPACE
##
toBeagle <- function(finalRep) {
## Function to turn final report file into an input BEAGLE file for the purpose of imputing
## Conditions for data in final report 'finalRep':
## 1. final report is a data frame with first 5 columns listed in the following order:
##     a. snp name, as factor, with equal number of snps per animal
##     b. animal id as factor
##     c. allele 1 (one character: A, C, G, or T)
##     d. allele 2 (one character: A, C, G, or T)
##     e. gencall
## 2. all snps of animal 1 are listed first followed by snps of animal 2, and so on.
## 3. snps are listed for each animal in the same order


snpNo <- nlevels(as.factor(finalRep[,1]))
animalNo <- nlevels(as.factor(finalRep[,2]))
k <- sample(2:(animalNo-1), 1)

## Check condition #3
if(animalNo > 1) {
  if(!all(finalRep[1:snpNo, 1] == finalRep[((k*snpNo+1):((k+1)*snpNo)), 1]))
    stop("Sequence of SNPs is not unique across animals")
}

## Check condition #1
if(sum(diff(table(finalRep[,2])))) stop("Variable number of SNPs per animal")
snps <- as.character(finalRep[1:snpNo,1]); animals <- as.character(unique(finalRep[,2]))

bf <- finalRep[1:snpNo, 3:4]
if(animalNo > 1) {
  for(i in 2:animalNo) bf <- cbind(bf, finalRep[((i-1)*snpNo+1):(i*snpNo), 3:4])
}

bf <- cbind(snps, bf)
bf <- cbind(rep("M", snpNo), bf)
colnames(bf) <- c("I", "id", rep(animals, each=2))
bf

} # End toBeagle
## x <- toBeagle(finalReport)
## write.table(x, file="./x.bgl", sep=" ", quote=F, row.names=F)


## Append two beagle-formatted arrays - add individuals 'horizontally'
## ###################################################################
biAppend = function(leftA, rightA) {

  ##if(all.equal(leftA[,2], rightA[,2])==TRUE)
  if(setequal(leftA[,2], rightA[,2]))
    leftA <- cbind(leftA, rightA[,-(1:2)])
  else {
    ## 1. Files that follow may contain 4275 SNPs without the 'ARS-' prefix
    ## The prefix exists in Illumina files such as allele report and final
    ## SNP list -- the prefix continues to version 2.
    u <- which((!is.element(rightA[,2], leftA[,2])))
    rightA[u,2] <- paste("ARS-", rightA[u,2], sep="")
    if(sum(is.element(rightA[,2], leftA[,2])) != dim(leftA)[1])
      stop("Something other than ARS-, chech and interfere manually!")

    ## 2. Order of SNPs is not the same
    leftA = cbind(leftA, rightA[leftA[,2], -(1:2)])
  }
  leftA
}

## Append two beagle-formatted arrays -- add markers 'vertically'
## ##############################################################
bmAppend = function(topA, bottomA) {

  ##if(all.equal(colnames(topA), colnames(bottomA))==TRUE)
  if(setequal(colnames(topA), colnames(bottomA)))
    topA <- rbind(topA, bottomA)
  else stop("To add markers, individuals in topA and bottomA must all be equal")
  topA
}

#
# Function to read phased/imputed beagle output file
# ##################################################
readBeagle <- function(bfile) {
  imputed <- read.table(file=bfile, header=T, check.names=F)
  imputed
}

#
# Function to transform beagle format to array format
# ###################################################
b2Array <- function(imputed) {
  row.names(imputed) <- imputed$id
  imputed$I <- imputed$id <- NULL

  Nsamples <- ncol(imputed)/2
  Nsnps <- nrow(imputed)

  sampleIDs <- colnames(imputed)[seq(1, 2*Nsamples, 2)]
  snpNames <- row.names(imputed)

  a <- matrix(0, nrow=Nsnps, ncol=Nsamples, dimnames=list(snpNames,sampleIDs))
  for(i in 1:Nsamples) a[,i] <- paste(imputed[,(2*i-1)], imputed[,(2*i)], sep="")

  t(a)
}

