##
## Functions related to findhap
##


##
## Function to read 'genotypes.filled'
#'@export
read.findhap <- function(Nanim, Nmark, file="./genotypes.filled") {

g <- array(0, (Nanim*Nmark) );
h <- array(0, (Nanim*3));

##
## The C Function, cReadGenotypes.c
##
out <- .C("cReadGenotypes", as.integer(Nanim), as.integer(Nmark),
  h=as.integer(h), g=as.integer(g), as.character(file), PACKAGE="gdmp")

g <- out$g
h <- out$h

cbind(matrix(out$h, nrow=Nanim, ncol=3, byrow=T),
  matrix(out$g, nrow=Nanim, ncol=Nmark, byrow=T))

}
