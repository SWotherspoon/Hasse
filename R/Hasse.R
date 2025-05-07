##' Create a Hasse Diagram.
##'
##' @title Hasse Diagram
##' @param formula A formula object.
##' @param random A character vector of terms considered random.
##' @param data A data frame.
##' @return DiagrammeR script for the Hasse diagram.
##' @examples
##' library(DiagrammeR)
##' df <- expand.grid(A=factor(1:2),B=factor(1:7),Rep=factor(1:2))
##' grViz(hasseDiagram(~ A*B,random="B",data=df))
##'
##' @references Oehlert, G. W. (2010). A first course in design and
##'   analysis of experiments.
##' @importFrom stats setNames
##' @export
hasseDiagram <- function(formula,random=NULL,data=NULL) {
  tmsobj <- terms(formula)
  hasse <- hasseGraph(tmsobj,random,data=data)
  fac <- attr(tmsobj,"factors")
  vrs <- rownames(fac)
  tms <- colnames(fac)
  ord <- setNames(attr(tmsobj,"order"),tms)

  scripts <- !is.null(hasse$neffs) && !is.null(hasse$ndfs) && !is.null(hasse$nobs)

  if(scripts) {
    nodes <- paste("    TM",
                   seq_along(hasse$random),
                   " [label=<<TABLE border='0' cellborder='0' cellspacing='-2'>",
                   "<TR><TD rowspan='2'>",
                   ifelse(hasse$random,"(",""),
                   names(hasse$random),
                   ifelse(hasse$random,")",""),
                   "</TD><TD><font point-size='8'>",
                   hasse$neffs,
                   "</font></TD></TR><TR><TD><font point-size='8'>",
                   hasse$ndfs,
                   "</font></TD></TR></TABLE>>]",
                   sep="")
  } else {
    nodes <- paste("    TM",
                   seq_along(hasse$random),
                   " [label=<",
                   ifelse(hasse$random,"(",""),
                   names(hasse$random),
                   ifelse(hasse$random,")",""),
                   ">]",
                   sep="")
  }

  gr <- paste(
    "graph Hasse {",
    "  node [shape=plaintext]",
    "  ## Mean",
    "  {rank=min;",
    if(scripts)
      paste("    M [label=<<TABLE border='0' cellborder='0' cellspacing='-2'>",
            "<TR><TD rowspan='2'>Mean</TD><TD><font point-size='8'>1</font></TD></TR>",
            "<TR><TD><font point-size='8'>1</font></TD></TR></TABLE>>]",
            sep="")
    else
      "    M [label=<Mean>]",
    "  }",sep="\n")

  for(k in unique(ord))
    gr <- paste(gr,"  {rank=same;",paste(nodes[ord==k],collapse="\n"),"  }",sep="\n")

  gr <- paste(gr,
              "  ## Error",
              "  {rank=max;",
              if(scripts)
                paste("    E [label=<<TABLE border='0' cellborder='0' cellspacing='-2'>",
                      "<TR><TD rowspan='2'>(Error)</TD><TD><font point-size='8'>",
                      hasse$nobs,
                      "</font></TD></TR><TR><TD><font point-size='8'>",
                      hasse$nobs-sum(hasse$ndfs)-1,
                      "</font></TD></TR></TABLE>>]",
                      sep="")
              else
                "    E [label=<(Error)>]",
              "  }",
              "  edge [style=solid; color=grey60]",sep="\n")

  gr <- paste(gr,paste("  M -- TM",which(ord==1),sep="",collapse="\n"),sep="\n")
  for(k in seq_along(ord))
    if(any(hasse$hasse[,k]==-1))
      gr <- paste(gr,paste("  TM",k," -- TM",which(hasse$hasse[,k]==-1),sep="",collapse="\n"),sep="\n")
    else
      gr <- paste(gr,paste("  TM",k," -- E",sep="",collapse="\n"),sep="\n")

  paste(gr,"}",sep="\n")
}

##' Create the Hasse graph
##'
##' @param tmsobj A terms object.
##' @param random A character vector of terms considered random.
##' @param restricted should the restricted or unrestricted model be calculated?
##' @param data A data frame.
##' @return A Hasse graph.
##' @importFrom stats setNames
hasseGraph <- function(tmsobj,random=NULL,restricted=TRUE,data=NULL) {
  ## Decompose terms object
  fac <- attr(tmsobj,"factors")
  vrs <- rownames(fac)
  tms <- colnames(fac)
  ord <- setNames(attr(tmsobj,"order"),tms)

  ## Variables declared random
  rvrs <- setNames(vrs %in% random,vrs)
  ## Terms explicitly declared random
  rdm <- setNames(tms %in% random,tms)
  ## Terms that are implied random
  rtms <- setNames(rdm | sapply(seq_along(tms),function(k) any(rvrs & fac[,k]>0)),tms)

  ## Adjacency matrix for the Hasse graph
  hasse <- matrix(0,ncol(fac),ncol(fac),dimnames=list(tms,tms))

  ## Eligble random terms for EMS
  ertms <- matrix(0,ncol(fac),ncol(fac),dimnames=list(tms,tms))

  ## Test denominators
  denom <- setNames(integer(ncol(fac)),tms)

  ## Search up Hasse graph to determine if term i is eligible to
  ## contribute to the ems for term j under the restricted model
  eligible <- function(i,j) {
    ## If explicitly declared random term is eligible
    if(rdm[i]) return(TRUE)

    if(rtms[i]) {
      ## Term is random - search terms above to see if ineligble
      for(k in which(hasse[,i]==1))
        if(!eligible(k,j))
          return(FALSE)
    } else {
      ## Term is fixed - not eligible if not marginal to term j
      if(!all(fac[,i]==0 | fac[,j]>0))
        return(FALSE)
    }
    return(TRUE)
  }

  ## Build the Hasse adjacency matrix
  for(j in seq_len(ncol(fac))) {
    ## Loop over terms of order-1
    for(i in which(ord==ord[j]-1)) {
      ## If i is (strictly) marginal to j
      if(all(fac[,i]==0 | fac[,j]>0)) {
        ## Mark up (decreasing order) and down (increasing order)
        hasse[i,j] <- 1
        hasse[j,i] <- -1
        if(rtms[i]) rtms[j] <- TRUE
      }
    }
    ## Declare nested random terms explicitly random
    if(rtms[j] && sum(hasse[,j]==1) < 2) rdm[j] <- TRUE
  }

  ## Find all random terms below each term
  if(any(rtms))
    for(j in seq_len(ncol(fac))) {
      rs <- integer(0)
      r <- which(rowSums(hasse[,j,drop=FALSE]==-1)>0)
      while(length(r)>0) {
        rs <- unique(c(rs,r))
        r <- which(rowSums(hasse[,r,drop=FALSE]==-1)>0)
      }
      ## Remove ineligible terms if restricted model
      if(restricted && length(rs)>0)
        rs <- rs[sapply(rs,eligible,j=j)]
      ertms[rs,j] <- ord[rs]
    }

  ## Denominator is (unique) highest order eligible term
  for(j in seq_len(ncol(fac))) {
    m <- ertms[,j]
    m <- c(m,max(m)+1)
    is <- which(m==min(m[m>0]))
    denom[j] <- if(length(is)>1) NA else is
  }

  if(!is.null(data)) {
    ## Number of levels for each var
    nlevs <- sapply(data[,vrs],nlevels)
    ## Number effects for each term
    neffs <- apply(ifelse(fac >0,nlevs,1),2,prod)
    ## Add up degrees of freedom by walking hasse graph
    ndfs <- integer(ncol(fac))
    for(j in seq_len(ncol(fac))) {
      dof <- 0
      ts <- which(rowSums(hasse[,j,drop=FALSE]==1)>0)
      while(length(ts)>0) {
        dof <- dof+sum(ndfs[ts])
        ts <- which(rowSums(hasse[,ts,drop=FALSE]==1)>0)
      }
      ndfs[j] <- neffs[j]-dof-1
    }
    neffs <- setNames(neffs,tms)
    ndfs <- setNames(ndfs,tms)
    nobs <- nrow(data)
  } else {
    neffs <- ndfs <- nobs <- NULL
  }

  list(hasse=hasse,ertms=ertms,denom=denom,random=rtms,order=ord,neffs=neffs,ndfs=ndfs,nobs=nobs)
}


##' Create the anova table for a mixed model anova.
##'
##' @title Mixed Model Anova Table
##' @param aov and aov object.
##' @param random A character vector of terms considered random.
##' @param restricted should the restricted or unrestricted model be calculated?
##' @return Aov object with recalculated anova table.
##' @examples
##' library(DiagrammeR)
##' df <- expand.grid(A=factor(1:2),B=factor(1:7),Rep=factor(1:2))
##' df$y <- c(12,14)[df$A]+rnorm(7,0,0.4)[df$B]+rnorm(nrow(df),0,0.3)
##' grViz(hasseDiagram(~ A*B,random="B",data=df))
##' ## A is tested against A:B
##' mixed(aov(y ~ A*B,data=df),random="B")
##'
##' @references Oehlert, G. W. (2010). A first course in design and
##'   analysis of experiments.
##' @importFrom stats pf formula replications terms
##' @export
mixed <- function(aov,random=NULL,restricted=TRUE) {
  ## Test for balance
  if(!is.null(aov$model) &&
     !is.list(replications(formula(aov),aov$model)))
    warning("Design may not be balanced")

  ## Recalculate anova table based on correct test denominator
  hasse <- hasseGraph(terms(aov),random,restricted=restricted)
  smry <- summary(aov)
  for(k in seq_along(smry)) {
    tab <- smry[[k]]
    tab <- cbind(Rdm=1,Den=nrow(tab),tab)
    denom <- pmin(hasse$denom,nrow(tab))
    denom <- ifelse(seq_along(denom)==denom,NA,denom)
    rs <- seq_along(denom)
    tab[rs,"F value"] <- tab[rs,"Mean Sq"]/tab[denom,"Mean Sq"]
    tab[rs,"Pr(>F)"] <- pf(tab[rs,"F value"],tab[rs,"Df"],tab[denom,"Df"],lower.tail = FALSE)
    tab[rs,"Den"] <- denom
    tab[rs,"Rdm"] <- as.integer(hasse$random)
    class(tab) <- class(smry[[k]])
    smry[[k]] <- tab
  }
  smry
}
