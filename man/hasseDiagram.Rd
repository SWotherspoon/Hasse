% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Hasse.R
\name{hasseDiagram}
\alias{hasseDiagram}
\title{Hasse Diagram}
\usage{
hasseDiagram(formula, random = NULL, data = NULL)
}
\arguments{
\item{formula}{A formula object.}

\item{random}{A character vector of terms considered random.}

\item{data}{A data frame.}
}
\value{
DiagrammeR script for the Hasse diagram.
}
\description{
Create a Hasse Diagram.
}
\examples{
library(DiagrammeR)
df <- expand.grid(A=factor(1:2),B=factor(1:7),Rep=factor(1:2))
grViz(hasseDiagram(~ A*B,random="B",data=df))

}
\references{
Oehlert, G. W. (2010). A first course in design and
analysis of experiments.
}
