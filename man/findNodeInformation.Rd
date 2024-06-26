% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findNodeInformation.R
\name{findNodeInformation}
\alias{findNodeInformation}
\title{Reports the mean inferential variance for that node and if the transcript
to gene mapping exists provide the genes mapped to the node. Also reports
this information for the descendant nodes if requested by the user.}
\usage{
findNodeInformation(tse, node, type = NULL, txpsAnn = NULL)
}
\arguments{
\item{tse}{TreeSumarizedExperiment obtained as the output of running
\code{buildTSE}}

\item{node}{numeric index of the node whose information is required}

\item{txpsAnn}{Optional, data.frame txp to gene mapping, with rownames set to
transcripts and gene column set to 'gene_id'}
}
\value{
data.frame
}
\description{
Reports the mean inferential variance for that node and if the transcript
to gene mapping exists provide the genes mapped to the node. Also reports
this information for the descendant nodes if requested by the user.
}
\examples{
example(buildTSE)
node <- 300
nodeInf <- findNodeInformation(tse, node = node, type = "children")
print(nodeInf)
nodeInf <- findNodeInformation(tse, node = node, type = NULL)
print(nodeInf)
}
