% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cdgnet-functions.R
\name{check_WT_biomarkers}
\alias{check_WT_biomarkers}
\title{function that checks whether any of the drugs have biomarkers that are wild-type
if any of the biomarkers are wild-type, then check the MP report
if they are not wild-type there, then remove them from the list}
\usage{
check_WT_biomarkers(drugs_bio, MP_df)
}
\arguments{
\item{drugs_bio}{data frame with drugs and biomarkers}

\item{MP_df}{data frame with molecular alterations}
}
\value{
data frame
}
\description{
function that checks whether any of the drugs have biomarkers that are wild-type
if any of the biomarkers are wild-type, then check the MP report
if they are not wild-type there, then remove them from the list
}
