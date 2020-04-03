rm(list=ls())

library(shiny, warn.conflicts = F)
library(DT,warn.conflicts = F)
library(ggplot2, warn.conflicts = F)
library(gsubfn, warn.conflicts = F)
library(reshape2, warn.conflicts = F)
library(gtools, warn.conflicts = F)
library(gplots, warn.conflicts = F)
library(RColorBrewer, warn.conflicts = F)
library(ggfortify, warn.conflicts = F)

# source("helper.R")

source('global.R', local = TRUE)

function(input, output, session) {
	source('server_dropdownmenu.R', local = TRUE)
	
	# source('server_transcriptomics.R', local = TRUE)
	
	source("server_fastqc.R", local = TRUE)
	
	source("server_alignment.R", local = TRUE)
	
	source("server_geneExp.R", local = TRUE)
	
	source("server_diffGeneExp.R", local = TRUE)
	
	# source("server_isoformExp.R", local = TRUE)
	# source("server_diffIsoformExp.R", local = TRUE)
}
