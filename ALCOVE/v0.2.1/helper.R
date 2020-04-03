
project_assign <- read.table("/local/projects/RNASEQ/karthik_rnaseq_gui/data/RNASEQ.pipelines.txt", stringsAsFactors = FALSE, header = T)
alignment = read.table(file = "/local/projects/RNASEQ/karthik_rnaseq_gui/test/02_Alignment/Summary.txt", header = T)
alignment <- data.frame(alignment)

options(digits = 5, device = "png")
par()
dev.off()

tweaks <-
   list(tags$head(tags$style(HTML("
                                 .multicol {
                                   height: 150px;
                                   -webkit-column-count:3; /* Chrome, Safari, Opera */
                                   -moz-column-count: 4;    /* Firefox */
                                   column-count: 3;

                                   -column-fill: balance;
                                   -webkit-column-fill: balance;

                                 }
                                 "))
   ))

fastqc_dir <- list.dirs("/local/projects/RNASEQ/karthik_rnaseq_gui/test/01_FastQC/10706343630_fastqc/i1/", full.names = T, recursive = F)
fastqc_dir <- mixedsort(fastqc_dir)

fastqc_dir_names <- list.dirs("/local/projects/RNASEQ/karthik_rnaseq_gui/test/01_FastQC/10706343630_fastqc/i1/", full.names = F, recursive = F)
fastqc_dir_names <- mixedsort(fastqc_dir_names)

fastqc_files <- NULL

for (i in 1:length(fastqc_dir)){
   fastqc_files[i] <- list.dirs(path = fastqc_dir[i], recursive = F, full.names = T)[1]
}

fastqc_filename <- vector()
for(i in 1:length(fastqc_dir)){
   fastqc_filename[i] <- list.dirs(path = fastqc_dir[i], recursive = F, full.names = F)[1]
}

project_options <- strapplyc(fastqc_filename, "(.*)_1_1", simplify = T)
# project_options <- sort(project_options)


names(fastqc_dir_names) <- project_options


c(lapply(1:length(project_options), function(i){ project_options[i] == fastqc_dir_names[i]}))

gene_Exp_dir <- list.dirs("/local/projects/RNASEQ/karthik_rnaseq_gui/test/03_Gene_Expressiion/10706343630_exon_counts/i1", full.names = T)

file_path <- list.files(path = gene_Exp_dir[1], pattern = "*.counts")

for (i in 2:length(gene_Exp_dir)){
   file_path[i]<- list.files(path = gene_Exp_dir[i], pattern = "*.counts", full.names = T)
}

file_name <- vector()
for (i in 2:length(gene_Exp_dir)){
   file_name[i]<- list.files(path = gene_Exp_dir[i], pattern = "*.counts", full.names = F)
}



datalist = lapply(file_path[-1], FUN = "read.table", header = F)

merged.data = Reduce(function(...) merge(..., by = "V1", all = T, sort = F), datalist)
merged.data = merged.data[1:(nrow(merged.data)-5),]

col_names <-strapplyc(file_name[-1], "(.*).accepted", simplify = T)

names(merged.data) = c("Gene_ID",col_names)
merged.data <- merged.data[,order(names(merged.data))]

merged.data.filtered <- subset(merged.data, rowSums((merged.data[,2:ncol(merged.data)])/(ncol(merged.data)-1)) >= 1)


#merged.data.filtered<- merged.data.filtered[rowSums(is.na(merged.data.filtered)) != (ncol(merged.data.filtered)-1), ]

merged.data_mat <- data.matrix(merged.data.filtered[,-1])


merged.data.norm <- (merged.data_mat %*% diag(1/colSums(merged.data_mat))) * 1000000
merged.data.norm <- apply(merged.data.norm,2, function(i) {round(i, digits = 4)})


merged.data.norm.df <- data.frame(merged.data.filtered[,1], merged.data.norm)
names(merged.data.norm.df) <- names(merged.data)
merged.data.norm.df[merged.data.norm.df == 0] <- as.numeric(min(merged.data.norm.df[merged.data.norm.df !=0 ]))
merged.data.norm.full <- merged.data.norm.df

m_merged.data <- melt(merged.data.filtered, id = "Gene_ID")
names(m_merged.data)[2:3] <- c("Sample_ID", "Read_Count")

m_merged.data <- m_merged.data[order(m_merged.data$Sample_ID),]

#merged.data.norm.full[merged.data.norm.full < 1] <- NA
merged.data.norm.full <- merged.data.norm.full[rowSums(merged.data.norm.full[,-1] < 1) != (ncol(merged.data.norm.full)-1),]
#merged.data.norm.full[is.na(merged.data.norm.full)] <- 0

m_merged.data.norm.full <- melt(merged.data.norm.full, id ="Gene_ID")
names(m_merged.data.norm.full)[2:3] <- c("Sample_ID", "CPM")

merged.data.filtered[,1] <- paste("<a href='","http://useast.ensembl.org/Multi/Search/Results?q=", merged.data.filtered[,1],";site=ensemblgene","'target='_blank'>",merged.data.filtered[,1] ,"</a>",sep = "")

merged.data.norm.df[,1] <- paste("<a href='","http://useast.ensembl.org/Multi/Search/Results?q=", merged.data.norm.df[,1],";site=ensemblgene","'target='_blank'>",merged.data.norm.df[,1] ,"</a>",sep = "")

diffeq_dir <- list.dirs("/local/projects/RNASEQ/karthik_rnaseq_gui/test/04_Differential_Gene_Expression/10706343630_differential_expression/i1/", full.names = T, recursive = F)
diffeq_dir <- mixedsort(diffeq_dir)

diffeq_dir_names <- list.dirs("/local/projects/RNASEQ/karthik_rnaseq_gui/test/04_Differential_Gene_Expression/10706343630_differential_expression/i1/", full.names = F, recursive = F)
diffeq_dir_names <- mixedsort(diffeq_dir_names)

diffeq_files <- NULL

for (i in 1:length(diffeq_dir)){
   diffeq_files[i] <- list.files(path = diffeq_dir[i],pattern = "*.de_genes.txt" ,recursive = F, full.names = T)
}

diffeq_filename <- vector()
for(i in 1:length(diffeq_dir)){
   diffeq_filename[i] <- list.files(path = diffeq_dir[i],pattern = "*.de_genes.txt" , recursive = F, full.names = F)
}

diffeq_files_count <- NULL

for (i in 1:length(diffeq_dir)){
   diffeq_files_count[i] <- list.files(path = diffeq_dir[i],pattern = "all_counts_noZero" ,recursive = F, full.names = T)
}

diffeq_filename <- unlist(strapplyc(diffeq_filename, "(.*).de_genes.txt"))

names(diffeq_dir)<- diffeq_filename


infinity_eliminator <- function(x){
   x[is.infinite(x) & x < 0] <- floor(min(x[!is.infinite(x)], 0))
   x[is.infinite(x) & x > 0] <- ceiling(max(x[!is.infinite(x)],0))
   x
}
