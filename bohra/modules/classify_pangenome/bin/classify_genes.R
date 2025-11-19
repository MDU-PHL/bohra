#! /usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(optparse)

option_list <- list( 
  make_option(c("-p", "--presence_absence"),  type="character", metavar = "FILE", 
              help="Required: presence_absence.Rtab from pan-genome anlaysis (Roary, Panaroo...)"),
  make_option(c("-g", "--grouping"),  type="character",metavar = "FILE", 
              help="Required: tab separated file, with the first column the genome name, second column its group"),
  make_option(c("-o", "--output_directory"), type="character", default="out/", 
              help="Directory for output files [default %default]",
              metavar="DIR"),
  make_option(c("-s", "--min_size"), type="integer", default=10, 
              help="Minimum number of genomes per group to be considered in analysis [default %default]",
              metavar="NUMBER"),
  make_option(c("-c", "--core_threshold"), type="double", default=0.95, 
              help="Threshold used to define a core gene within each group [default %default]",
              metavar="FLOAT"),
  make_option(c("-r", "--rare_threshold"), type="double", default=0.15, 
              help="Threshold used to define a rare gene within each group [default %default]",
              metavar="FLOAT")
)


opt = parse_args(OptionParser(option_list=option_list))

## input params
input_presence_absence_file = opt$presence_absence
input_grouping_file = opt$grouping
if (is.null(input_presence_absence_file) || is.null(input_grouping_file)) {
  stop("Must provide both a presence absence file with '-p' and grouping file with '-g'")
}
out = opt$output_directory
min_size = opt$min_size
core_thresold = opt$core_threshold
if (core_thresold > 1 || core_thresold < 0) { stop("Core threshold must be between 0 and 1")}
rare_threshold = opt$rare_threshold
if (rare_threshold > 1 || rare_threshold < 0) { stop("Rare threshold must be between 0 and 1")}
if (rare_threshold >= core_thresold) { stop("Rare threshold must be smaller than core threshold") }

## check if output directory exists, and create if not
if (!dir.exists(out)) { dir.create(out) }

output_frequnecy_matrix = file.path(out,"frequencies.csv")
output_classification_table = file.path(out,"classification.tab")
output_genes_per_isolate = file.path(out, "genes_per_isolate.tab")

## read in the roary/panaroo presence absence file (Rtab)
print('Reading presence absence file...')
complete_presence_absence = fread(input_presence_absence_file, sep = "\t", header = T, stringsAsFactors = F)
## read in the grouping file, tab separated, no header
grouping = read.table(input_grouping_file, sep =  "\t", comment.char = "", stringsAsFactors = F, header = F)
colnames(grouping) = c("ID","grouping")


## check if all the genomes in groups are in presence_absence
missing =  which(!grouping$ID %in% colnames(complete_presence_absence))


if (length(missing) > 0)  {
  names_missing = grouping$ID[missing]
  print(names_missing)
  grouping = grouping[-missing,]
  warning(paste("The IDs of",length(missing),"names in the grouping file do not match the names in the presence absence file! (Written above this warning.)
                \nIf all should be included, please correct these names to match between the two files."))
  ## remove the missing IDs from the grouping file
  
}


## The column with the genome names is called ID in my casgroups = unique(grouping$grouping)  ## names of groups
group_sizes = table(grouping$grouping)

groups_to_keep = names(which(group_sizes >= min_size)) # ignoring anything smaller than min_size
num_groups = length(groups_to_keep) ## count how many groups there are

# if (num_groups < 2) { stop("Not enough groups to compared with at least min_size isolates! Try decreasing min_size to compare smaller groups.") }

## extract the names of the genes, they should be in the first column
gene_names = as.character(unlist(complete_presence_absence[,1]))

## create a dataframe where each group has the frequency of each gene
group_freqs = data.frame(matrix(ncol = num_groups, nrow = length(gene_names)))
colnames(group_freqs) = groups_to_keep

## create a vector of frequencies for each group
for (group in groups_to_keep) {
  print(paste("Calculating frequency vector for group: ", group, "...", sep = ""))
  # print(group)

  curr_columns = which(colnames(complete_presence_absence) %in% grouping$ID[grouping$grouping == group])
  # print(curr_columns)
  curr_presence_absence = data.frame(complete_presence_absence[,..curr_columns])
  # group_freqs[is.na(group_freqs)] <- 0
  # print(head(curr_presence_absence))
  group_freqs[,which(colnames(group_freqs) == group)] = rowSums(curr_presence_absence)/dim(curr_presence_absence)[2]
}

## save the frequency table -> this is quite useful
rownames(group_freqs) = gene_names

write.table(x = data.frame(Gene = gene_names, group_freqs),
            file = output_frequnecy_matrix, sep = ',', col.names = T, row.names = F, quote = F)

# ## when debugging -> read in a table
# group_freqs = read.table(file = "~/cholera_club/process_clustering/out/frequencies.csv", sep = ",", comment.char = "", stringsAsFactors = F,
#            header = T)

# initiate a dataframe for the output
classification = data.frame(gene_name = gene_names,
                            core = rep(0, length(gene_names)),
                            inter = rep(0, length(gene_names)),
                            rare = rep(0, length(gene_names)),
                            total = rep(0, length(gene_names)),
                            details = rep("", length(gene_names)), stringsAsFactors = F   )


## loop to count for each gene -> this takes a while
for (i in 1:dim(classification)[1]) {
  print(paste("Calculating values for gene in index: ", i, " of ", dim(classification)[1],"...", sep = ""))
  curr_freqs = group_freqs[i,] ## get the frequencies of current gene
  ## count how many core/inter/rare
  core = colnames(curr_freqs)[which(curr_freqs >= core_thresold)]
  inter = colnames(curr_freqs)[which(curr_freqs < core_thresold & curr_freqs >= rare_threshold)]
  rare = colnames(curr_freqs)[which(curr_freqs < rare_threshold & curr_freqs > 0)]
  
  ## I've added this for you so you can see where it's core/rare/inter
  classification$details[i] = paste("Core:", paste(core, collapse = "+"), 
                                    "Inter:", paste(inter, collapse = "+"), 
                                    "Rare:", paste(rare, collapse = "+"), sep = " ")
  ## update the output dataframe
  classification$core[i] = length(core)
  classification$inter[i] = length(inter)
  classification$rare[i] = length(rare)
  classification$total[i] = length(c(core, rare, inter))
}

## now the genes can be classified based on their total presence
print(head(classification))
## First classify generall into varied, core, inter and rare (see my schematic)
classification$general_class = rep("Varied", dim(classification)[1]) ## initiate the default 
classification$general_class[which(classification$core == classification$total)] = "Core"  ##always core
classification$general_class[which(classification$inter == classification$total)] = "Intermediate" ## always inter
classification$general_class[which(classification$rare == classification$total)] = "Rare" ## always rare
classification$general_class[which(classification$total == 0)] = "Absent in labelled lineages"
# print(head(classification))
# ## Now sub-classify based on the precise combinations -> as you can see it's not very clever!
classification$specific_class = rep("Core, intermediate and rare", dim(classification)[1])
classification$specific_class[which(classification$general_class == "Core" & !(classification$total %in% c(1,num_groups)))] = "Multi-lineage core"
classification$specific_class[which(classification$general_class == "Core" & classification$total == num_groups)] = "Collection core"
classification$specific_class[which(classification$general_class == "Core" & classification$total == 1)] = "Lineage specific core"
classification$specific_class[which(classification$general_class == "Intermediate" & classification$total == num_groups)] = "Collection intermediate"
classification$specific_class[which(classification$general_class == "Intermediate" & classification$total == 1)] = "Lineage specific intermediate"
classification$specific_class[which(classification$general_class == "Intermediate" & classification$total != 1 & classification$total != num_groups)] = "Multi-lineage intermediate"
classification$specific_class[which(classification$general_class == "Rare" & classification$total == num_groups)] = "Collection rare"
classification$specific_class[which(classification$general_class == "Rare" & classification$total == 1)] = "Lineage specific rare"
classification$specific_class[which(classification$general_class == "Rare" & classification$total != 1 & classification$total != num_groups)] = "Multi-lineage rare"
classification$specific_class[which(classification$general_class == "Varied" & classification$rare == 0)] = "Core and intermediate"
classification$specific_class[which(classification$general_class == "Varied" & classification$inter == 0)] = "Core and rare"
classification$specific_class[which(classification$general_class == "Varied" & classification$core == 0)] = "Intermediate and rare"
classification$specific_class[which(classification$general_class == "Absent in labelled lineages")] = "Absent in labelled lineages"
print(head(classification))

# ## finally, save the file
write.table(classification, file = output_classification_table, sep = "\t", col.names = T, row.names = F, quote = F)

# # classification = read.table("~/cholera_club/process_clustering/out/classification.tab", sep = "\t", header = T,
# #                             comment.char = "", stringsAsFactors = F)

# ## Calculate typical values for the whole dataset, and for each lineage
# genes_per_isolate = data.frame(group = character(0),
#                              class = character(0),
#                              count= numeric(0), stringsAsFactors = F)

# for (curr in groups_to_keep) {
#   print(curr)
#   sample_names = grouping$ID[which(grouping$grouping == curr)]
#   samples = which(colnames(complete_presence_absence) %in% sample_names)
#   curr_presence_absence = data.frame(complete_presence_absence[,..samples])
#   for (gene_class in unique(classification$specific_class)) {
#     print(gene_class)
#     gene_class_presence_absence = curr_presence_absence[which(gene_names %in% classification$gene_name[which(classification$specific_class == gene_class)]),]
#     print(gene_class_presence_absence)
#     num_genes = colSums(gene_class_presence_absence)
#     genes_per_isolate = rbind(genes_per_isolate ,
#                             data.frame(cluster = rep(curr, length(num_genes)),
#                                        class = rep(gene_class, length(num_genes)),
#                                        count = num_genes, stringsAsFactors = F))
    
#   }
# }  

# write.table(genes_per_isolate, file = output_genes_per_isolate, sep = "\t", col.names = T, row.names = T, quote = F)
# # genes_per_isolate = read.table("~/cholera_club/process_clustering/out/genes_per_isolate.tab", sep = "\t",
# #                                comment.char = "", stringsAsFactors = F, header = T)

# ############# PLOTTING ############# 


### generate all the output plots and measurements
plots_out = file.path(out,"plots/")
if (!dir.exists(plots_out)) { dir.create(plots_out) }

colours = data.frame(
  Class = c( "Lineage specific core","Multi-lineage core", "Collection core","Lineage specific intermediate",
             "Multi-lineage intermediate","Collection intermediate","Lineage specific rare", "Multi-lineage rare" ,
             "Collection rare", "Intermediate and rare","Core, intermediate and rare","Core and rare", "Core and intermediate",
             "Absent in labelled lineages"),
  Colour = c("#542788","#8c96c6","#08519c","#fa9fb5","#c51b8a","#7d1158",
             "#fec44f","#d95f0e","#b1300b","#edf8e9","#bae4b3","#74c476","#238b45", "#d3d3d3"), stringsAsFactors = F
)

classification$label = paste(classification$total, "/", num_groups ,sep ="")
classification$label = factor(classification$label, paste(1:num_groups, "/",num_groups,sep="")) 
classification$specific_class = factor(classification$specific_class, colours$Class)

# ## how many from each class, stratified
# A = ggplot(classification, aes(x = label, fill = specific_class)) + geom_bar(color = "black", lwd = 0.2) +
#   facet_grid(general_class~., scales = "free") + 
#   xlab("") + theme_classic(base_size = 12) + ylab("Number of genes") +
#   scale_fill_manual(values = colours$Colour, drop = F, name = "Distribution\nclass")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   xlab("Number of lineages in which gene is present") 
# ggsave(A, filename = file.path(plots_out, "barplot_per_lineage_count.pdf"), height = 10, width= 10)
print("Plot B")
# ## total counts per class
B = ggplot(classification, aes(x = general_class, fill = specific_class)) + geom_bar( color = "black", lwd = 0.2) +
  scale_fill_manual(values = colours$Colour, "Distribution class", drop = F) + theme_classic(base_size = 12) +
  scale_y_continuous( expand = c(0,0,0.1,0)) +
  ylab("Genes") + xlab("Distribution class") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave(B, filename = file.path(plots_out, "total_class_counts.pdf"), height = 5, width= 6)


# ## hex plot showing plain

# ## first calculate means
# mean_without_zeros <- function(x) {
#   x = as.numeric(x)
#   zeros = which(x == 0)
#   if (length(zeros > 0 )) {
#     return (mean(x[-zeros]))
#   }
#   return(mean(x))
# }
# classification$means =  apply(X = group_freqs, 1, FUN = mean_without_zeros)
# classification$means[classification$total == 0] = NA
# classification$jitter_total = jitter(classification$total, amount = 0.1)
# C =  ggplot(classification, aes(y = means, x = jitter_total, fill = specific_class)) + geom_hex(bins = 70)+
#   scale_fill_manual(values = colours$Colour, name = "Distribution\nclass", drop = F) + theme_classic(base_size = 12) +
#   ylab("Mean frequency\nwhen present") + xlab("Number of lineages in\nwhich gene is present") + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_x_continuous(breaks = 1:num_groups)
# ggsave(C, filename = file.path(plots_out, "hex_plain.pdf"), height = 5, width= 7)


# #### typical E. coli
# mean_per_lineage = aggregate(x = genes_per_isolate$count, by = list(genes_per_isolate$cluster, genes_per_isolate$class), median)
# mean_all = aggregate(mean_per_lineage$x, by = list(mean_per_lineage$Group.2), median)
# colnames(mean_all) = c("Class","Count")
# mean_all$Count = round(mean_all$Count , digits = 0)
# write.table(mean_all, file = file.path(plots_out, "typical_genome_counts.tab"), sep = "\t", col.names = T, row.names = F, quote = F)


# mean_all$Class = factor(mean_all$Class, colours$Class)
# D = ggplot(mean_all,aes(fill = Class, x= "", y = Count))+ geom_bar(stat = "identity", color = "black", lwd = 0.1) + 
#   coord_polar("y", start=0) +
#   scale_fill_manual(values = colours$Colour,  drop = F, name = "Distribution\nclass") +  theme_minimal()+ 
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.border = element_blank(),
#     panel.grid=element_blank(),
#     axis.ticks = element_blank()) +
#   theme(axis.text.x=element_blank())
# ggsave(D, filename = file.path(plots_out, "typical_genome.pdf"), height = 5, width= 7)


# ## typical per lineage
# plots_out = file.path(out,"plots/typical_per_class/")
# if (!dir.exists(plots_out)) { dir.create(plots_out) }

# for (curr_class in unique(genes_per_isolate$class)) {
#   curr = genes_per_isolate[genes_per_isolate$class == curr_class,]
#   curr$cluster = factor(curr$cluster, names(sort(group_sizes, decreasing = T)))
#   med_all = median(aggregate(by = list(curr$cluster), x = curr$count, FUN = median)$x)
#   p = ggplot(curr, aes(x = cluster, y = count)) + geom_boxplot(fill = "#eeeeee") +
#     theme_classic(base_size = 12) + ylab(paste("Number of  '", curr_class, "'  genes\nper genome", sep = "")) +
#     xlab("Group") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     geom_hline(yintercept = med_all, col = "red")
#   out_file =  file.path(plots_out, paste(gsub(curr_class, pattern =  " ", replacement = "_"), ".pdf",sep = ""))
#   ggsave(plot = p, filename = out_file,
#          height = 5, width= 7)
# }

# ## PCA plots
# create_pca_plot <- function(class, comp = F){
#   if (class %in% c("Absent in labelled lineages","Collection core",
#                    "Lineage specific rare", "Lineage specific core" , "Lineage specific intermediate" )){
#       return()}
#   curr_freqs = group_freqs[which(rownames(group_freqs) %in% classification$gene_name[classification$specific_class == class]),]
#   for_pca = t(curr_freqs)
#   remove = c()
#   for (i in 1:dim(for_pca)[2]) {
#     if (length(unique(for_pca[,i])) == 1) { ## no variation in gene
#       remove = c(remove, i)
#     }
#   }
#   if (length(remove) > 0) {
#     for_pca = for_pca[,-remove]
#   }
#   ### PCA plot of the clusters -> what are the relationships between the clusters based on the frequencies of all genes
#   freqs.pca = prcomp(for_pca , center = T)
#   summary_pca = summary(freqs.pca)
#   importance = round(summary_pca$importance[2,1:2] * 100, digits = 2)
#   freqs.pca = data.frame(freqs.pca$x)
#   freqs.pca$label = rownames(freqs.pca)
#   p = ggplot(freqs.pca, aes(x = PC1, y = PC2, label = label)) + geom_text() +
#     theme_classic(base_size = 12) +
#     xlab(paste("PC1 (", importance[1], "%)", sep = "")) + 
#     ylab(paste("PC2 (", importance[2], "%)", sep = "")) 
#   out_file =  file.path(plots_out, paste(gsub(class, pattern =  " ", replacement = "_"), ".pdf",sep = ""))
#   ggsave(plot = p, filename = out_file,
#          height = 5.5, width= 6)
#   return() 
  
# }

# plots_out = file.path(out,"plots/pca_per_class/")
# if (!dir.exists(plots_out)) { dir.create(plots_out) }
# sapply(unique(genes_per_isolate$class), FUN = create_pca_plot)







