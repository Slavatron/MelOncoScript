# PREPARE BOXPLOTS COMPARING MUTATION LOAD FOR EVERY GENE

# DEFINE TARGET GENE LIST
#gene_list = c("BRAF", "NF1", "HRAS")

# READ ANNOTATIONS TABLE PRODUCED BY OncoScript_Sort.R
anno_table = read.table("Onco_Data.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# DEFINE VALUES THAT INDICATE A GENE HAS BEEN MUTATED
vals = c("HOT", "LOF", "MUT") # ignore copy number changes
# COMPILE LIST OF ALL GENES
all_names = names(anno_table)
gene_names = names(anno_table[7:length(all_names)])
#length(gene_names)
temp_mat = as.matrix(anno_table)

# MODIFY DATA FRAME TO SHAVE OFF ALL MUTATION DATA AFTER ";"
for (i in c(7:length(all_names))) {
  x = sapply(strsplit(temp_mat[,i], ";"), '[', 1)
  anno_table[i] = x
}

# RUN FISHER
# DEFINE CONTINGENCY MATRIX
#lb = anno_table[anno_table$RespType == "LB",]
#nb = anno_table[anno_table$RespType == "NB",]
#lbg = dim(lb[lb[,i] %in% vals,])[1]
#lbno = dim(lb[! lb[,i]  %in% vals,])[1]
#nbg = dim(nb[nb[,i] %in% vals,])[1]
#nbno = dim(nb[! nb[,i]  %in% vals,])[1]
#mat = matrix(c(lbg,lbno,nbg,nbno), nrow = 2, ncol = 2, dimnames = list(c("Mutation", "None"),c("LB", "NB")))
#mat
#f = fisher.test(mat)
#pval = f[1]

fishergenes = data.frame(Gene=character(), P_Value = numeric()) # LIST OF GENES THAT PASS FISHER TEST
#fishergenes[1] = c("Gene", "P-Value")
# FISHER LOOP FOR ALL GENES
for (i in gene_names) {
  #write(i, file = "Fisher_Data.txt", append = TRUE)
  sink("Fisher_Data.txt", append = TRUE)
  print(i)
  lb = anno_table[anno_table$RespType == "LB",] # DEFINE SUBTABLES BASED ON RESPONSE
  nb = anno_table[anno_table$RespType == "NB",]
  lbg = dim(lb[lb[,i] %in% vals,])[1]     # LB WITH GENE MUTATION
  lbno = dim(lb[! lb[,i]  %in% vals,])[1] # LB WITHOUT GENE MUTATION
  nbg = dim(nb[nb[,i] %in% vals,])[1]
  nbno = dim(nb[! nb[,i]  %in% vals,])[1]
  # COMPILE MATRIX FOR FISHER TEST
  mat = matrix(c(lbg,lbno,nbg,nbno), nrow = 2, ncol = 2, dimnames = list(c("Mutation", "None"),c("LB", "NB")))
  print(mat)
  #write.table(mat, file = "Fisher_Data.txt", row.names = TRUE, col.names = TRUE, append = TRUE)
  #  mat
  # RUN FISHER TEST
  f = fisher.test(mat)
  print(f)
  sink()
  #write(d, file = "Fisher_Data.txt", append = TRUE)
  #lapply(f, write, "Fisher_Data.txt", append = TRUE)
  pval = f[1] # DEFINE P-VALUE
  # APPEND OUTPUT DATAFRAME
  add_frame = data.frame(Gene=i, P_Value=pval)
  fishergenes = rbind(fishergenes, add_frame)
  #  print(paste(i, "has p-value:", pval))
  # IF P-VALUE IS SUFFICIENTLY LOW, RUN WILCOXON TEST AND PLOT BOXPLOT
    if ( pval < 0.05 ) { 
  f_pval = round(pval[[1]], digits = 3)
  my_table = anno_table[,c("Mut_Count", i)]
  Mut_Count = log10(my_table$Mut_Count)
  # COUNT SAMPLES WITH AND WITHOUT MUTATION
  counting_dat = subset(my_table, my_table[,i] %in% vals)
  n_yes = dim(counting_dat)[1]
  n_total = dim(my_table)[1]
  n_no = n_total - n_yes
  my_table$Mutated = "no"
  my_table[my_table[,i] %in% vals,]$Mutated = "yes"
  x = my_table[my_table$Mutated == "yes",]$Mut_Count
  y = my_table[my_table$Mutated == "no",]$Mut_Count
  if ( length(x) > 1 & length(y) > 1) {
    w_pval = (wilcox.test(x,y))[[3]]
    w_pval = round(w_pval, digits = 3)
  } else { w_pval = "NULL"}
  png(paste(i, "_boxplot.png"))
  plot(as.factor(my_table$Mutated),Mut_Count, main = paste("Mutation Load vs Mutation in", i, "\nFisher p-value = ", f_pval, "        Wilcoxon p-value = ", w_pval, "\nn no = ", n_no, "     n yes = ", n_yes), ylab = "Mutation Load", xlab = paste("Mutated ", i, "?"))
  dev.off()
  #fishergenes[length(fishergenes) + 1] = i
    }
  write("\n", file = "Fisher_Data.txt", append = TRUE)
}
write.table(fishergenes, file = "Fisher_P_Values.tsv", sep = "\t", col.names = TRUE, row.names =FALSE, quote = FALSE)
