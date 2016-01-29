# SCRIPT TO PRODUCE ONCOPRINTS
# INPUT FORMAT: PATIENT ID IN THE FIRST COLUMN, GENES AS HEADERS IN ALL SUBSEQUENT COLUMNS
#               EACH SAMPLE_GENE PAIR SHOULD BE LABELED AS "HOMDEL", "AMP","BOTH" OR "MUT"
# ADAPTED FROM https://www.bioconductor.org/packages/release/bioc/vignettes/ComplexHeatmap/inst/doc/s8.oncoprint.html
require(methods)
require(ComplexHeatmap)
require(circlize)
require(xlsx)

# DEFINE NEW ONCOPRINT FUNCTION THAT:
# - MAKES TOP-ANNOTATION BASED ON ha, VARIABLE DEFINED OUTSIDE OF THIS FUNCTION
# 
my_onco = function (mat, get_type = function(x) x, alter_fun_list, col, 
                    row_order = oncoprint_row_order(), column_order = oncoprint_column_order(), 
                    show_column_names = FALSE, pct_gp = gpar(), axis_gp = gpar(fontsize = 8), 
                    show_row_barplot = TRUE, row_barplot_width = unit(2, "cm"), 
                    show_column_barplot = TRUE, column_barplot_height = unit(2, 
                                                                             "cm"), remove_empty_columns = FALSE, heatmap_legend_param = list(title = "Alterations"), 
                    ...) 
{
  if (inherits(mat, "matrix")) {
    all_type = unique(unlist(lapply(mat, get_type)))
    all_type = all_type[!is.na(all_type)]
    all_type = all_type[grepl("\\S", all_type)]
    mat_list = lapply(all_type, function(type) {
      m = sapply(mat, function(x) type %in% get_type(x))
      dim(m) = dim(mat)
      dimnames(m) = dimnames(mat)
      m
    })
  }
  else if (inherits(mat, "list")) {
    mat_list = mat
    all_type = names(mat_list)
    mat_list = lapply(mat_list, function(x) {
      oattr = attributes(x)
      x = as.logical(x)
      attributes(x) = oattr
      x
    })
    if (length(unique(sapply(mat_list, nrow))) > 1) {
      stop("All matrix in 'mat_list' should have same number of rows.")
    }
    if (length(unique(sapply(mat_list, ncol))) > 1) {
      stop("All matrix in 'mat_list' should have same number of columns.")
    }
  }
  else {
    stop("Incorrect type of 'mat'")
  }
  if (missing(alter_fun_list) && missing(col)) {
    if (length(mat_list) == 1) {
      alter_fun_list = list(function(x, y, w, h) grid.rect(x, 
                                                           y, w * 0.9, h * 0.9, gp = gpar(fill = "red", 
                                                                                          col = NA)))
      col = "red"
    }
    else if (length(mat_list) == 2) {
      alter_fun_list = list(function(x, y, w, h) grid.rect(x, 
                                                           y, w * 0.9, h * 0.9, gp = gpar(fill = "red", 
                                                                                          col = NA)), function(x, y, w, h) grid.rect(x, 
                                                                                                                                     y, w * 0.9, h * 0.4, gp = gpar(fill = "blue", 
                                                                                                                                                                    col = NA)))
      col = c("red", "blue")
    }
    names(alter_fun_list) = names(mat_list)
    names(col) = names(mat_list)
  }
  arr = array(FALSE, dim = c(dim(mat_list[[1]]), length(all_type)), 
              dimnames = c(dimnames(mat_list[[1]]), list(all_type)))
  for (i in seq_along(all_type)) {
    arr[, , i] = mat_list[[i]]
  }
  oncoprint_row_order = function() {
    order(rowSums(count_matrix), decreasing = TRUE)
  }
  oncoprint_column_order = function() {
    scoreCol = function(x) {
      score = 0
      for (i in 1:length(x)) {
        if (x[i]) {
          score = score + 2^(length(x) - i * 1/x[i])
        }
      }
      return(score)
    }
    scores = apply(count_matrix[row_order, ], 2, scoreCol)
    order(scores, decreasing = TRUE)
  }
  count_matrix = apply(arr, c(1, 2), sum)
  if (is.null(row_order)) 
    row_order = seq_len(nrow(count_matrix))
  row_order = row_order
  if (is.character(column_order)) {
    column_order = structure(seq_len(dim(arr)[2]), names = dimnames(arr)[[2]])[column_order]
  }
  column_order = column_order
  names(column_order) = as.character(column_order)
  if (remove_empty_columns) {
    l = rowSums(apply(arr, c(2, 3), sum)) > 0
    arr = arr[, l, , drop = FALSE]
    column_order = structure(seq_len(sum(l)), names = which(l))[as.character(intersect(column_order, 
                                                                                       which(l)))]
  }
  if (is.null(alter_fun_list$background)) 
    alter_fun_list$background = function(x, y, w, h) grid.rect(x, 
                                                               y, w, h, gp = gpar(fill = "#CCCCCC", col = NA))
  sdf = setdiff(all_type, names(alter_fun_list))
  if (length(sdf) > 0) {
    stop(paste0("You should define shape function for: ", 
                paste(sdf, collapse = ", ")))
  }
  all_type = names(alter_fun_list)
  all_type = setdiff(all_type, "background")
  arr = arr[, , all_type, drop = FALSE]
  sdf = setdiff(all_type, names(col))
  if (length(sdf) > 0) {
    stop(paste0("You should define colors for:", paste(sdf, 
                                                       collapse = ", ")))
  }
  add_oncoprint = function(type, x, y, width, height) {
    alter_fun_list[[type]](x, y, width, height)
  }
  pct = rowSums(apply(arr, 1:2, any))/ncol(mat_list[[1]])
  pct = paste0(round(pct * 100), "%")
  ha_pct = rowAnnotation(pct = row_anno_text(pct, just = "right", 
                                             offset = unit(1, "npc"), gp = pct_gp), width = grobWidth(textGrob("100%", 
                                                                                                               gp = pct_gp)))
  anno_row_bar = function(index, k = NULL, N = NULL) {
    n = length(index)
    count = apply(arr, c(1, 3), sum)[index, , drop = FALSE]
    max_count = max(rowSums(count))
    pushViewport(viewport(xscale = c(0, max_count * 1.1), 
                          yscale = c(0.5, n + 0.5)))
    for (i in seq_len(nrow(count))) {
      if (any(count[i, ] > 0)) {
        x = count[i, ]
        x = x[x > 0]
        x2 = cumsum(x)
        type = all_type[count[i, ] > 0]
        grid.rect(x2, n - i + 1, width = x, height = 0.8, 
                  default.units = "native", just = "right", gp = gpar(col = NA, 
                                                                      fill = col[type]))
      }
    }
    breaks = grid.pretty(c(0, max_count))
    if (k == 1) {
      grid.xaxis(at = breaks, label = breaks, main = FALSE, 
                 gp = axis_gp)
    }
    upViewport()
  }
  ha_row_bar = rowAnnotation(row_bar = anno_row_bar, width = row_barplot_width)
  anno_column_bar = function(index) {
    n = length(index)
    count = apply(arr, c(2, 3), sum)[index, , drop = FALSE]
    max_count = max(rowSums(count))
    pushViewport(viewport(yscale = c(0, max_count * 1.1), 
                          xscale = c(0.5, n + 0.5)))
    for (i in seq_len(nrow(count))) {
      if (any(count[i, ] > 0)) {
        y = count[i, ]
        y = y[y > 0]
        y2 = cumsum(y)
        type = all_type[count[i, ] > 0]
        grid.rect(i, y2, height = y, width = 0.8, default.units = "native", 
                  just = "top", gp = gpar(col = NA, fill = col[type]))
      }
    }
    breaks = grid.pretty(c(0, max_count))
    grid.yaxis(at = breaks, label = breaks, gp = axis_gp)
    upViewport()
  }
  ha_column_bar = HeatmapAnnotation(column_bar = anno_column_bar, 
                                    which = "column", height = column_barplot_height)
  pheudo = c(all_type, rep(NA, nrow(arr) * ncol(arr) - length(all_type)))
  dim(pheudo) = dim(arr[, , 1])
  dimnames(pheudo) = dimnames(arr[, , 1])
  if (show_column_barplot) {
    ht = Heatmap(pheudo, col = col, rect_gp = gpar(type = "none"), 
                 cluster_rows = FALSE, cluster_columns = FALSE, row_order = row_order, 
                 column_order = column_order, cell_fun = function(j, 
                                                                  i, x, y, width, height, fill) {
                   z = arr[i, j, ]
                   add_oncoprint("background", x, y, width, height)
                   for (type in all_type[z]) {
                     add_oncoprint(type, x, y, width, height)
                   }
                 }, show_column_names = show_column_names, top_annotation = ha1,   
                 heatmap_legend_param = heatmap_legend_param, ...)
  }
  else {
    ht = Heatmap(pheudo, rect_gp = gpar(type = "none"), cluster_rows = FALSE, 
                 cluster_columns = FALSE, row_order = row_order, column_order = column_order, 
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   z = arr[i, j, ]
                   add_oncoprint("background", x, y, width, height)
                   for (type in all_type[z]) {
                     add_oncoprint(type, x, y, width, height)
                   }
                 }, show_column_names = show_column_names, ...)
  }
  if (show_row_barplot) {
    ht_list = ha_pct + ht + ha_row_bar
  }
  else {
    ht_list = ha_pct + ht
  }
  return(ht_list)
}
# COMMAND LINE ARGUMENTS: $1 = INPUTFILE, $2 = ANNOTATIONS FILE, $3 OUTPUT FILE, $4 GENE ORDER
args <- commandArgs(trailingOnly = TRUE)
inputfile=args[1]
input2=args[2]
outputfile=args[3]
print_file=args[4] # ORDERED LIST OF GENES, TO LIST THEM ON THE ONCOPRINT

# OPEN PNG GRAPHICS DEVICE, DEFINE SIZE
png(outputfile, height = 1000, width = 2000)

# READ AND PRE-PROCESS FILES
mat = read.table(inputfile, header = TRUE,stringsAsFactors=FALSE, sep = "\t")
#annofile = read.table(input2, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# SORT ANNOFILE BASED ON RESPTYPE THEN ON MUTATION COUNT
annofile = read.table(input2, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# DEFINE PRINT ORDER AND REMOVE GENES NOT SEEN IN RESULTS FILE
vals = as.vector(names(mat))
#vals = sub("\\.", "\\-", vals)
print_order = read.table(print_file, stringsAsFactors = FALSE)
my_row_order = as.vector(print_order$V1)
my_row_order = sub("\\-", "\\.", my_row_order)
my_row_order = my_row_order[my_row_order %in% vals == TRUE]

# CUT OUT LINES FROM ANNOFILE THAT DON'T EXIST IN MAT
samples = mat$Sample.ID
annofile = annofile[annofile$Sample %in% samples, ]

#######    SAMPLE   ORDER   DEFINED   HERE   #######
# DEFINE BINARY VARIABLE BASED ON RESPONSE TYPE
annofile$C = 1
annofile[annofile$RespType == "LB", ]$C= 2
# DEFINE NEW VARIABLE FOR SORTING
# VARIABLE EQUALS: BINARY VARIABLE  +  MUT_COUNT  +  MAXIMUM MUT_COUNT FOR HALF THE LIST
# ORIGINAL #annofile$sort = annofile$C + annofile$Mut_Count + (annofile$C * max(annofile$Mut_Count))
# CHANGE ORDER BY SUBTYPE THEN BY RESPONSE
annofile$sort = annofile$Mut_Count + (annofile$C * max(annofile$Mut_Count)) + (10 * annofile$Subtype * max(annofile$Mut_Count))

annofile = annofile[order(annofile$sort),]

# ENCODE TEXT VERSION OF SUBTYPE 
annofile$Type = "BRAF"
if ( any(annofile$Subtype == 2)) { annofile[annofile$Subtype == 2,]$Type = "RAS" }
if ( any(annofile$Subtype == 3)) {annofile[annofile$Subtype == 3,]$Type = "NF1" }
if ( any(annofile$Subtype == 4)) {annofile[annofile$Subtype == 4,]$Type = "WT" }

out_table = annofile
# DEFINE SAMPLE ORDER TO TO BE USED IN PLOT
sample_order = annofile$Sample
#######    SAMPLE   ORDER   DEFINED   HERE   #######


# SORT mat SO THAT IT HAS THE SAME ORDER AS ANNOFILE
mat = mat[order(match(mat[,1],annofile[,1])),]
#annofile = annofile[order(match(annofile[,1],mat[,1])),]
#annofile = annofile[,-1] # CUT SAMPLE NAMES 
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
#mat=  mat[, -ncol(mat)]
out_table = cbind(out_table, mat)
mat = t(as.matrix(mat))

# DEFINE ANNOTATIONS
Mut_Count = annofile$Mut_Count

# LOG TRANSFORM
Mut_Count = (log10(Mut_Count))

##############################################################
# USES NETHLA
#annofile = annofile[,c("RespType", "netHLA", "Subtype")]
#ha1 = HeatmapAnnotation(df = annofile, barplot = anno_barplot(Mut_Count), col = list(RespType = c("NB" = "black", "LB" = "green"), netHLA = colorRamp2(c(1,6), c("green", "red")), Subtype = colorRamp2(c(1,4), c("blue", "red"))), gap = unit(1, "mm"), annotation_height = c(1,1,1,5))
##########################################################
#ha2 = HeatmapAnnotation(bar = anno_barplot(Mut_Count), annotation_height = 4)

annofile = annofile[,c("RespType", "Type")]
# DEFINE COLOR SCHEME FOR "Type" DEPENDENT ON THE PRESENCE OF EACH TYPE
type_colors = c()
if ( any(annofile$Type == "BRAF") ) {
  code = c("BRAF" = "red3")
  type_colors = append(type_colors, code)
}
if ( any(annofile$Type == "RAS")) {
  code = c("RAS" = "orangered")
  type_colors = append(type_colors, code)
}
if ( any(annofile$Type == "NF1")) {
  code = c("NF1" = "orange")
  type_colors = append(type_colors, code)
}
if ( any(annofile$Type == "WT")) {
  code = c("WT" = "sienna")
  type_colors = append(type_colors, code)
}
#ha1 = HeatmapAnnotation(df = annofile, barplot = anno_barplot(Mut_Count), col = list(RespType = c("NB" = "black", "LB" = "green"), Type = c("BRAF" = "royalblue", "RAS" = "red2", "NF1" = "gold", "WT" = "chocolate")), gap = unit(1, "mm"), annotation_height = c(1,1,5))
ha1 = HeatmapAnnotation(df = annofile, barplot = anno_barplot(Mut_Count), col = list(RespType = c("NB" = "darkblue", "LB" = "cadetblue1"), Type = type_colors), gap = unit(1, "mm"), annotation_height = c(1,1,5))
#ha2 = HeatmapAnnotation(bar = anno_barplot(Mut_Count), annotation_height = 4)

#draw(ha, 1:64)

# DEFINE GRAPHICAL FUNCTIONS TO INCLUDE ONLY MUTATIONS
alter_fun_list = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  #  HOMDEL = function(x, y, w, h) {
  #    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
  #  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w*0.33, h*0.33, gp = gpar(fill = "#008000", col = NA))
  }
  
)
#col = c("MUT" = "#008000", "HOMDEL" = "blue")
col = c("MUT" = "#008000")
my_labels = "Mutation"  

##########################################################
# CHECK FOR "HOT" AND "LOF", APPEND my_labels AS NEEDED
if ( any(grepl("HOT", as.character(mat))) == TRUE ) {
  my_count = length(alter_fun_list) + 1
  HOT = function(x,y,w,h) {
    grid.rect(x,y, w*0.33, h*0.33, gp = gpar(fill = "blue", col = NA))
  }
  alter_fun_list[[my_count]]= HOT
  names(alter_fun_list)[my_count] = "HOT"
  my_count = length(col) + 1
  col[my_count] = c("blue")
  names(col)[my_count] = "HOT"
  my_labels = c(my_labels, "Hotspot")
}
if ( any(grepl("LOF", as.character(mat))) == TRUE ) {
  my_count = length(alter_fun_list) + 1
  LOF = function(x,y,w,h) {
    grid.rect(x,y, w*0.33, h*0.33, gp = gpar(fill = "black", col = NA))
  }
  alter_fun_list[[my_count]]= LOF
  names(alter_fun_list)[my_count] = "LOF"
  my_count = length(col) + 1
  col[my_count] = c("black")
  names(col)[my_count] = "LOF"
  my_labels = c(my_labels, "Loss of Function")
}
if ( any(grepl("DEL", as.character(mat))) == TRUE ) {
  my_count = length(alter_fun_list) + 1
  DEL = function(x,y,w,h) {
    grid.rect(x,y, w*0.8, h*0.8, gp = gpar(fill = "blue", col = NA))
  }
  alter_fun_list[[my_count]]= DEL
  names(alter_fun_list)[my_count] = "DEL"
  my_count = length(col) + 1
  col[my_count] = c("blue")
  names(col)[my_count] = "DEL"
  my_labels = c(my_labels, "Deletion")
}
if ( any(grepl("AMP", as.character(mat))) == TRUE ) {
  my_count = length(alter_fun_list) + 1
  AMP = function(x,y,w,h) {
    grid.rect(x,y, w*0.8, h*0.8, gp = gpar(fill = "red", col = NA))
  }
  alter_fun_list[[my_count]]= AMP
  names(alter_fun_list)[my_count] = "AMP"
  my_count = length(col) + 1
  col[my_count] = c("red")
  names(col)[my_count] = "AMP"
  my_labels = c(my_labels, "Amplification")
}
alter_fun_list = rev(alter_fun_list)
my_labels = rev(my_labels)
##########################################################
my_types = names(col)

x = my_onco(mat, get_type = function(x) strsplit(x, ";")[[1]],
            # USE GRAPHICAL PARAMTERS DEFINED ABOVE
            alter_fun_list = alter_fun_list, col = col,
            # DEFINE COLUMN NAME PRESENCE AND FONT SIZE
            show_column_names = TRUE,
            row_order = my_row_order,
            column_order = sample_order, top_annotation_height = unit(5, "cm"),
            column_names_gp = gpar(fontsize = 11),
            column_title = "\"All Genes\" in NEJM Data",
            heatmap_legend_param = list(title = "Alternations", at = my_types, 
                                        labels = my_labels))
print(x) # PLOT MUST BE PRINTED EXPLICITLY
dev.off()

drop = c("C", "sort")
out = out_table[,!(names(out_table) %in% drop)]
write.xlsx(out, "Onco_Data.xlsx")
write.table(out, file = "Onco_Data.tsv", sep = "\t", col.names = NA)
