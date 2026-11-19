# Native-R truth fixtures only; the installed BReMAP reader must not invoke R.
suppressPackageStartupMessages({library(SeuratObject); library(SingleCellExperiment); library(Matrix)})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
dest <- args[[1]]; dir.create(dest, recursive=TRUE, showWarnings=FALSE)
x <- Matrix(matrix(c(1,0,3, 0,2,2, 4,1,0, 2,2,1), nrow=3), sparse=TRUE)
rownames(x) <- c('CD3D','MS4A1','LYZ'); colnames(x) <- c('cell2','cell1','cell4','cell3')
meta <- data.frame(donor_id=c('donorB','donorA','donorB','donorC'),
                  state=factor(c('T','B','T','B')), numeric_note=c(2,NA,4,3),row.names=colnames(x))
logs <- x; logs@x <- log1p(logs@x*.25)
z <- matrix(seq_len(8),nrow=4,dimnames=list(colnames(x),c('PC_1','PC_2')))
make <- function(assay) {
  obj <- CreateSeuratObject(assay= 'RNA', counts=assay, meta.data=meta)
  obj[['pca']] <- CreateDimReducObject(z,key='PC_',assay='RNA')
  obj
}
v4 <- make(CreateAssayObject(counts=x)); v4[['RNA']]$data <- logs
saveRDS(v4,file.path(dest,'seurat_assay.rds'))
v5 <- make(CreateAssay5Object(counts=x)); v5[['RNA']]$data <- logs
saveRDS(v5,file.path(dest,'seurat_assay5.rds'))
only <- make(CreateAssay5Object(data=logs))
saveRDS(only,file.path(dest,'seurat_data_only.rds'))
split <- v5; split[['RNA']] <- split(split[['RNA']],f=c('batchB','batchA','batchB','batchA'))
saveRDS(split,file.path(dest,'seurat_split.rds'))
subset <- v5; subset[['RNA']]$small <- logs[1:2,1:2]
saveRDS(subset,file.path(dest,'seurat_subset_layer.rds'))
sce <- SingleCellExperiment(assays=list(counts=x,logcounts=logs),colData=meta)
reducedDim(sce,'PCA') <- z
saveRDS(sce,file.path(dest,'sce.rds'))
saveRDS(x,file.path(dest,'sparse_matrix.rds'))
saveRDS(list(not_a_single_cell_object=TRUE),file.path(dest,'unsupported.rds'))
# Human-readable native accessor truth for numerical cross-checks.
write.csv(as.matrix(t(x)),file.path(dest,'counts_truth.csv'),quote=TRUE)
write.csv(as.matrix(t(logs)),file.path(dest,'logcounts_truth.csv'),quote=TRUE)
write.csv(z,file.path(dest,'embedding_truth.csv'),quote=TRUE)
cat('Native fixtures created:',dest,'\n')
saveRDS(as.matrix(x),file.path(dest,'dense_matrix.rds'))
bad <- x; colnames(bad)[2] <- colnames(bad)[1]
saveRDS(bad,file.path(dest,'duplicate_cells.rds'))
unnamed <- as.matrix(x); dimnames(unnamed) <- NULL
saveRDS(unnamed,file.path(dest,'unnamed_matrix.rds'))

na_names <- x; rownames(na_names)[1] <- NA_character_
saveRDS(na_names,file.path(dest,'missing_gene.rds'))
saveRDS(DelayedArray::DelayedArray(x),file.path(dest,'delayed_matrix.rds'))
plain <- v5; plain@reductions <- list()
saveRDS(plain,file.path(dest,'seurat_no_reductions.rds'))
