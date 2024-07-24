args = commandArgs(trailingOnly=T)
strD <- args[1]

if(!dir.exists(strD)) stop("Missing qs folder: ",strD)

for(one in list.files(strD,"qs$",full.names=T)){
    message("working on ",basename(one))
    D <- qs::qread(one)
    D$gID <- rownames(D)
    data.table::fwrite(D,gsub("qs$","csv",one))
}


