args = commandArgs(trailingOnly=T)
if(length(args)<2){
  stop("Two EA objects are required!")
}
obj1rdata <- args[1]
obj2rdata <- args[2]

# check input
stopifnot(grepl("rdata$",obj1rdata,ignore.case=T))
stopifnot(grepl("rdata$",obj2rdata,ignore.case=T))
stopifnot(file.exists(obj1rdata))
stopifnot(file.exists(obj2rdata))

# load QuickOmics objects
ll = load(obj1rdata)
obj1 = list()
for (id in ll) {obj1[[id]] = get(id)}
ll = load(obj1rdata)
obj2 = list()
for (id in ll) {obj2[[id]] = get(id)}

# check the deg results
message("Checking the deg results ...")
rl1 = obj1[['results_long']]
rl2 = obj2[['results_long']]
rl1 = rl1[order(rl1$UniqueID,rl1$test),]
rl2 = rl2[order(rl2$UniqueID,rl2$test),]
cat('comparing results_long between the two models\n')
print(all.equal(rl1,rl2))
cat('***************************\n')

# check the expression
message("Checking the expression ...")
dw1 = obj1[['data_wide']]
dw2 = obj2[['data_wide']]
stopifnot(all.equal(sort(colnames(dw1)),sort(colnames(dw2))))
stopifnot(all.equal(sort(row.names(dw1)),sort(row.names(dw2))))
dw1 = dw1[row.names(dw2),colnames(dw2)]
cat('comparing data_wide between the two models\n')
print(all.equal(as.data.frame(dw1),as.data.frame(dw2)))
cat('***************************\n')

for (id in setdiff(ll,c('results_long','data_wide')))
{
  cat('comparing ',id,' between the two models\n')
  print(all.equal(obj1[[id]],obj2[[id]]))
  cat('***************************\n')
}


