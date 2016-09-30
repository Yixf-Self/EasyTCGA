library(glmnet)
n = 100 
p = 200

y = rep("blue", n)
y[(n/2):n] = "red"

y = factor(y) # label to "factor"
y_names = seq(n,1,-1)

x = matrix(rnorm(n*p), nrow=n, ncol=p)
x_names = seq(1,n,1)plot

idx = match(y_names, x_names)

y_names = y_names[idx]
y = y[idx]

# glmnet
obj = glmnet(x, y, family="binomial")
plot(obj)

c = coef(obj)
idx = which(abs(c[-1, 10])>0)

# print the names of miRNA corresponding to the indices

nc = ncol(acc.miRSeq_reshaped)
na.ind = sapply(1:nc, function(i) { any(is.na(acc.miRSeq_reshaped[,i])) })

X = acc.miRSeq_reshaped[,!na.ind]
barcode = rownames(X)
y = acc.label$vital_status
names(y) = acc.label$tcga_participant_barcode

idx = match(barcode, names(y))
y = y[idx,drop=FALSE]


y = factor(y)

out = cv.glmnet(X, y, family="binomial", intercept=FALSE)
c = as.numeric(coef(out)[-1])
names(c) = colnames(X)
idx = order(abs(c), decreasing=TRUE)
c = c[idx,drop=FALSE]
idx = abs(c) > 0
data.frame(coef=c[idx])



https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4742203/
http://link.springer.com/article/10.1007/s13277-015-4630-5


