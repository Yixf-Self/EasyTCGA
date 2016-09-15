library(glmnet)

n = 100 
p = 200

y = rep("blue", n)
y[(n/2):n] = "red"

y = factor(y) # label to "factor"
y_names = seq(n,1,-1)

# data = dn_miRSeq(miRNA_ID[sample.int(2588,100)], "ACC", "", page.Size=500, sort_by="tcga_participant_barcode")
x = matrix(rnorm(n*p), nrow=n, ncol=p)
x_names = seq(1,n,1)

idx = match(y_names, x_names)

y_names = y_names[idx]
y = y[idx]

# glmnet
obj = glmnet(x, y, family="binomial")
plot(obj)

c = coef(obj)
idx = which(abs(c[-1, 10])>0)

# print the names of miRNA corresponding to the indices