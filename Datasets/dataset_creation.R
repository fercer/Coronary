set.seed(3333)

### 10
test.set <- sample(1:60, 10)
test.set.dir <- paste("../dataset/",test.set, ".png", sep="")
test.set.dir <- c(length(test.set.dir), test.set.dir)
test.gt.dir <- paste("../dataset/",test.set, "_gt.png", sep="")
write.table(test.set.dir, "data_base_10.dir", quote=F, col.names=F, row.names=F)
write.table(test.set.dir, "data_gt_10.dir", quote=F, col.names=F, row.names=F)

set.seed(7777)
### 20
test.set <- sample(1:60, 20)
test.set.dir <- paste("../dataset/",test.set, ".png", sep="")
test.set.dir <- c(length(test.set.dir), test.set.dir)
test.gt.dir <- paste("../dataset/",test.set, "_gt.png", sep="")
write.table(test.set.dir, "data_base_20.dir", quote=F, col.names=F, row.names=F)
write.table(test.set.dir, "data_gt_20.dir", quote=F, col.names=F, row.names=F)


set.seed(1111)
### 40
test.set <- sample(1:60, 40)
test.set.dir <- paste("../dataset/",test.set, ".png", sep="")
test.set.dir <- c(length(test.set.dir), test.set.dir)
test.gt.dir <- paste("../dataset/",test.set, "_gt.png", sep="")
write.table(test.set.dir, "data_base_40.dir", quote=F, col.names=F, row.names=F)
write.table(test.set.dir, "data_gt_40.dir", quote=F, col.names=F, row.names=F)


### 60
test.set <- 1:60
test.set.dir <- paste("../dataset/",test.set, ".png", sep="")
test.set.dir <- c(length(test.set.dir), test.set.dir)
test.gt.dir <- paste("../dataset/",test.set, "_gt.png", sep="")
write.table(test.set.dir, "data_base_60.dir", quote=F, col.names=F, row.names=F)
write.table(test.set.dir, "data_gt_60.dir", quote=F, col.names=F, row.names=F)
