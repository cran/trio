### R code from vignette source 'trio.Rnw'

###################################################
### code chunk number 1: trio.Rnw:143-147
###################################################
library(trio)
data(trio.data)
str(trio.ped1)
trio.ped1[1:10,1:12]


###################################################
### code chunk number 2: trio.Rnw:190-195
###################################################
trio.ped1[,2] <- paste(trio.ped1[,1], trio.ped1[,2], sep="_")
ids <- trio.ped1[,3] != 0
trio.ped1[ids,3] <- paste(trio.ped1[ids,1], trio.ped1[ids,3], sep="_")
trio.ped1[ids,4] <- paste(trio.ped1[ids,1], trio.ped1[ids,4], sep="_")
trio.ped1[1:5, 1:4]


###################################################
### code chunk number 3: trio.Rnw:200-202
###################################################
geno <- ped2geno(trio.ped1)
geno[1:5,]


###################################################
### code chunk number 4: trio.Rnw:209-211
###################################################
trio.gen1[1:5, 3:12]
table(trio.gen1[,3:12] == geno)


###################################################
### code chunk number 5: trio.Rnw:222-223
###################################################
tdt(mat.test[,1])


###################################################
### code chunk number 6: trio.Rnw:233-234
###################################################
tdt(mat.test[,1], model="dominant")


###################################################
### code chunk number 7: trio.Rnw:238-239
###################################################
tdt(mat.test[,1], model="recessive")


###################################################
### code chunk number 8: trio.Rnw:247-248
###################################################
tdtGxG(mat.test[,1], mat.test[,2])


###################################################
### code chunk number 9: trio.Rnw:270-271
###################################################
tdtGxG(mat.test[,1], mat.test[,2], test="screen")


###################################################
### code chunk number 10: trio.Rnw:279-281
###################################################
tdt.out <- colTDT(mat.test)
tdt.out


###################################################
### code chunk number 11: trio.Rnw:287-288
###################################################
print(tdt.out, 3)


###################################################
### code chunk number 12: trio.Rnw:294-295
###################################################
print(tdt.out, 10)


###################################################
### code chunk number 13: trio.Rnw:303-305
###################################################
max.stat <- colTDTmaxStat(mat.test)
max.stat


###################################################
### code chunk number 14: trio.Rnw:313-314
###################################################
max.out <- colTDTmaxTest(mat.test, perm=1000)


###################################################
### code chunk number 15: trio.Rnw:319-320
###################################################
max.out


###################################################
### code chunk number 16: trio.Rnw:340-342
###################################################
genes <- paste("G", rep(1:2, c(2,4)), sep="")
genes


###################################################
### code chunk number 17: trio.Rnw:348-350
###################################################
tdt2.out <- colGxG(mat.test, genes=genes)
tdt2.out


###################################################
### code chunk number 18: trio.Rnw:354-355
###################################################
print(tdt2.out, 8)


###################################################
### code chunk number 19: trio.Rnw:366-367
###################################################
sex <- rep(0:1, e=50)


###################################################
### code chunk number 20: trio.Rnw:372-374
###################################################
gxe.out <- colGxE(mat.test, sex)
gxe.out


###################################################
### code chunk number 21: trio.Rnw:394-395
###################################################
print(gxe.out, onlyGxE=TRUE)


###################################################
### code chunk number 22: trio.Rnw:403-405
###################################################
dat.top3 <- getGxEstats(gxe.out, top=3, sortBy="lrt2df")
dat.top3


###################################################
### code chunk number 23: trio.Rnw:415-417
###################################################
a.out <- allelicTDT(mat.test)
a.out


###################################################
### code chunk number 24: trio.Rnw:433-435
###################################################
a.out2 <- allelicTDT(mat.test, correct=TRUE)
a.out2


###################################################
### code chunk number 25: trio.Rnw:446-448
###################################################
s.out <- scoreTDT(mat.test)
s.out


###################################################
### code chunk number 26: trio.Rnw:500-503
###################################################
trio.tmp <- trio.check(dat=trio.ped1)
str(trio.tmp, max=1)
trio.tmp$trio[1:6,]


###################################################
### code chunk number 27: trio.Rnw:540-545
###################################################
set.seed(123456)
trio.bin <- trio.prepare(trio.dat=trio.tmp, blocks=c(1,4,2,3))
str(trio.bin, max=1)

trio.bin$bin[1:8,]


###################################################
### code chunk number 28: trio.Rnw:557-564
###################################################
str(trio.gen1)
trio.gen1[1:10,1:12]
trio.tmp <- trio.check(dat=trio.gen1, is.linkage=FALSE)
set.seed(123456)
trio.bin2 <- trio.prepare(trio.dat=trio.tmp, blocks=c(1,4,2,3))

trio.bin2$bin[1:8,]


###################################################
### code chunk number 29: trio.Rnw:574-578
###################################################
str(trio.ped2)

trio.tmp <- trio.check(dat=trio.ped2)
trio.tmp$trio[1:6,]


###################################################
### code chunk number 30: trio.Rnw:588-591
###################################################
set.seed(123456)
trio.bin3 <- trio.prepare(trio.dat=trio.tmp, blocks=c(1,4,2,3))
trio.bin3$bin[1:8,]


###################################################
### code chunk number 31: trio.Rnw:597-602
###################################################
str(trio.gen2)
trio.tmp <- trio.check(dat=trio.gen2, is.linkage=FALSE)
set.seed(123456)
trio.bin4 <- trio.prepare(trio.dat=trio.tmp, blocks=c(1,4,2,3))
trio.bin4$bin[1:8,]


###################################################
### code chunk number 32: trio.Rnw:612-619
###################################################
trio.tmp <- trio.check(dat=trio.gen2, is.linkage=FALSE)
set.seed(123456)
trio.imp <- trio.prepare(trio.dat=trio.tmp, blocks=c(1,4,2,3), logic=FALSE)
str(trio.imp, max=1)
trio.imp$miss[c(1:6),]
trio.gen2[1:6,]
trio.imp$trio[1:6,]


###################################################
### code chunk number 33: trio.Rnw:624-628
###################################################
trio.tmp <- trio.check(dat=trio.ped2)
set.seed(123456)
trio.imp2 <- trio.prepare(trio.dat=trio.tmp, blocks=c(1,4,2,3), logic=FALSE)
trio.imp$trio[1:6,]


###################################################
### code chunk number 34: trio.Rnw:639-642
###################################################
trio.tmp <- trio.check(dat=trio.ped.err)
str(trio.tmp, max=1)
trio.tmp$errors


###################################################
### code chunk number 35: trio.Rnw:647-649
###################################################
trio.tmp$trio.err[1:3, c(1,2, 11:12)]
trio.ped.err[1:3,c(1:2, 23:26)]


###################################################
### code chunk number 36: trio.Rnw:662-665
###################################################
trio.rep <- trio.check(dat=trio.ped.err, replace=TRUE)
str(trio.rep, max=1)
trio.rep$trio[1:3,11:12]


###################################################
### code chunk number 37: trio.Rnw:671-676
###################################################
trio.tmp <- trio.check(dat=trio.gen.err, is.linkage=FALSE)
trio.tmp$errors
trio.tmp$trio.err[1:6, c(1,2,7), drop=F]
trio.rep <- trio.check(dat=trio.gen.err, is.linkage=FALSE, replace=TRUE)
trio.rep$trio[1:6,c(1,2,7)]


###################################################
### code chunk number 38: trio.Rnw:694-696
###################################################
str(freq.hap)
freq.hap[1:6,]


###################################################
### code chunk number 39: trio.Rnw:700-706
###################################################
trio.tmp <- trio.check(dat=trio.gen2, is.linkage=FALSE)
set.seed(123456)
trio.imp3 <- trio.prepare(trio.dat=trio.tmp, freq=freq.hap, logic=FALSE)
str(trio.imp3, max=1)
trio.gen2[1:6,]
trio.imp3$trio[1:6,]


###################################################
### code chunk number 40: trio.Rnw:726-727
###################################################
my.control <- lrControl(start=1, end=-3, iter=1000, output=-4)


###################################################
### code chunk number 41: trio.Rnw:751-753
###################################################
lr.out <- trioLR(trio.bin, control=my.control, rand=9876543)
lr.out


###################################################
### code chunk number 42: trio.Rnw:760-761
###################################################
print(lr.out, asDNF=TRUE)


###################################################
### code chunk number 43: trio.Rnw:776-778
###################################################
n.trios <- 100 
y <- rep(c(3, 0, 0, 0), n.trios)


###################################################
### code chunk number 44: trio.Rnw:785-786
###################################################
x <- trio.bin$bin[,-1]


###################################################
### code chunk number 45: trio.Rnw:791-793
###################################################
lr.out2 <- trioLR(x, y, control=my.control, rand=9876543)
lr.out2


###################################################
### code chunk number 46: trio.Rnw:833-835
###################################################
lr.out3 <- trioLR(trio.bin, nleaves=c(3,5), control=my.control, rand=9876543)
lr.out3


###################################################
### code chunk number 47: trio.Rnw:855-857
###################################################
lr.out4 <- trioLR(trio.bin, search="greedy", rand=9876543)
lr.out4


###################################################
### code chunk number 48: trio.Rnw:870-872
###################################################
lr.out5 <- trioLR(trio.bin, search="mcmc", control=my.control, rand=9876543)
lr.out5


###################################################
### code chunk number 49: trio.Rnw:907-909
###################################################
fs.out <- trioFS(trio.bin, B=5, control=my.control, rand=9876543)
fs.out


###################################################
### code chunk number 50: trio.Rnw:946-947
###################################################
ld.out <- getLD(LDdata, asMatrix=TRUE)


###################################################
### code chunk number 51: trio.Rnw:954-955
###################################################
round(ld.out$Dprime[1:10,1:10], 2)


###################################################
### code chunk number 52: trio.Rnw:958-959
###################################################
round(ld.out$rSquare[1:10,1:10], 2)


###################################################
### code chunk number 53: trio.Rnw:984-986
###################################################
blocks <- findLDblocks(LDdata)
blocks


###################################################
### code chunk number 54: trio.Rnw:992-995
###################################################
ld.out2 <- getLD(LDdata, addVarN=TRUE)
blocks2 <- findLDblocks(ld.out2)
blocks2


###################################################
### code chunk number 55: trio.Rnw:1031-1033
###################################################
hap <- as.vector(table(blocks$blocks))
hap


###################################################
### code chunk number 56: trio.Rnw:1071-1077
###################################################
str(simuBkMap)
simuBkMap[1:7,]
sim <- trio.sim(freq=simuBkMap, interaction="1R and 5R", prev=.001, OR=2, 
n=20, rep=1)
str(sim)
sim[[1]][1:6, 1:12]


###################################################
### code chunk number 57: trio.Rnw:1089-1099
###################################################
trio.tmp <- trio.check(dat=trio.gen1, is.linkage=FALSE)
trio.impu <- trio.prepare(trio.dat=trio.tmp, blocks=c(1,4,2,3), logic=TRUE)

str(trio.impu, max=2)
trio.impu$freq[1:7,]
sim <- trio.sim(freq=trio.impu$freq, interaction="1R and 5R", prev=.001, OR=2, 
n=20, rep=1)

str(sim)
sim[[1]][1:6, ]


###################################################
### code chunk number 58: trio.Rnw:1108-1112
###################################################
sim <- trio.sim(freq=freq.hap, interaction="1R or 4D", prev=.001, OR=2, 
n=10, rep=2)
str(sim)
sim[[1]][1:6,]


###################################################
### code chunk number 59: trio.Rnw:1140-1144
###################################################
sim <- trio.sim(freq=freq.hap, interaction="1R or (6R and 10D)", prev=.001, 
OR=2, n=10, rep=1)
str(sim)
sim[[1]][1:6,]


###################################################
### code chunk number 60: trio.Rnw:1147-1151
###################################################
sim <- trio.sim(freq=freq.hap, interaction="1R or (6R and 10D)", prev=.001, 
OR=3, n=10, rep=1, step.save="step3way")
str(sim, max=1)
sim[[1]][1:6,]


