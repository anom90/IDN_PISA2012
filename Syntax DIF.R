
# Package
library(haven)
library(dplyr)
library(mirt)
library(difR)


# Data
data_idn <- read_sav("IDN_PISA2012.sav")
summary(data_idn)
table(data_idn$st04q01)
data_idn[data_idn=="Not reached"] <- NA
data_idn <- data_idn[apply(data_idn[,-c(1,2)], 1, function(x) !all(is.na(x))), ]

data_idn[, nama_var] <- ifelse(data_idn[, nama_var] == 'Score 1', 1, 0)
data_idn[is.na(data_idn)] <- 0

#Data dan Q-Matrix
data <- data_idn[,-c(1:2)]
Q <- matrix(c(0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0,
            1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0,
            0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
            0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
            0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0,
            0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1,
            0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,
            0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1,
            0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1),
            nrow = 12, ncol = 11, byrow = T)
Q <- as.data.frame(Q)
colnames(Q) <- c('N1', 'N2', 'N3', 'N4', 'P1', 'P2', 'P3', 'C1', 'C2', 'C3', 'C4')

#Model Fit
model <- '
       F1 = 1-12
       CONSTRAIN = (1-12, a1)
       '

mod_rasch <- mirt(data = data, model = 1,  itemtype = 'Rasch', SE =T)
mod_1pl <- mirt(data = data, model = model, SE = T)
mod_2pl <- mirt(data = data, model = 1,  itemtype = '2PL', SE = T)
mod_3pl <- mirt(data = data, model = 1,  itemtype = '3PL', SE = T, TOL = 0.001)
mod_4pl <- mirt(data = data, model = 1,  itemtype = '4PL', SE = T, TOL = 0.001)

anova(mod_2pl, mod_4pl)
summary(mod_2pl)


#Asumsi IRT
#Asumsi Unidimensi

#Asumsi Invariansi Parameter
#5.3.1 Ivariansi Parameter Butir
#Peserta Nomor Ganjil
th_ganjil <- seq(1, nrow(data), by=2)
mod_2pl_th_ganjil <- mirt(data[th_ganjil,], model=1, itemtype = "2PL")
koef_2pl_ganjil <- coef(mod_2pl_th_ganjil, IRTpars=T, simplify=T)
koef_2pl_ganjil <- as.data.frame(koef_2pl_ganjil)
koef_2pl_ganjil$Grup <- 1
#Peserta Nomor Genap
th_genap <- seq(2, nrow(data),by=2)
mod_2pl_th_genap <- mirt(data[th_genap,], model=1, itemtype = "2PL")
koef_2pl_genap <- coef(mod_2pl_th_genap, IRTpars=T, simplify=T)
koef_2pl_genap <- as.data.frame(koef_2pl_genap)
koef_2pl_genap$Grup <- 2
#Plot Invariansi Paramter Diskriminan
plot(koef_2pl_ganjil[,1], koef_2pl_genap[,1], xlab = "Diskriminan Data Ganjil",
     ylab = "Diskriminan Data Genap")
abline(0,1)
cor(koef_2pl_ganjil[,1], koef_2pl_genap[,1])

#Plot Invariansi Paramter Tingkat Kesulitan
plot(koef_2pl_ganjil[,2], koef_2pl_genap[,2], xlab = "Tingkat Kesulitan Data Ganjil",
     ylab = "Tingkat Kesulitan Data Genap")
abline(0,1)
cor(koef_2pl_ganjil[,2], koef_2pl_genap[,2])

# Model fit
tol <- 0.001
fitGDINA <- GDINA(dat = data_dif[,-1], Q = Q, model = 'GDINA', mono.constraint = T) 
fitDINA <- GDINA(dat = data_dif[,-1], Q = Q, model = 'DINA', mono.constraint = T)
fitDINO <- GDINA(dat = data_dif[,-1], Q = Q, model = 'DINO', mono.constraint = T)
fitACDM <- GDINA(dat = data_dif[,-1], Q = Q, model = 'ACDM', mono.constraint = T)
fitRRUM <- GDINA(dat = data_dif[,-1], Q = Q, model = 'RRUM', mono.constraint = T)
fitLLM <- GDINA(dat = data_dif[,-1], Q = Q, model = 'LLM', mono.constraint = T)


qval1 <- Qval (fitGDINA, method = "PVAF", eps = 0.95)
PVAFvalues1 <- extract (qval1, what = "PVAF")
PVAFvalues1
suggQ1 <- extract (qval1, what = "sug.Q")
suggQ1
Q
a <- GDINA::CA(fitGDINA, what = "EAP")
summary(fitGDINA)
GDINA::extract(fitDINA)

#Reliabilitas
#CTT
reliability(data_dif[,-1])
#IRT
mirt::marginal_rxx(mod_2pl)
#CDM
rel_cdm <- GDINA::CA(fitGDINA)
rel_cdm$tau_k
## Relative Fit indices
anova(fitGDINA, fitDINA, fitDINO, fitACDM, fitRRUM, fitLLM)
summary(fitLLM)
summary(fitDINA)

GDINA::modelfit(fitLLM)

#Local Dependence
LD <- CDM::gdina(dat = data_dif[,-1], q.matrix = Q, group = data_dif[,1], rule = "LLM")
summary(LD)
summary(CDM::modelfit.cor.din (LD))

#Coef
coef(fitLLM, what = "gs", withSE = TRUE)
coef(fitDINA, what = "delta")

extract(fitLLM, what = "discrim")


#DIF
data_dif <- data.frame(SEX = data_idn[,2], data)
data_dif$SEX <- as.factor(data_dif$SEX)


# DIF-CTT
#MH
mh <- difMH(data_dif[,-1], group = data_dif[,1], focal.name = "Female")
mh
plot(mh)

#LR
logreg <- difLogReg(data_dif[,-1], group = data_dif[,1], focal.name = "Female")
logreg
plot(logreg)


#DIF-IRT
#Raju
raju <- difRaju(data_dif[,-1], group = data_dif[,1], focal.name = "Female",  model = "2PL")
raju
plot(raju)

#Lord
lord <- difLord(data_dif[,-1], group = data_dif[,1], focal.name = "Female",  model = "2PL")
lord
plot(lord)

#DIF-CDM
wald <- dif(data_dif[,-1], Q = Q, group = data_dif[,1], method = "wald", 
            mono.constraint = T, SE.type = 3, model = "RRUM")
wald 


# DIF using LR test
LR <- dif(data_dif[,-1], Q = Q, group = data_dif[,1], method = "LR", 
          mono.constraint = T, SE.type = 3, model = "GDINA")
LR

options(scipen = 99)

#Atribut by Gender
fitGDINA_gruop <- GDINA(dat = data_dif[,-1], Q = Q, model = 'GDINA', 
                        mono.constraint = T, group = data_dif[,1]) 


att <- extract(fitGDINA_gruop,what = "prevalence")
att <- cbind(Female = att$Female[,2], Male = att$Male[,2])
att <- round(att, 2)
x <- extract(fitGDINA_gruop,what = "posterior.prob")

y <- round(t(x),2)
y.f <- data.frame(Female = y[,1])

y.0 <- round(y[lc_name,],2)


pdf('5. Plot MH.pdf', paper = 'a4')
plot(mh)
dev.off()

pdf('6. Plot LogReg.pdf', paper = 'a4')
plot(logreg)
dev.off()

pdf('7. Plot Raju.pdf', paper = 'a4')
plot(raju)
dev.off()

pdf('8. Plot Lord.pdf', paper = 'a4')
plot(lord)
dev.off()





