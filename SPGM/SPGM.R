library(XMRF)

data("brcadat")
tpgm.fit <- XMRF(brcadat, method = 'SPGM', stability='bootstrap', lambda=0.1, parallel=FALSE)

saveRDS(tpgm.fit, file = "/Users/yanyan/Documents/Duke University/Courses/STA 841 Categotrical Data Analysis/Poisson-Graphical-Models/TPGM/TPGM_fit.rds")