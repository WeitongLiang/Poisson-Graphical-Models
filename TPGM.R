library(XMRF)

data("brcadat")
tpgm.fit <- XMRF(brcadat, method = 'TPGM', stability='bootstrap', lambda=0.1, parallel=FALSE)
