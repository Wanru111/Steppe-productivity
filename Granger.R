
pkgs <- c("readxl","tseries","dynlm","lmtest")
for(pkg in pkgs){
  if(!require(pkg,character.only=TRUE)){
    install.packages(pkg,dependencies=TRUE)
    library(pkg,character.only=TRUE)
  }
}


ndvi_df <- readxl::read_excel(
  "F:/1-Dzungar/scocial history/GCA/ndvi.xlsx"
)
conf_df <- readxl::read_excel(
  "F:/1-Dzungar/scocial history/GCA/Conflict.xlsx"
)
# 修正可能的列名拼写
if("NDVI"    %in% names(ndvi_df))  names(ndvi_df)[names(ndvi_df)=="NDVI"] <- "NDVI"
if("Conflict"%in% names(conf_df)) names(conf_df)[names(conf_df)=="Conflict"] <- "Conflict"

df <- merge(ndvi_df, conf_df, by="Year")
df <- df[order(df$Year), ]
df <- df[complete.cases(df[,c("NDVI","Conflict")]), ]
cat("样本长度：", nrow(df), "\n")

adf_safe <- function(x){
  if(length(unique(x))<2) return(list(p.value=0))
  r <- tryCatch(tseries::adf.test(x), error=function(e) list(p.value=0))
  list(p.value = r$p.value[1])
}
r1 <- adf_safe(df$NDVI)
r2 <- adf_safe(df$Conflict)
cat("ADF p -- NDVI:",r1$p.value,"Conflict:",r2$p.value,"\n")
if(r1$p.value>0.05 || r2$p.value>0.05){
  df <- data.frame(
    Year     = df$Year[-1],
    NDVI     = diff(df$NDVI),
    Conflict = diff(df$Conflict)
  )
  df <- df[complete.cases(df), ]
  cat("差分后样本长度：", nrow(df), "\n")
}


library(dynlm)
library(lmtest)

windows <- c(15,20)       
n       <- nrow(df)
cols    <- c("blue","green","purple")


plot(NA,NA,
     xlim=range(df$Year), ylim=c(0,1),
     xlab="Window End Year", ylab="p-value",
     main="Rolling GCA: NDVI → Conflict\n(auto lag in 1..3)")
abline(h=0.05, col="red", lty=2)

for(j in seq_along(windows)){
  W <- windows[j]
  yrs <- numeric(n - W + 1)
  pvs <- numeric(n - W + 1)
  
  for(i in 1:(n - W + 1)){
    sub <- df[i:(i+W-1), ]
    yrs[i] <- sub$Year[W]
    
    
    if(var(sub$NDVI)==0 || var(sub$Conflict)==0){
      pvs[i] <- NA
      next
    }
    
   
    ts_sub <- ts(sub[,c("Conflict","NDVI")],
                 start=sub$Year[1], frequency=1)
    
    
    candidate_p <- 1:3
    aics <- sapply(candidate_p, function(p){
      m <- dynlm(Conflict ~ 
                   L(Conflict, 1:p) +
                   L(NDVI,      1:p),
                 data = ts_sub)
      AIC(m)
    })
    best_p <- candidate_p[ which.min(aics) ]
    
    
    gt <- tryCatch(
      grangertest(Conflict ~ NDVI,
                  order = best_p,
                  data  = sub),
      error = function(e) NULL
    )
    pvs[i] <- if(is.null(gt)) NA else gt$`Pr(>F)`[2]
  }
  
  lines(yrs, pvs, type="b", pch=20, col=cols[j])
}

legend("topright",
       legend=paste0(windows,"-yr"),
       col=cols, pch=20, lty=1, bg="white")


yrs15 <- numeric(n - 15 + 1)
p15   <- numeric(n - 15 + 1)
for (i in seq_len(n - 15 + 1)) {
  sub      <- df[i:(i + 15 - 1), ]
  yrs15[i] <- sub$Year[15]
  if (var(sub$NDVI) == 0 || var(sub$Conflict) == 0) {
    p15[i] <- NA
  } else {
    ts_sub <- ts(sub[, c("Conflict","NDVI")],
                 start = sub$Year[1], frequency = 1)
    
    aics   <- sapply(1:3, function(p) {
      AIC(dynlm(Conflict ~ L(Conflict,1:p) + L(NDVI,1:p),
                data = ts_sub))
    })
    best_p <- which.min(aics)
    gt     <- tryCatch(
      grangertest(Conflict ~ NDVI, order = best_p, data = sub),
      error = function(e) NULL
    )
    p15[i] <- if (is.null(gt)) NA else gt$`Pr(>F)`[2]
  }
}
df15 <- data.frame(year = yrs15, p15yr = p15)


yrs20 <- numeric(n - 20 + 1)
p20   <- numeric(n - 20 + 1)
for (i in seq_len(n - 20 + 1)) {
  sub      <- df[i:(i + 20 - 1), ]
  yrs20[i] <- sub$Year[20]
  if (var(sub$NDVI) == 0 || var(sub$Conflict) == 0) {
    p20[i] <- NA
  } else {
    ts_sub <- ts(sub[, c("Conflict","NDVI")],
                 start = sub$Year[1], frequency = 1)
    aics   <- sapply(1:3, function(p) {
      AIC(dynlm(Conflict ~ L(Conflict,1:p) + L(NDVI,1:p),
                data = ts_sub))
    })
    best_p <- which.min(aics)
    gt     <- tryCatch(
      grangertest(Conflict ~ NDVI, order = best_p, data = sub),
      error = function(e) NULL
    )
    p20[i] <- if (is.null(gt)) NA else gt$`Pr(>F)`[2]
  }
}
df20 <- data.frame(year = yrs20, p20yr = p20)


res <- merge(df15, df20, by = "year", all = TRUE)
print(res)

write.csv(
  res,
  file         = "F:/1-Dzungar/scocial history/GCA/rolling_pvalues_NDVI-Conflict.csv",
  row.names    = FALSE,
  fileEncoding = "UTF-8"
)

