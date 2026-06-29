# Sasha Kramer | HPLC Pigment Clustering — R
# Translation of Kramer_cluster.m | Kramer et al., 2022
# https://github.com/sashajane19/HPLCcluster

library(R.matlab)

# --- Column names ---
pigment_cols <- c(
  "Tchla","Tchlb","Tchlc","ABcaro","ButFuco",
  "HexFuco","Allo","Diadino","Diato","Fuco",
  "Perid","Zea","MVChla","DVchla","Chllide",
  "MVChlb","DVchlb","Chlc12","Chlc3","Lut",
  "Neo","Viola","Phytin","Phide","Pras")

# --- Load data ---
mat_data     <- readMat("Global_HPLC_all.mat")
Global_RHPLC <- as.data.frame(mat_data$Global.RHPLC)
colnames(Global_RHPLC) <- pigment_cols

# --- QC ---
detection_limits <- list(
  Phide=0.003, Perid=0.003, MVChlb=0.003, Phytin=0.003,
  ButFuco=0.002, Fuco=0.002, Neo=0.002,
  Pras=0.002, HexFuco=0.002, DVchla=0.002)

for (pigment in names(detection_limits)) {
  limit <- detection_limits[[pigment]]
  mask  <- Global_RHPLC[[pigment]] != 0 & Global_RHPLC[[pigment]] <= limit
  Global_RHPLC[mask, pigment] <- 0
  cat(sprintf("  %8s: %3d value(s) <= %.3f set to 0\n", pigment, sum(mask), limit))
}

# --- Check % below detection ---
n_samples <- nrow(Global_RHPLC)
below_det <- sapply(pigment_cols, function(col)
  100 * sum(Global_RHPLC[[col]] <= 0.001, na.rm=TRUE) / n_samples)
below_det_df <- data.frame(Pigment=names(below_det),
                           Pct=round(below_det,1), row.names=NULL)
print(below_det_df[order(-below_det_df$Pct), ])

# --- Remove degradation products ---
Rpigcluster1 <- Global_RHPLC[, !colnames(Global_RHPLC) %in%
                               c("Chllide","Phytin","Phide")]

# --- Remove redundant / below-detection pigments ---
redundant_pigments <- c("Tchlb","Tchlc","ABcaro","Diadino","Diato",
                        "MVChla","DVchlb","Lut","Pras")
Rpigcluster2 <- Rpigcluster1[, !colnames(Rpigcluster1) %in% redundant_pigments]
label2 <- colnames(Rpigcluster2)

# --- Cluster pigments (absolute) ---
D2 <- as.dist(1 - cor(Rpigcluster2, use="pairwise.complete.obs"))
Z2 <- hclust(D2, method="ward.D2")
plot(Z2, labels=label2, main="Dendrogram — Absolute Pigment Values",
     xlab="", ylab="Linkage Distance", sub="", hang=-1)

# --- Normalize to Tchla ---
normlabel  <- label2[label2 != "Tchla"]
normchl_df <- Rpigcluster2[, normlabel] / Global_RHPLC[["Tchla"]]

# --- Cluster pigments (normalized) ---
D3 <- as.dist(1 - cor(normchl_df, use="pairwise.complete.obs"))
Z3 <- hclust(D3, method="ward.D2")
plot(Z3, labels=normlabel, main="Dendrogram — Tchla-Normalised Pigment Ratios",
     xlab="", ylab="Linkage Distance", sub="", hang=-1)

# --- Cophenetic correlation coefficient ---
coph_dist <- cophenetic(Z3)
c_val     <- cor(D3, coph_dist)
coph_test <- cor.test(as.vector(as.dist(D3)), as.vector(coph_dist), method="pearson")
cat(sprintf("c = %.4f  |  rho = %.4f  |  p = %.4e\n",
            c_val, coph_test$estimate, coph_test$p.value))

# --- Cluster samples ---
maxclust  <- 4  # update based on dendrogram
D_samples <- as.dist(1 - cor(t(normchl_df), use="pairwise.complete.obs"))
Z_samples <- hclust(D_samples, method="ward.D2")
C         <- cutree(Z_samples, k=maxclust)
print(table(C))
results_df <- Global_RHPLC
results_df$Cluster <- C