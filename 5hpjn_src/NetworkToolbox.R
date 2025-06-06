###NetworkToolbox: Methods and Measures for Brain,
###Cognitive, and Psychometric Network Analysis
###Alexander P. Christensen, University of North Carolina at Greensboro
###R Journal

#Load NetworkToolbox
library(NetworkToolbox)

#Load NEO-PI-3 data
data("neoOpen")

####Network Construction####
#Construct TMFG network
tmfg <- TMFG(neoOpen)$A

#Construct LoGo network
logo <- LoGo(neoOpen, partial = TRUE)

#Load qgraph
library(qgraph)

#NEO-PI-3 defined facets
facets <- c(rep("actions", 8), rep("aesthetics", 8), rep("fantasy", 8),
            rep("feelings", 8), rep("ideas", 8), rep("values", 8))

#Visualize TMFG
A <- qgraph(tmfg, groups = facets, palette = "ggplot2")
#Visualize LoGo
B <- qgraph(logo, groups = facets, palette = "ggplot2")

#Visualize TMFG and LoGo side-by-side
layout(t(1:2))
Layout <- averageLayout(A,B)
qgraph(A, layout = Layout, esize = 20, title = "TMFG")
qgraph(B, layout = Layout, esize = 20, title = "LoGo")

#Dependency matrix
depmat <- depend(neoOpen)

#Construct dependency TMFG network
deptmfg <- TMFG(depmat, depend = TRUE)$A

#Visualize dependency TMFG network
layout(c(1,1))
qgraph(deptmfg,layout="spring", groups = facets, palette = "ggplot2",
       label.prop = 1, vsize = 4, title = "Dependency TMFG")

####Local Network Characteristics####
#Randomized shortest paths betweenness centrality
bc <- betweenness(tmfg)
rbc <- rspbc(tmfg, beta = .01)

#Plot for comparison
plot(cbind(log(bc+1), log(rbc+1)), xlab="Standard BC (log)", ylab="RSPBC (log)",
     main="RSPBC on Standard BC", pch=16)
text(6, 2, labels = paste("r = ", round(cor(bc, rbc, method = "kendall"), 2)))

####Meso-scale Network Characteristics####
#Unique facets
uniq <- unique(facets)

#Initialize matrix
corr <- matrix(0, nrow = length(uniq), ncol = 2)

#Name rows
row.names(corr) <- uniq

#Name columns
colnames(corr) <- c("Scale2Inv", "CommCent")

#Compute scale-to-inventory correlations
for(i in 1:length(uniq))
{
    #Identify facet members
    target <- which(facets == uniq[i])
    corr[i, 1] <- cor(rowMeans(neoOpen[, -target]), rowMeans(neoOpen[, target]))
}

#Compute community closeness centrality
corr[, 2] <- comm.close(tmfg, comm = facets)

#Compute correlation
round(cor(corr, method = "kendall"), 2)

#Communities via Exploratory Graph Analysis
#Load EGA
library(EGA)

#Estimate dimensions
ega <- EGA(neoOpen, model = "TMFG")
ega.var <- as.character(ega$dim.variables$items)

#Order by item
ega.ord <- ega$dim.variables[match(colnames(neoOpen),ega.var),2]

#Visualize theoretical factors
A <- qgraph(tmfg, groups = as.factor(facets), palette = "ggplot2")

#Visualize walktrap factors
B <- qgraph(tmfg, groups = as.factor(ega.ord), palette = "ggplot2")

#Compare theoretical and walktrap factors
layout(t(1:2))
Layout <- averageLayout(A,B)
qgraph(A, layout = Layout, esize = 20, title = "Theoretical")
qgraph(B, layout = Layout, esize = 20, title = "EGA")

#Unique facets
uniq <- unique(ega.ord)

#Initialize matrix
corr <- matrix(0, nrow = length(uniq), ncol = 2)

#Name rows
row.names(corr) <- uniq

#Name columns
colnames(corr) <- c("Scale2Inv", "CommCent")

#Compute scale-to-inventory correlations
for(i in 1:length(uniq))
{
    #Identify facet members
    target <- which(ega.ord == uniq[i])
    corr[i, 1] <- cor(rowMeans(neoOpen[, -target]),rowMeans(neoOpen[, target]))
}

#Compute community closeness centrality
corr[, 2] <- comm.close(tmfg, comm = ega.ord, weighted = FALSE)

#Compute correlation
round(cor(corr, method = "kendall"), 2)

####Network Adjusted Means/Sums####
#Latent variable model
#Load lavaan
library(lavaan)

#Build NEO model
neo.model <- 'actions =~ Act1 + Act2 + Act3 + Act4 + Act5 + Act6 + Act7 + Act8
aesthetics =~ Aes1 + Aes2 + Aes3 + Aes4 + Aes5 + Aes6 + Aes7 + Aes8
fantasy =~ Fan1 + Fan2 + Fan3 + Fan4 + Fan5 + Fan6 + Fan7 + Fan8
feelings =~ Fee1 + Fee2 + Fee3 + Fee4 + Fee5 + Fee6 + Fee7 + Fee8
ideas =~ Ide1 + Ide2 + Ide3 + Ide4 + Ide5 + Ide6 + Ide7 + Ide8
values =~ Val1 + Val2 + Val3 + Val4 + Val5 + Val6 + Val7 + Val8
open =~ actions + aesthetics + fantasy + feelings + ideas + values'

#Fit CFA model
fit <- cfa(neo.model, data = neoOpen, estimator = "WLSMV")

#Compute latent variable scores
cfaScores <- lavPredict(fit)
cfaLV <- scale(cfaScores[,7])

#Compute network latent variable scores
netScores <- nams(neoOpen, tmfg, comm = facets)
netLV <- netScores$Standardized$overall

#Plot network latent variable on CFA latent variable
layout(t(c(1,2)))
plot(cfaLV, netLV,
     main = "Network Latent Variable\non CFA Latent Variable",
     ylab = "Network Latent Variable",
     xlab = "CFA Latent Variable")

#Compute correlation
netCor <- cor(cfaLV, netLV)

#Compute root mean square error
netRoot <- rmse(cfaLV, netLV)

#Add text to plot
text(-2, 2, labels = paste("r = ", round(netCor, 2),
                           "\nrmse = ", round(netRoot, 3)))

#Compute standardized participant means
pmeans <- scale(rowMeans(neoOpen))

#Plot participant means on CFA latent variable
plot(cfaLV, pmeans,
     main = "Participant Means on\nCFA Latent Variable",
     ylab = "Participant Means",
     xlab = "CFA Latent Variable")

#Compute correlation
meansCor <- cor(cfaLV, pmeans, method = "spearman")

#Compute root mean square error
meansRoot <- rmse(cfaLV, pmeans)

#Add text to plot
text(-2, 2, labels = paste("r = ", round(meansCor, 2), 
                           "\nrmse = ", round(meansRoot, 3)))


#Examine facet-level latent abilities
#Network latent facet scores
netFacet <- matrix(0, nrow = 6, ncol = 2)

for(i in 1:6)
{
    netFacet[i,1] <- cor(cfaScores[,i], netScores$Standardized[,i])
    netFacet[i,2] <- rmse(scale(cfaScores[,i]), netScores$Standardized[,i])
}

#Identify unique facets
uniq <- unique(facets)

#Participant facet means
meanFacet <- matrix(0, nrow = 6, ncol = 2)

for(i in 1:6)
{
    meanFacet[i, 1] <- cor(cfaScores[, i],rowMeans(neoOpen[, which(facets == uniq[i])]))
    meanFacet[i, 2] <- rmse(scale(cfaScores[, i]),scale(rowMeans(neoOpen[, which(facets == uniq[i])])))
}

#Compare network latent facet scores and participant facet means
comp <- cbind(netFacet, meanFacet)
row.names(comp) <- c("actions", "aesthetics",
                   "fantasy", "feelings",
                   "ideas", "values")
colnames(comp) <- c("netCor", "netRMSE", "meanCor", "meanRMSE")
comp
