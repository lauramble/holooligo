library(dplyr)
library(networkD3)

FUT1 = read.csv("../TEST2-FUT1/majhapind.csv", row.names = 1)
FUT2 = read.csv("../TEST2-FUT2/majhapind.csv", row.names = 1)
FUT3 = read.csv("../TEST2-FUT3/majhapind.csv", row.names = 1)

pop = lapply(strsplit(colnames(FUT1), split=".", fixed = T), `[`, 1)
pop = unlist(pop)
pop.hap = rep(pop, each = 2)

mat13 = matrix(nrow = 0, ncol = 2)
mat23 = mat13

for (i in 1:nrow(FUT1)) {
  f1 = rownames(FUT1)[i]
  het = colnames(FUT1)[which(FUT1[i,] == 1)]
  homo = colnames(FUT1)[which(FUT1[i,] == 2)]
  comb = c(het, rep(homo, each = 2))
  combhet = rowSums(FUT3[,comb]==1)
  combhomo = rowSums(FUT3[,comb]==2)*2
  comb = combhet + combhomo
  c = c()
  for (i in 1:length(comb)) {
    c = c(c, rep(names(comb)[i], comb[i]))
  }
  temp = cbind(rep(f1, length(c)), c)
  mat13 = rbind(mat13, temp)
}

mat13[,1] = paste0("FUT1_", mat13[,1])
mat13[,2] = paste0("FUT3_", mat13[,2])
df = data.frame(mat13)
colnames(df) = c("FUT1", "FUT3")
df = df %>% count(FUT1, FUT3, name = "n")

links = df

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes = data.frame(
  name=c(as.character(links$FUT1), 
         as.character(links$FUT3)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$FUT1, nodes$name)-1 
links$IDtarget <- match(links$FUT3, nodes$name)-1

# Make the Network
p1 <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "n", NodeID = "name", fontSize = 14,
                   sinksRight=FALSE)
p1

for (i in 1:nrow(FUT2)) {
  f1 = rownames(FUT2)[i]
  het = colnames(FUT2)[which(FUT2[i,] == 1)]
  homo = colnames(FUT2)[which(FUT2[i,] == 2)]
  comb = c(het, rep(homo, each = 2))
  combhet = rowSums(FUT3[,comb]==1)
  combhomo = rowSums(FUT3[,comb]==2)*2
  comb = combhet + combhomo
  c = c()
  for (i in 1:length(comb)) {
    c = c(c, rep(names(comb)[i], comb[i]))
  }
  temp = cbind(rep(f1, length(c)), c)
  mat23 = rbind(mat23, temp)
}

mat23[,1] = paste0("FUT2_", mat23[,1])
mat23[,2] = paste0("FUT3_", mat23[,2])
df = data.frame(mat23)
colnames(df) = c("FUT2", "FUT3")
df = df %>% count(FUT2, FUT3, name = "n")

links = df

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes = data.frame(
  name=c(as.character(links$FUT2), 
         as.character(links$FUT3)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$FUT2, nodes$name)-1 
links$IDtarget <- match(links$FUT3, nodes$name)-1

# Make the Network
p2 <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "n", NodeID = "name", fontSize = 14,
                   sinksRight=FALSE)
p2

