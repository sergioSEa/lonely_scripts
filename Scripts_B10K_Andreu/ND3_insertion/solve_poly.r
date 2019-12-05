#Fist: remove internal node names: nw_topology -I TREE
#Second: solve polytomy
TREE="Total_tree2.tree"
library(ape)
t = read.tree(TREE)
tb = multi2di(t)
write.tree(tb,file=TREE)
