library("treeio")
library("ggtree")
library("ggrepel")

args<-commandArgs(T)

tree_file <- args[1]
table_file <- args[2]
outdir <- args[3]


tree <- read.tree(tree_file)

# Add group inofrmation to tree
group_file <- read.delim(table_file,header = T,row.names = 1)
groupInfo <- split(row.names(group_file), group_file$crispr_type)
tree <- groupOTU(tree,groupInfo)
p <- ggtree(tree, layout="daylight", ladderize = FALSE, branch.length = "none",aes(color=group))
p <- p + theme(legend.position="right",legend.text=element_text(size=12))
p <- p + scale_color_manual(values=c('Cyan','red','Aquamarine3','coral','orange','purple','Magenta4','#F38974','#F7B793','#9BF09B','#43F29B','#19A15F','green','blue','steelblue','yellow','goldenrod','#9B6212','#A2E0F2'))
p <- p + scale_fill_manual(values=c('Cyan','red','Aquamarine3','coral','orange','purple','Magenta4','#F38974','#F7B793','#9BF09B','#43F29B','#19A15F','green','blue','steelblue','yellow','goldenrod','#9B6212','#A2E0F2'))
p <- p + geom_tiplab(aes(subset=group=='Novel class 2 effector protein',label=group))

save_img_file <- paste(outdir,'/effector_tree.png',sep='')
png(save_img_file,height = 800,width = 1000)
p
dev.off()