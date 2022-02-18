library(ggplot2)
library(gridExtra)
library(scales)

comp<-c(rep("   Corer vs. Panaroo  ",3), rep("Corer vs. SibeliaZ",3),
rep("Panaroo vs. SibeliaZ",3))
tool <- factor(rep(c("Tool 1", "both", "Tool 2"),3), levels=c("Tool 1","both","Tool 2"))

values1<-c(516,23273,2135,99,24514,894,2007,22606,1183)
d1<-data.frame(comp,tool,values1)

p1<-ggplot(d1, aes(fill=tool, y=values1, x=comp)) +
geom_bar(position="stack", stat="identity") + xlab("") + 
ylab("Number of core genes") + guides(fill="none") + 
theme(axis.text.x = element_text(angle = 90))  + ggtitle("B. animalis") + 
theme(plot.title = element_text(face="italic")) +
scale_y_continuous(labels = comma)


values2<-c(2193,174514,15691,44,184076,6129,13957,170163,6544)
d2<-data.frame(comp,tool,values2)


p2<-ggplot(d2, aes(fill=tool, y=values2, x=comp)) +
geom_bar(position="stack", stat="identity") + xlab("") + ylab("") +
guides(fill="none") + 
theme(axis.text.x = element_text(angle = 90)) + 
ggtitle("Y. pestis")+ theme(plot.title = element_text(face="italic")) +
scale_y_continuous(labels = comma)


values3<-c(6668,280586,35225,57,54869,260942,4999,49927,237327)
d3<-data.frame(comp,tool,values3)




p3<-ggplot(d3, aes(fill=tool, y=values3, x=comp)) +
geom_bar(position="stack", stat="identity") + xlab("") + ylab("") +
guides(fill="none") + 
theme(axis.text.x = element_text(angle = 90)) + 
ggtitle("E. faecium")+ theme(plot.title = element_text(face="italic")) +
scale_y_continuous(labels = comma)


comp4<-c(rep("   Corer vs. Panaroo  ",3))
tool4 <- factor(rep(c("Tool 1", "both", "Tool 2"),1), levels=c("Tool 1","both","Tool 2"))

values4<-c(26894,632905,50507)
d4<-data.frame(comp4,tool4,values4)


p4<-ggplot(d4, aes(fill=tool4, y=values4, x=comp4)) +
geom_bar(position="stack", stat="identity") + xlab("") + ylab("") +
guides(fill=guide_legend(title="")) + 
theme(axis.text.x = element_text(angle = 90)) + 
ggtitle("L. monocytogenes")+ theme(plot.title = element_text(face="italic")) +
scale_y_continuous(labels = comma)

pdf(file="tab2.pdf",width=10, height=4)
grid.arrange(p1, p2, p3, p4, nrow = 1)
dev.off()
