## R script for creating graphs and running moderation analyses
## QSM Cortical Columns Project, 2023
## By Jenna Merenstein


## 1. LOAD DATA & LIBRARIES ##
library(readxl)
library(ggplot2)
library(scales)
library(ggprism)
library(RColorBrewer)
library(ggbreak)
library(readr)


setwd("/Users/jenna/Library/Mobile Documents/com~apple~CloudDocs/Jenna/duke")
graphs_curv <- read_excel("qsm_columns/data/graphs_curv.xlsx")
graphs_depth <- read_excel("qsm_columns/data/graphs_depth.xlsx")
depth <- read_excel("qsm_columns/data/global_depth.xlsx")
data <- read_csv("qsm_columns/data/subject_data.csv")


## 2. QSM by depth graphs (Figure 2) ##
## AD = black, controls = blue
ggplot(graphs_depth, aes(x = depth, y = lh_pos,color=dx)) + geom_line(size=1) + geom_errorbar(aes(ymin = lh_pos - lh_pos_se, ymax = lh_pos + lh_pos_se)) + labs(y = "LH Positive Susceptibility", x = NULL) + theme_prism(base_size = 25,base_fontface = "plain") + theme(legend.position = "none",text=element_text(family="Arial"),axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(),axis.line = element_line(size = 1)) + scale_color_manual(values = c("#000000","#004C99")) + coord_cartesian(ylim = c(0.008, 0.016))
ggsave("qsm_columns/brain_A0/figures/lh_pos_depth.tiff", units="in", width = 7.8, height = 6.2, device='tiff', dpi=300)


ggplot(graphs_depth, aes(x = depth, y = rh_pos,color=dx)) + geom_line(size=1) + geom_errorbar(aes(ymin = rh_pos - rh_pos_se, ymax = rh_pos + rh_pos_se)) + labs(y = "RH Positive Susceptibility", x = NULL) + theme_prism(base_size = 25,base_fontface = "plain") + theme(legend.position = "none",text=element_text(family="Arial"),axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(),axis.line = element_line(size = 1)) + scale_color_manual(values = c("#000000","#004C99")) + coord_cartesian(ylim = c(0.008, 0.016))
ggsave("qsm_columns/brain_A0/figures/rh_pos_depth.tiff", units="in", width = 7.8, height = 6.2, device='tiff', dpi=300)

ggplot(graphs_depth, aes(x = depth, y = lh_neg,color=dx)) + geom_line(size=1) + geom_errorbar(aes(ymin = lh_neg - lh_neg_se, ymax = lh_neg + lh_neg_se)) + labs(y = "LH Negative Susceptibility", x = NULL) + theme_prism(base_size = 25,base_fontface = "plain") + theme(legend.position = "none",text=element_text(family="Arial"),axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(),axis.line = element_line(size = 1)) + scale_color_manual(values = c("#000000","#004C99")) + coord_cartesian(ylim = c(-0.016,-0.008))
ggsave("qsm_columns/brain_A0/figures/lh_neg_depth.tiff", units="in", width = 7.8, height = 6.2, device='tiff', dpi=300)

ggplot(graphs_depth, aes(x = depth, y = rh_neg,color=dx)) + geom_line(size=1) + geom_errorbar(aes(ymin = rh_neg - rh_neg_se, ymax = rh_neg + rh_neg_se)) + labs(y = "RH Negative Susceptibility", x = NULL) + theme_prism(base_size = 25,base_fontface = "plain") + theme(legend.position = "none",text=element_text(family="Arial"),axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(),axis.line = element_line(size = 1)) + scale_color_manual(values = c("#000000","#004C99")) + coord_cartesian(ylim = c(-0.016,-0.008))
ggsave("qsm_columns/brain_A0/figures/rh_neg_depth.tiff", units="in", width = 7.8, height = 6.2, device='tiff', dpi=300)


## 3. DEPTH x GROUP MODERATIONS ##
# Must first run the PROCESS macro code
# Controls coded as 2, AD as 1
process(data = depth, y = "lh_pos", x = "depth", w = "group",model = 1,boot=5000,decimal=10.7)
process(data = depth, y = "rh_pos", x = "depth", w = "group",model = 1,boot=5000,decimal=10.7)
process(data = depth, y = "lh_neg", x = "depth", w = "group",model = 1,boot=5000,decimal=10.7)
process(data = depth, y = "rh_neg", x = "depth", w = "group",model = 1,boot=5000,decimal=10.7)


## 4. QSM by curvature graphs (Figure 4) ##
## AD = black, controls = blue
ggplot(graphs_curv, aes(x = curvature, y = lh_pos,fill=dx)) + geom_bar(stat="identity",position="dodge") + geom_errorbar(aes(ymin = lh_pos - lh_pos_se, ymax = lh_pos + lh_pos_se,width=0.5),position=position_dodge(0.9)) + labs(y = "LH Positive Susceptibility", x = NULL) + theme_prism(base_size = 12,base_fontface = "plain") + theme(legend.position = "none",text=element_text(family="Arial"),axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(),axis.line = element_line(size = 0.5)) + scale_fill_manual(values = c("#000000","#004C99")) + coord_cartesian(ylim = c(0, 0.02)) 
ggsave("qsm_columns/brain_A0/figures/lh_pos_curv.tiff", units="in", width = 3.9, height = 3.1, device='tiff', dpi=300)

ggplot(graphs_curv, aes(x = curvature, y = rh_pos,fill=dx)) + geom_bar(stat="identity",position="dodge") + geom_errorbar(aes(ymin = rh_pos - rh_pos_se, ymax = rh_pos + rh_pos_se, width=0.5),position=position_dodge(0.9)) + labs(y = "RH Positive Susceptibility", x = NULL) + theme_prism(base_size = 12,base_fontface = "plain") + theme(legend.position = "none",text=element_text(family="Arial"),axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(),axis.line = element_line(size = 0.5)) + scale_fill_manual(values = c("#000000","#004C99")) + coord_cartesian(ylim = c(0, 0.02))
ggsave("qsm_columns/brain_A0/figures/rh_pos_curv.tiff", units="in", width = 3.9, height = 3.1, device='tiff', dpi=300)

ggplot(graphs_curv, aes(x = curvature, y = lh_neg,fill=dx)) + geom_bar(stat="identity",position="dodge") + geom_errorbar(aes(ymin = lh_neg - lh_neg_se, ymax = lh_neg + lh_neg_se,width=0.5),position=position_dodge(0.9)) + labs(y = "LH Negative Susceptibility", x = NULL) + theme_prism(base_size = 12,base_fontface = "plain") + theme(legend.position = "none",text=element_text(family="Arial"),axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(),axis.line = element_line(size = 0.5)) + scale_fill_manual(values = c("#000000","#004C99")) + coord_cartesian(ylim = c(-.02, 0))
ggsave("qsm_columns/brain_A0/figures/lh_neg_curv.tiff", units="in", width = 3.9, height = 3.1, device='tiff', dpi=300)

ggplot(graphs_curv, aes(x = curvature, y = rh_neg,fill=dx)) + geom_bar(stat="identity",position="dodge") + geom_errorbar(aes(ymin = rh_neg - rh_neg_se, ymax = rh_neg + rh_neg_se,width=0.5),position=position_dodge(0.9)) + labs(y = "RH Negative Susceptibility", x = NULL) + theme_prism(base_size = 12,base_fontface = "plain") + theme(legend.position = "none",text=element_text(family="Arial"),axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(),axis.line = element_line(size = 0.5)) + scale_fill_manual(values = c("#000000","#004C99")) + coord_cartesian(ylim = c(-.02, 0))
ggsave("qsm_columns/brain_A0/figures/rh_neg_curv.tiff", units="in", width = 3.9, height = 3.1, device='tiff', dpi=300)
