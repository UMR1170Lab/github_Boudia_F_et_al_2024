

pdf("/path/to/dotplot_pour_review_v3.pdf",width = 25,height = 20)


DotPlot(selected_integration, features = c("KIT","HDC","AZU1","MPO","LYZ","RETN","PLAC8","CD163","KLF1","GYPA","GYPB","TOP2A","MKI67","GATA1","MEIS1","ITGA2B","TUBB1","PF4"),dot.scale = 15,cols = "RdBu") +theme(axis.text = element_text(size =20,face = "bold"),legend.key.height = unit(1.8, "cm"),legend.key.width = unit(1, "cm"))+
 scale_y_discrete(limits =c("13" = "13","1" = "1", "4" = "4", "6" = "6" , "0" = "0", "5" = "5", "3" = "3", "8" = "8","12" = "12","2" = "2","10" = "10","7" = "7","14" = "14","11" = "11","15" = "15","9" = "9","17" = "17"))

dev.off()
