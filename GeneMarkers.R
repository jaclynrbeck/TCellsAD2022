# List of gene markers for various types of T cells. Most markers come from 
# Andreatta, et. al 2021 (https://www.nature.com/articles/s41467-021-23324-4)

markers.core = list(
  "Naive" = c("Sell", "Tcf7", "Lef1", "Il7r", "Ccr7"),
  
  "CM" = c("Sell", "Lef1", "Ccr7", "Tcf7", "Il7r", "Cd44", "Id3", "Cd27", 
           "Eomes", "Stat3", "Bcl6", "Klf2"),
  
  "EM" = c("Gzmb", "Ifng", "Klrg1", "Cx3cr1", "Il7r", "Tbx21", "Prdm1", "Eomes",
           "Cd44", "Id2"),
  
  "RM" = c("Itga1", "Cxcr6", "Itgae", "Gzmb", "Cd69", "Zfp683", "Il7r"),
  
  "Activated" = c("Klrg1", "Klrk1", "Ccl4", "Ifng", "Tnf", "Ccl5", "Cd69"),
  
  "Exhausted" = c("Layn", "Ctla4", "Pdcd1", "Lag3", "Havcr2", "Tigit", "Cd244a",
                  "Cd160"),
  
  "Effector" = c("Klrg1", "Gzmb", "Ccl5", "Prf1", "Cx3cr1", "Il2ra"),
  
  "MPEC" = c("Il7r", "Tcf7", "Eomes", "Id3", "Bcl6", "Stat3", "Foxo1", 
             "Serpina3g"),
  
  "Anergic" = c("Lag3", "Pdcd1", "Rela", "Nfkb1", "Ikzf1", "Egr1", "Egr2"),
  
  "IFNResp" = c("Ifit3", "Irf7", "Mx1", "Bst2", "Irf1", "Irf2", "Stat1", 
                "Stat2"),
  
  "Cytotoxic" = c("Gzma", "Gzmk", "Nkg7", "Prf1"),
  
  "Th1" = c("Cd4", "Ifng", "Tnf", "Ccr5", "Klrd1", "Il18r1"),
  
  "Th2" = c("Cd4", "Il4", "Il5", "Il13", "Il1r1", "Gata3"),
  
  "Th9" = c("Cd4", "Il9", "Gata3", "Irf4", "Stat6"), 
  
  "Th17" = c("Cd4", "Il17a", "Il21", "Il22", "Il25", "Il17f", "Il1r1"),
  
  "Th22" = c("Cd4", "Il22", "Ccr6", "Ccr4", "Ahr", "Ccr10"),
  
  "Tfh" = c("Cd4", "Il21", "Cd40lg", "Icos", "Bcl6", "Slamf1"), 
  
  "Treg" = c("Il2ra", "Ctla4", "Tnfrsf4", "Foxp3", "Tgfb1", "Il10"),
  
  "GD" = c("Tcrg-V1", "Trgv2", "Tcrg-V3", "Tcrg-V4", "Tcrg-V5", "Tcrg-V6", 
           "Tcrg-V7", "Tcrg-C1", "Tcrg-C2", "Tcrg-C3", "Tcrg-C4", "Trdv1",
           "Trdv2-1", "Trdv2-2", "Trdv3", "Trdv4", "Trdv5", "Trdc", "Sox13", 
           "Blk"),
  
  "MAIT" = c("Traj33", "Trav1", "Zbtb16", "Trbv13-1", "Trbv13-2", "Trbv13-3", 
             "Trbv19"),
  
  "iNKT" = c("Traj18", "Zbtb16", "Il17rb", "Trav11", "Trbv13-1", "Trbv13-2", 
             "Trbv13-3", "Trav11d", "Trbv29"),
  
  # Gamma-delta, MAIT, iNKT subclasses
  
  "Type1" = c("Ifng", "Tnf", "Tbx21", "Il2rb", "Cxcr3", "Cd27", "Nkg7", "Ccl5"),
  "Type2" = c("Il4", "Il5", "Il13", "Zbtb16", "Cd27", "Il17rb", "Gata3"),
  "Type17" = c("Il17a", "Il17f", "Il22", "Il17rb", "Rorc", "Ccr6", "Cd44")
)


markers.extended = list(
  "Naive" = c("Sell", "Tcf7", "Lef1", "Il7r", "Ccr7", "S100a10", "Klf6", "Crem",
              "Gpr183", "Ltb", "Satb1", "Klf2", "S1pr1", "Fos", "Fosb", "Junb"),
  
  "CM" = c("Sell", "Lef1", "Ccr7", "Tcf7", "Il7r", "Cd44", "Id3", "Cd27", 
           "Eomes", "Stat3", "Bcl6", "Klf2", "S100a10", "Gpr183", "Ltb", 
           "Satb1", "S1pr1", "Il2ra", "Il15ra", "Il2rb", "Ctnnb1", "Itgal", 
           "Fos", "Fosb", "Junb", "Id2", "Fasl", "Hopx", "Rora", "Atf2", 
           "Bhlhe40", "Klf4", "Runx2", "Cebpb"),
  
  "EM" = c("Gzmb", "Ifng", "Klrg1", "Cx3cr1", "Il7r", "Tbx21", "Prdm1", "Eomes",
           "Cd44", "Id2", "Ccl4", "Gzmk", "Ccl5", "S1pr5", "Prf1", "Tnf", 
           "Bcl6", "Stat4", "Rora", "Fos", "Fosb", "Itgal", "Atf2", "Junb", 
           "Bhlhe40", "Klf4", "Runx2", "Cebpb", "Klf2"),
  
  "RM" = c("Itga1", "Cxcr6", "Itgae", "Gzmb", "Cd69", "Zfp683", "Gzmk", "Ifng", 
           "Ccl4", "Prdm1", "Il7r", "Cxcr3"),
  
  "Activated" = c("Klrg1", "Klrk1", "Ccl4", "Ifng", "Tnf", "Ccl5", "Il2", 
                  "Il4ra", "Nme1", "Gzmb", "Cx3cr1", "Il2ra", "Pdcd1", "Cd69"),
  
  "Exhausted" = c("Layn", "Ctla4", "Pdcd1", "Lag3", "Havcr2", "Tigit", "Cd244a",
                  "Cd160", "Cxcl13", "Irf4", "Stat1", "Hspb1", "Gimap6", 
                  "Hsph1", "Crem", "Cxcr6", "Tnfrsf9", "Batf", "Btla", "Eomes",
                  "Prdm1", "Fas", "Foxp1", "Nkg7", "Ctsw", "Gzmb", "Prf1", 
                  "Ccl4", "Ccl3", "Ifng"),
  
  "Effector" = c("Klrg1", "Gzmb", "Ccl5", "Prf1", "Cx3cr1", "Il2ra", "Ifng", 
                 "Tnf", "Itga4", "Runx1", "S1pr5", "Eomes", "Cd69", "Id2", 
                 "Tnfrsf9", "Prdm1", "Tbx21", "Stat4", "Spn", "Fas", "Runx3", 
                 "Stat5a"),
  
  "MPEC" = c("Il7r", "Tcf7", "Eomes", "Id3", "Bcl6", "Stat3", "Foxo1", 
             "Serpina3g"),
  
  "Anergic" = c("Lag3", "Pdcd1", "Rela", "Nfkb1", "Ikzf1", "Egr1", "Egr2"), 
  
  "IFNResp" = c("Ifit3", "Irf7", "Mx1", "Bst2", "Irf1", "Irf2", "Stat1", 
                "Stat2"),
  
  "Cytotoxic" = c("Gzma", "Gzmk", "Nkg7", "Prf1"),
  
  "Th1" = c("Cd4", "Ifng", "Tnf", "Ccr5", "Klrd1", "Cxcr3", "Havcr2", "Tbx21", 
            "Csf2", "Il2", "Dpp4", "Ifngr1", "Ccr1", "Tnfsf11", "Ltbr", "Lta", 
            "Il18r1"),
  
  "Th2" = c("Cd4", "Il4", "Il5", "Il13", "Il1r1", "Gata3", "Ccr7", "Icos", 
            "Maf", "Cxcr4", "Il10", "Csf2", "Il6", "Ccr8", "Ccr3", "Ccr4", 
            "Ptgdr2", "Havcr1"), 
  
  "Th9" = c("Cd4", "Il9", "Gata3", "Irf4", "Stat6"),
  
  "Th17" = c("Cd4", "Il17a", "Il21", "Il22", "Il25", "Il17f", "Il1r1", "Rora", 
             "Ccr6", "Stat3", "Tgfb1", "Ccr4", "Rorc", "Cd38", "Klrb1c"),
  
  "Th22" = c("Cd4", "Il22", "Ccr6", "Ccr4", "Ahr", "Ccr10"),
  
  "Tfh" = c("Cd4", "Il21", "Cd40lg", "Icos", "Bcl6", "Slamf1", "Pdcd1", 
            "Tnfsf4", "Stat3", "Cd84", "Il6ra", "Cxcr5"),
  
  "Treg" = c("Il2ra", "Ctla4", "Tnfrsf4", "Foxp3", "Tgfb1", "Il10", "Sell", 
             "Lag3", "Itgae", "Entpd1", "Nt5e", "Ccr4", "Izumo1r", "Lrrc32", 
             "Tnfrsf18", "Ikzf2", "Irf4", "Batf")
  
)


markers.immune <- list(
  "Cytokines" = c("Csf2", "Ifng", "Il1a", "Il2", "Il4", "Il7", "Il10", "Il13", 
                  "Il15", "Il16", "Il17a", "Il17f", "Il18", "Il21", "Il22", 
                  "Lif", "Csf1", "Tnfsf11", "Tnf", "Tnfsf4", "Tnfsf10", "Tnfsf8",
                  "Tnfsf13b", "Tnfsf12", "Tnfsf11", "Tnfsf9", "Tnfsf14", "Lta",
                  "Ltb", "Fasl", "Cd70", "Cd40lg"),
  "Chemokines" = c("Cxcl3", "Cxcl2", "Cxcl10", "Ccl3", "Ccl6", "Ccl9", "Ccl1",
                   "Ccl4", "Ccl5", "Ccl25", "Ccl27a", "Cxcl16", "Xcl1"),
  "GrowthFactors" = c("Btc", "Vegfa", "Tgfb1", "Tgfb3", "Tgfa"),
  "Receptors" = c("Il1r2", "Il1r1", "Il1rl2", "Il1rl1", "Il18r1", "Il18rap", 
                  "Il2ra", "Il15ra", "Il6ra", "Il11ra1", "Il12rb2", 
                  "Il23r", "Il17re", "Il17rc", "Il17ra", "Il4ra", "Il21r", 
                  "Il12rb1", "Il27ra", "Il10ra", "Il20rb", "Il20ra", "Il9r",
                  "Il3ra", "Il17rd", "Il17rb", "Il7r", "Il2rb", "Il1rap", 
                  "Il10rb", "Il2rg", "Il18bp",
                  "Cxcr2", "Cxcr4", "Cxcr5", "Cxcr6", "Cxcr3", "Cx3cr1",
                  "Ifngr1", "Ifngr2", "Ifnar1", "Ifnar2",
                  "Tgfbr1", "Tgfbr3", "Tgfbr2", "Tgfbrap1",
                  "Ltbr", "Tnfrsf11a", "Tnfrsf1b", "Tnfrsf8", "Tnfrsf9",
                  "Tnfrsf25", "Tnfrsf14", "Tnfrsf4", "Tnfrsf18", "Tnfrsf1a",
                  "Tnfrsf26", "Tnfrsf22", "Tnfrsf23", "Tnfrsf13b", "Tnfrsf10b",
                  "Tnfrsf10b", "Tnfrsf13c", "Tnfrsf12a", "Tnfrsf21", "Fas",
                  "Pglyrp1", "Cd40"),
  "P2r" = c("P2ry14", "P2ry12", "P2ry1", "P2rx7", "P2rx4", "P2ry10", "P2ry10b"),
  "TxFactors" = c("Tbx21", "Gata3", "Stat4", "Stat6", "Stat3", "Rorc", "Rora", 
                  "Foxp3", "Spi1", "Irf4", "Ahr", "Bcl6", "Rela", "Nfkb1", 
                  "Stat1", "Stat2", "Stat5b", "Stat5a", "Irf6", "Irf5", "Irf3",
                  "Irf7", "Irf2", "Irf8", "Irf1", "Irf9")
)



