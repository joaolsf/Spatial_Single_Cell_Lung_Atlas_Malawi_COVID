ortho <- read.csv("consolidated_edited.csv")
cluster <- read.csv("DC_cells_final_v2.csv")

idx <- match(ortho$ROI, cluster$ROI)
ortho$ImageNumber <- cluster$ImageNumber[ idx ]
ortho$ObjectNumber <- cluster$ObjectNumber[ idx ]
ortho$ROI <- cluster$ROI[ idx ]
ortho$TMAID <- cluster$TMAID[ idx ]
ortho$CaseID <- cluster$CaseID[ idx ]
ortho$Region <- cluster$Region[ idx ]
ortho$Patient <- cluster$Patient[ idx ]
ortho$Group <- cluster$Group[ idx ]
ortho$Diagnosis <- cluster$Diagnosis[ idx ]

write.csv(ortho, file="Denoise_DC_cells_final.csv")
