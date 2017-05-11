library(BoolNet)
library(colorRamps)

hgf2 <- loadNetwork("EGFR_new.txt")
# knockout cjun
#hgf2 <- fixGenes(hgf2, "cjun", 0)
out2 <- getPathToAttractor(hgf2,c(1,rep(0,65)))

strt <- out2[nrow(out2),]  

# turn on Pai-1
idx <- grep("NPY","SERPINE1",names(strt))
strt[idx] <- 1

out_late <- getPathToAttractor(hgf2,as.numeric(strt))
if(nrow(out_late) == 0) {
	out_late <- strt
}

# Turn off Met
strt_met_off <- out_late[nrow(out_late),]
strt_met_off[1:2] <- 0
out_met_off <- getPathToAttractor(hgf2,as.numeric(strt_met_off))

if(nrow(out_met_off) == 0) {
	out_met_off <- strt_met_off
}
# Turn off Pai-1
strt_met_off_serpin_off <- strt_met_off
idx <- grep("PAI1",names(strt_met_off_serpin_off))
strt_met_off_serpin_off[idx] <- 0
out_met_off_serpin_off <- getPathToAttractor(hgf2,as.numeric(strt_met_off_serpin_off))

# Turn off EGFR
strt_met_off_egfr_off <- strt_met_off
idx <- grep("EGFR",names(strt_met_off_egfr_off))
# knockout egfr
hgfe <- fixGenes(hgf2, "EGFR", 0)
strt_met_off_egfr_off[idx] <- 0
out_met_off_egfr_off <- getPathToAttractor(hgfe,as.numeric(strt_met_off_egfr_off))


#myorder <- c(1:9,26:29,49,34,41,55,42,24:25,19:21,43:44,10:15,59:61,63,62,16,22:23,46,47,17:18,37:38,35:36,50:54,56:58,45,30,64,39,31:33,40,48)

myorder <- c(1:9,26:29,49,34,41,55,42,24:25,19:21,43:44,10:15,59:61,63,62,16,22:23,46,47,17:18,37:38,35:36,50:54,56:58,45,30,64,39,31:33,40,48)

myorder <- c(1:9,26:29,48,34,41,54,42,24:25,19:21,43:44,10:15,58:60,61,62,16,
              22:23,46,47,17:18,37:38,35:36,49:53,55:57,45,30,64,39,31:33,40,63)

# Create the heatmap
out_disp1 <- out2[,myorder]
                   for(i in 1:5){
                #   	   out_disp1 <- rbind(out_disp1,out_disp1[nrow(out_disp1),])
                   }
                   	   
out_disp2 <- out_late[,myorder]
                   for(i in 1:5){
                 #  	   out_disp2 <- rbind(out_disp2,out_disp2[nrow(out_disp2),])
                   }
                   
out_disp3 <- out_met_off[,myorder] 
                   for(i in 1:5){
                  # 	   out_disp3 <- rbind(out_disp3,out_disp3[nrow(out_disp3),])
                   }
                   
out_disp4 <- out_met_off_serpin_off[,myorder]                    
 
out_disp5 <- out_met_off_egfr_off[,myorder]                    

                   out <- cbind(t(out_disp1),t(out_disp2),t(out_disp2),t(out_disp4))

                   mycols <- c(rep(1,2),rep(2,7),rep(3,5),rep(4,6),rep(5,5),rep(6,12),rep(7,3),
                   	 	rep(8,7),rep(9,8),rep(10,8),rep(11,1))
 mycol <- matlab.like2(12)
  mycol[1] <- "#FFFFFF"
  
for(i in 1:ncol(out)){
	out[,i] <- out[,i]*mycols
}
#heatmap.2(out,Rowv=F,Colv=F,trace=c("none"),dendrogram = c("none"),key=F,colsep=c(9,15),rowsep=c(9,15))
#pheatmap(out,cluster_cols=F,cluster_rows=F,cellwidth=10,color=mycol)

pheatmap(t(out_disp1)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol)
pheatmap(t(out_disp2)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol)
# Met_off
out_disp3 <- out_disp2[3:6,]
out_disp3[1:3,] <- out_disp3[4,]
out_disp3[,1:2] <- 0
pheatmap(t(out_disp3)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol)
# serin off
colnames(out_disp4) <- toupper(colnames(out_disp4))

pheatmap(t(out_disp4)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,fontsize=8)

pheatmap(t(out_disp5)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,fontsize=8)
