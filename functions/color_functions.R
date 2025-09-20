library(viridis)
#martincolorscale=c("#6ef914ff","#060cf1ff","#f106e3ff","#1f1c25ff","#176462ff","#f3fa54ff","#54f0faff","#fa6d54ff","#da32daff","#fbf2f9ff","#fa54a6ff","#54fac4ff","#602646ff","#a96350ff","#d1720cff","#e4eac1ff","#deee82ff","#187695ff","#203655ff","#989865ff","#f2e7f7ff");

martincolourscale=c("green4","royalblue4","magenta4","yellow3","turquoise4","orangered4","yellow4","hotpink4","olivedrab4","#f3fa54ff","#1f1c25ff","deeppink4","#176462ff","#54f0faff","#fa6d54ff","#da32daff","#54fac4ff","gray4","#fbf2f9ff","#602646ff","#a96350ff","orange2","#e4eac1ff","#deee82ff","blue4","#203655ff","#989865ff","#f2e7f7ff","#6ef914ff","#f106e3ff","#060cf1ff","#1f1c25ff");


heatmap_cols=c(rgb(0,1,0),rgb(1,0,1))

add_alpha <- function(alpha=1,cols=heatmap_cols){
	  apply(sapply(cols, col2rgb)/255, 2, function(x)  rgb(x[1], x[2], x[3], alpha=alpha))  
}


schematic_cols_alpha<-function(alpha){
	schematic_cols=c(rgb(0,1,0,alpha),rgb(1,0,1,alpha),rgb(0,0,1,alpha))
	return(schematic_cols)
}

#generates a palette of colours for a heatmap from a vector (called internally)
colPalette<-function(vec,alpha_val=1,cols=heatmap_cols,col_depth=100,Viridis=F){
	if(Viridis){
		cols<-viridis(col_depth)[cut(vec,breaks=col_depth,labels=F)]
	}else{
		colPal<-colorRampPalette(colors=add_alpha(alpha_val,cols=cols),alpha=T)
		cols<-colPal(col_depth)[cut(vec,breaks=col_depth,labels=F)]
	}
	return(cols)
}

