options(stringsAsFactors=F)

Args <- commandArgs(T)
input <- Args[1]
output <- Args[2]

dat <- read.table(input,header=F)

chr_vec <- unique(dat[,1])
chr_len_vec <- chr_vec
## format pos ##
for (i in 1:length(chr_vec)) {
	tmp_chr_dat <- dat[dat[,1] == chr_vec[i],]
	chr_len_vec[i] <- as.numeric(tail(tmp_chr_dat,1)[2])
}
chr_len_vec <- c(0,chr_len_vec)

for (i in 1:nrow(dat)) {
	tmp_chr <- dat[i,1]
	tmp_pos <- dat[i,2]
	format_pos <- tmp_pos + sum(chr_len_vec[1:tmp_chr])
	dat[i,2] <- format_pos
}

#tar_pos <- 87864065 + sum(chr_len_vec[1:11])

colnames(dat) <- c("chr","pos","hapIndex","SNPIndex","ED4","index001","index005","ED001","ED005")

plot_dat <- dat
plot_dat[,3] <- abs(plot_dat[,3])
plot_dat[,4] <- abs(plot_dat[,4])
ymax <- max(plot_dat[,3:4],na.rm=T)
## start plot ##
col1 <- rgb(37,67,54,alpha=130,maxColorValue=255)
col2 <- rgb(247,197,102,alpha=130,maxColorValue=255)

pdf(output,width=14,height=3.5)
plot(x=plot_dat[,2],y=plot_dat[,3],xlab="chromosome",ylab="SNP-index / hap-index",bty="l",pch=20, xaxt="n",xaxs="i",yaxs="i",col=col1,ylim=c(0,ymax * 1.1))
# SNP index #
points(x=plot_dat[,2],y=plot_dat[,4],col=col2,pch=20)

sig005col <- rgb(209,187,158,alpha=140,maxColorValue=255)
sig001col <- rgb(167,146,119,alpha=140,maxColorValue=255)

lines(x=plot_dat[,2],y=plot_dat[,6],col=sig001col)
lines(x=plot_dat[,2],y=plot_dat[,7],col=sig005col)


# get chr sep pos #
chr_sep <- NULL
for (i in 2:(length(chr_len_vec)-1)) {
	chr_sep <- c(chr_sep,sum(chr_len_vec[1:i])) 
}
abline(v=chr_sep,lty=2,col=grey(0.8))

# get chr text pos #
chr_text_pos <- NULL
for (i in 2:length(chr_len_vec)) {
	chr_text_pos <- c( chr_text_pos, (sum(chr_len_vec[1:i-1]) + sum(chr_len_vec[1:i]))/2 )
}

text(x=chr_text_pos, y=-ymax*0.09, labels=paste("chr",chr_vec), xpd=NA)

# true QTL #
# points(x=tar_pos,y=ymax*1.15,pch=25,cex=2,col="darkred",xpd=NA)

# legend #

dev.off()
