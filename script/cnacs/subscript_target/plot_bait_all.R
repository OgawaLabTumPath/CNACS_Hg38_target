inputPath <- commandArgs()[5];
centromerePath <- commandArgs()[6];
outputPath <- commandArgs()[7];

tmp <- scan(inputPath, sep=",")
number <- length(tmp) / 3
input <- matrix(tmp, 3, number)

idx_all <- ( input[2,] == 1 )
pos_all <- input[1,idx_all]

idx_snp <- ( input[3,] == 1 )
pos_snp <- input[1,idx_snp]


tmp2 <- scan(centromerePath)
number2 <- length(tmp2)
centromere <- matrix(tmp2, 1, number2)


if (Sys.getenv("HUMAN_GENOME_VERSION") == '' || Sys.getenv("HUMAN_GENOME_VERSION") == 'hg38') {
	chr_pos <- c(0, 249, 491, 689, 880, 1061, 1232, 1391, 1536, 1675, 1809, 1944, 2077, 2191, 2298, 2400, 2491, 2574, 2654, 2713, 2777, 2824, 2875, 3031)
} else {
	chr_pos <- c(0, 249, 492, 690, 882, 1063, 1234, 1393, 1539, 1680, 1816, 1951, 2085, 2200, 2307, 2410, 2500, 2581, 2659, 2719, 2782, 2830, 2881, 3036)
}
label_pos <- numeric(23)
for ( i in 1:23 ) {
	label_pos[i] <- ( chr_pos[i] + chr_pos[i+1] ) * 0.5
}
chr_label <- 1:22
chr_label <- c(chr_label, 'X')

dummy_x <- c()
dummy_y <- c()

pdf(file=outputPath, onefile=FALSE, width=960/72, height=320/72)
plot(dummy_x, dummy_y, xlim=c(0, chr_pos[24]), ylim=c(0,3), xlab="", ylab="", xaxt="n", yaxt="n")
pa <- par("usr")
rect(pa[1],pa[3],pa[2],pa[4], col=rgb(1,1,1), border=NA)
par(new=T)
for (i in 1:11) {
	rect(chr_pos[2 * i], pa[3], chr_pos[2 * i + 1], pa[4], col=rgb(0.97, 0.97, 0.97), border=NA)
	par(new=T)
}
segments(pos_all, 1.7, pos_all, 2.3, lwd=0.05)
segments(pos_snp, 0.7, pos_snp, 1.3, lwd=0.05)
segments(0, 1, chr_pos[24], 1)
segments(0, 2, chr_pos[24], 2)

for (i in 1:24) {
	abline(v=chr_pos[i], lty=2, col=rgb(0.3,0.3,0.3))
}

abline(v=centromere, lty=3, col=rgb(0.7,0.7,0.7))
axis(1, at = chr_pos, labels=F, las = 0, lwd.ticks=1)
mtext(chr_label, side=1, line=1, at=label_pos, cex=0.8)
dev.off()
