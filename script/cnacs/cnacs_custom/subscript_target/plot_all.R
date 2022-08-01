# settings 
baf_bottom_up <- 0.6
snp_bottom_up <- -0.2

# pseudo-chromosome end/start position
# this values should be arguments/
chr_pos <- c(0, 249, 491, 689, 880, 1061, 1232, 1391, 1536, 1675, 1809, 1944, 2077, 2191, 2298, 2400, 2491, 2574, 2654, 2713, 2777, 2824, 2875, 3031)


# arguments
inputPath <- commandArgs()[5];
segmentPath <- commandArgs()[6];
centromerePath <- commandArgs()[7];

upper_lod <- commandArgs()[8];
# Because signal.txt.tmp has "4" added copy number value to show same axis figure with alelric information
# So, we needs addition of 4 to LOD.
upper_lod <- as.numeric(upper_lod)+4

lower_lod <- commandArgs()[9];
lower_lod <- as.numeric(lower_lod)+4

outputPath <- commandArgs()[10];


# load signal information like below (pseudo-pos, CN+4, SNP-freq)
# copy number value in this matrix were modified. +4
# 3.624127,7.11002857388,NA
# 3.638575,6.47261433652,0.920989117772
tmp <- scan(inputPath, sep=",")
number <- length(tmp) / 3
input <- matrix(tmp, 3, number)

# Pseudo-position taking into account all chromosome lengths
pos <- input[1,]
# Pseudo ploidy value. ploidy + 4
ploidy <- input[2,]
# actually, these are NOT BAF. BAF x 2 value
baf <- input[3,]
# prepare y position of a-allele and b-allele frequency
aaf_mod <- 2 - baf + baf_bottom_up
baf_mod <- baf + baf_bottom_up


# load segment infromation like below (pseudo-pos1, pseudo-pos2, CN+4, SNP-freq)
# copy number value in this matrix were modified. +4
# 692.26802,878.187024,6.18514907151,0.876011229868
# 1816.440519,1934.235184,6.84950653757,0.694642221311
tmp2 <- scan(segmentPath, sep=",")
number2 <- length(tmp2) / 4
segment <- matrix(tmp2, 4, number2)

# load pseudo-centromere position
tmp3 <- scan(centromerePath)
number3 <- length(tmp3)
centromere <- matrix(tmp3, 1, number3)


# high amp means CN > 4 in this test. because CN in matrix were added "4", we must extract CN > 8
high_idx1 <- ( ploidy > 8 )
amp_tmp <- pos[high_idx1]
high_idx2 <- !is.na(amp_tmp)
# prepare x-axis position of marking 
amp_x1 <- amp_tmp[high_idx2]
amp_x2 <- amp_x1 - 10
amp_x3 <- amp_x1 + 10
amp_num <- length(amp_x1)

hetero_ind <- !is.na(input[3,])
hetero_pos <- pos[hetero_ind]
homo_ind <- is.na(input[3,])
homo_pos <- pos[homo_ind]

x1 <- segment[1,]
x2 <- segment[2,]
y <- segment[3,]
y1 <- segment[4,]
y2 <- 2 - y1


label_pos <- numeric(23)
for ( i in 1:23 ) {
	label_pos[i] <- ( chr_pos[i] + chr_pos[i+1] ) * 0.5
}
chr_label <- 1:22
chr_label <- c(chr_label, 'X')

dummy_x <- c()
dummy_y <- c()

# png(file=outputPath, width=960/72, height=320/72)
png(file=outputPath, width=1920, height=640)
par(bty="n")
plot(dummy_x, dummy_y, xlim=c(0,chr_pos[24]), ylim=c(snp_bottom_up,8), xlab="", ylab="", xaxt="n", yaxt="n")
pa <- par("usr")
rect(pa[1],pa[3],pa[2],pa[4], col=rgb(1,1,1), border=NA)

# Pale yellow back
par(new=T)
rect(chr_pos[2],pa[3],chr_pos[3],pa[4], col=rgb(1,1,0.95), border=NA)
par(new=T)
rect(chr_pos[4],pa[3],chr_pos[5],pa[4], col=rgb(1,1,0.95), border=NA)
par(new=T)
rect(chr_pos[6],pa[3],chr_pos[7],pa[4], col=rgb(1,1,0.95), border=NA)
par(new=T)
rect(chr_pos[8],pa[3],chr_pos[9],pa[4], col=rgb(1,1,0.95), border=NA)
par(new=T)
rect(chr_pos[10],pa[3],chr_pos[11],pa[4], col=rgb(1,1,0.95), border=NA)
par(new=T)
rect(chr_pos[12],pa[3],chr_pos[13],pa[4], col=rgb(1,1,0.95), border=NA)
par(new=T)
rect(chr_pos[14],pa[3],chr_pos[15],pa[4], col=rgb(1,1,0.95), border=NA)
par(new=T)
rect(chr_pos[16],pa[3],chr_pos[17],pa[4], col=rgb(1,1,0.95), border=NA)
par(new=T)
rect(chr_pos[18],pa[3],chr_pos[19],pa[4], col=rgb(1,1,0.95), border=NA)
par(new=T)
rect(chr_pos[20],pa[3],chr_pos[21],pa[4], col=rgb(1,1,0.95), border=NA)
par(new=T)
rect(chr_pos[22],pa[3],chr_pos[23],pa[4], col=rgb(1,1,0.95), border=NA)

# lod_region gray back
par(new=T)
rect(0,upper_lod,chr_pos[23],lower_lod, col=rgb(0.85,0.85,0.85), border=NA)


# copy number plot
par(new=T)
plot(pos, ploidy, type = "p", col=rgb(0.65,0.75,1), pch=20, cex=1.5, xlim=c(0, chr_pos[24]), ylim=c(0, 8), axes=F, xlab="", ylab="")
# allele freq plot, use baf_mod and aaf_mod to set position
par(new=T)
plot(pos, baf_mod, type = "p", col=rgb(0.6,0.9,0.5), pch=20, cex=1.5, xlim=c(0, chr_pos[24]), ylim=c(0, 8), axes=F, xlab="", ylab="")
par(new=T)
plot(pos, aaf_mod, type = "p", col=rgb(0.9,0.5,1), pch=20, cex=1.5, xlim=c(0, chr_pos[24]), ylim=c(0, 8), axes=F, xlab="", ylab="")

# copy number alteration pos line
par(new=T)
segments(x1, y, x2, y, col=rgb(0.85,0,0), lwd=5)





# abnormal allele freq. position line
# baf
segments(x1, y1+baf_bottom_up, x2, y1+baf_bottom_up, col=rgb(0.2,0.3,0), lwd=5)
# aaf
segments(x1, y2+baf_bottom_up, x2, y2+baf_bottom_up, col=rgb(1,0,0), lwd=5)

# hetero SNP position
if ( length(homo_pos) > 0 ) {
	segments(homo_pos, 0+snp_bottom_up, homo_pos, 0.2+snp_bottom_up, col=rgb(0.85,0.85,0.85), lwd=1)
}
if ( length(hetero_pos) > 0 ) {
	segments(hetero_pos, 0+snp_bottom_up, hetero_pos, 0.2+snp_bottom_up, col=rgb(0,0.55,0), lwd=1)
}


# baselines "copy number" 0 to 4 and "allele freq." 0 to 2
# black line, CN=2 and BAF=0.5(allele ratio=1.0)
segments(0,1+baf_bottom_up,chr_pos[24],1+baf_bottom_up)
segments(0,6,chr_pos[24],6)
# gray line, CN=0, 1, 3, 4 and BAF=0, 1
segments(0,baf_bottom_up,chr_pos[24],baf_bottom_up,col=rgb(0.8,0.8,0.8))
segments(0,2+baf_bottom_up,chr_pos[24],2+baf_bottom_up,col=rgb(0.8,0.8,0.8))
segments(0,5,chr_pos[24],5,col=rgb(0.8,0.8,0.8))
segments(0,4,chr_pos[24],4,col=rgb(0.8,0.8,0.8))
segments(0,7,chr_pos[24],7,col=rgb(0.8,0.8,0.8))
segments(0,8,chr_pos[24],8,col=rgb(0.8,0.8,0.8))

# chromosome borderline
abline(v=chr_pos[1], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[2], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[3], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[4], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[5], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[6], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[7], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[8], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[9], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[10], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[11], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[12], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[13], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[14], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[15], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[16], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[17], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[18], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[19], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[20], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[21], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[22], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[23], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=chr_pos[24], lty=2, col=rgb(0.3,0.3,0.3))
abline(v=centromere, lty=3, col=rgb(0.4,0.4,0.5))

# white masking at center (chromosome number area)
rect(-10,2.05+baf_bottom_up,chr_pos[24]+10,3.95, col=rgb(1,1,1), border=NA)


# label
mtext("Copy number", side=2, line=1.0, at=6, cex=1.8)
mtext("Allele freq. × 2", side=2, line=1.0, at=1.6, cex=1.8)
mtext("Hetero SNPs", side=1, line=-1.3, at=-110, cex=1.5)

text(label_pos, (2+baf_bottom_up+4)/2, chr_label, cex=1.8)
mtext("Chromosome", side=1, line=-15.1, at=-110, cex=1.5)


# axis
# par(mar=c(2, 5.5, 2, 2))
par(mar=c(0, 5.5, 0, 0))
axis(2, at = c(0+baf_bottom_up, 1+baf_bottom_up, 2+baf_bottom_up), labels = c(0, 1, 2), las = 0, lwd.ticks=1, cex.axis=1.8, pos = -30)
axis(2, at = c(4, 5, 6, 7, 8), labels = c(0, 1, 2, 3, 4), las = 0, lwd.ticks=1, cex.axis=1.8, pos = -30)


# white masking and mark high amp
rect(-10,8.05,chr_pos[24]+10,10.5, col=rgb(1,1,1), border=NA)
par(new=T)
if ( amp_num > 0 ) {
    for ( i in 1:amp_num ) {
        if ( i > 1 ) {
            diff <- amp_x1[i] - amp_x1[i-1]
            if ( diff < 20 ) {
                next
            }
        }
        polygon(c(amp_x1[i], amp_x2[i], amp_x3[i]), c(8.45, 8.05, 8.05), col="tomato1")
        # lines(c(amp_x2[i], amp_x3[i]), c(8, 8), col="tomato1")
        # points(amp_x1[i], 8, pch=18, cex=1, col="indianred1")
    }
}

# white masking and mark high amp
rect(-10,snp_bottom_up,chr_pos[24]+10,-0.2+snp_bottom_up, col=rgb(1,1,1), border=NA)

dev.off()
