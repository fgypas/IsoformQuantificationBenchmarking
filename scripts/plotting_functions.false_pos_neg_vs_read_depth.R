## Alexander Kanitz
## 24-JAN-2015

#<--- FALSE POSITIVE RATE VS READ DEPTH --->#
plot_false_pos_depth <- function(df, prefix, min_y, max_y) {
        colors <- c("#B58900", "#CB4B16", "#DC322f", "#D33682", "#A300A3", "#6C71c4", "#268BD2", "#2AA198", "#073642", "#859900", "#009914", "#000000", "#252525", "#737373", "#bdbdbd", "#f0f0f0")
        characters <- 0:15
	reads <- c(1,3,10,30,100)
	xlab <- "Library size (million reads)"
	ylab <- "False positive rate"

        ## Load required package
        if ( suppressWarnings(suppressPackageStartupMessages(require("KernSmooth"))) == FALSE ) { stop("Package 'KernSmooth' required!\nExecution aborted.") }

	## Data
        pdf(paste(prefix, "pdf", sep="."))
        par(mar=c(5.1,5.1,1.1,1.1))
        plot(0, type="n", xlim=c(0,2), ylim=c(min_y,max_y), xlab=xlab, ylab=ylab, axes=FALSE, cex.lab=1.5)
        axis(1,at=log10(reads), labels=reads)
        axis(2, at=seq(min_y,max_y,0.1))
        lapply(1:nrow(df), function(nrow) {
                lines(log10(reads), df[nrow,],  type='o', pch=characters[nrow], col=colors[nrow])
        })
        dev.off()

        png(paste(prefix, "png", sep="."), pointsize=20)
        par(mar=c(5.1,5.1,1.1,1.1))
        plot(0, type="n", xlim=c(0,2), ylim=c(min_y,max_y), xlab=xlab, ylab=ylab, axes=FALSE, cex.lab=1.5)
        axis(1,at=log10(reads), labels=reads)
        axis(2, at=seq(min_y,max_y,0.1))
        lapply(1:nrow(df), function(nrow) {
                lines(log10(reads), df[nrow,],  type='o', pch=characters[nrow], col=colors[nrow])
        })
        dev.off()

        svg(paste(prefix, "svg", sep="."))
        par(mar=c(5.1,5.1,1.1,1.1))
        plot(0, type="n", xlim=c(0,2), ylim=c(min_y,max_y), xlab=xlab, ylab=ylab, axes=FALSE, cex.lab=1.5)
        axis(1,at=log10(reads), labels=reads)
        axis(2, at=seq(min_y,max_y,0.1))
        lapply(1:nrow(df), function(nrow) {
                lines(log10(reads), df[nrow,],  type='o', pch=characters[nrow], col=colors[nrow])
        })
        dev.off()

	## Legends
        pdf(paste(prefix, "legend", "pdf", sep="."))
        plot(0, type="n", axes=FALSE, cex.sub=1.5, xlab="", ylab="", main="")
        legend("center", rownames(df), lty=1, bty="n", col=colors[1:length(rownames(df))], pch=characters[1:length(rownames(df))])
        dev.off()

        png(paste(prefix, "legend", "png", sep="."), pointsize=20)
        plot(0, type="n", axes=FALSE, cex.sub=1.5, xlab="", ylab="", main="")
        legend("center", rownames(df), lty=1, bty="n", col=colors[1:length(rownames(df))], pch=characters[1:length(rownames(df))])
        dev.off()

        svg(paste(prefix, "legend", "svg", sep="."))
        plot(0, type="n", axes=FALSE, cex.sub=1.5, xlab="", ylab="", main="")
        legend("center", rownames(df), lty=1, bty="n", col=colors[1:length(rownames(df))], pch=characters[1:length(rownames(df))])
        dev.off()

}

#<--- FALSE NEGATIVE RATE VS READ DEPTH --->#
plot_false_neg_depth <- function(df, prefix, min_y, max_y) {
        colors <- c("#B58900", "#CB4B16", "#DC322f", "#D33682", "#A300A3", "#6C71c4", "#268BD2", "#2AA198", "#073642", "#859900", "#009914", "#000000", "#252525", "#737373", "#bdbdbd", "#f0f0f0")
        characters <- 0:15
        reads <- c(1,3,10,30,100)
        xlab <- "Library size (million reads)"
        ylab <- "False negative rate"

        ## Data
        pdf(paste(prefix, "pdf", sep="."))
        par(mar=c(5.1,5.1,1.1,1.1))
        plot(0, type="n", xlim=c(0,2), ylim=c(min_y,max_y), xlab=xlab, ylab=ylab, axes=FALSE, cex.lab=1.5)
        axis(1,at=log10(reads), labels=reads)
        axis(2, at=seq(min_y,max_y,0.1))
        lapply(1:nrow(df), function(nrow) {
                lines(log10(reads), df[nrow,],  type='o', pch=characters[nrow], col=colors[nrow])
        })
        dev.off()

        png(paste(prefix, "png", sep="."), pointsize=20)
        par(mar=c(5.1,5.1,1.1,1.1))
        plot(0, type="n", xlim=c(0,2), ylim=c(min_y,max_y), xlab=xlab, ylab=ylab, axes=FALSE, cex.lab=1.5)
        axis(1,at=log10(reads), labels=reads)
        axis(2, at=seq(min_y,max_y,0.1))
        lapply(1:nrow(df), function(nrow) {
                lines(log10(reads), df[nrow,],  type='o', pch=characters[nrow], col=colors[nrow])
        })
        dev.off()

        svg(paste(prefix, "svg", sep="."))
        par(mar=c(5.1,5.1,1.1,1.1))
        plot(0, type="n", xlim=c(0,2), ylim=c(min_y,max_y), xlab=xlab, ylab=ylab, axes=FALSE, cex.lab=1.5)
        axis(1,at=log10(reads), labels=reads)
        axis(2, at=seq(min_y,max_y,0.1))
        lapply(1:nrow(df), function(nrow) {
                lines(log10(reads), df[nrow,],  type='o', pch=characters[nrow], col=colors[nrow])
        })
        dev.off()

        ## Legends
        pdf(paste(prefix, "legend", "pdf", sep="."))
        plot(0, type="n", axes=FALSE, cex.sub=1.5, xlab="", ylab="", main="")
        legend("center", rownames(df), lty=1, bty="n", col=colors[1:length(rownames(df))], pch=characters[1:length(rownames(df))])
        dev.off()

        png(paste(prefix, "legend", "png", sep="."), pointsize=20)
        plot(0, type="n", axes=FALSE, cex.sub=1.5, xlab="", ylab="", main="")
        legend("center", rownames(df), lty=1, bty="n", col=colors[1:length(rownames(df))], pch=characters[1:length(rownames(df))])
        dev.off()

        svg(paste(prefix, "legend", "svg", sep="."))
        plot(0, type="n", axes=FALSE, cex.sub=1.5, xlab="", ylab="", main="")
        legend("center", rownames(df), lty=1, bty="n", col=colors[1:length(rownames(df))], pch=characters[1:length(rownames(df))])
        dev.off()

}

#<--- TRUE POSITIVE RATE VS READ DEPTH --->#
plot_true_pos_depth <- function(df, prefix, min_y, max_y) {
        colors <- c("#B58900", "#CB4B16", "#DC322f", "#D33682", "#A300A3", "#6C71c4", "#268BD2", "#2AA198", "#073642", "#859900", "#009914", "#000000", "#252525", "#737373", "#bdbdbd", "#f0f0f0")
        characters <- 0:15
        reads <- c(1,3,10,30,100)
        xlab <- "Library size (million reads)"
        ylab <- "True positive rate"

        ## Data
        pdf(paste(prefix, "pdf", sep="."))
        par(mar=c(5.1,5.1,1.1,1.1))
        plot(0, type="n", xlim=c(0,2), ylim=c(min_y,max_y), xlab=xlab, ylab=ylab, axes=FALSE, cex.lab=1.5)
        axis(1,at=log10(reads), labels=reads)
        axis(2, at=seq(min_y,max_y,0.1))
        lapply(1:nrow(df), function(nrow) {
                lines(log10(reads), df[nrow,],  type='o', pch=characters[nrow], col=colors[nrow])
        })
        dev.off()

        png(paste(prefix, "png", sep="."), pointsize=20)
        par(mar=c(5.1,5.1,1.1,1.1))
        plot(0, type="n", xlim=c(0,2), ylim=c(min_y,max_y), xlab=xlab, ylab=ylab, axes=FALSE, cex.lab=1.5)
        axis(1,at=log10(reads), labels=reads)
        axis(2, at=seq(min_y,max_y,0.1))
        lapply(1:nrow(df), function(nrow) {
                lines(log10(reads), df[nrow,],  type='o', pch=characters[nrow], col=colors[nrow])
        })
        dev.off()

        svg(paste(prefix, "svg", sep="."))
        par(mar=c(5.1,5.1,1.1,1.1))
        plot(0, type="n", xlim=c(0,2), ylim=c(min_y,max_y), xlab=xlab, ylab=ylab, axes=FALSE, cex.lab=1.5)
        axis(1,at=log10(reads), labels=reads)
        axis(2, at=seq(min_y,max_y,0.1))
        lapply(1:nrow(df), function(nrow) {
                lines(log10(reads), df[nrow,],  type='o', pch=characters[nrow], col=colors[nrow])
        })
        dev.off()

        ## Legends
        pdf(paste(prefix, "legend", "pdf", sep="."))
        plot(0, type="n", axes=FALSE, cex.sub=1.5, xlab="", ylab="", main="")
        legend("center", rownames(df), lty=1, bty="n", col=colors[1:length(rownames(df))], pch=characters[1:length(rownames(df))])
        dev.off()

        png(paste(prefix, "legend", "png", sep="."), pointsize=20)
        plot(0, type="n", axes=FALSE, cex.sub=1.5, xlab="", ylab="", main="")
        legend("center", rownames(df), lty=1, bty="n", col=colors[1:length(rownames(df))], pch=characters[1:length(rownames(df))])
        dev.off()

        svg(paste(prefix, "legend", "svg", sep="."))
        plot(0, type="n", axes=FALSE, cex.sub=1.5, xlab="", ylab="", main="")
        legend("center", rownames(df), lty=1, bty="n", col=colors[1:length(rownames(df))], pch=characters[1:length(rownames(df))])
        dev.off()

}

#<--- TRUE NEGATIVE RATE VS READ DEPTH --->#
plot_true_neg_depth <- function(df, prefix, min_y, max_y) {
        colors <- c("#B58900", "#CB4B16", "#DC322f", "#D33682", "#A300A3", "#6C71c4", "#268BD2", "#2AA198", "#073642", "#859900", "#009914", "#000000", "#252525", "#737373", "#bdbdbd", "#f0f0f0")
        characters <- 0:15
        reads <- c(1,3,10,30,100)
        xlab <- "Library size (million reads)"
        ylab <- "True negative rate"

        ## Data
        pdf(paste(prefix, "pdf", sep="."))
        par(mar=c(5.1,5.1,1.1,1.1))
        plot(0, type="n", xlim=c(0,2), ylim=c(min_y,max_y), xlab=xlab, ylab=ylab, axes=FALSE, cex.lab=1.5)
        axis(1,at=log10(reads), labels=reads)
        axis(2, at=seq(min_y,max_y,0.1))
        lapply(1:nrow(df), function(nrow) {
                lines(log10(reads), df[nrow,],  type='o', pch=characters[nrow], col=colors[nrow])
        })
        dev.off()

        png(paste(prefix, "png", sep="."), pointsize=20)
        par(mar=c(5.1,5.1,1.1,1.1))
        plot(0, type="n", xlim=c(0,2), ylim=c(min_y,max_y), xlab=xlab, ylab=ylab, axes=FALSE, cex.lab=1.5)
        axis(1,at=log10(reads), labels=reads)
        axis(2, at=seq(min_y,max_y,0.1))
        lapply(1:nrow(df), function(nrow) {
                lines(log10(reads), df[nrow,],  type='o', pch=characters[nrow], col=colors[nrow])
        })
        dev.off()

        svg(paste(prefix, "svg", sep="."))
        par(mar=c(5.1,5.1,1.1,1.1))
        plot(0, type="n", xlim=c(0,2), ylim=c(min_y,max_y), xlab=xlab, ylab=ylab, axes=FALSE, cex.lab=1.5)
        axis(1,at=log10(reads), labels=reads)
        axis(2, at=seq(min_y,max_y,0.1))
        lapply(1:nrow(df), function(nrow) {
                lines(log10(reads), df[nrow,],  type='o', pch=characters[nrow], col=colors[nrow])
        })
        dev.off()

        ## Legends
        pdf(paste(prefix, "legend", "pdf", sep="."))
        plot(0, type="n", axes=FALSE, cex.sub=1.5, xlab="", ylab="", main="")
        legend("center", rownames(df), lty=1, bty="n", col=colors[1:length(rownames(df))], pch=characters[1:length(rownames(df))])
        dev.off()

        png(paste(prefix, "legend", "png", sep="."), pointsize=20)
        plot(0, type="n", axes=FALSE, cex.sub=1.5, xlab="", ylab="", main="")
        legend("center", rownames(df), lty=1, bty="n", col=colors[1:length(rownames(df))], pch=characters[1:length(rownames(df))])
        dev.off()

        svg(paste(prefix, "legend", "svg", sep="."))
        plot(0, type="n", axes=FALSE, cex.sub=1.5, xlab="", ylab="", main="")
        legend("center", rownames(df), lty=1, bty="n", col=colors[1:length(rownames(df))], pch=characters[1:length(rownames(df))])
        dev.off()

}
