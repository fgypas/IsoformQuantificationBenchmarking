## Foivos Gypas & Alexander Kanitz
## 24-JAN-2015

plot_accuracy_per_subset <- function(data, prefix, left_margin=NA, colors=NULL, chars=NULL, min_x=NULL, max_x=NULL, xlab=NULL) {

        ## Load required packages
        if ( suppressWarnings(suppressPackageStartupMessages(require("Hmisc"))) == FALSE ) { stop("Package 'Hmisc' required!\nExecution aborted.") } # for converting coordinates
        if ( suppressWarnings(suppressPackageStartupMessages(require("KernSmooth"))) == FALSE ) { stop("Package 'KernSmooth' required!\nExecution aborted.") }

        ## Parameters
        if (is.null(xlab)) xlab <- expression("Accuracy (r"[s]*")")
        yticks <- rev( seq(1, (nrow(data) * 2 - 1), 2) / (2 * nrow(data)))
        if (is.null(colors)) colors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999", "#000000")
        if (is.null(chars)) chars <- c(0, 1, 2, 5, 6, 3, 4, 7, 10, 12)
        margins <- c(1.02, left_margin, 0.22, 0.22)
        line_color <- "#555555"

        ## Check if dataframe dimensions are in allowed range
        if(ncol(data) > 10) stop("Too many subsets (columns) in dataframe. The maximum number is 10.")

        ## Get/set x axis limits
        if (is.null(min_x)) min_x <- floor(min(data))
        if (is.null(max_x)) max_x <- ceiling(max(data))

	## Get width of largest label
	pdf("/dev/null")
	max_width <- max(strwidth(rownames(data), units="inches")) + 0.44
	dev.off()

	## Set left margin if not fixed
	if (is.na(margins[[2]])) margins[[2]] <- max_width

        #---> DATA <---#

        ## PDF
        pdf(paste(prefix, "pdf", sep="."), width=max_width + 5, height=1.24 + nrow(data) * 0.36)
        par(mai=margins)
        plot(0, type="n", xlim=c(min_x, max_x), ylim=c(0, 1), main="", xlab=xlab, ylab="", cex.lab=1.5, axes=FALSE)
	yticks <- cnvrt.coords(NA, yticks, input="plt")$usr$y
        axis(1, at=seq(min_x, max_x, 0.1))
        axis(2, at=yticks, labels=rownames(data), par(las=1))
        for (index in 1:nrow(data) ) {
                lines(c(min_x, max_x), rep(yticks[index], 2), col=line_color)
                points(data[index, ], rep(yticks[index], ncol(data)), type="p", pch=chars[1:ncol(data)], col=colors[1:ncol(data)])
        }
        dev.off()

        ## PNG
        png(paste(prefix, "png", sep="."), width=max_width + 5, height=1.24 + nrow(data) * 0.36, res=90, units="in", pointsize=20)
        par(mai=margins)
        plot(0, type="n", xlim=c(min_x, max_x), ylim=c(0, 1), main="", xlab=xlab, ylab="", cex.lab=1.5, axes=FALSE)
        axis(1, at=seq(min_x, max_x, 0.1))
        axis(2, at=yticks, labels=rownames(data), par(las=1))
        for (index in 1:nrow(data) ) {
                lines(c(min_x, max_x), rep(yticks[index], 2), col=line_color)
                points(data[index, ], rep(yticks[index], ncol(data)), type="p", pch=chars[1:ncol(data)], col=colors[1:ncol(data)])
        }
        dev.off()

        ## SVG
        svg(paste(prefix, "svg", sep="."), width=max_width + 5, height=1.04 + nrow(data) * 0.36)
        par(mai=margins)
        plot(0, type="n", xlim=c(min_x, max_x), ylim=c(0, 1), main="", xlab=xlab, ylab="", cex.lab=1.5, axes=FALSE)
        axis(1, at=seq(min_x, max_x, 0.1))
        axis(2, at=yticks, labels=rownames(data), par(las=1))
        for (index in 1:nrow(data) ) {
                lines(c(min_x, max_x), rep(yticks[index], 2), col=line_color)
                points(data[index, ], rep(yticks[index], ncol(data)), type="p", pch=chars[1:ncol(data)], col=colors[1:ncol(data)])
        }
        dev.off()

	#---> LEGEND <---#

	## PDF
        pdf(paste(prefix, "legend", "pdf", sep="."))
        plot(0, type="n", axes=FALSE, cex.sub=1.5, xlab="", ylab="", main="")
        legend("center", colnames(data), bty="n", col=colors[1:length(colnames(data))], pch=chars[1:length(colnames(data))])
        dev.off()

	## PNG
        png(paste(prefix, "legend", "png", sep="."), pointsize=20)
        plot(0, type="n", axes=FALSE, cex.sub=1.5, xlab="", ylab="", main="")
        legend("center", colnames(data), bty="n", col=colors[1:length(colnames(data))], pch=chars[1:length(colnames(data))])
        dev.off()

	## SVG
        svg(paste(prefix, "legend", "svg", sep="."))
        plot(0, type="n", axes=FALSE, cex.sub=1.5, xlab="", ylab="", main="")
        legend("center", colnames(data), bty="n", col=colors[1:length(colnames(data))], pch=chars[1:length(colnames(data))])
        dev.off()

}
