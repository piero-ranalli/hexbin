## the following code comes from the hydrosanity package
## and is (c) 2007 Felix Andrews <felix@nfrac.org>
## https://github.com/cran/hydrosanity/blob/1428e15f158b6489314184ef0cc2aa2d90f43c33/R/timeblob_plots.R
##
##
##  hydrosanity is released under a GPL>=2 license,
##  while hexbin with a GPL==2 license, so we can include
##  hydrosanity's code in hexbin


### AXIS STUFF ###

# logLim should be given in log10 scale
logAxisComponents <- function(logLim, label=T) {
	if (length(logLim) != 2) { stop("'logLim' must be of length 2") }
	lim <- logLim
	# integer log powers should always be labelled
	at <- seq(ceiling(min(lim)), max(lim))
	hasLabel <- rep(T, length(at))
	mags <- seq(floor(min(lim)), max(lim))
	if (diff(range(lim)) < 3) {
		# make it a linear sequence log-transformed (gradient effect)
		newAt <- log10(
			sequence(rep(9,length(mags))) * 
				rep(10^mags, each=9) )
		newAt <- newAt[(min(lim) <= newAt) & (newAt <= max(lim))]
		# draw labels for the previous sequence (integer log powers)
		hasLabel <- rep(F, length(newAt))
		hasLabel[match(at, newAt)] <- T
		at <- newAt
	}
	if ((diff(range(lim)) < 1.1)) {
		# label all ticks
		hasLabel[] <- T
	}
	if ((diff(range(lim)) < 0.5)) {
		# go into more detail: ticks between sub-orders of magnitude
		newAt <- log10(
			0.1*(9+sequence(rep(90,length(mags)))) * 
				rep(10^mags, each=90) )
		newAt <- newAt[(min(lim) <= newAt) & (newAt <= max(lim))]
		# keep labels only for the previous sequence
		hasLabel <- rep(F, length(newAt))
		hasLabel[match(at, newAt)] <- T
		at <- newAt
	}
	if ((diff(range(lim)) < 0.35)) {
		# label all ticks
		hasLabel[] <- T
	}
	if ((diff(range(lim)) < 0.05)) {
		# go into more detail: ticks between sub-sub-orders of magnitude
		newAt <- log10(
			0.01*(99+sequence(rep(900,length(mags)))) * 
				rep(10^mags, each=900) )
		newAt <- newAt[(min(lim) <= newAt) & (newAt <= max(lim))]
		# keep labels only for the previous sequence
		hasLabel <- rep(F, length(newAt))
		hasLabel[match(at, newAt)] <- T
		at <- newAt
	}
	if ((diff(range(lim)) < 0.035)) {
		# label all ticks
		hasLabel[] <- T
	}
	makeLabels <- function(x) {
		# use sapply rather than format directly, so not common format
		newLabels <- sapply(10^x, format, digits=2, scientific=1)
		# large powers of 10 are better presented with an exponent
		simplePowers <- (x >= 5) & (x == round(x))
		newLabels[simplePowers] <- paste("10^",x[simplePowers],sep='')
		newLabels
	}
	# blank labels
	atLabels <- rep("", length(at))
	if (label) {
		atLabels[hasLabel] <- makeLabels(at[hasLabel])
	}
	return(list(at=at, label=atLabels))
}

grid.yaxis.log <- function(logLim=as.numeric(convertY(unit(c(0,1), "npc"), "native")), label=T, draw=T, name=NULL, ...) {
	axisStuff <- logAxisComponents(logLim, label=label)
	if (label==F) { axisStuff$label <- F }
	tmp <- yaxisGrob(at=axisStuff$at, label=axisStuff$label, name=name, ...)
	if (label) {
            tmp <- editGrob(tmp, gPath("labels"), check.overlap=F)
        }
        ## put ticks inside, see here
        ## http://stackoverflow.com/questions/28949001/mirroring-axis-ticks-in-ggplot2
        tmp$children[[2]]$x1 = tmp$children[[2]]$x1 - unit(-0.5, "cm")
	if (draw) { grid.draw(tmp) }
	tmp
}

grid.xaxis.log <- function(logLim=as.numeric(convertX(unit(c(0,1), "npc"), "native")), label=T, draw=T, name=NULL, ...) {
	axisStuff <- logAxisComponents(logLim, label=label)
	if (label==F) { axisStuff$label <- F }
	tmp <- xaxisGrob(at=axisStuff$at, label=axisStuff$label, name=name, ...)
	if (label) {
            tmp <- editGrob(tmp, gPath("labels"), check.overlap=F)
        }
        ## put ticks inside
        tmp$children[[2]]$y1 = tmp$children[[2]]$y1 - unit(-0.5, "cm")
	if (draw) { grid.draw(tmp) }
	tmp
}

lattice.y.sqrt <- function(lim, ...) {
	arglist <- list(...)
	tmp <- yscale.components.default(lim, ...)
	##tmp$left$labels$labels <- format(tmp$left$labels$at ^ 2)
	tmp$left$ticks$at <- sqrt(pretty(pmax(0,tmp$num.limit) ^ 2))
	tmp$left$labels$at <- tmp$left$ticks$at
	tmp$left$labels$labels <- format(tmp$left$labels$at ^ 2)
	return(tmp)
}

lattice.x.sqrt <- function(lim, ...) {
	arglist <- list(...)
	tmp <- xscale.components.default(lim, ...)
	##tmp$bottom$labels$labels <- format(tmp$bottom$labels$at ^ 2)
	tmp$bottom$ticks$at <- sqrt(pretty(pmax(0,tmp$num.limit) ^ 2))
	tmp$bottom$labels$at <- tmp$bottom$ticks$at
	tmp$bottom$labels$labels <- format(tmp$bottom$labels$at ^ 2)
	return(tmp)
}


lattice.y.prettylog <- function(lim, ...) {
	arglist <- list(...)
	have.log <- (!is.null(arglist$logsc)) && (!identical(arglist$logsc, F))
	tmp <- yscale.components.default(lim, ...)
	if (have.log) {
		axisStuff <- logAxisComponents(lim)
		tmp$left$ticks$at <- axisStuff$at
		tmp$left$labels$at <- axisStuff$at
		tmp$left$labels$labels <- axisStuff$label
		tmp$left$labels$check.overlap <- F
	}
	return(tmp)
}

lattice.x.prettylog <- function(lim, ...) {
	arglist <- list(...)
	have.log <- (!is.null(arglist$logsc)) && (!identical(arglist$logsc, F))
	tmp <- xscale.components.default(lim, ...)
	if (have.log) {
		axisStuff <- logAxisComponents(lim)
		tmp$bottom$ticks$at <- axisStuff$at
		tmp$bottom$labels$at <- axisStuff$at
		tmp$bottom$labels$labels <- axisStuff$label
		tmp$bottom$labels$check.overlap <- F
	}
	return(tmp)
}

