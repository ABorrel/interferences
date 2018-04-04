#!/usr/bin/env Rscript

library(ggplot2)


radial.plot = function (lengths, radial.pos = NULL, labels = NA, label.pos = NULL, 
                        radlab = TRUE, start = 0, clockwise = FALSE, rp.type = "r", 
                        label.prop = 1.4, main = "", xlab = "", ylab = "", line.col = par("fg"), 
                        lty = par("lty"), lwd = 4, mar = c(10, 10, 10, 10), 
                        show.grid = TRUE, show.grid.labels = 3, show.radial.grid = TRUE, 
                        grid.col = "#333333", grid.bg = "transparent", grid.left = FALSE, 
                        grid.unit = NULL, point.symbols = NULL, point.col = NULL, 
                        show.centroid = FALSE, radial.lim = NULL, radial.labels = NULL, 
                        boxed.radial = TRUE, poly.col = NULL, add = FALSE, ...) 
{
  
  #modif environement
  par(cex.axis = 3)
  par(cex.lab = 2)	
  #print (par())
  
  if (is.null(radial.lim)) 
    radial.lim <- range(lengths)
  length.dim <- dim(lengths)
  if (is.null(length.dim)) {
    npoints <- length(lengths)
    nsets <- 1
    lengths <- matrix(lengths, nrow = 1)
  }
  else {
    npoints <- length.dim[2]
    nsets <- length.dim[1]
    lengths <- as.matrix(lengths)
  }
  lengths <- lengths - radial.lim[1]
  lengths[lengths < 0] <- NA
  if (is.null(radial.pos[1])) 
    radial.pos <- seq(0, pi * (2 - 2/npoints), length.out = npoints)
  radial.pos.dim <- dim(radial.pos)
  if (is.null(radial.pos.dim)) 
    radial.pos <- matrix(rep(radial.pos, nsets), nrow = nsets, 
                         byrow = TRUE)
  else radial.pos <- as.matrix(radial.pos)
  if (clockwise) 
    radial.pos <- -radial.pos
  if (start) 
    radial.pos <- radial.pos + start
  if (show.grid) {
    if (length(radial.lim) < 3) 
      grid.pos <- pretty(radial.lim)
    else grid.pos <- radial.lim
    if (grid.pos[1] < radial.lim[1]) 
      grid.pos <- grid.pos[-1]
    maxlength <- max(grid.pos - radial.lim[1])
    angles <- seq(0, 1.96 * pi, by = 0.04 * pi)
  }
  else {
    grid.pos <- NA
    maxlength <- diff(radial.lim)
  }
  oldpar <- par("xpd", "mar", "pty")
  if (!add) {
    par(mar = mar, pty = "s")
    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), 
         type = "n", axes = FALSE, main = main, xlab = xlab, 
         ylab = ylab)
    if (show.grid) {
      for (i in seq(length(grid.pos), 1, by = -1)) {
        xpos <- cos(angles) * (grid.pos[i] - radial.lim[1])
        ypos <- sin(angles) * (grid.pos[i] - radial.lim[1])
        polygon(xpos, ypos, border = grid.col, col = grid.bg)
      }
    }
  }
  par(xpd = TRUE)
  if (length(line.col) < nsets) 
    line.col <- 1:nsets
  if (length(rp.type) < nsets) 
    rp.type <- rep(rp.type, length.out = nsets)
  if (length(point.symbols) < nsets) 
    point.symbols <- rep(point.symbols, length.out = nsets)
  if (length(point.col) < nsets) 
    point.col <- rep(point.col, length.out = nsets)
  if (length(poly.col) < nsets) 
    poly.col <- rep(poly.col, length.out = nsets)
  if (length(lty) < nsets) 
    lty <- rep(lty, length.out = nsets)
  if (length(lwd) < nsets) 
    lwd <- rep(lwd, length.out = nsets)
  for (i in 1:nsets) {
    if (nsets > 1) {
      linecol <- line.col[i]
      polycol <- poly.col[i]
      pointcol <- point.col[i]
      pointsymbols <- point.symbols[i]
      ltype <- lty[i]
      lwidth <- lwd[i]
    }
    else {
      linecol <- line.col
      polycol <- poly.col
      pointcol <- point.col
      pointsymbols <- point.symbols
      ltype <- lty
      lwidth <- lwd
    }
    rptype <- unlist(strsplit(rp.type[i], ""))
    if (match("s", rptype, 0)) {
      if (is.null(pointsymbols)) 
        pointsymbols <- i
      if (is.null(pointcol)) 
        pointcol <- i
    }
    xpos <- cos(radial.pos[i, ]) * lengths[i, ]
    ypos <- sin(radial.pos[i, ]) * lengths[i, ]
    if (match("r", rptype, 0)) 
      segments(0, 0, xpos, ypos, col = linecol, lty = ltype, 
               lwd = lwidth, ...)
    if (match("p", rptype, 0)){ 
      polygon(xpos, ypos, border = linecol, col = polycol, 
              lty = ltype, lwd = lwidth, ...)
      print(ypos[1])
    }
    if (match("s", rptype, 0)) 
      points(xpos, ypos, pch = pointsymbols, col = pointcol, 
             ...)
    if (show.centroid) 
      if (match("p", rptype, 0)) {
        nvertices <- length(xpos)
        polygonarea <- xpos[nvertices] * ypos[1] - xpos[1] * 
          ypos[nvertices]
        for (vertex in 1:(nvertices - 1)) polygonarea <- polygonarea + 
          xpos[vertex] * ypos[vertex + 1] - xpos[vertex + 
                                                   1] * ypos[vertex]
        polygonarea <- polygonarea/2
        centroidx <- (xpos[nvertices] + xpos[1]) * (xpos[nvertices] * 
                                                      ypos[1] - xpos[1] * ypos[nvertices])
        centroidy <- (ypos[nvertices] + ypos[1]) * (xpos[nvertices] * 
                                                      ypos[1] - xpos[1] * ypos[nvertices])
        for (vertex in 1:(nvertices - 1)) {
          centroidx <- centroidx + (xpos[vertex] + xpos[vertex + 
                                                          1]) * (xpos[vertex] * ypos[vertex + 1] - 
                                                                   xpos[vertex + 1] * ypos[vertex])
          centroidy <- centroidy + (ypos[vertex] + ypos[vertex + 
                                                          1]) * (xpos[vertex] * ypos[vertex + 1] - 
                                                                   xpos[vertex + 1] * ypos[vertex])
        }
        points(centroidx/(6 * polygonarea), centroidy/(6 * 
                                                         polygonarea), col = point.col[i], pch = point.symbols[i], 
               cex = 2, ...)
      }
    else points(mean(xpos), mean(ypos), col = pointcol, 
                pch = pointsymbols, cex = 2, ...)
  }
  if (!add) {
    if (is.na(labels[1])) {
      label.pos <- seq(0, 1.8 * pi, length = 9)
      labels <- as.character(round(label.pos, 2))
    }
    if (is.null(label.pos[1])) {
      lablen <- length(labels)
      label.pos <- seq(0, pi * (2 - 2/lablen), length.out = lablen)
    }
    if (clockwise) 
      label.pos <- -label.pos
    if (start) 
      label.pos <- label.pos + start
    xpos <- cos(label.pos) * maxlength
    ypos <- sin(label.pos) * maxlength
    if (show.radial.grid) 
      segments(0, 0, xpos, ypos, col = grid.col)
    xpos <- cos(label.pos) * maxlength * label.prop
    ypos <- sin(label.pos) * maxlength * label.prop
    if (radlab) {
      for (label in 1:length(labels)) {
        labelsrt <- (180 * label.pos[label]/pi) + 180 * 
          (label.pos[label] > pi/2 && label.pos[label] < 
             3 * pi/2)
        
        text(xpos[label], ypos[label], labels[label], 
             cex = par("cex.axis"), srt = labelsrt)
      }
    }
    else{
      ypos[as.integer(length (ypos)/4) + 2]=ypos[as.integer(length (ypos)/4) +2] + 0.01
      ypos[as.integer(length (ypos)/4)+1 ]=ypos[as.integer(length (ypos)/4)+1 ] + 0.025
      ypos[as.integer(length (ypos)/4)+3 ]=ypos[as.integer(length (ypos)/4)+3 ] + 0.01
      
      
      ypos[as.integer(length (ypos)/1.3) ]=ypos[as.integer(length (ypos)/1.3) ] - 0.02
      ypos[as.integer(length (ypos)/1.3) +1 ]=ypos[as.integer(length (ypos)/1.3) ] + 0.015
      ypos[as.integer(length (ypos)/1.3) -1 ]=ypos[as.integer(length (ypos)/1.3) ] + 0.015
      
      xpos[as.integer(length (xpos)/1.3) +1 ]=xpos[as.integer(length (xpos)/1.3) ] + 0.2
      xpos[as.integer(length (xpos)/1.3) -1 ]=xpos[as.integer(length (xpos)/1.3) ]- 0.2
      
      
      boxed.labels(xpos, ypos, labels, ypad = 0.7, border = FALSE, 
                   cex = par("cex.axis"))
    }
    if (show.grid.labels) {
      if (show.grid.labels%%2) {
        ypos <- grid.pos - radial.lim[1]
        xpos <- rep(0, length(grid.pos))
        if (show.grid.labels == 1) 
          ypos <- -ypos
      }
      else {
        xpos <- grid.pos - radial.lim[1]
        ypos <- rep(0, length(grid.pos))
        if (show.grid.labels == 2) 
          xpos <- -xpos
      }
      if (is.null(radial.labels)) 
        radial.labels = as.character(grid.pos)
      if (!is.null(grid.unit)) 
        radial.labels[length(grid.pos)] <- paste(radial.labels[length(grid.pos)], 
                                                 grid.unit)
      if (boxed.radial) 
        boxed.labels(xpos, ypos, radial.labels, border = FALSE, 
                     cex = par("cex.lab"))
      else text(xpos, ypos, radial.labels, cex = par("cex.lab"))
    }
  }
  invisible(oldpar)
}




###########
#  MAIN   #
###########

args = commandArgs(TRUE)
pMIC = args[1]

pMIC = "/home/aborrel/fluoroquinolones/MIC_currated.csv"
prout = "/home/aborrel/fluoroquinolones/results/MIC/"

dMIC = read.csv(pMIC, header = TRUE)
print (dMIC)
dMIC = dMIC[order(dMIC[,2]),]


svg(paste(prout, "MIC.svg", sep = ""), 35, 5)
plot(seq(1,dim(dMIC)[1]), dMIC[,2], ylim=c(min(dMIC[,-1]),max(dMIC[,-1])), type = "l", ylab="MIC", xlab="", las = 2, axis = "n")
axis(1, seq(1, length(dMIC[,1])), labels = dMIC[,1], las = 2,  cex.axis = 0.3)
lines(dMIC[,2],type="l",col="black")
lines(dMIC[,3],type="l",col="red")
lines(dMIC[,4],type="l",col="blue")
lines(dMIC[,5],type="l",col="green")
legend("topright", legend = colnames(dMIC[,-1]), col = c("black", "red", "blue", "green"), pch = 19)
grid()
dev.off()

svg(paste(prout, "log10.svg", sep = ""), 35, 5)
plot(seq(1,dim(dMIC)[1]), -log10(dMIC[,2]), ylim=c(min(-log10(dMIC[,-1])),max(-log10(dMIC[,-1]))), type = "n", ylab="MIC", xlab="", axis = "n")
axis(1, seq(1, length(dMIC[,1])), labels = dMIC[,1], las = 2, cex.axis = 0.3)
lines(-log10(dMIC[,2]),type="l",col="black")
lines(-log10(dMIC[,3]),type="l",col="red")
lines(-log10(dMIC[,4]),type="l",col="blue")
lines(-log10(dMIC[,5]),type="l",col="green")
legend("topright", legend = colnames(dMIC[,-1]), col = c("black", "red", "blue", "green"), pch = 19)
grid()
dev.off()

dlog = dMIC
dlog[,c(2,3,4,5)] = -log10(dlog[,c(2,3,4,5)])
rownames(dlog) = dMIC[,1]

dM = apply(dlog[,c(2,3,4,5)], 1, mean)

for(orga in colnames(dMIC)[-1]){
  svg (paste (prout, orga, "_radar.svg", sep = ""), 30, 30)
  par(mar=c(0,0,0,0))
  dplot = abs(dM-dlog[,orga])
  names(dplot) = rownames(dlog)
  dplot = dplot[which(dplot >= 1.5)]
  radial.plot(dplot, labels=names(dplot),rp.type="p",main="", line.col="black", mar=c(25,25,25,25), cex.lab = 0.5, radial.lim=c(0,3))
  
  dev.off()
}

