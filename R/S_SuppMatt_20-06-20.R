### Richness map for supplementary material ###

cell.r <- cellFromXY(r, coords[row.names(W),])


values_cell <- rep(NA, ncell(r))
names(values_cell) <- 1:ncell(r)
W.cells <- 1:ncell(r) %in% cell.r
values_cell[W.cells] <- rowSums(W)

r.W <- raster::setValues(r, values = values_cell)

windows(5,5)
plot(r.W, xlab = "Longitude", ylab = "Latitude")
plot(costline, add=T)

#### model plot for supplementary material ####
####Models for supplementary material 
curve.nlslrc = nlsLM(NRI ~ a+b/(latitude),
                     start=list(a=(max(NRI)-min(NRI)),
                                b=-min(NRI)))
