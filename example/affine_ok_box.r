# read data file as character vector
data   <- scan(file="affine_ok_box.out", what="character")
#
# set some line types
no_line  <- 0
solid    <- 1
dashed   <- 2
dotted   <- 3
dotdash  <- 4
longdash <- 5
twodast  <- 6
#
# minimum allowable value for x2(t)
x2_min <- -1
#
# convert to a character matrix
len    <- length(data)
nc     <- 5
#
# separate from the header line
header <- data[1 : nc]
data   <- data[(nc+1) : len]
#
# convert data to double precision
data   <- as.double(data)
#
# convert to a matrix
nr     <- (len - 1)  / nc
data   <- matrix(data, nrow=nr, ncol=nc, byrow=TRUE)
#
# plot the data
circle <- 21
t      <- data[, "t" == header]
x2_z   <- data[, "z1" == header]
plot(
	x = t, 
	y = x2_z, 
	type = "p", 
	main = "Affine Kalman-Bucy Smoother",
	xlab = "t", 
	ylab = "x2", 
	lty=no_line
)
#
# plot true values with dotted line
x2_true <- data[, "x2_true" == header]
lines(x = t, y = x2_true, type = "l", lty=dotted)
#
# plot lower and upper bounds with a straight line
x1_bnd   <- c(0, 2*pi)
x2_low   <- c(x2_min, x2_min)
x2_up    <- c(2 + x2_min, 2 + x2_min)
lines(x = x1_bnd, y = x2_low, type = "l", lty=solid)
lines(x = x1_bnd, y = x2_up,  type = "l", lty=solid)
#
# plot constrained estimate with longdash line
x2_con   <- data[, "x2_con" == header]
lines(x = t, y = x2_con, type = "l", lty=longdash)
#
# plot unconstrained estimate with dashed line
x2_free <- data[, "x2_free" == header]
lines(x = t, y = x2_free, type = "l", lty=dashed)
#
# add legend to plot
yleg    <- 2.4
xleg    <- .1
no_sym  <- -1
legend (
	x = xleg , 
	y = yleg, 
	legend = c( "meas",  "true",     "con",    "free",   "bound"), 
	lty    = c(no_line,   dotted,  longdash,  dashed,     solid),
	pch    = c( circle,  no_sym,    no_sym,    no_sym,    no_sym)
)
#
# save the plot in encapsulated postscript
savePlot(filename = "affine_ok_box", type = "eps");
