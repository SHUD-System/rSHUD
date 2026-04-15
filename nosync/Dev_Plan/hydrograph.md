In R, write a function plotHG()to plot a hydrograph with ggplot2. 
The dimension of input data x is Nr * Nc, Nr is the number of rows, and Nc is the number of columns. The x is the type of xts, so the index(x) is time format. The first column of x is precipitation P, and the rest columns are discharge Q. 
requirement of plotHG():
1. The axis is time format. 
2. Plot the P with tile. Plot the Qs as lines.
3. the P uses the right y-axis, and the Q uses the left y-axis. 
4. The right y-axis is up-down flipped.  Please show the right Y-axis reversely.

The test R code is following:

nt=100
x = as.xts(matrix(rnorm(3*nt), nrow=nt, ncol=3), as.Date('2000-01-01')+1:nt)
plotHG(x)





I have a two-column data (Rainfall P and discharge Q). please plot the data in R with ggplot2.
1. plot rainfall with geom_col, along with the right side Y-axis.  The right side y-axis is reverse. 
2. plot discharge Q with geom_line, along with the left side Y-axis.

See the data. 
Rainfall_P = runif(10, 0, 1500)
Discharge_Q = runif(10, 0, 500)
data=cbind(Rainfall_P, Discharge_Q)

Please reverse the RIGHT SIDE y-axis.  Please make sure the geom_col() also plot reversely. 