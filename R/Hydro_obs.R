#' Linefit the observation and simulation
#' \code{LineFit} 
#' @param x Observation data or matrix/data.frame consist of observation and simulation
#' @param y Simulation data
#' @param xlab xlab, default is 'Observation'
#' @param ylab ylab, default is 'Simulation'
#' @param ... more options for stat_poly_eq()
#' @return Triangle mesh of model domain.
#' @export
#' @examples 
#' obs <- rnorm(100)
#' sim <- obs+rnorm(100)/2
#' p <- LineFit(cbind(obs,sim))
#' p
#' if(require(ggpmisc)){
#' formula = y ~ x
#' p+ ggpmisc::stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#' label.x.npc = "right", label.y.npc = 0.15,
#' formula = formula, parse = TRUE, size = 3)
#' }
LineFit <-function(x, y=NULL,
                   xlab='Observation', ylab='Simulation', 
                   ...){
  colnames(x) = c(xlab, ylab)
  if(is.null(y)){
    df = as.data.frame(x)
  }else{
    df = data.frame(x,y)
  } 
  # lim=apply(apply(df, 2, range, na.rm=TRUE), 2, diff)
  # vmin=apply(df, 2, min, na.rm=TRUE)
  # loc= c(vmin[1]+lim[1] *eq.loc,
  #        vmin[2]+lim[2] *eq.loc[2])
  
  # lm_eqn <- function(df){
  #   colnames(df)=c('x','y')
  #   m <- stats::lm(x ~ y, df);
  #   a=as.numeric(round(stats::coef(m)[1], 2))
  #   b=as.numeric(round(stats::coef(m)[2], 2))
  #   r2=as.numeric(round(summary(m)$r.squared, 3))
  #   eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
  #                    list(a = a, 
  #                         b = b, 
  #                         r2 = r2))
  #   as.character(as.expression(eq));                 
  # }
  lim=range(df, na.rm = TRUE)
  p=ggplot2::ggplot(data=df, 
                    ggplot2::aes_string(x=xlab, y=ylab)) + 
    ggplot2::geom_smooth(method = "lm", se = T, col=2) +
    ggplot2::geom_point() + 
    ggplot2::geom_abline() + 
    ggplot2::coord_fixed(ratio=1, expand = TRUE, xlim=lim, ylim=lim)
  
  # p=ggplot2::ggplot(data=df, 
  #                   ggplot2::aes_string(x=xlab, y=ylab)) + 
  #   ggplot2::geom_smooth(method = "lm", se = T, col=2) +
  #   ggplot2::geom_point() + 
  #   ggplot2::geom_abline() + 
  #   ggplot2::coord_fixed(ratio=1)+
  #   ggplot2::geom_text(x = loc[1],
  #           y = loc[2], label = lm_eqn(df),
  #           parse = TRUE)
  p
}
#' Default Comparison of the observation and simulation
#' \code{QvsO} 
#' @param qq Matrix or data.frame. Column 1 is observation, while colum 2 is simulation.
#' @param sim Simulation time-seris data;
#' @param obs Observation time-seris data;
#' @param ... more options for plot()
#' @export
#' @examples 
#' obs=rnorm(100)
#' sim=obs+rnorm(100)/2
#' QvsO(cbind(obs,sim))
QvsO <-function(qq,sim=qq[,2],obs=qq[,1], ...){
  hydroGOF::ggof(sim=sim,obs=obs,col=c( 'blue','red'),...)
}



#' Default Comparison of the observation and simulation
#' \code{fdc} 
#' @param val Matrix or data.frame, number of column=2
#' @return ggplot handler
#' @export
#' @examples 
#' obs=rnorm(100)
#' sim=obs+rnorm(100)/2
#' fdc(cbind(obs, sim))
fdc <- function(val,
                   xlab='Exceedance (%)',
                   ylab='value'){
  vname=deparse(substitute(val))
  if(xts::is.xts(val) | zoo::is.zoo(val)){
    val = zoo::coredata(val)
  }
  if(is.matrix(val) | is.data.frame(val)){
    nx=nrow(val)
    x = 1:nx /(nx-1) *100
    y = apply(val, 2, sort,  decreasing = TRUE)
    if(is.null(colnames(y))){
      colnames(y)=vname 
    }
    md0=data.frame(x,y)
    colnames(md0)=c('x', colnames(y))
  }else{
    nx=length(val)
    x = 1:nx /(nx-1) *100
    y = sort(val, decreasing = TRUE)
    md0=data.frame(x=x,y=y)
  }
  md = reshape2::melt(md0, id=1)
  p=ggplot2::ggplot(md, ggplot2::aes(x=x, y=value, color=variable))+
    ggplot2::geom_line()+
    ggplot2::xlab(xlab) + 
    ggplot2::ylab(ylab)+ 
    ggplot2::theme(legend.position=c(.8, .9),
          legend.direction = 'horizontal',
          legend.title = ggplot2::element_blank())
  p
}