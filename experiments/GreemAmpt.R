#' Calculate infiltration with Green Ampt method
#' \code{GreenAmpt} 
#' @param theta  Soil moisture content
#' @param thetas  Saturated moisture content
#' @param ksat  Saturated hydraulic conductivity
#' @param Hfront Depth of wetting front.
#' @param y0 Depth of ponding water.
#' @param Time Time.
#' @return Infiltration rate 
#' @references Green, W.H. and G. Ampt. 1911. Studies of soil physics, part I – the flow of air and water through soils. J. Ag. Sci. 4:1-24.
#' @export
GreenAmpt <- function(theta, thetas, ksat, Hfront, y0=0.01, Time=0){
  nt=length(Time)
  if(nt <= 1){
    dt = thetas - theta
    dh = Hfront # dh = Infiltration Depth
    dh= max(dh, .0010) # Minimum D
    # if(y0 > ksat) y0 = ksat
    if(dt <= 0){
      q = 0
    }else{
      q = ksat * (y0 + dh) / (dh * dt)
      # d0 = ksat * dh / dt
      # q = dt * sqrt(d0 / 2 / 1)
    }
    if(q > y0){
      q = y0
    }
    ret = rbind(cbind(Time, q, theta, 
                  thetas, y0, ksat, Hfront) )
  }else{
    dtime= c(0, diff(Time))
    ret = NULL
    hf=Hfront
    for(it in 1:nt){
      tt = Time[it]
      x=GreenAmpt(theta = theta, thetas=thetas, 
                  ksat=ksat, y0=y0,  Hfront = hf, Time=tt)
      hf = hf + x[2]
      # print(hf)
      ret= rbind(ret, x)
    }
  }
  colnames(ret) = c('Time', 'Infiltration', 'Theta', 
                    'ThetaS', 'Y0', 'Ksat', 'WetFront')
  return (ret)
}