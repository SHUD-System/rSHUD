Atmosphericpressure <- function(z){
  # Atmospheric pressure (P)
  # http://www.fao.org/docrep/X0490E/x0490e07.htm
  P = 101.3 * ( (293 - 0.0065 * z) / 293 ) ^ (5.26)
  return(P)
}
LatentHeatVaporization<- function(K, pressure, r_a){
  # Latent heat of vaporization (l)
  # http://www.fao.org/docrep/X0490E/x0490e07.htm
  CC = 2.45e6 # 2.45MJ/kg
  return(CC)
}
SlopeSatVaporPressure <- function(){
  # Slope of saturation vapour pressure curve (D )
  # http://www.fao.org/docrep/X0490E/x0490e07.htm
  # D slope of saturation vapour pressure curve at air temperature T [kPa °C-1],
  # T air temperature [°C],
  # exp[..] 2.7183 (base of natural logarithm) raised to the power [..].
  delta = 4098 * ( 0.6108* exp(17.27 * T / (T + 237.3))  ) /  ( (T + 237.3)^2 )
  return(delta)
}


AerodynamicResistance <- function(Uz, hc, Z_u, Z_e=2.0, VON_KARMAN= 0.41){
  # /* Allen, R. G., S, P. L., Raes, D., & Martin, S. (1998).
  # Crop evapotranspiration : Guidelines for computing crop water requirements
  # by Richard G. Allen ... [et al.].
  # FAO irrigation and drainage paper: 56.
  # Eq(4) in reference
  # ra aerodynamic resistance [s m-1],
  # zm height of wind measurements [m],
  # zh height of humidity measurements [m],
  # d zero plane displacement height [m],
  # zom roughness length governing momentum transfer [m],
  # zoh roughness length governing transfer of heat and vapour [m],
  # k von Karman's constant, 0.41 [-],
  #       uz wind speed at height z [m s-1].
  #    Or:
  #       Eq 4.2.25 in David R Maidment, Handbook of Hydrology
  #    */
  #   //    r_a = 12 * 4.72 * log(Ele[i].windH / rl) / (0.54 * Vel / UNIT_C / 60 + 1) / UNIT_C / 60;    return r_a;
  d = 0.67 * hc;
  Z_om = 0.123 * hc;
  Z_ov = 0.0123 * hc;
  r_a = log( (Z_u - d) / Z_om ) * log( (Z_e - d) / (Z_ov)) /
    (VON_KARMAN * VON_KARMAN* Uz); 
  return(r_a);
}
BulkSurfaceResistance <-function(R_ref, lai){
  # // R_ref     bulk stomatal resistance of the well-illuminated leaf [min m-1],
  return (R_ref * 2. / lai)
  # /* Allen(1998) */
}
PressureElevation <- function( z){
  return (101.325 *((293. - 0.0065 * z) / 293)^5.26)
  # /*Pressure based on Elevation*/
  #   //    P atmospheric pressure [kPa],
  # //    z elevation above sea level [m],
  # /* Allen(1998) Eq(7) */
}

VaporPressure_Sat = function(T_in_C){
  # /* Eq 4.2.2 in David R Maidment, Handbook of Hydrology*/
    return (0.6108 * exp( 17.27 * T_in_C / (T_in_C + 237.3)) )
}
AirDensity = function( P,  T){
  # /* Eq 4.2.4 in David R Maidment, Handbook of Hydrology*/
  #   /*  P  in [kP]
  # T in [C]
  # rho  -- Density of Air. kg/m3
  # */
    return (3.486 * P / (275 + T))
}
SlopeSatVaporPressure <- function( T,  ES){
   tt = (T + 237.3);
   delta = 4098. * ES /  ( tt * tt);
  return (delta)
}
CanopyResistance <- function( rmin,  lai){
  rc = rmin * 2. / lai;
  return (rc)
}

# v=1:100
# plot(v, AerodynamicResistance(Uz=v, hc=0.01,  Z_u=2))
