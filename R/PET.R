
#' Latent Heat, Eq 4.2.1 in David R Maidment, Handbook of Hydrology
#' \code{shud.ic} 
#' @param Temp Air temperature, [C]
#' @return Latent heat, lambda in [MJ/kg]
#' @export
LatentHeat <- function(Temp){
  # Eq 4.2.1 in David R Maidment, Handbook of Hydrology */
  # Temp in [C], lambda in [MJ/kg]*/
  # Latent heat of vaporization */
  return (2.501 - 0.002361 * Temp)
}

#' Bulk Surface Resistance, Eq 4.2.22 in David R Maidment, Handbook of Hydrology
#' \code{BulkSurfaceResistance} 
#' @param lai Leaf Area Index [m2 m-2]
#' @return Bulk Surface Resistance [s m-1]
#' @export
BulkSurfaceResistance <- function( lai){
  # # Eq 4.2.22 in David R Maidment, Handbook of Hydrology */
  # #     return 200. / 60. / lai;  # min/m*/
  # # [s m-1] */
  return (200. / lai)
}

#' Air pressure at the elevation, Allen(1998) Eq(7)
#' \code{PressureElevation} 
#' @param z Elevation [m]
#' @return Air pressure at the elevation [kPa]
#' @export
PressureElevation <- function( z){
  return (101.325 * ((293. - 0.0065 * z) / 293)^5.26 )
  # #Pressure based on Elevation*/
  #   #     P atmospheric pressure [kPa],
  # #     z elevation above sea level [m],
  # # Allen(1998) Eq(7) */
}

#' Psychrometric Constant, Eq 4.2.1 in David R Maidment, Handbook of Hydrology
#' \code{PsychrometricConstant} 
#' @param Pressure Air pressure [kPa]
#' @param lambda Latent heat of vaporized water. [MJ/kg]
#' @return  γ psychrometric constant [kPa °C-1]
#' @export
PsychrometricConstant <- function(Pressure, lambda){
  # # Eq 4.2.1 in David R Maidment, Handbook of Hydrology */
  #   # Presser in [kPa], lambda in [MJ/kg]. */
  return (0.0016286 * Pressure / lambda)
  # #     γ psychrometric constant [kPa °C-1],
}

#' Saturated vapor pressure, Eq 4.2.2 in David R Maidment, Handbook of Hydrology
#' \code{VaporPressure_Sat} 
#' @param Temp Air temperature, [C]
#' @return Saturated vapor pressure [kPa]
#' @export
VaporPressure_Sat <- function(Temp){
  # # Eq 4.2.2 in David R Maidment, Handbook of Hydrology*/
  return (0.6108 * exp( 17.27 * Temp / (Temp + 237.3)))
}

#' Actual vapor pressure 
#' \code{VaporPressure_Act} 
#' @param esat Saturated vapor pressure
#' @param rh Relative Humidity [-] 0~1
#' @return 
#' @export
VaporPressure_Act  <- function(esat,  rh){
  return (esat * rh)
}

#' Aerodynamic Resistance. Eq 4.2.25 in David R Maidment, Handbook of Hydrology
#' \code{AerodynamicResistance} 
#' @param Uz Wind speed at height z [m s-1]
#' @param hc Height of crop [m]
#' @param Z_u height of wind measurements [m]
#' @param Z_e height of humidity measurements [m]
#' @param VON_KARMAN von Karman's constant. Default = 0.41 [-]
#' @return Aerodynamic resistance [s m-1]
#' @export 
#' #Page 4.13, Aerodynamic Resistance. in David R Maidment, Handbook of Hydrology
#' # Grass(Hc=0.1m, ra = 45 s/m), 
#' # Agriculture crop(Hc=1m, ra = 18 s/m) and
#' # tree (Hc=10m, ra = 6.5 s/m)
#' hc=cbind(0.1, 1, 10)
#' uz=5
#' zu=2
#' ze = 2
#' AerodynamicResistance(uz, hc, zu, ze)
AerodynamicResistance <- function( Uz,  hc,  Z_u,  Z_e, 
                                   VON_KARMAN = 0.4){
  # Allen, R. G., S, P. L., Raes, D., & Martin, S. (1998).
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
  # uz wind speed at height z [m s-1].
  #  Or:
  #     Eq 4.2.25 in David R Maidment, Handbook of Hydrology
  # #     r_a = 12 * 4.72 * log(Ele[i].WindHeight / rl) / (0.54 * Vel / UNIT_C / 60 + 1) / UNIT_C / 60;    return r_a;
  d = 0.67 * hc
  Z_om = 0.123 * hc;
  Z_ov = 0.0123 * hc;
  r_a = log( abs(Z_u - d) / Z_om ) * log( abs(Z_e - d) / (Z_ov)) / (VON_KARMAN * VON_KARMAN* Uz);
  return(r_a)
}
# hc=cbind(0, 0.1, 1, 10)
# uz=5
# zu=2
# ze = 2
# lai = seq(0.01, 10, 0.1)
# hc = LAI2hc(lai)
# ra =AerodynamicResistance(uz, hc, zu, ze)
# rs =BulkSurfaceResistance(lai)
# par(mfrow=c(2, 1))
# plot(lai,hc)
# matplot(type='l', hc, cbind(ra,rs), log='y')
# stop()

#' Air Density, Eq 4.2.4 in David R Maidment, Handbook of Hydrology
#' \code{shud.ic} 
#' @param P Air pressure [kP]
#' @param Temp Air temperature [C]
#' @return Density of Air. [kg m-3]
#' @export
AirDensity  <- function( P,  Temp){
  #  Eq 4.2.4 in David R Maidment, Handbook of Hydrology
  #  P  in [kP]
  # Temp in [C]
  # rho  -- Density of Air. kg/m3
  return (3.486*P / (275. + Temp))
}

#' Slope of saturated vapor pressure, Eq 4.2.3 in David R Maidment, Handbook of Hydrology
#' \code{SlopeSatVaporPressure} 
#' @param Temp Tempearture [C]
#' @param ES 
#' @return 
#' @export
SlopeSatVaporPressure <- function( Temp,  ES){
  tt = (Temp + 237.3);
  delta = 4098. * ES /  ( tt * tt);
  return (delta)
}


#' Canopy resistance
#' \code{CanopyResistance} 
#' @param rmin Minimum Resistance [s m-1]
#' @return lai Leaf Area Index (LAI) [m2 m-2]
#' @export
CanopyResistance <- function( rmin,  lai){
  rc = rmin * 2. / lai;
  return (rc)
}


#' Calculate vegetation height from Leaf Area Index (LAI), Eq 4.2.23 in David R Maidment, Handbook of Hydrology
#' \code{LAI2hc} 
#' @param lai Leaf Area Index (LAI) [m2 m-2]
#' @return Height of crop [m]
#' @export
#' @examples 
#' lai = seq(0, 10, length.out=100)
#' hc = LAI2hc(lai)
#' # lai=5.5+log(hc)*1.5
#' plot(lai, hc)
LAI2hc <- function(lai){
  hc = exp( (lai - 5.5)/1.5)
  # idx = which(lai < 1)
  # hc[idx] = lai[idx]
  hc
}
#' Potential Evapotranspiration with Pennmann-Monteith equation.
#' \code{PET_PM} 
#' @param Wind Windspeed [m s-1]
#' @param Temp Air temperature, [C]
#' @param RH Relative Humidity [-] 0~1.
#' @param RadNet Net Solar Radiation [w m-2]
#' @param Press Pressure [kPa]
#' @param WindHeight Height that wind was measured [m]
#' @param Veg_height Height of vegetation [m]
#' @param albedo Albedo [-]
#' @param Res_surf Aerodynamic surface resistance [s m-1]
#' @param Elevation Elevation [m]
#' @return PET at [mm m-2 s-1] or [kg m-2 s-1]
#' @export
PET_PM <- function(Wind, Temp, RH, RadNet, Press, 
                   WindHeight=10, 
                   Veg_height=0.12, 
                   albedo= 0.23, 
                   Res_surf=70, 
                   Elevation = NULL){
  # FA056:
  # Veg_height = 0.12 By defining the reference crop as a hypothetical crop 
  # with an assumed height of 0.12 m having a surface resistance of 70 s m-1 and
  # an albedo of 0.23, closely resembling the evaporation of an extension surface 
  # of green grass of uniform height, actively growing and adequately watered, 
  # the FAO Penman-Monteith method was developed.
  eqn <- function(Press, Rad, rho,
                  ed, Delta, r_a, r_s,
                  Gamma, lambda){
    # http:#www.fao.org/docrep/X0490E/x0490e06.htm#penman%20monteith%20equation
    # Rn net radiation at the crop surface [MJ m-2 min-1],
    # G soil heat flux density [MJ m-2 min-1],
    # es saturation vapour pressure [kPa],
    # ea actual vapour pressure [kPa],
    # es - ea saturation vapour pressure deficit [kPa],
    # ∆ slope vapour pressure curve [kPa °C-1],
    # γ psychrometric constant [kPa °C-1].
    #       Gamma;  # psychrometric constant
    #       Delta;  # the slope of the saturation vapour pressure temperature relationship
    U2 = WindSpeed2m(Wind, WindHeight)
    # r_a = 208/U2
    Cp = 1.013e-3    # cp specific heat at constant pressure, 1.013E-3 [MJ kg-1 °C-1] Allen(1998) eq(8) 
    gm = Gamma * (1+0.33 * U2)
    Rad_MJ =  Rad*1e-6 * (1-albedo)
    method1 <- function(){
      ETp = (Delta * Rad  / lambda + rho * Cp * ed / r_a / lambda) /
        (Delta + Gamma * (1 + r_s / r_a)); # eq 4.2.27 [ ]
      x = ETp / lambda; # [ kg/s = 0.001 m3/m2/s=0.001 m/s ]
      message('\n x1: ')
      message('Delta * Rad  / lambda\t', Delta * Rad  / lambda)
      message('rho * Cp * ed / r_a  / lambda\t' , rho * Cp * ed / r_a/ lambda)
      message('Delta + Gamma * (1 + r_s / r_a)  \t' ,Delta + Gamma * (1 + r_s / r_a) )
      
      x
    }
    method2 <- function(){
      tsec = 24*3600
      Rad = RadNet * 24*3600 * 1e-6  # [W m-2] to [MJ m-2 day-1]  24hour daylight
      ETp1 = Delta*(Rad)  / lambda / (Delta + Gamma)
      ETp2 = Gamma * (6.43 * (1+0.536 * U2) * ed) /(Delta + Gamma)
      x = (ETp1 + ETp2)
      # x2 = x2* tsec #4.2.30
      message('\n x2: ')
      message('Etp1\t',ETp1, '\tETP2', ETp2)
      message(' Delta*(Rad)  / lambda  \t' ,  Delta*(Rad)  / lambda )
      message('(Delta + Gamma) \t' , (Delta + Gamma))
      x
    }

    method3 <- function(){
      x = (Delta*Rad / lambda + Gamma * (900 / (Temp + 275)) * ed * U2 ) / (Delta + gm)  # 4.2.31
      # x3=x3 * tsec
      message('\n x3: ')
      message(' x3 (Delta*Rad / lambda + Gamma * (900 / (Temp + 275)) * ed * U2 ) \t' ,
              (Delta*Rad / lambda + Gamma * (900 / (Temp + 275)) * ed * U2 )  )
      message(' (Delta + gm)\t' , (Delta + gm))
      x
    }

    method4 <-function(){
      x = (0.408*Delta*Rad + Gamma*900/(T+273)*U2 * ed ) / (Delta + gm)  # FAO (6)
      # x4 = x4 * tsec
      message('\n x4: ')
      message('(0.408*Delta*Rad /lambda + Gamma*900/(T+273)*U2 * ed )  \t' ,
              (0.408*Delta*Rad /lambda + Gamma*900/(T+273)*U2 * ed ) )
      message(' (Delta + gm)\t' , (Delta + gm))
      x
    }
    x1=method1()
    x2=method2()
    x3=method3()
    x4=method4()
    print(c(x1, x2, x3, x4))
    
    return (ETp)
  }
  CONST_RC = 0.01
  CONST_HC = 0.01
  #define CONST_RH 0.01  #0.01 is the minimum value for Relative Humidity. [m]
  #define CONST_HC 0.01  #0.01 is the minimum Height of CROP. [m]
  if(is.null(Press)){
    Press = PressureElevation(Elevation)
  }
  lambda = LatentHeat(Temp);                      #  eq 4.2.1  [MJ/kg]
  Gamma = PsychrometricConstant(Press, lambda); #  eq 4.2.28  [kPa C-1]
  es = VaporPressure_Sat(Temp);                   #  eq 4.2.2 [kpa]
  ea = es * RH;   #  [kPa]
  ed = es - ea ;  #  [kPa]
  Delta = SlopeSatVaporPressure(Temp, es);        #  eq 4.2.3 [kPa C-1]
  rho = AirDensity(Press, Temp);    #  eq 4.2.4 [kg m-3]
  rs = AerodynamicResistance(Wind, CONST_HC, WindHeight, 2.);     #  eq 4.2.25  [s m-1]
  ra = AerodynamicResistance(Wind, Veg_height, WindHeight, 2.);     #  eq 4.2.25  [s m-1]
  RG = RadNet * 0.9;  # R - G in the PM equation.*/
  etp = eqn(Press, RG, rho, ed, Delta, ra, rs, Gamma, lambda);  # [m/s]
  # xx=cbind(Press, RG, rho, ed, Delta, ra, rs, Gamma, lambda, etp);
  # colnames(xx) = c('Press', 'RG', 'rho', 'ed', 'Delta', 'ra', 'rs', 'Gamma', 'lambda', 'ETP_mm-day');
  # View(xx)
  # print(xx[1, ])
  return(etp)
}

#' Calculate windspeed at Zx heights.
#' \code{WindProfile} 
#' @param zx Desired height [m]
#' @param Um Windspeed at measured height [m s-1]
#' @param zm Measured height [m s-1]
#' @param d Zero plane [m]
#' @param roughness roughness [m]
#' @return windspeed at any desired heights.
#' @export
WindProfile <- function( zx, Um,zm = 10,  
                         d = 0.12,  roughness=0.012){
  # /* Equation 2.3 in [@abtew2012evaporation] */
  if(any((zx - d) < 0) ){
    message('Warning: ZERO plane (', d, ' m) is below desired height')
    print(zx)
  }
  if(any((zx - roughness) < 0 )){
    message('Warning: roughness (', roughness, ' m) is below desired height')
    print(roughness)
  }
  if(any((d - roughness) < 0 )){
    message('Warning: ZERO plane (', d, ' m) is below roughness (', roughness, ' m)')
    print(roughness)
  }
    return (Um * log( (zx - d) / roughness ) / log ((zm - d) /  roughness))
}
# zm = seq(0.2, 20, 0.01)
# U2 = WindProfile(2, Um = 2, zm = zm, d=0.12, roughness = 0.06)
# plot(U2, zm)
# head(cbind(U2, zm))
# # PET_PM(Wind=1, Temp = 25, RH = 0.5, RadNet = 400, Press = pp, WindHeight = 10)*86400
# 
# WindProfile(2, Um = 1.8, zm = 10, d=0, roughness = 0.06)

x=seq(0, 11, 0.01)
hc = LAI2hc(x)
plot(x, hc)
