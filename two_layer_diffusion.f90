program two_layer_diffusion
!
!------------------------------------------------------------------------- 
!           Practice 2-layer Sream Temperature  DIFFUSION model
!                 Ryan Niemeyer, January 2016
!     based on Strzepek et al. (2015) equations 6 & 7, which is based
!      off of Chapra 1997
!
!-------------------------------------------------------------------------

implicit none

!-------------------------------------------------------------------------
!    define variables
!-------------------------------------------------------------------------

! 1-dim. real array for 10 years
real, dimension(365*10) :: flow_in,flow_out,flow_eout,flow_hout, flow_Tin
real, dimension(365*10) ::  temp_epil,temp_hypo, temp_out_tot
real, dimension(365*10) :: temp_change_ep, temp_change_hyp, energy
real, dimension(365*10) :: energy_tot, diffusion_tot, T_in_tot, T_out_tot

REAL, PARAMETER :: Pi = 3.1415927, prcnt_flow_epil = 0.2, prcnt_flow_hypo=0.8
real, parameter :: v_t = 2.1 !diffusion coefficient (m/day)
!    diffusion coefficient - based on Snodgrass, 1974

CHARACTER*49 :: path='/raid3/rniemeyr/practice/practice_fortran/output/'
real  :: flow_constant
integer  :: i
real :: x1, x2, x3

real  :: depth_total, depth_e, depth_h, width, length, volume_e_x, outflow_x
real :: energy_x, volume_h_x, area, density, heat_c, temp_change, delta_t
real  :: flow_in_hyp_x, flow_in_epi_x, flow_out_epi_x, flow_out_hyp_x
real  :: epix, hypox, dif_epi_x, dif_hyp_x, x, flow_epi_x, flow_hyp_x

! CHARACTER(*), PARAMETER :: fileplace = "/raid3/rniemeyr/practice/practice_fortran/output/" !output file

open(unit=10, file="temp_change_epilim.dat")

depth_total = 30
width = 200
length = 17000
area = width*length
delta_t = 1 ! time is days,  assumes all units in equations are in days

!-------------------------------------------------------------------------
!     generate flow and energey temporal variables
!-------------------------------------------------------------------------

! --------------------------- flow -------------------------

! --------- constant flow paramaeter -------
! constant  to change day of
!  flow to go "up and down" with sin wave
! with 0 and 365 being "0" point
flow_constant = 365/(2*Pi)

! generates flow eacy day as a sin wave with the peak flow on April 1
! at 90000 cfs, and lowest point is 30000 cfs on October 1
! days are calendar year (day = 1 = January 1)


do  i=0,10*365

   ! ------ get flow in to the reservoir for the year ------
   flow_in(i) =  sin( i/flow_constant)
   flow_in(i) = (flow_in(i) + 2)*30000   ! gets flow to vary from 30000 to 90000 cfs
   flow_in(i) = (flow_in(i)/35.315)*60*60*24  ! converts ft3/sec to m3/day
   ! ------ get flow in to the reservoir for the year ------
   flow_out(i) =  sin( i/flow_constant)
   flow_out(i) = (flow_out(i) + 2)*30000   ! gets flow to vary from 30000 to 90000
   flow_out(i) = (flow_out(i)/35.315)*60*60*24   ! converts ft3/sec to m3/sec

   ! NOTE: don't need these two flows - tried to vary these individual, but
   ! decided not to.

   ! ------ get flow out of reservoir for the year ------
   flow_hout(i) =  sin( i/flow_constant)
   flow_hout(i) = (flow_hout(i) + 2)*24000   ! gets flow to vary from 24000 to 72000
   flow_hout(i) = (flow_hout(i)/35.315)*60*60*24   ! converts cfs to m3/sec

   ! ------ get flow out of reservoir for the year ------
   flow_eout(i) =  sin( i/flow_constant)
   flow_eout(i) = (flow_eout(i) + 2)*6000  ! gets outflow to vary from 6000 to 18000
   flow_eout(i) = (flow_eout(i)/35.315)*60*60*24   ! converts cfs to m3/sec

end  do

! print *, flow(366:456)

! --------------------------- flow - temperature   -------------------------

do  i=1,10*365
   flow_Tin(i)  =  cos((i/flow_constant)+ Pi)
   flow_Tin(i) = (flow_Tin(i) + 1.5)*10   ! gets temperature to vary form 5 to 25C
end  do

! print *, "Flow temperature from day 1 to 200"
! print *, flow_Tin(1:200)

! print *, "Flow temperature from day 3000 to 3200"
! print *, flow_Tin(3000:3200)

! --------------------------- energy  -------------------------

! --------- constant flow paramaeter -------
!  energy in W/m2
do  i=0,10*365
   energy(i) =  cos((i/flow_constant)+ Pi)
  ! gets net energy (positive downward) to vary from:
  ! " +0.7)*120":  -36 to 204 W/m2
  ! " +0.5)*120":  -60 to 180 W/m2
  ! units are  W/m2 or Joules/m2 * sec
   energy(i) =( (energy(i) + 0.5)*120)*60*60*24 !converts to Joules/m2 * day 
end  do

! print *, "energy from 2000 to 2100"
! print *, energy(2000:2100)

!-------------------------------------------------------------------------
!              Calculate two-layer stream temperature   
!-------------------------------------------------------------------------

density = 1000 !  density of water in kg / m3
heat_c = 4180  !  heat capacity of water in joules/ kg * C

! initial variables
depth_total = 20
depth_e = 5
depth_h = 15
volume_e_x = area*depth_e
volume_h_x = area*depth_h
temp_epil(1) = 5 ! starting epilimnion temperature at 5 C
temp_hypo(1) = 5 ! starting hypolimnion temperature at 5 C

! start at 2, because need to have an epil and hypo temperature to start out
do  i=2,10*365

  ! calculate incoming net energy (J/m2*day) to epilimnion
  energy_x = energy(i)*area

  ! divide incoming flow into reservoir to epilimnion and hypolimnion
  flow_in_hyp_x = flow_in(i)*prcnt_flow_epil
  flow_in_epi_x = flow_in(i)*prcnt_flow_hypo

  flow_out_hyp_x = flow_out(i)*prcnt_flow_epil
  flow_out_epi_x = flow_out(i)*prcnt_flow_hypo

  ! calculate temperature change due to diffusion
  ! NOTE: don't need to multiply by heat capacity or density of water because
  !       energy component is divided by those two
  dif_epi_x  = v_t * area * density * heat_c * (temp_hypo(i-1) - temp_epil(i-1))
  dif_hyp_x  = v_t * area * density * heat_c * (temp_epil(i-1) - temp_hypo(i-1))

  ! calculate change in EPILIMNION  temperature (celsius)
  flow_epi_x = flow_in_epi_x * density * heat_c * flow_Tin(i)
  flow_epi_x = flow_epi_x - (flow_out_epi_x * density * heat_c *  temp_epil(i-1))
  temp_change_ep(i) = flow_epi_x + energy_x + dif_epi_x
  temp_change_ep(i) = temp_change_ep(i)/(volume_e_x * density * heat_c)
  temp_change_ep(i) = temp_change_ep(i) * delta_t

  ! save each temperature change component
  energy_tot(i) = energy_x
  diffusion_tot(i) = dif_epi_x
  T_in_tot(i) = flow_in_epi_x*flow_Tin(i)/volume_e_x
  T_out_tot(i) = flow_out_epi_x*temp_epil(i-1)/volume_e_x
  ! update epilimnion volume for next time step
  volume_e_x = volume_e_x + (flow_in_epi_x - flow_out_epi_x)
  temp_epil(i) = temp_epil(i-1) +  temp_change_ep(i)

  ! calculate change in HYPOLIMNION  temperature (celsius)
  flow_hyp_x  = flow_in_hyp_x * density * heat_c * flow_Tin(i)
  flow_hyp_x = flow_hyp_x - (flow_out_hyp_x * density * heat_c * temp_hypo(i-1))
  temp_change_hyp(i) = flow_hyp_x  +  dif_hyp_x  !add horizontal advection and  diffusion 
  temp_change_hyp(i) = temp_change_hyp(i)/(volume_h_x * density * heat_c)
  temp_change_hyp(i) = temp_change_hyp(i) * delta_t

  ! update hypolimnion volume for next time step
  volume_h_x = volume_h_x + (flow_in_hyp_x - flow_out_hyp_x)
  temp_hypo(i) = temp_hypo(i-1) +  temp_change_hyp(i)

  ! calculate combined (hypo. and epil.) temperature of outflow
  outflow_x = flow_out_epi_x + flow_out_hyp_x
  epix = temp_epil(i)*(flow_out_epi_x/outflow_x)  ! portion of temperature from epilim. 
  hypox= temp_hypo(i)*(flow_out_hyp_x/outflow_x)  ! portion of temperature from hypol.
  temp_out_tot(i) = epix + hypox

 if (i==2 .or. i==10  .or. i==20 .or. i==30 .or. i==40 .or. i==3649) then
  print *, "run: ", i
!  print *,"energy of incoming radiation " , energy(i)/(60*60*24)
!  print *, "energy - joules/day*m2", energy(i)
!  print *, "area:                 ", area
!  print *, "numeriator  of energy equation: ", x1
!  print *, "heat capacity  ", heat_c
!  print *, "density   ", density
!  print *, "volume of epilimnion: ", volume_e_x
!  print *, "denominator of energy equation: ", x2
!  print *, "volume of epilimnion: ", volume_e_x
!   print *, "change in volume - epilim.: ", flow_in_epi_x -  flow_out_epi_x
!   print *, "change in volume - hypolim.: ", flow_in_hyp_x - flow_out_hyp_x
!   print *, "depth of epilimnion: ", volume_e_x/area
!   print *, "depth of hypolimnion: ", volume_h_x/area
  print *, "energy temp change is : ",energy_x/(volume_e_x * density * heat_c)
  print *, "diffus temp change  is: ", dif_epi_x/(volume_e_x * density *heat_c)
  print *, "Qin temp change: ", (flow_in_epi_x * flow_Tin(i))/volume_e_x
!  print *, "Qout temp change: ", (flow_out_epi_x * temp_epil(i-1))/volume_e_x
!  print *, "flow epilim. temp change: ", flow_epi_x/(volume_e_x*density* heat_c)
!   print *, "flow hypolim.  temp change is: ", flow_hyp_x
!  print *, "depth of epilimnion: ", volume_e_x/area 
  print *, "temperature change in epilimnion:  ", temp_change_ep(i)
   print *, "outflow temperature from epilimnion is: ", temp_epil(i)
!   print *, "depth of hypolimnion: ", volume_h_x/area
!  print *, "temperature change in hypolimnion:  ", temp_change_hyp(i)
!  print *, "outflow temperature from hypolimnion is: ", temp_hypo(i)
!   print *, "outflow (combined)  temperature is: ", temp_out_tot(i)
  print *, "  "
 end if

end  do
 
!-------------------------------------------------------------------------
!         Calculations to print at the end 
!-------------------------------------------------------------------------

 open(unit=30, file=path//"hypo_temp_change.txt",action="write",status = "replace")
 open(unit=31, file=path//"epil_temp_change.txt",action="write",status ="replace")
 open(unit=40, file=path//"temperature_epil.txt",action="write",status ="replace")
 open(unit=41, file=path//"temperature_hypo.txt",action="write",status ="replace")
 
 write(30,*),temp_change_hyp
 write(31,*),temp_change_ep
 write(40,*),temp_epil
 write(41,*),temp_hypo

 close(30)
 close(31)
 close(40)
 close(41)

 print *,"total simulation temperature change in epilim.: ", sum(temp_change_ep(1:3650))
 print *,"total simulation  temperature change in hypolim.: ",  sum(temp_change_hyp(1:3650))

! print *, "sum of 1:720 temp change due to flow in: ", sum(T_in_tot(1:720))
! print *, "sum of 1:720 temp change due to flow out: ", sum(T_out_tot(1:720))
! print *, "sum of 1:720 temp change due to energy: ", sum(energy_tot(1:720))
! print *, "sum of 1:720 temp change due to diffusion: ", sum(diffusion_tot(1:720))


end program two_layer_diffusion
