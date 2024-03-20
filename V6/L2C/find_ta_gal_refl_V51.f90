!   V5.1
!   TM 02/25/2023
!   Added the residual DTA gal refl for V5.1

    subroutine find_ta_gal_refl(ta_gal_refl_tab,gallat,gallon,wind, ta_gal_refl)
    use external_files_l2c_module
    
    implicit none
      
      
    real(4)         ::  ta_gal_refl_tab(3,5),gallat,gallon,wind,  ta_gal_refl(2)
    real(4), save   ::  dgalta_smooth(1440,720,2,4)
    real(4), save   ::  dgalta_V51   (1440,720,2)   ! IQ basis. no wind speed dependence.
 	
  	integer(4),save ::  istart=1
  	integer(4)      ::  iwin1,iwin2,jwin1,jwin2,ilatg,ilong,ipol
	real(4)         ::  brief,a1,a2,b1,b2,ta1,ta2,dta1,dta2,xta_res,y_ta_gal_refl

      
    if(istart.eq.1) then 
    istart=0
    call openbig(3,gal_adj_map_file,'old')
   	read(3) dgalta_smooth
 	close(3)
 	endif
 	
 	call openbig(3,gal_adj_map_file_2,'old') ! new adjustment file for V6.0
   	read(3) dgalta_V51
 	close(3)
 	endif

    brief=(wind+2)/5.  !2m/s add to wind
	if(brief.gt.3.999) brief=3.999
	iwin1=1+brief
	iwin2=iwin1+1
	a1=iwin1-brief
	a2=1-a1


	if(wind.lt.5) then
    brief=(wind-1.5)/3.5
    if(brief.lt.0) brief=0
	jwin1=1+brief  
	jwin2=jwin1+1
	b1=jwin1-brief  
	b2=1-b1
	else
    brief=wind/5.  
	if(brief.gt.2.999) brief=2.999
	jwin1=1+brief
	jwin2=jwin1+1
	b1=jwin1-brief
	b2=1-b1
	endif
	
    ilatg=1+int(4*(gallat+90.))
	ilong=1+int(4* gallon)
	if(ilatg.lt.1 .or. ilatg.gt. 720) then
	write(*,*) gallat,ilatg
	stop 'error in ilatg in '
    endif
	
	if(ilong.lt.1 .or. ilong.gt.1440) then
	write(*,*) gallon,ilong
	stop 'error in ilong'
	endif

    do ipol=1,2
      
	if(ipol.eq.1) then
	ta1= 0.5*(ta_gal_refl_tab(1,iwin1) + ta_gal_refl_tab(2,iwin1))
	ta2= 0.5*(ta_gal_refl_tab(1,iwin2) + ta_gal_refl_tab(2,iwin2))
	xta_res = 0.5*(dgalta_V51(ilong,ilatg,1) + dgalta_V51(ilong,ilatg,2))   ! (I+Q)/2
	else
	ta1= 0.5*(ta_gal_refl_tab(1,iwin1) - ta_gal_refl_tab(2,iwin1))
	ta2= 0.5*(ta_gal_refl_tab(1,iwin2) - ta_gal_refl_tab(2,iwin2))
	xta_res = 0.5*(dgalta_V51(ilong,ilatg,1) - dgalta_V51(ilong,ilatg,2))	! (I-Q)/2
	endif
	
    dta1=dgalta_smooth(ilong,ilatg,ipol,jwin1)
    dta2=dgalta_smooth(ilong,ilatg,ipol,jwin2)

    y_ta_gal_refl=a1*ta1 + a2*ta2 + b1*dta1 + b2*dta2 + xta_res
    if (y_ta_gal_refl <0.0) y_ta_gal_refl=0.0 ! this avoids that all the correction add up to a total ta gal ref that is < 0
    ta_gal_refl(ipol) = y_ta_gal_refl
    enddo  !ipol
      
    return
    end subroutine find_ta_gal_refl
 