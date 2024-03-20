  ! Richard Lindsley + Thomas Meissner 
  ! RSS
  ! sun glint flag for SMAP V5.1
   
  ! Sunglint QC flag
  !
  ! sunglt: the sun glint angle in degrees
  ! winspd: the wind speed in m/s
  ! alpha:  look angle in deg. 
  ! 0 = forward, +90 = left of forward (inot sun), 
  ! +180 = aft, +270 = right of forward (Aquarius configuration)
  !
  ! output = 1 if the QC flag for sun intrusion should be set   = 0 otherwise
  
  
subroutine compute_sun_qc(sunglt, winspd, alpha,   iflag_sun)
implicit none
    
    real(4), intent(in)         :: sunglt, winspd, alpha
    integer(4), intent(out)     :: iflag_sun

    real(4)                     :: thresh_wind
    
    integer(4), parameter       :: nalpha=18
    real(4), parameter          :: balpha=20.0
    ! tie points
    integer(4), parameter, dimension(0:nalpha) :: g0_tab = &
    (/130., 115., 100., 90., 80., 80., 90., 100., 115., 130., 140., 155., 160., 160., 160., 160., 155., 140., 130./)
    integer(4), parameter, dimension(0:nalpha) :: h0_tab = &
    (/50., 40., 30., 20., 20., 20., 20., 20., 30., 50., 70., 80., 80., 90., 90., 90., 80., 70., 50./)
    integer(4) :: ialpha, i1, i2
    real(4)    :: g1, g2, h1, h2, brief, gmax, hmin, xalpha, disc1, disc2
    ! elliptical shapes
    real(4), parameter :: ax=25.0, ay=4.0, bx=40.0, by=9.0


    ! default:
    iflag_sun=0
    if (sunglt<0) then   ! sun eclipsed
        return
    endif

    ! 1. Richard's V4.0 flag. usid in V5.0
    ! Create a customized sunglint QC flag: it's a polynomial
    ! function of sun glint angle between 30 and 50 degrees, where
    ! if the wind speed threshold is larger than the function, the
    ! mask applies. Greater than 50 degrees, it's never masked, and
    ! between 0 and 30 degrees it's always masked.

    if (sunglt < 30.0) then
       iflag_sun = 1
       return
    endif

    ! When the sun glint is between 30 and 50 degrees, compare the
    ! wind speed against the threshold level
    thresh_wind = 1.0 / 8000.0 * (sunglt - 30.0)**4

    if (winspd > thresh_wind) then
        iflag_sun = 1
        return
    endif

    ! 2. Tighten for low wind speeds
    ! This seems necessary as there is some sunglint sneaking in at high latitudes 
    if (winspd <= 5.0 .and. sunglt <= 45.) then
        iflag_sun = 1
    endif 

    ! 3. Low winds and large glint angles. sun enters sidelobes. Cut out elliptical shapes in (sunglt,wspd) plane
    xalpha=alpha

    if (xalpha .lt. 0.) return ! sun ray does not intersect the Earth
    if (xalpha .gt. 359.999) xalpha=359.999

    ialpha= floor(xalpha/balpha)

    if (ialpha < 0 .or. ialpha >= nalpha) then 
        write(*,*),alpha,xalpha,ialpha, ' ialpha OOB. stop'
        stop
    endif

    i1 = ialpha
    i2 = i1 + 1

    ! center of the ellipse 
    g1 = g0_tab(i1)
    g2 = g0_tab(i2)

    brief = (xalpha - i1*balpha)/balpha
    gmax = (1.0-brief)*g1 + brief*g2          

    ! Ellipse around (gmax,0) with half axes 25,4
    disc1 = ((sunglt-gmax)/ax)**2 + (winspd/ay)**2
    if (disc1 .le. 1.0) then ! lies wihtin ellipse
        iflag_sun = 1
        return
    endif
    
    ! 4. high wind speeds. possible sun scatter. Cut out elliptical shapes in (sunglt,wspd) plane
    ! center of the ellipse 
    h1 = h0_tab(i1)
    h2 = h0_tab(i2)

    brief = (xalpha - i1*balpha)/balpha
    hmin = (1.0-brief)*h1 + brief*h2          

    ! Ellipse around (hmin,20) with half axes 40,9
    disc2 = ((sunglt-hmin)/bx)**2 + ((winspd-20.)/by)**2
    if (winspd .le. 20. .and. disc2 .le. 1.0) then ! lies wihtin ellipse
        iflag_sun = 1
        return
    endif
    if (winspd .gt. 20. .and. sunglt .lt. hmin+bx) then 
    ! discard everything with sun angle below hmin+40 if wind speed is above 20 m/s
        iflag_sun = 1
        return
    endif
    
    
return
end subroutine compute_sun_qc

