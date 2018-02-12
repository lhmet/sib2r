!
! Routine 'sib2_offline_lib' to be used as dynamical library from R package
!
! To compile:
!
! =======================
! THIS WAY DOES NOT WORK!
! =======================
! export PKG_FCFLAGS='-freal-4-real-8'
! R CMD SHLIB -o sib2_offline_lib.so comsibc_h.f90 pardif_h.f90 cal2jul.f90 
! derive_trans_new.f aero.f90 momopt_stockli.f90 SIB2sub.f90 sib2_offline_lib.f90
!
! ==============
! THIS WAY WORK!
! ==============
! MAKEFLAGS="FCFLAGS=-ffixed-line-length-none\ -freal-4-real-8\ -O2 FFLAGS=-ffixed-line-length-none 
! -freal-4-real-8\ -O2 FLIBS= LIBR0= LIBR1= LIBR=" 
! R CMD SHLIB -o sib2_offline_lib.so !comsibc_h.f90 pardif_h.f90 cal2jul.f90 derive_trans_new.f aero.f90 momopt_stockli.f90 !SIB2sub.f90 sib2_offline_lib.f90
!
! ==============
! THIS WAY WORK!
! ==============
! gfortran -ffixed-line-length-none -freal-4-real-8 -O2 -Wall -fPIC -c comsibc_h.f90 pardif_h.f90 cal2jul.f90 derive_trans_new.f aero.f90 momopt_stockli.f90 SIB2sub.f90 sib2_offline_lib.f90
! gfortran -shared -Wl,-z,relro -o sib2_offline_mylib.so comsibc_h.o pardif_h.o cal2jul.o derive_trans_new.o aero.o momopt_stockli.o SIB2sub.o sib2_offline_lib.o
!
! Last revision 20180208 by nelsonvn
!

!subroutine sib2_offline_lib(&
!      !SIB_RUN
!      n1, infile2, n2, dir_out, n3, id, &                            !06
!      date1, date2, dtt, itrunk, ilw, nlai, aero_parms_method, &     !13
!      !SIB_INVARS
!      ivtype, istype, isnow, ipbl, idirr, decay, zlat, &             !20
!      poros0, phsat0, satco0, bee0, slope0, slpp0, tc0, tg0, td0, &  !29
!      capac0, snoww0, www0, gwdep, app, bpp, cpp, &                  !36
!      !SIB_CAEROD
!      zwind0, zmet0, &                                               !38
!      !SIB_VDERIV
!      gmudmu, &                                                      !39
!      !MORPHO1
!      zs, dsfc, g10, ztz00, &                                        !43
!      !PHYSIO1
!      shti0, slti0, trda0, trdm0, trop0, btheta0, &                  !49
!      !SIB_PARMS
!      z10, z20, zc, chil0, leafw, leafl, vcover0, rootd0, phc0, &    !58
!      tranlv, tranln, trandv, trandn, reflv, refln, refdv, refdn, &  !66
!      sorefv, sorefn, effcon0, gradm0, binter0, respcp0, &           !72
!      atheta0, hlti0, hhti0, vmax0, sodep0)                          !77

subroutine sib2_offline_lib(n1, infile2, n2, dir_out, n3, id, &
      vrun, vini, vsoil, vmorpho, vphysio, vparms)

   use const   !pi, g, asnow, clai, cpair, cw, epsfac, kappa, rhoair, 
               !snomel, stefan, tf, vkc, rcp, timcon, psy, hlat, snofac
   use vdyijt   !zlt, green, fparc
   use caerod, only : ha, g2, g3, corb1, corb2, zwind, zmet, g1, ztz0
   use vderiv, only : z0d, dd, cc1, cc2
   use soils, only : poros, phsat, satco, bee, slope, slpp
   use stepv   !tc, tg, td, capac, snoww, www
   use atmos, only : swdown, rnetm, em, tm, um, sunang
   use soilij   !sodep, soref
   use vstate   !z1, z2, vcover, chil, tran, ref, rootd, phc, effcon, gradm, 
                !binter, respcp, atheta, btheta, trda, trdm, trop, 
                !slti, hlti, shti, hhti

   implicit none

!   !Dummy arguments
!   integer, intent(in) :: n1, n2, n3   !'SIB_RUN'
!   character (len = n1) :: infile2
!   character (len = n2) :: dir_out
!   character (len = n3) :: id
!   integer, intent(in) :: date1, date2, itrunk, ilw, nlai, aero_parms_method
!   real, intent(in) :: dtt
!   integer, intent(in) :: ivtype, istype, isnow, ipbl, idirr   !'SIB_INVARS'
!   real, intent(in) :: decay, zlat, poros0, phsat0, satco0, bee0, slope0
!   real, intent(in) :: slpp0, tc0, tg0, td0, capac0(2), snoww0(2), www0(3)
!   real, intent(in) :: gwdep, app, bpp, cpp
!
!   real, intent(in) :: zwind0, zmet0   !SIB_CAEROD
!   real, intent(in) :: gmudmu   !SIB_VDERIV
!   real, intent(in) :: zs, dsfc, g10, ztz00   !MORPHO1
!   real, intent(in) :: shti0, slti0, trda0, trdm0, trop0, btheta0   !PHYSIO1
!
!   real, intent(in) :: z10, z20, zc, chil0, leafw, leafl, vcover0, rootd0, phc0   !SIB_PARMS
!   real, intent(in) :: tranlv, tranln, trandv, trandn, reflv, refln, refdv, refdn
!   real, intent(in) :: sorefv, sorefn, effcon0, gradm0, binter0, respcp0
!   real, intent(in) :: atheta0, hlti0, hhti0, vmax0, sodep0

   !Dummy arguments
   integer, intent(in) :: n1, n2, n3
   character (len = n1) :: infile2
   character (len = n2) :: dir_out
   character (len = n3) :: id
   integer, intent(in) :: vrun(11)
   real, intent(in) :: vini(19)
   real, intent(in) :: vsoil(6)
   real, intent(in) :: vmorpho(6)
   real, intent(in) :: vphysio(6)
   real, intent(in) :: vparms(27)

   !Local variables (NEW)
   integer :: date1, date2, itrunk, ilw, nlai, aero_parms_method
   real :: dtt
   integer :: ivtype, istype, isnow, ipbl, idirr
   real :: decay, zlat
   real :: gwdep, app, bpp, cpp
   real :: gmudmu
   real :: zs, dsfc
   real :: zc, leafw, leafl
   real :: tranlv, tranln, trandv, trandn, reflv, refln, refdv, refdn
   real :: sorefv, sorefn, vmax0



   !Local parameters (OLD)
   integer, parameter :: npar = 9
   integer, parameter :: in2 = 10
   integer, parameter :: out1 = 20
   integer, parameter :: ichmet = 7
   integer, parameter :: iu = 8
   integer, parameter :: nchar_max = 200
   integer, parameter :: icho1 = 21
   integer, parameter :: icho2 = 22
   integer, parameter :: iout1 = 38
   real, parameter :: decmax = 23.5 * (pi / 180.0)



   !Local variables (OLD)
   !'govern' module
   integer :: year, month, day, hour
   real :: sols, realday, time, season, dec, sindec, cosdec, tandec, coshr
   real :: totwb   !'checkbal' module
   real :: rsnow, tsnow   !'snow' module
   real :: salb(2,2)   !'site' module
   real :: tprec, zlwd   !'readin' module
   
   integer :: nymd, jul, dayloc, sta_id
   integer :: iter, i, ierr

   real :: lat, sinlat, coslat

   real, allocatable :: aero_parms(:,:)   ! Taken off from 'sib2river_inc' module
   integer, external :: cal2jul



   namelist /sib_run/ infile2, dir_out, id, date1, date2, nlai, &
      aero_parms_method, isnow, ipbl, ilw, itrunk, ivtype, istype, idirr
   namelist /sib_ini/ zlat, zwind, zmet, dtt, tc, tg, td, gwdep, &
      gmudmu, app, bpp, cpp, capac, snoww, www
   namelist /sib_soil/ bee, decay, poros, phsat, satco, slpp
   namelist /sib_morpho/ zs, dsfc, g1, ztz0, sodep, slope
   namelist /sib_physio/ shti, slti, trda, trdm, trop, btheta
   namelist /sib_parms/ z1, z2, zc, chil, leafw, leafl, vcover, rootd, phc, &
      tranlv, tranln, trandv, trandn, reflv, refln, refdv, refdn, &
      sorefv, sorefn, effcon, gradm, binter, respcp, atheta, hlti, hhti, vmax0



!   !Copy arguments values to modules variables (OLD)
!   !SIB_INVARS
!   poros = poros0
!   phsat = phsat0
!   satco = satco0
!   bee = bee0
!   slope = slope0
!   slpp = slpp0
!   tc = tc0
!   tg = tg0
!   td = td0
!   capac = capac0
!   snoww = snoww0
!   www = www0
!
!   !SIB_CAEROD
!   zwind = zwind0
!   zmet = zmet0
!
!   !MORPHO1
!   g1 = g10
!   ztz0 = ztz00
!
!   !PHYSIO1
!   shti = shti0
!   slti = slti0
!   trda = trda0
!   trdm = trdm0
!   trop = trop0
!   btheta = btheta0
!
!   !SIB_PARMS
!   z1 = z10
!   z2 = z20
!   chil = chil0
!   vcover = vcover0
!   rootd = rootd0
!   phc = phc0
!   effcon = effcon0
!   gradm = gradm0
!   binter = binter0
!   respcp = respcp0
!   atheta = atheta0
!   hlti = hlti0
!   hhti = hhti0
!   sodep = sodep0

   !Copy arguments values to modules variables (NEW)
   date1 = vrun(1)
   date2 = vrun(2)
   nlai = vrun(3)
   aero_parms_method = vrun(4)
   isnow = vrun(5)
   ipbl = vrun(6)
   ilw = vrun(7)
   itrunk = vrun(8)
   ivtype = vrun(9)
   istype = vrun(10)
   idirr = vrun(11)

   zlat = vini(1)
   zwind = vini(2)
   zmet = vini(3)
   dtt = vini(4)
   tc = vini(5)
   tg = vini(6)
   td = vini(7)
   gwdep = vini(8)
   gmudmu = vini(9)
   app = vini(10)
   bpp = vini(11)
   cpp = vini(12)
   capac = vini(13:14)
   snoww = vini(15:16)
   www = vini(17:19)

   bee = vsoil(1)
   decay = vsoil(2)
   poros = vsoil(3)
   phsat = vsoil(4)
   satco = vsoil(5)
   slpp = vsoil(6)

   zs = vmorpho(1)
   dsfc = vmorpho(2)
   g1 = vmorpho(3)
   ztz0 = vmorpho(4)
   sodep = vmorpho(5)
   slope = vmorpho(6)

   shti = vphysio(1)
   slti = vphysio(2)
   trda = vphysio(3)
   trdm = vphysio(4)
   trop = vphysio(5)
   btheta = vphysio(6)

   z1 = vparms(1)
   z2 = vparms(2)
   zc = vparms(3)
   chil = vparms(4)
   leafw = vparms(5)
   leafl = vparms(6)
   vcover = vparms(7)
   rootd = vparms(8)
   phc = vparms(9)
   tranlv = vparms(10)
   tranln = vparms(11)
   trandv = vparms(12)
   trandn = vparms(13)
   reflv = vparms(14)
   refln = vparms(15)
   refdv = vparms(16)
   refdn = vparms(17)
   sorefv = vparms(18)
   sorefn = vparms(19)
   effcon = vparms(20)
   gradm = vparms(21)
   binter = vparms(22)
   respcp = vparms(23)
   atheta = vparms(24)
   hlti = vparms(25)
   hhti = vparms(26)
   vmax0 = vparms(27)



   !Write input data on file
  ! open(unit=51, file='config_file_sib2.txt', status='unknown')
  ! write(51,nml=sib_run)
  ! write(51,nml=sib_ini)
  ! write(51,nml=sib_soil)
  ! write(51,nml=sib_morpho)
  ! write(51,nml=sib_physio)
  ! write(51, nml=sib_parms)
  ! close(unit=51, status='keep')



   !Set values from input data (OLD)
   lat = zlat * (pi / 180.0)   !To avoid RADC2_2 routine. 20180105
   sinlat = sin(lat)
   coslat = cos(lat)
   tran(1,1) = tranlv
   tran(2,1) = tranln
   tran(1,2) = trandv
   tran(2,2) = trandn
   ref(1,1) = reflv
   ref(2,1) = refln
   ref(1,2) = refdv
   ref(2,2) = refdn
   soref(1) = sorefv
   soref(2) = sorefn



   allocate( aero_parms(npar,nlai) )

   !Computes SiB2 aeordynamical parameters (ZLT means LAI)
   do i = 1, nlai
      zlt = 0.1 * real(i)
      if (aero_parms_method == 1) then   !NVN. 20180201
         call derive_trans_new(z1, z2, zc, chil, leafw, leafl, vcover, zlt, &
            zs, g1, ztz0, ha, z0d, dd, g2, g3, cc1, cc2, corb1, corb2)
      else
         call momopt_stockli(z1, z2, zc, chil, leafw, leafl, vcover, zlt, &
            zs, zmet, zwind, ha, z0d, dd, g2, g3, cc1, cc2, corb1, corb2)
      end if

      aero_parms(1,i) = ha
      aero_parms(2,i) = z0d
      aero_parms(3,i) = dd
      aero_parms(4,i) = g2
      aero_parms(5,i) = g3
      aero_parms(6,i) = cc1
      aero_parms(7,i) = cc2
      aero_parms(8,i) = corb1
      aero_parms(9,i) = corb2
   end do

   ! Computes the location on the year of the initial date
   year = date1 / 1000000
   month = mod(date1, 1000000) / 10000
   day = mod(date1, 10000) / 100
   jul = cal2jul(year, month, day)
   dayloc = jul - cal2jul(year, 1, 1) + 1
   sols = (4141.0 / 24.0) + 0.25 * real(mod(year + 3, 4))

   call varcal(dsfc, sodep, rootd, decay, satco)   !NVN. Moved from main loop on 20180201.

   open(unit=in2, file=trim(infile2), status='old')
   read(in2,*)
   open(unit=iout1, file=trim(dir_out)//'sib2diag'//trim(id)//'.txt', status='unknown')

   !
   ! Main loop
   !
   iter = 0
   do
!      read(in2,*,iostat=ierr) sta_id, nymd, swdown, rnetm, em, tm, um, tprec, zlt, green, zlwd
      read(in2,*,iostat=ierr) sta_id, nymd, swdown, rnetm, em, tm, um, tprec!, zlt, green, zlwd
      if (ierr /= 0 .or. nymd > date2) exit

      if (nymd >= date1) then
         iter = iter + 1

         ! Extract YEAR, MONTH, DAY and HOUR from NYMD variable
         year = nymd / 1000000
         month = mod(nymd, 1000000) / 10000
         day = mod(nymd, 10000) / 100

         !!! TIRAR o '-1' DA LINHA ABAIXO POIS O ARQUIVO VAI DE 1 A 24 E DEVERIA IR DE 0 A 23 
         hour = mod(nymd, 100) - 1

         hour = mod(hour, 24)   ! Useful only if file have daily data

         ! Yearly update
         if (month == 1 .and. day == 1 .and. hour == 0) then
            dayloc = 1
            sols = (4141.0 / 24.0) + 0.25 * real(mod(year + 3, 4))
         end if

         realday = real(dayloc) + real(hour) / 24.0
         season = (realday - sols) / 365.2
         dec = decmax * cos(twopi * season)
         sindec = sin(dec)
         cosdec = cos(dec)
         tandec = tan(dec)
         time = real(hour) + 0.5
         coshr = cos(-pi + time / 24.0 * twopi)   ! Melhorar isto
         sunang = sinlat * sindec + coslat * cosdec * coshr   !To avoid RADC2_2 routine. 201800205
         sunang = max(0.01, sunang)

         ! The following three variables must be read from data file
         zlt = 3.54
         green = 0.8
         zlwd = 4.903E-3 / 24.0 * (tm ** 4) * (0.66 + 0.039 * sqrt(em / 10.0)) / 3600.0

         i = nint(10.0 * zlt)
         i = min(nlai, max(i, 1))
         ha      = aero_parms(1,i)   !NVN. Taken from VARCAL. 20180201
         z0d     = aero_parms(2,i)
         dd      = aero_parms(3,i)
         g2      = aero_parms(4,i)
         g3      = aero_parms(5,i)
         cc1     = aero_parms(6,i)
         cc2     = aero_parms(7,i)
         corb1   = aero_parms(8,i)
         corb2   = aero_parms(9,i)


         call driver(tprec, zlwd, isnow, iter)
         call balan(1, tprec, zlwd, totwb, icho2, nymd, dtt, iter, ivtype, istype)
         call inter2(idirr, app, bpp, cpp, dtt)
         call rada2(ilw, salb)
         call begtem(rsnow, tsnow)
         call endtem(ipbl, ztz0, itrunk, dtt, vmax0, gmudmu, rsnow, tsnow)
         call updat2(idirr, gwdep, dtt)
         call balan(2, tprec, zlwd, totwb, icho2, nymd, dtt, iter, ivtype, istype)
         call outer(iout1, nymd, tprec, dtt, iter, vmax0, salb, ivtype, istype)

         if (hour == 23) dayloc = dayloc + 1
      end if   ! (nymd>=date1)
   end do   ! main loop

   deallocate( aero_parms )
   close(unit=in2, status='keep')
   close(unit=iout1, status='keep')

   end subroutine sib2_offline_lib

