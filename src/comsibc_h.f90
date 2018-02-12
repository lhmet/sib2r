!
! This module replaces 'COMSIBC.H' file
!

module stepv
   implicit none
   real :: tc
   real :: tg
   real :: td
   real :: capac(2)
   real :: snoww(2)
   real :: www(3)
end module stepv



module const
   implicit none
   real, parameter :: pi = 3.14159265358979323846
   real, parameter :: twopi = 6.28318530717958647692
   real, parameter :: g = 9.81
   real, parameter :: asnow = 13.2
   real, parameter :: clai = 4.2 * 1000.0 * 0.2
   real, parameter :: cpair = 1010.0
   real, parameter :: cw = 4.2 * 1000.0 * 1000.0
   real, parameter :: epsfac = 0.622
   real, parameter :: kappa = 0.286
   real, parameter :: rhoair = 1.225
   real, parameter :: snomel = 370518.5 * 1000.0
   real, parameter :: stefan = 5.669 * 10e-9
   real, parameter :: tf = 273.16
   real, parameter :: vkc = 0.41
   real, parameter :: rcp = rhoair * cpair
   real, parameter :: timcon = pi / 86400.0
   real :: psy
   real :: hlat
   real :: snofac
end module 



module  atchem
   implicit none
   real, parameter :: po2m = 20900.0
   real, parameter :: pco2m = 34.0
end module 



!THIS MODULE IS NOT MORE NECESSARY
!module gridij
!   implicit none
!   integer :: ivtype 
!   integer :: istype 
!end module 



module vstate
   implicit none
   real :: z2
   real :: z1
   real :: vcover
   real :: chil
   real :: tran(2,2)
   real :: ref(2,2)
!!!cxx      rootd, ph1, ph2, &
   real :: rootd
   real :: phc
   real :: effcon
   real :: gradm
   real :: binter
   real :: respcp
   real :: atheta
   real :: btheta
!!!cxx      trda, trdm, trop, tpsa, tpsb, tpsc
   real :: trda 
   real :: trdm
   real :: trop
   real :: slti
   real :: hlti
   real :: shti
   real :: hhti
end module 



module vdyijt
   implicit none
   real :: zlt
   real :: green
   real :: fparc
end module 



module vderiv
   implicit none
   real :: z0d
   real :: dd
   real :: cc1
   real :: cc2
   real :: vmax0
   real :: gmudmu
end module 



module soilij
   implicit none
   real :: sodep
   real :: soref(2)
end module 



module soils
   implicit none
   real :: bee
   real :: phsat
   real :: poros
   real :: satco
   real :: slope
   real :: zdepth(3)   ! tatsch, 6 jun 2011
   real :: slpp        ! tatsch, 6 jun 2011
   real :: satcoz(3)   ! tatsch 15 set 2011
   real :: satcoz_a    ! tatsch 15 set 2011
   real :: zlay(3)     ! tatsch 15 set 2011
   real :: zlay_a      ! tatsch 15 set 2011
   real :: decay       ! tatsch 15 set 2011
!!!c
!!!c     satcoz = satco * exp(-decay*zlay)
!!!c     satcoz_a = sum(satcoz(i, i =1,3))/3
!!!c     zlay ! Mid-level of soil layer (m)
!!!c     zlay_a! Mid-level of soil layer (m)
end module 



module atmos
   implicit none
   real, parameter :: bps = 1.0
   real, parameter :: psur = 1000.0
   real :: em
   real :: tm
   real :: um
   real :: zm
   real :: ppc
   real :: ppl
   real :: radn(3,2)
   real :: sunang
   real :: swdown
   real :: rnetm
   real :: cloud
end module 



module caerod
   implicit none
   real :: corb1
   real :: corb2
   real :: ha
   real :: g1
   real :: g2
   real :: g3
   real :: ztz0
   real :: zwind
   real :: zmet
   real :: ht
end module 



!THIS MODULE IS NOT MORE NECESSARY
!module site
!   implicit none
!   real :: zlong
!   real :: zlat
!   real :: salb(2,2)
!   real :: rab(2,3,2)
!end module 



!THIS MODULE IS NOT MORE NECESSARY
!module steps
!   implicit none
!   integer :: itrunk 
!   integer :: ilw 
!   integer :: niter 
!   integer :: iter 
!!   integer :: ispare   ! NOT USED
!   real :: dtt
!end module 



!THIS MODULE IS NOT MORE NECESSARY
!!! Note that, by default, TIME, YEAR, DAY and HOUR are REAL and MONTH is INTEGER
!!! Variables MONTH and HOUR are NEVER USED
!module govern
!   implicit none
!   integer :: year
!   integer :: month 
!   integer :: day
!   integer :: hour
!!   real :: year
!!   real :: day
!!   real :: hour
!   real :: time
!   real :: realday
!   real :: sols
!   real :: season
!   real :: dec
!   real :: sindec
!   real :: cosdec
!   real :: tandec
!   real :: coshr
!end module 



module donor
   implicit none
   real :: etmass
   real :: hflux
   real :: roff
   real :: zlwup
   real :: drag
end module 



module rause
   implicit none
   real :: z0
   real :: d
   real :: rbc
   real :: rdc
end module 



module aerorx
   implicit none
   real :: ra
   real :: rb
   real :: rd
   real :: rbbest              ! adicionado por tatsch
end module 



module grads
   implicit none
   real :: tgs
   real :: ta
   real :: ea
   real :: etc
   real :: etgs
   real :: getc
   real :: getgs
   real :: u2
   real :: ustar
end module 



module radabs
   implicit none
   real :: albedo(2,2,2)
   real :: radfac(2,2,2)
   real :: radt(2)
   real :: thermk
   real :: exrain
   real :: tgeff
end module 



module surfrs
   implicit none
   real :: rst
   real :: rstfac(4)
   real :: rsoil
   real :: cog1
   real :: cog2
   real :: hr
   real :: fc
   real :: fg
end module 



module hydrol
   implicit none
   real :: satcap(2)
   real :: wc
   real :: wg
   real :: canex
   real :: areas
end module 



module stores
   implicit none
   real :: ccx
   real :: cg
   real :: csoil
end module 



module delts
   implicit none
   real :: dtc
   real :: dtg
   real :: dtd
   real :: dth
   real :: dqm
end module 



module carbio
   implicit none
   real :: assimn
   real :: respc
   real :: respg
   real :: pco2i
   real :: gsh2o
end module 



module flux
   implicit none
   real :: ec
   real :: eg
   real :: hc
   real :: hg
   real :: chf
   real :: shf
   real :: gflux   ! adicionado por tatsch
   real :: ect
   real :: eci
   real :: egi
   real :: egs
   real :: ecmass
   real :: egmass
   real :: heaten
end module 



!THIS MODULE IS NOT MORE NECESSARY
!module snow
!   implicit none
!   real :: tsnow
!   real :: rsnow
!end module 



!THIS MODULE IS NOT MORE NECESSARY
!module initialv
!   implicit none
!   real :: tc_ini
!   real :: tg_ini
!   real :: td_ini
!   real :: www_ini(3)
!end module 



!THIS MODULE IS NOT MORE NECESSARY
!!! Are these variables integer or real ones?
!module readin
!   implicit none
!   integer :: iqcaws   !NOT USED
!   integer :: iqchyd   !NOT USED
!   integer :: mevap   !NOT USED
!   integer :: msensh   !NOT USED
!   integer :: mustar
!   real :: zlwd
!   real :: tprec
!end module 



!THIS MODULE IS NOT MORE NECESSARY
!module checkbal
!   implicit none
!   real :: totwb
!end module 



module output
   implicit none
   real :: roff1
   real :: roff2
   real :: roff3
   real :: roff4
   real :: gwsoil
   real :: otest1
   real :: otest2
   real :: otest3
   real :: roffGA
   real :: roffq3g
   real :: roffq3g2   ! tatsch 11 set 2011
end module 



!THIS MODULE IS NOT MORE NECESSARY
!module temp
!   implicit none
!   integer :: item01
!   integer :: item02
!   real :: GWdep
!end module 



!THIS MODULE IS NOT MORE NECESSARY
!module suroff
!   implicit none
!   real :: surdep
!   real :: finfil
!end module 



!THIS MODULE IS NOT MORE NECESSARY
!module preccoff
!   implicit none
!   real :: app
!   real :: bpp
!   real :: cpp
!end module 



!THIS MODULE IS NOT MORE NECESSARY
!module irrgrid
!   implicit none
!   integer :: idirr
!end module 



