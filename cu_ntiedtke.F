!=================================================================================================================
 module cu_ntiedtke_common
 use ccpp_kinds,only: kind_phys
 use mpas_log


 implicit none
 save


 real(kind=kind_phys),parameter:: t13   = 1.0/3.0
 real(kind=kind_phys),parameter:: tmelt = 273.16
 real(kind=kind_phys),parameter:: c1es  = 610.78
 real(kind=kind_phys),parameter:: c3les = 17.2693882
 real(kind=kind_phys),parameter:: c3ies = 21.875
 real(kind=kind_phys),parameter:: c4les = 35.86
 real(kind=kind_phys),parameter:: c4ies = 7.66

 real(kind=kind_phys),parameter:: rtwat = tmelt
 real(kind=kind_phys),parameter:: rtber = tmelt-5.
 real(kind=kind_phys),parameter:: rtice = tmelt-23.

 integer,parameter:: momtrans = 2
 real(kind=kind_phys),parameter:: entrdd  = 2.0e-4
 real(kind=kind_phys),parameter:: cmfcmax = 1.0
 real(kind=kind_phys),parameter:: cmfcmin = 1.e-10
 real(kind=kind_phys),parameter:: cmfdeps = 0.30
 real(kind=kind_phys),parameter:: zdnoprc = 2.0e4
 real(kind=kind_phys),parameter:: cprcon  = 1.4e-3
 real(kind=kind_phys),parameter:: pgcoef  = 0.7

 real(kind=kind_phys):: rcpd
 real(kind=kind_phys):: c2es
 real(kind=kind_phys):: c5les
 real(kind=kind_phys):: c5ies
 real(kind=kind_phys):: r5alvcp
 real(kind=kind_phys):: r5alscp
 real(kind=kind_phys):: ralvdcp
 real(kind=kind_phys):: ralsdcp
 real(kind=kind_phys):: ralfdcp
 real(kind=kind_phys):: vtmpc1
 real(kind=kind_phys):: zrg

 logical,parameter:: nonequil = .true.
 logical,parameter:: lmfpen   = .true.
 logical,parameter:: lmfmid   = .true.
 logical,parameter:: lmfscv   = .true.
 logical,parameter:: lmfdd    = .true.
 logical,parameter:: lmfdudv  = .true.


!=================================================================================================================
 end module cu_ntiedtke_common
!=================================================================================================================

 module cu_ntiedtke
 use ccpp_kinds,only: kind_phys
 use cu_ntiedtke_common


 implicit none
 private
 public:: cu_ntiedtke_init,  &
          cu_ntiedtke_final


 contains


!=================================================================================================================
 subroutine cu_ntiedtke_init(con_cp,con_rd,con_rv,con_xlv,con_xls,con_xlf,con_grav)
!=================================================================================================================

!input arguments:
 real(kind=kind_phys),intent(in):: &
    con_cp,  &
    con_rd,  &
    con_rv,  &
    con_xlv, &
    con_xls, &
    con_xlf, &
    con_grav

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine cu_ntiedtke_init:')

 rcpd    = 1.0/con_cp
 c2es    = c1es*con_rd/con_rv
 c5les   = c3les*(tmelt-c4les)
 c5ies   = c3ies*(tmelt-c4ies)
 r5alvcp = c5les*con_xlv*rcpd
 r5alscp = c5ies*con_xls*rcpd
 ralvdcp = con_xlv*rcpd
 ralsdcp = con_xls*rcpd
 ralfdcp = con_xlf*rcpd
 vtmpc1  = con_rv/con_rd-1.0
 zrg     = 1.0/con_grav

 call mpas_log_write('--- end subroutine cu_ntiedtke_init:')

 end subroutine cu_ntiedtke_init

!=================================================================================================================
 subroutine cu_ntiedtke_final(errmsg,errflg)
!=================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg

!-----------------------------------------------------------------------------------------------------------------

 errmsg = 'cu_ntiedtke_final OK'
 errflg = 0

 end subroutine cu_ntiedtke_final

!=================================================================================================================
 end module cu_ntiedtke
!=================================================================================================================

