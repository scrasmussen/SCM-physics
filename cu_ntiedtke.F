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
 public:: cu_ntiedtke_final


 contains


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

