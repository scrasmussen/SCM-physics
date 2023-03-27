!=================================================================================================================
 module cu_ntiedtke_common
 use ccpp_kinds,only: kind_phys
 use mpas_log


 implicit none
 save

 real(kind=kind_phys):: alf
 real(kind=kind_phys):: als
 real(kind=kind_phys):: alv
 real(kind=kind_phys):: cpd
 real(kind=kind_phys):: g
 real(kind=kind_phys):: rd
 real(kind=kind_phys):: rv

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
 public:: cu_ntiedtke_init,          &
          cu_ntiedtke_final,         &
          cu_ntiedtke_timestep_init, &
          cu_ntiedtke_timestep_final


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

 alf = con_xlf
 als = con_xls
 alv = con_xlv
 cpd = con_cp
 g   = con_grav
 rd  = con_rd
 rv  = con_rv

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
 subroutine cu_ntiedtke_timestep_init(its,ite,kts,kte,itimestep,grav,dz,pres,presi,t,rho,qv,qc,qi,u,v,w,qvften, &
                        thften,qvftenz,thftenz,prsl,ghtl,tf,qvf,qcf,qif,uf,vf,prsi,ghti,omg,errmsg,errflg)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,kts,kte
 integer,intent(in):: itimestep

 real(kind=kind_phys),intent(in):: grav
 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: dz,pres,t,rho,qv,qc,qi,u,v
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte):: qvften,thften
 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte+1):: presi,w

!--- inout arguments:
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte):: tf,qvf,qcf,qif,uf,vf
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte):: ghtl,omg,prsl
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte):: qvftenz,thftenz
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte+1):: ghti,prsi

!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg

!--- local variables and arrays:
 integer:: i,k,pp,zz

 real(kind=kind_phys),dimension(its:ite,kts:kte):: zl,dot
 real(kind=kind_phys),dimension(its:ite,kts:kte+1):: zi

!-----------------------------------------------------------------------------------------------------------------

 k = kts
 do i = its,ite
    zi(i,k) = 0.
 enddo
 do k = kts,kte
    do i = its,ite
       zi(i,k+1) = zi(i,k)+dz(i,k)
    enddo
 enddo
 do k = kts,kte
    do i = its,ite
       zl(i,k)  = 0.5*(zi(i,k)+zi(i,k+1))
       dot(i,k) = -0.5*grav*rho(i,k)*(w(i,k)+w(i,k+1))
    enddo
 enddo

 pp = 0
 do k = kts,kte+1
    zz = kte + 1 - pp
    do i = its,ite
       ghti(i,zz) = zi(i,k)
       prsi(i,zz) = presi(i,k)
    enddo
    pp = pp + 1
 enddo
 pp = 0
 do k = kts,kte
    zz = kte-pp
    do i = its,ite
       ghtl(i,zz) = zl(i,k)
       omg(i,zz)  = dot(i,k)
       prsl(i,zz) = pres(i,k)
    enddo
    pp = pp + 1
 enddo

 pp = 0
 do k = kts,kte
    zz = kte-pp
    do i = its,ite
       tf(i,zz)  = t(i,k)
       qvf(i,zz) = qv(i,k)
       qcf(i,zz) = qc(i,k)
       qif(i,zz) = qi(i,k)
       uf(i,zz)  = u(i,k)
       vf(i,zz)  = v(i,k)
    enddo
    pp = pp + 1
 enddo

 if(itimestep == 1) then
    do k = kts,kte
       do i = its,ite
          qvftenz(i,k) = 0.
          thftenz(i,k) = 0.
       enddo
    enddo
 else
    pp = 0
    do k = kts,kte
       zz = kte-pp
       do i = its,ite
          qvftenz(i,zz) = qvften(i,k)
          thftenz(i,zz) = thften(i,k)
       enddo
       pp = pp + 1
    enddo
 endif

 errmsg = 'cu_ntiedtke_timestep_init OK'
 errflg = 0

 end subroutine cu_ntiedtke_timestep_init

!=================================================================================================================
 subroutine cu_ntiedtke_timestep_final(its,ite,kts,kte,stepcu,dt,exner,qv,qc,qi,t,u,v,qvf,qcf,qif,tf,uf,vf,rn, &
                              raincv,pratec,rthcuten,rqvcuten,rqccuten,rqicuten,rucuten,rvcuten,errmsg,errflg)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,kts,kte
 integer,intent(in):: stepcu

 real(kind=kind_phys),intent(in):: dt
 real(kind=kind_phys),intent(in),dimension(its:ite):: rn
 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: exner,qv,qc,qi,t,u,v,qvf,qcf,qif,tf,uf,vf

!--- inout arguments:
 real(kind=kind_phys),intent(inout),dimension(its:ite):: raincv,pratec
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte):: rqvcuten,rqccuten,rqicuten
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte):: rthcuten,rucuten,rvcuten

!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg

!--- local variables and arrays:
 integer:: i,k,pp,zz

 real(kind=kind_phys):: delt,rdelt

!-----------------------------------------------------------------------------------------------------------------

 delt  = dt*stepcu
 rdelt = 1./delt

 do i = its,ite
    raincv(i) = rn(i)/stepcu
    pratec(i) = rn(i)/(stepcu*dt)
 enddo

 pp = 0
 do k = kts,kte
    zz = kte - pp
    do i = its,ite
       rthcuten(i,k) = (tf(i,zz)-t(i,k))/exner(i,k)*rdelt
       rqvcuten(i,k) = (qvf(i,zz)-qv(i,k))*rdelt
       rqccuten(i,k) = (qcf(i,zz)-qc(i,k))*rdelt
       rqicuten(i,k) = (qif(i,zz)-qi(i,k))*rdelt
       rucuten(i,k)  = (uf(i,zz)-u(i,k))*rdelt
       rvcuten(i,k)  = (vf(i,zz)-v(i,k))*rdelt
    enddo
    pp = pp + 1
 enddo

 errmsg = 'cu_ntiedtke_timestep_final OK'
 errflg = 0

 end subroutine cu_ntiedtke_timestep_final

!=================================================================================================================
 end module cu_ntiedtke
!=================================================================================================================

