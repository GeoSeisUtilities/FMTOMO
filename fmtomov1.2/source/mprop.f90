!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: MODULE
! CODE: FORTRAN 90
! This module declares variable for global use, that is, for
! USE in any subroutine or function or other module.
! Variables whose values are SAVEd can have their most
! recent values reused in any routine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE globalp
IMPLICIT NONE
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10), DIMENSION(:,:,:), ALLOCATABLE :: idep
! idep = interface node depth
END MODULE globalp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program is designed to calculate an estimate of model 
! roughness and variance by dicing up a cubic B-spline grid.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM misfit
USE globalp
IMPLICIT NONE
INTEGER :: i,j,k,l,m,n,i1,j1,k1,itl,nlay,ii
INTEGER :: nvr,nvt,nvp,conr,cont,conp
INTEGER :: ndr,ndt,ndp,nnr,nnt,nnp,str,stt,stp
INTEGER :: nin,nit,nip,icnt,id1,id2,id3
INTEGER :: checkstat
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ivel
REAL(KIND=i10) :: gor,got,gop,gsr,gst,gsp,u,v,mvar,mrough
REAL(KIND=i10) :: rgsr,rgst,rgsp,rdm,sumi,sumj,sumk,earth
REAL(KIND=i10) :: rdm1,sumi1,sumj1,sumk1,dp,dt,dr,ri,risti
REAL(KIND=i10) :: gsti, gspi,goti,gopi,xt,zp,dept,depb
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: ui,vi,wi
REAL(KIND=i10), DIMENSION(:,:,:), ALLOCATABLE :: velsm,velrm,velr
REAL(KIND=i10), PARAMETER :: pi=3.14159265359
CHARACTER (LEN=25) :: rmfile,smfile,ifile
!
! nvr = number of vertices in radial direction
! nvt = number of vertices in theta (N-S) direction
! nvp = number of vertices in phi (E-W) direction
! gor = grid origin in radius
! got = grid origin in theta (N-S)
! gop = grid origin in phi (E-W)
! gsr = grid separation in radius
! gst = grid separation in theta (N-S)
! gsp = grid separation in phi (E-W)
! smfile = solution model file 
! rmfile = reference model file
! velr = velocity at refined grid point
! velsm = velocity at spline vertex for solution model
! velrm = velocity at spline vertex for reference model
! ndr,ndt,ndp = node dicing level in r,theta,phi
! nnr,nnt,nnp = number of diced nodes in r,theta,phi
! u,v = Cubic spline independent variables
! ui,vi,wi = Cubic spline basis functions
! rgsr,rgst,rgsp = Refined node spacing in r,theta,phi
! sumi,sumj,sumk = Summation variables for constructing spline
! conr,cont,conp = Counters for refining grid
! str,stt,stp = Refined grid location in r,theta,phi
! itl = target layer
! mvar = Model variance
! mrough = Model roughness
! sumi1,sumj1,sumk1= Summation variables for reference spline
! dp,dt,dr = Contributions to roughness Laplacian
! ri,risti = Denomenators for difference operator
! earth = Earth radius
! ifile = interface file
! nlay = number of layers
! ivel = switch for layer nodes (0 = outside; 1 = inside layer)
! nin = number of interfaces
! nit,nip = number of interface nodes in latitude and longitude
! gsti, gspi = interface node separation in latitude and longitude
! goti,gopi = grid origin in latitude and longitude
! icnt = counting integer
! xt,zp = diced grid location
! dept,depb = top and bottom bounds of layer
!
OPEN(UNIT=10,FILE='mprop.in',STATUS='old')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a25)')smfile
READ(10,'(a25)')rmfile
READ(10,'(a25)')ifile
READ(10,*)itl
READ(10,*)ndr,ndt,ndp
READ(10,*)earth
CLOSE(10)
!
! Read in B-spline grid of solution model
!
OPEN(UNIT=20,FILE=smfile,STATUS='old')
READ(20,*)nlay
IF(itl.lt.1.or.itl.gt.nlay)then
   write(6,*)'You have requested layer ',itl
   write(6,*)'But there are only ', nlay, 'layers in the model!!!'
   write(6,*)'Terminating prematurely!!!!'
   stop
ENDIF
DO ii=1,nlay
   READ(20,*)nvr,nvt,nvp
   READ(20,*)gsr,gst,gsp
   READ(20,*)gor,got,gop
   nvr=nvr-2
   nvt=nvt-2
   nvp=nvp-2
   IF(ii.eq.itl)THEN
      gor=gor+nvr*gsr-earth
      got=got+gst*(nvt-1)
      got=(pi/2-got)
      gop=gop+gsp
      ALLOCATE(velsm(0:nvr+1,0:nvt+1,0:nvp+1), STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL velv'
      ENDIF
      DO i=0,nvr+1
         DO j=0,nvt+1
            DO k=0,nvp+1
               READ(20,*)velsm(nvr+1-i,nvt+1-j,k)
            ENDDO
         ENDDO
      ENDDO
      EXIT
   ELSE
       DO i=0,nvr+1
         DO j=0,nvt+1
            DO k=0,nvp+1
               READ(20,*)rdm1
            ENDDO
         ENDDO
      ENDDO
   ENDIF
ENDDO
CLOSE(20)
!
! Read in B-spline grid of reference model
!
OPEN(UNIT=20,FILE=rmfile,STATUS='old')
READ(20,*)
DO ii=1,nlay
   READ(20,*)nvr,nvt,nvp
   READ(20,*)gsr,gst,gsp
   READ(20,*)gor,got,gop
   nvr=nvr-2
   nvt=nvt-2
   nvp=nvp-2
   IF(ii.eq.itl)THEN
      gor=gor+nvr*gsr-earth
      got=got+gst*(nvt-1)
      got=(pi/2-got)
      gop=gop+gsp
      ALLOCATE(velrm(0:nvr+1,0:nvt+1,0:nvp+1), STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL velv'
      ENDIF
      DO i=0,nvr+1
         DO j=0,nvt+1
            DO k=0,nvp+1
               READ(20,*)velrm(nvr+1-i,nvt+1-j,k)
            ENDDO
         ENDDO
      ENDDO
      EXIT
   ELSE
       DO i=0,nvr+1
         DO j=0,nvt+1
            DO k=0,nvp+1
               READ(20,*)rdm1
            ENDDO
         ENDDO
      ENDDO
   ENDIF
ENDDO
close(20)
! Calculate total numer of refined nodes in r,theta,phi
! and the refined grid spacing.
!
nnr=(nvr-1)*ndr+1
nnt=(nvt-1)*ndt+1
nnp=(nvp-1)*ndp+1
rgsr=gsr/ndr
rgst=gst/ndt
rgsp=gsp/ndp
ALLOCATE(velr(nnr,nnt,nnp), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL velr'
ENDIF
!
! Now create an array to define nodes that lie within a layer 
!
ALLOCATE(ivel(nnr,nnt,nnp))
ivel=0
!
! Now read in the bounding interfaces
!
open(unit=30,file=ifile,status='old')
read(30,*)nin
read(30,*)nit,nip
read(30,*)gsti,gspi
read(30,*)goti,gopi
nit=nit-2
nip=nip-2
goti=goti+gsti*(nit-1)
goti=(pi/2-goti)
gopi=gopi+gspi
ALLOCATE(idep(0:nit+1,0:nip+1,2))
do i=1,nin
   if(i.eq.itl.or.i.eq.itl+1)then
      do j=0,nit+1
         do k=0,nip+1
            read(30,*)idep(nit+1-j,k,i-itl+1)
         enddo
      enddo
      if(i.eq.itl+1)EXIT
   else
      do j=0,nit+1
         do k=0,nip+1
            read(30,*)rdm1
         enddo
      enddo
   endif
enddo
close(30)
!
! Now populate the matrix ivel
!
DO i=1,nnt
   DO j=1,nnp
!
!     For each latitude and longitude grid point, work out the
!     two layer bounding nodes in depth.
!     Work out which interface grid cell this point lies in and 
!     the  corresponding u and v values
!
      xt=got+rgst*(i-1)
      zp=gop+rgsp*(j-1)
      id1=INT((xt-goti)/gsti)+1
      id2=INT((zp-gopi)/gspi)+1
      u=(xt-goti-(id1-1)*gsti)/(gsti)
      v=(zp-gopi-(id2-1)*gspi)/(gspi)
      if(id1.eq.nit)then
         id1=id1-1
         u=1.0
      endif
      if(id2.eq.nip)then
         id2=id2-1
         v=1.0
      endif
!
! Return value for top interface
!
      id3=1
      call ibspline(id1,id2,u,v,dept,id3)
      dept=dept-earth
!
! Return value for bottom interface
!
      id3=2
      call ibspline(id1,id2,u,v,depb,id3)
      depb=depb-earth
!
! Now determine which nodes are bound by these points
!
     id1=int(gor-dept)/rgsr+2
     id2=int(gor-depb)/rgsr+1
     do k=id1,id2
        ivel(k,i,j)=1
     enddo
  enddo
enddo
!
! Calculate the values of the basis functions
!
ALLOCATE(ui(nvr+1,4), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL ui'
ENDIF
DO i=1,ndr+1
   u=ndr
   u=(i-1)/u
   ui(i,1)=(1.0-u)**3/6.0
   ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   ui(i,4)=u**3/6.0
ENDDO
ALLOCATE(vi(nvt+1,4), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL vi'
ENDIF
DO i=1,ndt+1
   u=ndt
   u=(i-1)/u
   vi(i,1)=(1.0-u)**3/6.0
   vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   vi(i,4)=u**3/6.0
ENDDO
ALLOCATE(wi(nvp+1,4), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL wi'
ENDIF
DO i=1,ndp+1
   u=ndp
   u=(i-1)/u
   wi(i,1)=(1.0-u)**3/6.0
   wi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   wi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   wi(i,4)=u**3/6.0
ENDDO
!
! Calculate velocity values on refined grid
!
mvar=0.0
icnt=0
DO i=1,nvp-1
   conp=ndp
   IF(i==nvp-1)conp=ndp+1
   DO j=1,nvt-1
      cont=ndt
      IF(j==nvt-1)cont=ndt+1
      DO k=1,nvr-1
         conr=ndr
         IF(k==nvr-1)conr=ndr+1
         DO l=1,conp
            stp=ndp*(i-1)+l
            DO m=1,cont
               stt=ndt*(j-1)+m
               DO n=1,conr
                  str=ndr*(k-1)+n
                  sumi=0.0
                  sumi1=0.0
                  DO i1=1,4
                     sumj=0.0
                     sumj1=0.0
                     DO j1=1,4
                        sumk=0.0
                        sumk1=0.0
                        DO k1=1,4
                           rdm=ui(n,k1)*velsm(k-2+k1,j-2+j1,i-2+i1)
                           rdm1=ui(n,k1)*velrm(k-2+k1,j-2+j1,i-2+i1)
                           sumk=sumk+rdm
                           sumk1=sumk1+rdm1
                        ENDDO
                        sumj=sumj+vi(m,j1)*sumk
                        sumj1=sumj1+vi(m,j1)*sumk1
                     ENDDO
                     sumi=sumi+wi(l,i1)*sumj
                     sumi1=sumi1+wi(l,i1)*sumj1
                  ENDDO
                  velr(str,stt,stp)=sumi-sumi1
                  if(ivel(str,stt,stp).eq.1)then
                     mvar=mvar+(sumi-sumi1)**2
                     icnt=icnt+1
                  endif
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
mvar=mvar/REAL(icnt)
!
! Now estimate model roughness
!
mrough=0.0
icnt=0
DO i=2,nnp-1
   DO j=2,nnt-1
      DO k=2,nnr-1
         ri=gor-(k-1)*rgsr+earth
         risti=ri*sin(got+(j-1)*rgst)
         dp=velr(k,j,i+1)-2.0*velr(k,j,i)+velr(k,j,i-1)
         dp=dp/((risti*rgsp)**2)
         dt=velr(k,j+1,i)-2.0*velr(k,j,i)+velr(k,j-1,i)
         dt=dt/((ri*rgst)**2)
         dr=velr(k+1,j,i)-2.0*velr(k,j,i)+velr(k-1,j,i)
         dr=dr/(rgsr**2)
         if(ivel(k,j,i).eq.1)then
            mrough=mrough+ABS(dp)+ABS(dt)+ABS(dr)
            icnt=icnt+1
         endif
      ENDDO
   ENDDO
ENDDO
mrough=mrough/REAL(icnt)
WRITE(6,*)'Model variance in (km/s)**2 is ',mvar
WRITE(6,*)'Model roughness in (kms)**(-1) is ',mrough
DEALLOCATE(ui,vi,wi, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM cbsvel: REAL ui,vi,wi'
ENDIF
DEALLOCATE(velsm,velrm,velr,ivel, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM cbsvel: REAL velsm,velrm,velr,ivel'
ENDIF
STOP
END PROGRAM misfit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ibspline(gt,gp,u,v,sumi,iid)
USE globalp
IMPLICIT NONE
INTEGER :: i1,j1,stp,stt,iid,gt,gp
REAL(KIND=i10) :: u,v,sumi,sumj
REAL(KIND=i10), DIMENSION(4) :: ui,vi
!
! iid = Interface id.
! u,v = independent surface parameters
! ui,vi = bspline basis functions
! sumi,sumj = Summation variables for spline
! gt,gp =  i,j coordinate of current node cell
! dep = depth of point
!
! Compute the values of the basis functions
!
ui(1)=(1.0-u)**3/6.0
ui(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
ui(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
ui(4)=u**3/6.0
vi(1)=(1.0-v)**3/6.0
vi(2)=(4.0-6.0*v**2+3.0*v**3)/6.0
vi(3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
vi(4)=v**3/6.0
sumi=0.0
DO i1=1,4
   sumj=0.0
   DO j1=1,4
      sumj=sumj+ui(j1)*idep(gt-2+j1,gp-2+i1,iid)
   ENDDO
   sumi=sumi+vi(i1)*sumj
ENDDO
END SUBROUTINE ibspline

