!****************************************************************************************************************************
! 3D channel flow simulation with direct numerical simulation 
!
! space : 2nd order central
!
! time :
!        2nd order Adams-Bashforth(advection term, diffusion term)
!        1st order backward Euler method(pressure gradient, continuity eq)
!
! grid :
!        staggered (unequal interval for y direction)
!
! boundary condition : 
!        x,z : periodic
!        y : wall (non-slip) 
!******************************************************************************************************************************
      program channelflowDNS
      parameter (NX=64,NY=64,NZ=64)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*3  CT(0:100)
      character*4  file1
      double precision ::Ret,dt
      common /file1/file1 /RMS/URMSX,URMSC,VRMSC,WRMSC
      common /STEP/ISTEP /TSTEP/TSTEP /START/ISTART
      common /DIV/DIVX /COU/COMX /ENE/ENE /UME/UME,UMX
      common /ERRP/POIERR /ITRP/ITRP 
      data CT/'000', &      
      '001','002','003','004','005','006','007','008','009','010',  &
      '011','012','013','014','015','016','017','018','019','020',  &
      '021','022','023','024','025','026','027','028','029','030',  &
      '031','032','033','034','035','036','037','038','039','040',  &
      '041','042','043','044','045','046','047','048','049','050',  &
      '051','052','053','054','055','056','057','058','059','060',  &
      '061','062','063','064','065','066','067','068','069','070',  &
      '071','072','073','074','075','076','077','078','079','080',  &
      '081','082','083','084','085','086','087','088','089','090',  &
      '091','092','093','094','095','096','097','098','099','100'/

         file1='dns2'
         Ret = 300.d0   
         dt=  1.0d-4
         iskip = 40            
         isave = 100
         istop=iskip*isave   
         icount=0
         ICT=0
      call SBRMSH
      call SBRDTR(CT(ICT))    
      write( 6,1000)
      call SBRUMR
      call SBRST1
      call SBRCHK
      write( 6,1100) ISTART,TSTEP,ITRP,POIERR,DIVX,COMX  &
         ,UME,UMX,INT(Ret*UME),INT(Ret*UMX),ENE,URMSX,URMSC,VRMSC,WRMSC
 1000 FORMAT(3X,4HSTEP,7X,1HT,1X,4HITRP,4X,6HPOIerr,4X,6HDIVmax  &
      ,2X,6HCOUmax,3X,5HUmean,4X,4HUmax,4X,3HRem,4X,3HRex   &
      ,4X,6HEnergy,3X,5HU'max,5X,3HU'c,5X,3HV'c,5X,3HW'c)
 1100 FORMAT(I7,F8.4,I5,2E10.2,3F8.3,2I7,F10.5,4F8.4)

      !!! main loop  
      do
         call SBRCON
         call SBRVIS
         if(icount >= 1) then                      
            call SBRPRE(1.5d0,-0.5d0)               
         else
            call SBRPRE(1.0d0, 0.0d0)
         end if
         call SBRRHP                            
         ipstop=20
         call SBRSOR(ipstop)        
         call SBRCOR                 

         icount = icount + 1                                   
         istep  = istart + icount
         tstep  = tstep  + dt

         if(mod(istep,10) == 0 .OR. icount <= 10) then
            call SBRUMR
            call SBRST1
            call SBRCHK
            write( 6,1100) istep,tstep,ITRP,POIERR,DIVX,COMX  &
               ,UME,UMX,INT(Ret*UME),INT(Ret*UMX),ENE,URMSX,URMSC,VRMSC,WRMSC
         end if

         if(mod(istep,iskip) == 0) then    
            ICT=ICT+1
            call SBRDTS(CT(ICT))
         end if

         if(icount >= istop) exit

         if(istep >= 10.AND.DIVX >= 1.d+2) then 
            write(6,*) 'Stopped  because of the numerical instability'
         end if
      end do

      end program


! *********************************************************************
! *     SBR. MSH  :  MESH PARAMETERS                                  *
! *********************************************************************

      SUBROUTINE SBRMSH
      PARAMETER (NX=64,NY=64,NZ=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision ::dx,dz,Lx,Lz,Ret,dt

      common /Y/YV(0:NY),YP(0:NY1)
      common /DY/DY(NY) 

      dt=  1.0d-4

      call MESHY
      call MESHXZ

      WRITE(6,1000) NX,LX,dx,Ret*LX,Ret*dx,NZ,LZ,dz,Ret*LZ,Ret*dz  &
      ,NY,1.,DY(1),DY(NY/2),Ret*DY(1),Ret*DY(NY/2),YP(1),Ret*YP(1)  &
      ,INT(Ret),dt
 1000 FORMAT(/10(1H*),"  IN SBR.MSH  ",10(1H*)/  &
      /5X,"NX=",I4,5X,"LX=",F7.3,3X,"dx=",F7.5  &
      ,5X,"LX+=",F7.1,3X,"dx+=",F7.3  &
      /5X,"NZ=",I4,5X,"LZ=",F7.3,3X,"dz=",F7.5  &
      ,5X,"LZ+=",F7.1,3X,"dz+=",F7.3  &
      /5X,"NY=",I4,5X,"HY=",F7.3,3X,"DYmin=",F7.5,3X,"DYmax=",F7.5  &
      /30X,"DYmin+=",F6.3,3X,"DYmax+=",F6.2  &
      /5X,"YP(1)=",F7.5,3X,"YP(1)+=",F6.2  &
      /5X,"Ret=",I6/5X,"dt=",F9.5)

      call MDQ1VN
      call MDQ1VD
      call MDQ1P
      call MDQ2V
      call MDQ2PV
      call MDQ2PP
      return
      end
!---------------------------------------------------------------------
!-----  Mesh generation .....  for X(mainstream direction), Z directions (across the main stream)  ------------------
      SUBROUTINE MESHXZ
      PARAMETER (NX=64,NY=64,NZ=64)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision ::dx,dz,Lx,Lz,Ret                 
      common /IP/IP(NX) /IM/IM(NX) /KP/KP(NZ) /KM/KM(NZ)

      Ret = 300.d0
      dx = 18.d0/Ret
      dz = 9.d0/Ret
      LX = dble(NX)*dx
      LZ = dble(NZ)*dz

!  -----  List vectors to identify the neighbouring point   ---------
      do I=1,NX
         IP(I)=I+1
         IM(I)=I-1
         if(IP(I) > NX) IP(I)=IP(I)-NX
         if(IM(I) < 1) IM(I)=IM(I)+NX
      end do
      do k = 1, NZ
         KP(k) = k + 1
         KM(k) = k - 1
         if(KP(k) > NZ) KP(k) = KP(k) - NZ
         if(KM(k) < 1) KM(k) = KM(k) + NZ
      end do

      return
      end
!---------------------------------------------------------------------
!-----  Mesh generation .....  for Y direction (wall-wall) ----------------------
      SUBROUTINE MESHY
      PARAMETER (NX=64,NY=64,NZ=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common /Y/YV(0:NY),YP(0:NY1) /DY/DY(NY) /DYC/DYC(NY-1)      
      ALG=0.95d0           
      AT=dlog((1.d0+ALG)/(1.d0-ALG))/2.d0
      YV(0)=0.
      do J=1,NY-1
         ETA=AT*(-1.d0 + 2.d0*dble(j)/dble(NY))      
         YV(J)=(dtanh(ETA)/ALG+1.d0)/2.d0
      end do
      YV(NY)=1.D0
      do j=1,NY
         ETA=AT*(-1.d0+2.d0*(dble(j)-0.5d0)/dble(NY))
         YP(J)=(dtanh(ETA)/ALG+1.d0)/2.d0
      end do
! ... Outer points (half mesh)
      YP(0)=2.d0*YV(0)-YP(1)
      YP(NY1)=2.d0*YV(NY)-YP(NY)

      do j = 1,NY
         DY(J)=-YV(J-1)+YV(J)
      end do
      do j = 1,NY-1
         DYC(J)=-YP(J)+YP(J+1)
      end do
      return
      end
!-----------------------------------------------------------------------
!---- FDM operators at V points using data on P points and INSIDE walls 
      SUBROUTINE MDQ1VN
      PARAMETER (NY=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common /Y/YV(0:NY),YP(0:NY1) /D1VN/D1VNM(0:NY),D1VNP(0:NY)
      do j = 1,NY-1
         DYV=YP(J+1)-YP(J)
         D1VNM(J)=-1.d0/DYV
         D1VNP(J)= 1.D0/DYV
      end do
!    ----- For the Neumann B.C. at the wall (DP/DY=0) for Pressure
      D1VNM(0)= 0.d0
      D1VNP(0)= 0.d0
      D1VNM(NY)= 0.d0
      D1VNP(NY)= 0.d0
      return
      end

!-----------------------------------------------------------------------
!---- FDM operators at V points using data on P points and AT walls ----
      SUBROUTINE MDQ1VD
      PARAMETER (NY=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                      
      COMMON /Y/YV(0:NY),YP(0:NY1) 
      COMMON /D0V/D0VM(NY-1),D0VP(NY-1) /D1V/D1VM(0:NY),D1VP(0:NY)
      do j=1,NY-1
         DYV=YP(J+1)-YP(J)
         D0VM(J)=(YP(J+1)-YV(J))/DYV
         D0VP(J)=(YV(J)-YP(J))/DYV
         D1VM(J)=-1.d0/DYV
         D1VP(J)= 1.d0/DYV
      end do

!!!!Caution!  Stencils are shifted at the wall point
      DYP0=YP(2)-YP(1)
      D1VM(0)= YP(2)/YP(1)/DYP0
      D1VP(0)=-YP(1)/YP(2)/DYP0
      DYPN=YP(NY)-YP(NY-1)
      D1VM(NY)= (1.-YP(NY))/(1.-YP(NY-1))/DYPN
      D1VP(NY)=-(1.-YP(NY-1))/(1.-YP(NY))/DYPN
      return
      end
!-----------------------------------------------------------------------
!---  FDM operators for velocity at P points using data on V points ----
      SUBROUTINE MDQ1P
      PARAMETER (NY=64,NY1=NY+1,NYH1=NY-1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /Y/YV(0:NY),YP(0:NY1) /DY/DY(NY)             
      COMMON /D0P/D0PM(NY),D0PP(NY) /D1P/D1PM(NY),D1PP(NY)
      do j = 1,NY
         D0PM(J)=(YV(J)-YP(J))/DY(J)
         D0PP(J)=(YP(J)-YV(J-1))/DY(J)
         D1PM(J)=-1.d0/DY(J)
         D1PP(J)= 1.d0/DY(J)
      end do
      return
      end
!-----------------------------------------------------------------------
!---  FDM operators for viscous terms at V points  ---------------------
      SUBROUTINE MDQ2V
      PARAMETER (NY=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common /D2V/D2VM(NY-1),D2V0(NY-1),D2VP(NY-1) 
      common /D1V/D1VM(0:NY),D1VP(0:NY) /D1P/D1PM(NY),D1PP(NY)
      do j=1,Ny-1
         D2VM(J)=D1VM(J)*D1PM(J)
         D2V0(J)=D1VM(J)*D1PP(J)+D1VP(J)*D1PM(J+1)
         D2VP(J)=D1VP(J)*D1PP(J+1)
      end do
      return
      end
!-----------------------------------------------------------------------
!---  FDM operators for viscous terms at P points  ---------------------
      SUBROUTINE MDQ2PV
      PARAMETER (NY=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common /D2P/D2PM(2:NY),D2P0(NY),D2PP(NY-1) 
      common /D1V/D1VM(0:NY),D1VP(0:NY) /D1P/D1PM(NY),D1PP(NY)
      j = 1
      D2P0(J)=D1PM(J)*D1VM(J-1)+D1PP(J)*D1VM(J)
      D2PP(J)=D1PM(J)*D1VP(J-1)+D1PP(J)*D1VP(J)
      do j = 2,NY-1
         D2PM(J)=D1PM(J)*D1VM(J-1)
         D2P0(J)=D1PM(J)*D1VP(J-1)+D1PP(J)*D1VM(J)
         D2PP(J)=                  D1PP(J)*D1VP(J)
      end do
      J=NY
      D2PM(J)=D1PM(J)*D1VM(J-1)+D1PP(J)*D1VM(J)
      D2P0(J)=D1PM(J)*D1VP(J-1)+D1PP(J)*D1VP(J)
      return
      end

!-----------------------------------------------------------------------
!---  FDM operators for pressure's Poisson equation  -----------------
      SUBROUTINE MDQ2PP
      PARAMETER (NY=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DPP/DPPM(2:NY),DPP0(NY),DPPP(NY-1)
      COMMON /D1VN/D1VNM(0:NY),D1VNP(0:NY) /D1P/D1PM(NY),D1PP(NY)
      do j = 1,NY
         if (J > 1) DPPM(J)=D1PM(J)*D1VNM(J-1)
         DPP0(J)=D1PM(J)*D1VNP(J-1)+D1PP(J)*D1VNM(J)
         if (J < NY) DPPP(J)=D1PP(J)*D1VNP(J)
      end do
      return
      end

! *********************************************************************
! *     SBR. DTR  :  DATA READ                                        *
! *********************************************************************

      SUBROUTINE SBRDTR(CT)
      PARAMETER (NX=64,NY=64,NZ=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 CT
      CHARACTER*4 file1
      common /file1/file1 /start/istart /tstep/tstep
      common /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      common /P/P(NZ,NX,NY) 
      open(10,file=file1//CT//'u.d',status='old',form='formatted')     
         read(10,1000) istart,tstep           
         read(10,2000) U
      close(10)
      open(10,file=file1//CT//'v.d',status='old',form='formatted')
         read(10,1000) istart,tstep
         read(10,2000) V
      close(10)
      open(10,file=file1//CT//'w.d',status='old',form='formatted')
         read(10,1000) istart,tstep
         read(10,2000) W
      close(10)
      open(10,file=FILE1//CT//'p.d',status='old',form='formatted')
         read(10,1000) istart,tstep
         read(10,3000) P
      close(10)
 1000 format(I10,F12.6)
 2000 format(8F10.6)   
 3000 format(8F10.5)
      write(6,1200) istart,tstep
 1200 format(/'Data READ',10X,'step=',I7,10X,'t=',F10.3)
      return
      end

! *********************************************************************
! *     SBR. CON  :  NONLINEAR TERM                                   *
! *********************************************************************

      SUBROUTINE SBRCON
      call SBRNLU
      call SBRNLV
      call SBRNLW
      return
      end
! ----------------------------------------------------------------------
      SUBROUTINE SBRNLU     !!!! calculate first component of convection term ( = UF)
      PARAMETER (NX=64,NY=64,NZ=64)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision :: dx,dz,Ret
      common /UF/UF(NZ,NX,NY)                      
      common /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      common /W1/CX(NZ,NX,NY) /W2/CY(NZ,NX,0:NY) /W3/CZ(NZ,NX,NY)
      common /D1V/D1VM(0:NY),D1VP(0:NY) /D0P/D0PM(NY),D0PP(NY)
      common /IP/IP(NX) /IM/IM(NX) /KP/KP(NZ) /KM/KM(NZ)

      Ret=300.d0
      dx = 18.d0/Ret      
      dz = 9.d0/Ret
! ... CX:-U^XU_X  at P  (U^X ... average of U in x direction  U_X ... differentiation of U in x direction)
! ... CZ:-W^XU_Z  at P 
      do j = 1, Ny
         do i = 1, Nx
            do K=1, Nz
               CX(K,I,J)=-(U(K,IM(I),J)+U(K,I,J))*(-U(K,IM(I),J)+U(K,I,J))/(2.d0*dx)  ! u(du/dx)
               CZ(K,I,J)=-(W(K,I,J)+W(K,IP(I),J))*(-U(K,I,J)+U(KP(K),I,J))/(2.d0*dz)  ! w(du/dz)
            end do
         end do
      end do

! ... CY:-V^XU_Y  at UV  
      do j = 0, Ny
         if(J > 0.AND.J < NY) then
         do  i = 1, Nx
            do  k = 1, Nz
               CY(K,I,J)=-(D1VM(J)*U(K,I,J)+D1VP(J)*U(K,I,J+1))*(V(K,I,J)+V(K,IP(I),J))/2.d0  
            end do
         end do
         else
            do i = 1, Nx
               do  k = 1, Nz
                  CY(K,I,J)=0.d0
               end do
            end do
         end if
      end do
! ... UF:-(U^XU_X)^X-(V^XU_Y)^Y-(W^XU_Z)^Z  at U
      do j = 1, Ny
         do i = 1, Nx
            do k = 1, Nz
               UF(K,I,J)=(CX(K,I,J)+CX(K,IP(I),J)+CZ(KM(K),I,J)+CZ(K,I,J))/2.d0+D0PM(J)*CY(K,I,J-1)+D0PP(J)*CY(K,I,J)
            end do
         end do
      end do
      return
      end
! ----------------------------------------------------------------------
      SUBROUTINE SBRNLV      !!!! calculate second component of convection term ( = VF)
      PARAMETER (NX=64,NY=64,NZ=64)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision::dx,dz,Ret
      common /VF/VF(NZ,NX,0:NY)                    
      common /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      common /W1/CX(NZ,NX,NY) /W2/CY(NZ,NX,0:NY) /W3/CZ(NZ,NX,NY)
      common /D0V/D0VM(NY-1),D0VP(NY-1) 
      common /D0P/D0PM(NY),D0PP(NY) /D1P/D1PM(NY),D1PP(NY)
      common /IP/IP(NX) /IM/IM(NX) /KP/KP(NZ) /KM/KM(NZ)

      Ret = 300.d0
      dz = 9.d0/Ret
      dx = 18.d0/Ret
! ... CY:-(V^Y*V_Y)^Y at P
      do j = 1, NY
         do  i = 1, NX
            do  K= 1, NZ
               CY(K,I,J)=-(D0PM(J)*V(K,I,J-1)+D0PP(J)*V(K,I,J))*(D1PM(J)*V(K,I,J-1)+D1PP(J)*V(K,I,J))
            end do
         end do
      end do

! ... CX:-(U^Y*V_X) at UV, CZ:-(W^Y*V_Z)  at VW
      do j = 1, Ny-1
         do i = 1, Nx
            do k = 1, Nz
               CX(K,I,J)=-(D0VM(J)*U(K,I,J)+D0VP(J)*U(K,I,J+1))*(-V(K,I,J)+V(K,IP(I),J)) / dx
               CZ(K,I,J)=-(D0VM(J)*W(K,I,J)+D0VP(J)*W(K,I,J+1))*(-V(K,I,J)+V(KP(K),I,J)) / dz
            end do
         end do
      end do

! ... VF:-(U^Y*V_X)^X-(V^Y*V_Y)^Y-(W^Y*V_Z)^Z  at V
      do j = 1, Ny-1
         do i = 1, Nx
            do k = 1, Nz
               VF(K,I,J)=(CX(K,IM(I),J)+CX(K,I,J)+CZ(KM(K),I,J)+CZ(K,I,J))/2.D0  &
                           +D0VM(J)*CY(K,I,J)+D0VP(J)*CY(K,I,J+1)
            end do
         end do
      end do
      return
      end
! ----------------------------------------------------------------------
      SUBROUTINE SBRNLW     !!! calculate third component of convection term ( = WF)
      PARAMETER (NX=64,NY=64,NZ=64)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision :: dx,dz,Ret
      common /WF/WF(NZ,NX,NY)
      common /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      common /W1/CX(NZ,NX,NY) /W2/CY(NZ,NX,0:NY) /W3/CZ(NZ,NX,NY)
      common /D1V/D1VM(0:NY),D1VP(0:NY) /D0P/D0PM(NY),D0PP(NY)
      common /IP/IP(NX) /IM/IM(NX) /KP/KP(NZ) /KM/KM(NZ)

      Ret=300.d0
      dx = 18.d0/Ret   
      dz = 9.d0/Ret
! ... CX:-U^ZW_X  at UW,  CZ:-W^ZW_Z  at P
      do j = 1, Ny
         do i = 1, Nx
            do k = 1, Nz
               CX(K,I,J)=-(U(K,I,J)+U(KP(K),I,J))*(-W(K,I,J)+W(K,IP(I),J))/(2.d0*dx)
               CZ(K,I,J)=-(W(KM(K),I,J)+W(K,I,J))*(-W(KM(K),I,J)+W(K,I,J))/(2.d0*dz)
            end do
         end do
      end do

! ... CY:-V^ZW_Y  at VW
      do j = 0, Ny
         if(J > 0.AND.J < NY) then
            do i = 1,Nx
               do k = 1, Nz
                  CY(K,I,J)=-(D1VM(J)*W(K,I,J)+D1VP(J)*W(K,I,J+1))*(V(K,I,J)+V(KP(K),I,J))/2.D0
               end do
            end do
         else
            do i = 1, Nx
               do k = 1, Nz
                  CY(K,I,J)=0.d0
               end do
            end do
         end if
      end do

! ... WF:-(U^ZW_X)^X-(V^ZW_Y)^Y-(W^ZW_Z)^Z  at W
      do j = 1, Ny
         do i = 1,Nx
            do k = 1, Nz
               WF(K,I,J)=(CX(K,IM(I),J)+CX(K,I,J)+CZ(K,I,J)+CZ(KP(K),I,J))/2.d0  &
                        +D0PM(J)*CY(K,I,J-1)+D0PP(J)*CY(K,I,J)
            end do
         end do
      end do
      return
      end

! *********************************************************************
! *     SBR. VIS  :  VISCOUS TERM                                     *
! *********************************************************************

      SUBROUTINE SBRVIS
      PARAMETER (NX=64,NY=64,NZ=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision::Ret,dx,dz,C00
      common /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY) 
      common /UF/UF(NZ,NX,NY) /VF/VF(NZ,NX,0:NY) /WF/WF(NZ,NX,NY)
      common /D2V/D2VM(NY-1),D2V0(NY-1),D2VP(NY-1) 
      common /D2P/D2PM(2:NY),D2P0(NY),D2PP(NY-1) 
      common /IP/IP(NX) /IM/IM(NX) 
      common /KP/KP(NZ) /KM/KM(NZ) 

      Ret = 300.d0
      dx = 18.d0/Ret
      dz = 9.d0/Ret

      do j = 1, Ny
         C00 = -2.d0/(dx**2) - 2.d0/(dz**2) + D2P0(J)
         if(J == 1) then
            do i = 1,Nx
               do k = 1,Nz
                  UF(K,I,J)=UF(K,I,J)  &                                   ! here, UF = (convection term) + (viscous term) 
                     +(U(K,IM(I),J)/(dx**2) + U(KM(K),I,J)/(dz**2) + C00*U(K,I,J)   &
                     +U(KP(K),I,J)/(dz**2) + U(K,IP(I),J)/(dx**2) + D2PP(J)*U(K,I,J+1))/Ret   
                  WF(K,I,J)=WF(K,I,J)  &
                     +(W(K,IM(I),J)/(dx**2) + W(KM(K),I,J)/(dz**2) + C00*W(K,I,J)   &
                     +W(KP(K),I,J)/(dz**2) + W(K,IP(I),J)/(dx**2) + D2PP(J)*W(K,I,J+1))/Ret
               end do
            end do
         else if(J < NY) then
            do i = 1, Nx
               do k = 1, Nz
                  UF(K,I,J)=UF(K,I,J)+(D2PM(J)*U(K,I,J-1) + U(K,IM(I),J)/(dx**2) + U(KM(K),I,J)/(dz**2)  &
                        +C00*U(K,I,J)+U(KP(K),I,J)/(dz**2) + U(K,IP(I),J)/(dx**2) + D2PP(J)*U(K,I,J+1))/Ret
                  WF(K,I,J)=WF(K,I,J)+(D2PM(J)*W(K,I,J-1) + W(K,IM(I),J)/(dx**2) + W(KM(K),I,J)/(dz**2) &
                        +C00*W(K,I,J)+W(KP(K),I,J)/(dz**2) + W(K,IP(I),J)/(dx**2) + D2PP(J)*W(K,I,J+1))/Ret
               end do
            end do
         else
            do i = 1, NX
               do k = 1, NZ
                  UF(K,I,J)=UF(K,I,J)+(D2PM(J)*U(K,I,J-1) + U(K,IM(I),J)/(dx**2) + U(KM(K),I,J)/(dz**2) &
                        +C00*U(K,I,J)+U(KP(K),I,J)/(dz**2) + U(K,IP(I),J)/(dx**2) )/Ret
                  WF(K,I,J)=WF(K,I,J)+(D2PM(J)*W(K,I,J-1) + W(K,IM(I),J)/(dx**2) + W(KM(K),I,J)/(dz**2)  &
                        +C00*W(K,I,J)+W(KP(K),I,J)/(dz**2) + W(K,IP(I),J)/(dx**2))/Ret
               end do
            end do
         end if
      end do

      do j = 1,Ny-1
         C00 = -2.d0/(dx**2) - 2.d0/(dz**2) + D2V0(J)
         do i = 1, Nx
            do k = 1, Nz
               VF(K,I,J)=VF(K,I,J)+(D2VM(J)*V(K,I,J-1) + V(K,IM(I),J)/(dx**2) + V(KM(K),I,J)/(dz**2)  &
                     +C00*V(K,I,J)+V(KP(K),I,J)/(dz**2) + V(K,IP(I),J)/(dx**2) + D2VP(J)*V(K,I,J+1))/Ret
            end do
         end do
      end do
      return
      end

! *********************************************************************
! *     SBR. PRE  :  PREDICTION STEP                                  *
! *********************************************************************

      SUBROUTINE SBRPRE(AB,BB)
      PARAMETER (NX=64,NY=64,NZ=64)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision::dt
      common /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      common /UF/UF(NZ,NX,NY) /VF/VF(NZ,NX,0:NY) /WF/WF(NZ,NX,NY)
      common /UB/UB(NZ,NX,NY) /VB/VB(NZ,NX,NY-1) /WB/WB(NZ,NX,NY)

      dt=  1.0d-4

      DAB = AB*dt
      DBB = BB*dt
      do j = 1,Ny
         do i = 1,Nx
            do k = 1, Nz
               UF(K,I,J)=UF(K,I,J)+2.d0                                 !! fractional step method !!
               U(K,I,J)=U(K,I,J)+DAB*UF(K,I,J)+DBB*UB(K,I,J)            !! 2nd order Adams Bashforth
               W(K,I,J)=W(K,I,J)+DAB*WF(K,I,J)+DBB*WB(K,I,J)          
               UB(K,I,J)=UF(K,I,J)
               WB(K,I,J)=WF(K,I,J)
            end do
         end do
      end do

      do j = 1, Ny-1
         do i = 1,Nx
            do k = 1,Nz
               V(K,I,J)=V(K,I,J)+DAB*VF(K,I,J)+DBB*VB(K,I,J)
               VB(K,I,J)=VF(K,I,J)
            end do
         end do
      end do
      return
      end

! *********************************************************************
! *     SBR. RHP  :  R.H.S. OF POISSON EQ.                            *
! *********************************************************************

! --- R.H.S. FOR SCALER POTANTIAL -------------------------------------
      SUBROUTINE SBRRHP                   !!!calculate right hand side of poisson Eq
      PARAMETER (NX=64,NY=64,NZ=64)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision ::dx,dz,Ret,dt
      common /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      common /W3/Q(NZ,NX,NY)  
      common /D1P/D1PM(NY),D1PP(NY) /IM/IM(NX) /KM/KM(NZ)

      dt=  1.0d-4
      Ret=300.d0
      dx = 18.d0/Ret
      dz = 9.d0/Ret

      do j = 1,Ny
         do i = 1,Nx
            do k = 1,Nz
               Q(K,I,J)=((-U(K,IM(I),J)+U(K,I,J))/dx+D1PM(J)*V(K,I,J-1)+D1PP(J)*V(K,I,J)  &
                        +(-W(KM(K),I,J)+W(K,I,J))/dz )/dt                   !!!right hand side of poisson Eq
            end do
         end do
      end do
      return
      end

! *********************************************************************
! *     SBR. SOR  :  S.O.R. SCHEME FOR POISSON EQ.                    *
! *********************************************************************

      SUBROUTINE SBRSOR(ISOR)                                     
      PARAMETER (NX=64,NY=64,NZ=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision::Ret,dz,dx,c00
      common /DPP/DPPM(2:NY),DPP0(NY),DPPP(NY-1)
      common /P/P(NZ,NX,NY) /W3/Q(NZ,NX,NY) /ERRP/POIERR /ITRP/ITRP
      common /IP/IP(NX) /IM/IM(NX) 
      common /KP/KP(NZ) /KM/KM(NZ) 

      Ret=300.d0
      dz = 9.d0/Ret
      dx = 18.d0/Ret

      ITRP=0
      SUMS=0.D0

      do j = 1,Ny
         do i = 1, Nx
            do k= 1, Nz
               SUMS=SUMS+Q(K,I,J)**2
            end do
         end do
      end do

      do
         ITRP=ITRP+1
         SUMR=0.
         do j = 1,Ny
            C00 = -2.d0/(dx**2) - 2.d0/(dz**2) + DPP0(J)
            DPC=-1.5D0/C00
            if(J == 1) then
               do i = 1,Nx
                  do k = 1, Nz
                     RESI = P(K,IM(I),J)/(dx**2) + P(KM(K),I,J)/(dz**2) + C00*P(K,I,J)  &
                           +P(KP(K),I,J)/(dz**2) + P(K,IP(I),J)/(dx**2) +DPPP(J)*P(K,I,J+1)  &
                           -Q(K,I,J)
                     P(K,I,J)=P(K,I,J)+RESI*DPC
                     SUMR=SUMR+RESI**2
                  end do
               end do
            else if(J < NY) then
               do i = 1, Nx
                  do k = 1, Nz
                     RESI=DPPM(J)*P(K,I,J-1) + P(K,IM(I),J)/(dx**2) + P(KM(K),I,J)/(dz**2) + C00*P(K,I,J) &
                           +P(KP(K),I,J)/(dz**2) + P(K,IP(I),J)/(dx**2) +DPPP(J)*P(K,I,J+1)-Q(K,I,J)
                     P(K,I,J)=P(K,I,J)+RESI*DPC
                     SUMR=SUMR+RESI**2
                  end do
               end do
            else
               do i = 1,Nx
                  do k = 1,Nz
                     RESI=DPPM(J)*P(K,I,J-1) + P(K,IM(I),J)/(dx**2) + P(KM(K),I,J)/(dz**2)  &
                           +C00*P(K,I,J)+P(KP(K),I,J)/(dz**2) + P(K,IP(I),J)/(dx**2)  &
                           -Q(K,I,J)
                     P(K,I,J)=P(K,I,J)+RESI*DPC
                     SUMR=SUMR+RESI**2
                  end do
               end do
            end if
         end do

         POIERR=DSQRT(SUMR/SUMS)
         !WRITE(6,*) 'ITRP=',ITRP,'   ERR=',POIERR
         if(ITRP >= ISOR.OR.POIERR < 1.D-5) exit
      end do

      return
      end

! *********************************************************************
! *     SBR. COR  :  CORRECTION STEP                                  *
! *********************************************************************  

      SUBROUTINE SBRCOR
      PARAMETER (NX=64,NY=64,NZ=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision::dx,dz,Ret,dt
      common /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      common /P/P(NZ,NX,NY) /DY/DY(NY)
      common /IP/IP(NX) /KP/KP(NZ) /D1VN/D1VNM(0:NY),D1VNP(0:NY)

      dt=  1.0d-4
      Ret=300.d0
      dx = 18.d0/Ret
      dz = 9.d0/Ret

      TX=dt/dx
      TZ=dt/dz
      do j = 1,Ny
         do i=1,Nx
            do k = 1,Nz
               U(K,I,J)=U(K,I,J)-TX*(-P(K,I,J)+P(K,IP(I),J))
               W(K,I,J)=W(K,I,J)-TZ*(-P(K,I,J)+P(KP(K),I,J))
            end do
         end do
      end do

      do j=1,Ny-1
         do i=1,Nx
            do k=1,Nz
               V(K,I,J)=V(K,I,J)-dt*(D1VNM(J)*P(K,I,J)+D1VNP(J)*P(K,I,J+1))               
            end do
         end do
      end do
      return
      end

! *********************************************************************
! *     SBR. UMR  :  MEAN & RMS VALUES                                *
! *********************************************************************

      SUBROUTINE SBRUMR                   
      PARAMETER (NX=64,NY=64,NZ=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      common /UM/UM(NY),VM(0:NY),WM(NY) /S12/S12L(0:NY),S12V(0:NY)
      common /IP/IP(NX) /IM/IM(NX) /D0V/D0VM(NY-1),D0VP(NY-1)
      common /BVU/BVU2(NY),BVU3(NY),BVU4(NY)
      common /BVV/BVV2(0:NY),BVV3(0:NY),BVV4(0:NY)
      common /BVW/BVW2(NY),BVW3(NY),BVW4(NY)
      ARXZ=1./DBLE(NX*NZ)
      do j = 1, NY
         SUMU1=0.D0
         SUMU2=0.D0
         SUMW1=0.D0
         SUMW2=0.D0
         do i = 1, NX
            do k = 1, NZ
               SUMU1=SUMU1+U(K,I,J)
               SUMU2=SUMU2+U(K,I,J)**2
               SUMW1=SUMW1+W(K,I,J)
               SUMW2=SUMW2+W(K,I,J)**2
            end do
         end do
         UM(J)=ARXZ*SUMU1
         BVU2(J)=ARXZ*SUMU2
         WM(J)=ARXZ*SUMW1
         BVW2(J)=ARXZ*SUMW2
      end do
                        
      do j = 1,Ny-1
         SUMV1=0.D0
         SUMV2=0.D0
         SUMUV=0.D0
         do i = 1,Nx
            do k = 1,Nz
               SUMV1=SUMV1+V(K,I,J)
               SUMV2=SUMV2+V(K,I,J)**2
               SUMUV=SUMUV+(D0VM(J)*U(K,I,J)+D0VP(J)*U(K,I,J+1))*(V(K,I,J)+V(K,IP(I),J))  
            end do
         end do
         VM(J)=ARXZ*SUMV1
         BVV2(J)=ARXZ*SUMV2
         S12L(J)=ARXZ*SUMUV/2.D0
      end do
      return
      end

! *********************************************************************
! *     SBR. STT  :  TURBULENCE STATISTICS                            *
! *********************************************************************

      SUBROUTINE SBRST1           
      PARAMETER (NX=64,NY=64,NZ=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common /UM/UM(NY),VM(0:NY),WM(NY) /PM/PM(NY)
      common /UR/UR(NY),VR(0:NY),WR(NY) /PR/PR(NY)
      common /ENE/ENE /UME/UME,UMX /RMS/URMSX,URMSC,VRMSC,WRMSC
      common /DY/DY(NY) /DYC/DYC(NY-1)
      common /BVU/BVU2(NY),BVU3(NY),BVU4(NY)
      common /BVV/BVV2(0:NY),BVV3(0:NY),BVV4(0:NY)
      common /BVW/BVW2(NY),BVW3(NY),BVW4(NY)
      common /BVP/BVP2(NY)
      ENE=0.d0
      UME=0.d0
      UMX=0.d0
      URMSX=0.d0
      do j = 1,Ny
         UR(J)=dsqrt(BVU2(J)-UM(J)**2)
         WR(J)=dsqrt(BVW2(J)-WM(J)**2)
         PR(J)=dsqrt(BVP2(J)-PM(J)**2)
         ENE=ENE+DY(J)*(UR(J)**2+WR(J)**2)
         UME=UME+DY(J)*UM(J)
         UMX=DMAX1(UMX,UM(J))
         URMSX=DMAX1(URMSX,UR(J))
      end do
      do j = 1,Ny-1
         VR(J)=dsqrt(BVV2(J)-VM(J)**2)
         ENE=ENE+DYC(J)*VR(J)**2
      end do
      URMSC=(UR(NY/2)+UR(NY/2+1))/2.D0
      VRMSC=VR(NY/2)
      WRMSC=(WR(NY/2)+WR(NY/2+1))/2.D0
      ENE=ENE/2.D0
      return
      end

! *********************************************************************
! *     SBR. CHK  :  CHECK of DIVERGENCE and COURANT-NUMBER           *
! *********************************************************************

      SUBROUTINE SBRCHK
      PARAMETER (NX=64,NY=64,NZ=64)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision::dx,dz,Ret,dt
      common /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)      
      common /DY/DY(NY) /DIV/DIVX /COU/COMX
      common /D0P/D0PM(NY),D0PP(NY) /D1P/D1PM(NY),D1PP(NY)
      common /UR/UR(NY),VR(0:NY),WR(NY)
      common /IM/IM(NX) /KM/KM(NZ)

      dt=  1.0d-4
      Ret=300.d0
      dx = 18.d0/Ret
      dz = 9.d0/Ret

      DIVX=0.d0
      COMX=0.d0
      do j = 1,Ny
         DNOMAL=(dx*DY(J)*dz)**(1.D0/3.D0)/UR(J)
         do i = 1,Nx
            do k = 1,Nz
            DIV=(-U(K,IM(I),J)+U(K,I,J))/dx+D1PM(J)*V(K,I,J-1)+D1PP(J)*V(K,I,J)  &
               +(-W(KM(K),I,J)+W(K,I,J))/dz
            UCP=(U(K,IM(I),J)+U(K,I,J))/2.D0
            WCP=(W(KM(K),I,J)+W(K,I,J))/2.D0
            VCP=D0PM(J)*V(K,I,J-1)+D0PP(J)*V(K,I,J)
            DIVX=DMAX1(DIVX,DNOMAL*DIV)
            COU = dt*(DABS(UCP)/dx+DABS(VCP)/DY(J)+DABS(WCP)/dz)
            COMX=DMAX1(COMX,COU)
            end do
         end do
      end do
      return
      end

! *********************************************************************
! *     SBR. DTS  :  DATA SAVE         and calculate Q_criterion      *
! *********************************************************************

      SUBROUTINE SBRDTS(CT)
      PARAMETER (NX=64,NY=64,NZ=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*3  CT
      character*4  FILE1
      character(9) :: FileName='output/t_'                     
      integer :: i,j,k             
      double precision::dx,dz,Ret  
      double precision::Q_c(1:NZ,1:NX,1:NY)        
      double precision::dudx(1:NZ,1:NX,1:NY),dudy(1:NZ,1:NX,1:NY),dudz(1:NZ,1:NX,1:NY) 
      double precision::dvdx(1:NZ,1:NX,1:NY),dvdy(1:NZ,1:NX,1:NY),dvdz(1:NZ,1:NX,1:NY) 
      double precision::dwdx(1:NZ,1:NX,1:NY),dwdy(1:NZ,1:NX,1:NY),dwdz(1:NZ,1:NX,1:NY)                          
      common /FILE1/FILE1 /STEP/ISTEP /TSTEP/TSTEP
      common /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      common /P/P(NZ,NX,NY)                                 
      common /Y/YV(0:NY),YP(0:NY1)    
      common /IP/IP(NX) /IM/IM(NX) /KP/KP(NZ) /KM/KM(NZ)  
      common /DY/DY(NY)                

      Ret=300.d0
      dx = 18.d0/Ret
      dz = 9.d0/Ret

!********************************************************************************************************
!*       calculate Q-criterion of velocity gradient (Second invariant of the velocity gradient tensor)  * 
!********************************************************************************************************

      do j = 1,NY
         do i = 1,NX
            do k = 1,NZ
               dudx(k,i,j) = (-u(k,im(i),j)+u(k,i,j))/dx 
               dvdy(k,i,j) = (v(k,i,j)-v(k,i,j-1))/dy(j)
               dwdz(k,i,j) = (-w(km(k),i,j)+w(k,i,j))/dz
               ! calculate dudy
               if (1 < j .AND. j < NY ) then
                  dudy(k,i,j) = (-(u(k,i,j-1)+u(k,i,j))/2.d0 + (u(k,i,j+1)+u(k,i,j))/2.d0)/dy(j)
                  dwdy(k,i,j) = (-(w(k,i,j-1)+w(k,i,j))/2.d0 + (w(k,i,j)+w(k,i,j+1))/2.d0)/dy(j)
               else if (j == 1) then 
                  dudy(k,i,j) = (-u(k,i,j)/2.d0 + (u(k,i,j+1)+u(k,i,j))/2.d0)/dy(j)
                  dwdy(k,i,j) = (-w(k,i,j)/2.d0 + (w(k,i,j)+w(k,i,j+1))/2.d0)/dy(j)
               else 
                  dudy(k,i,j) = (-(u(k,i,j-1)+u(k,i,j))/2.d0 + u(k,i,j)/2.d0)/dy(j)
                  dwdy(k,i,j) = (-(w(k,i,j-1)+w(k,i,j))/2.d0 + w(k,i,j)/2.d0)/dy(j)
               end if
               dudz(k,i,j) = (-(u(km(k),i,j)+u(k,i,j))/2.d0 + (u(kp(k),i,j)+u(k,i,j))/2.d0)/dz
               dvdx(k,i,j) = (-(v(k,im(i),j)+v(k,i,j))/2.d0 + (v(k,i,j)+v(k,ip(i),j))/2.d0)/dx 
               dvdz(k,i,j) = (-(v(km(k),i,j)+v(k,i,j))/2.d0 + (v(kp(k),i,j)+v(k,i,j))/2.d0)/dz
               dwdx(k,i,j) = (-(w(k,im(i),j)+w(k,i,j))/2.d0 + (w(k,ip(i),j)+w(k,i,j))/2.d0)/dx
            end do
         end do
      end do

      do j=1,Ny
         do i=1,Nx
            do k=1,NZ
               Q_c(k,i,j) = - 0.5d0*(dudx(k,i,j)**2 + dvdy(k,i,j)**2 + dwdz(k,i,j)**2  &
                           + 2.d0*dudy(k,i,j)*dvdx(k,i,j) + 2.d0*dudz(k,i,j)*dwdx(k,i,j) & 
                           + 2.d0*dvdx(k,i,j)*dwdy(k,i,j))
            end do
         end do
      end do

      open(20, file=Filename//CT//'.vtk', status='unknown')

      write( 20, '(a)' ) '# vtk DataFile Version 1.0'
      write( 20, '(a)' ) trim( FileName )
      write( 20, '(a)' ) 'ASCII'
!    -------------- output grid data --------
      write( 20, '(a)' ) 'DATASET STRUCTURED_GRID'
      write( 20, '(a,3(a,i3))' ) 'DIMENSIONS', ' ', NX, ' ', NY, ' ', NZ    
      write( 20, '(a,i6,a)' ) 'POINTS ', NX*NY*NZ, ' float'             
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               write( 20, '(3(f8.5,a))' ) dx*dble(i), ' ', yp(j), ' ', dz*dble(k), ' '      
            end do
         end do
      end do
!       ------------- output u,v,w,p ---------------
      write( 20, '(a,i8)' ) 'POINT_DATA ', NX*NY*NZ

      write( 20, '(a)' ) 'SCALARS u float'   
      write( 20, '(a)' ) 'LOOKUP_TABLE default'
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX                                        
               write( 20, '(f8.5)' ) U(k,i,j)
            end do
         end do
      end do
   
      write( 20, '(a)' ) 'SCALARS v float'
      write( 20, '(a)' ) 'LOOKUP_TABLE default'
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               write( 20, '(f8.5)' ) V(k,i,j)
            end do
         end do
      end do

      write( 20, '(a)' ) 'SCALARS w float'
      write( 20, '(a)' ) 'LOOKUP_TABLE default'
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               write( 20, '(f8.5)' ) W(k,i,j)
            end do
         end do
      end do

      write( 20, '(a)' ) 'SCALARS p float'
      write( 20, '(a)' ) 'LOOKUP_TABLE default'
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               write( 20, '(f8.5)' ) P(k,i,j)
            end do
         end do
      end do

      !------------- output Q_c ---------------
      write( 20, '(a)' ) 'SCALARS QC float'   
      write( 20, '(a)' ) 'LOOKUP_TABLE default'
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX                                        
               write( 20, '(f8.5)' ) Q_c(k,i,j)
            end do
         end do
      end do

      write(6,1200) ISTEP,TSTEP
 1200 format(/'Data SAVE',10X,'STEP=',I7,10X,'T=',F10.3/)

      close(20)
      return
      end