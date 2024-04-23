      program dragcoefficient 
      implicit none 

      INTEGER :: l,ll,lmax,zone,zmax,emax,ipts,ipts2
      DOUBLE PRECISION :: dum,Ubulk,fForcing,pForcing,pLift
      DOUBLE PRECISION :: pref,mu,dx,dy,dz,AreaFront,fDrag,pDrag,AreaTop
      INTEGER,pointer,dimension(:) :: lzone
      DOUBLE PRECISION,pointer,dimension(:,:) :: x,y,z
      DOUBLE PRECISION,pointer,dimension(:,:) :: p,u,v,w,umag,taux,tauy,tauz,xE,yE,zE
      DOUBLE PRECISION,pointer,dimension(:,:) :: area,Px_area,Py_area,Pz_area
      DOUBLE PRECISION,pointer,dimension(:,:) :: Cpe,CfE,FpE,FfE
      DOUBLE PRECISION,pointer,dimension(:,:) :: ClE,FlE
      DOUBLE PRECISION,pointer,dimension(:,:) :: uv,uw,vw,ppref,pprefN
      DOUBLE PRECISION,pointer,dimension(:) :: Cp,Cf,Fp,Ff,Fl,Cl
      CHARACTER(300),pointer,dimension(:) :: zoneT
      CHARACTER(30),pointer,dimension(:) :: nzone
      CHARACTER(200) :: fin1,fin2,fout1,fout2
      

      mu=19.83*1E-6 ! N.m-2.s
      
      ! USER INPUT:
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Ubulk=0.5  ! (m.s-1) Inlet Bulk Velocity
      
      ! Projected Area - For Drag and Lift Coefficient
      AreaFront=0.008*0.01498 ! (m2) Front Projected Area (YZ)    c
      AreaTop=0.008*0.0965   ! (m2) Top Projected Area (XY)
      write(6,*) 'AreaFront:',AreaFront
      
      ! Grid Resolution:
      dx=0.0005  ! (m) 
      dy=0.001  ! (m)
      dz=0.0005  ! (m)
      
      ! Input file:
      fin1='CloudPoints_NACA-0012-2.5mm_FEM.dat' ! CloudPoints (output from FEM_For_Tecplot.f90)
      fin2='NACA-0012-2.5mm_FEM_Elements.dat'      ! FEM_Elements (output from FEM_For_Tecplot.f90)
      fout1='NACA-0012-2.5mm_FEM.dat'           
      fout2='NACA-0012-2.5mm_CD.dat'

      zmax=1 ! Number of zone FEM in fin1,fin2
      allocate(lzone(zmax),nzone(zmax))
      lzone(1)=672
!      lzone(2)=29980
!      lzone(3)=39462

      nzone(1)=' NACA-0012'     ! Give the zone name from fout1
!      nzone(2)=' Transit_Back'
!      nzone(3)=' Transit_Left'
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      lmax=maxval(lzone)
      emax=sum(lzone)
      
      allocate(x(lmax,zmax),y(lmax,zmax),z(lmax,zmax))
      allocate(area(lmax,zmax),Px_area(lmax,zmax),Py_area(lmax,zmax),Pz_area(lmax,zmax))
      allocate(u(lmax,zmax),v(lmax,zmax),w(lmax,zmax),p(lmax,zmax),umag(lmax,zmax))
      allocate(uw(lmax,zmax),uv(lmax,zmax),vw(lmax,zmax))
      allocate(ppref(lmax,zmax),pprefN(lmax,zmax))
      allocate(taux(lmax,zmax),tauy(lmax,zmax),tauz(lmax,zmax))
      allocate(FpE(lmax,zmax),FfE(lmax,zmax),Fp(zmax),Ff(zmax))
      allocate(CpE(lmax,zmax),CfE(lmax,zmax),Cp(zmax),Cf(zmax))
      allocate(zoneT(emax))
      allocate(xE(emax,3),yE(emax,3),zE(emax,3))

      allocate(ClE(lmax,zmax),FlE(lmax,zmax),Cl(zmax),Fl(zmax))

      WRITE(6,*) 
      WRITE(6,'(a)') ' Reading '//trim(adjustl(fin1))//' ...'
      open(200,file=fin1)
        read(200,*) 
        do zone=1,zmax
          read(200,*) 
          do l=1,lzone(zone)
            read(200,*) x(l,zone),y(l,zone),z(l,zone)   &
     &                 ,area(l,zone),Px_area(l,zone),Py_area(l,zone),Pz_area(l,zone)   &
     &                 ,p(l,zone),u(l,zone),v(l,zone),w(l,zone),pref,pprefN(l,zone)
                        umag(l,zone)=sqrt(u(l,zone)**2+v(l,zone)**2+w(l,zone)**2)
          end do
        end do
      close(200)
       
      do zone=1,zmax 
        if(zone.eq.2) then
          do l=1,lzone(zone)
            p(l,zone)=pprefN(l,zone)
          end do 
        else
          do l=1,lzone(zone)
            p(l,zone)=pprefN(l,zone)
          end do 
        end if
      end do
      p(:,:)=pprefN(:,:)!p(:,:)-pref
      
      ! Calculate AreaFront
!      AreaFront=0.d0
!      do l=1,lzone(1)
!        AreaFront = AreaFront + Px_area(l,1)
!        write(6,*) Px_Area(l,1),AreaFront
!      end do

      ! ELEMENTS Pressure Force + Pressure Coefficient CALCULATION:
      Cp(:)=0.d0 ; CpE(:,:)=0.d0 
      do zone=1,zmax ; do l=1,lzone(zone)
        FpE(l,zone) = (p(l,zone) * Px_area(l,zone))
        CpE(l,zone) = FpE(l,zone) / (0.5 * (Ubulk**2) * AreaFront)
      end do; end do
      
      ! FACE Pressure Force + Pressure Coefficient CALCULATION:
      Fp=0.d0 ; Cp=0.d0
      do zone=1,zmax ; do l=1,lzone(zone)
        Fp(zone) = Fp(zone) + FpE(l,zone)
        Cp(zone) = Cp(zone) + CpE(l,zone)
      end do ; end do

      ! TOTAL Friction CALCULATION:
      pDrag=0.d0 ; pForcing=0.d0
      do zone=1,zmax
        pForcing=pForcing+abs(Fp(zone))
        pDrag=pDrag+abs(Cp(zone))
      end do

      WRITE(6,*) 
      WRITE(6,'(a,E14.4,F20.1)') ' AreaFront, Pref:',AreaFront,Pref
      WRITE(6,*) 
      WRITE(6,*) 'PRESSURE DRAG:'
      WRITE(6,*) '=============='
      do zone=1,zmax
!        if(zone.le.3) then
          WRITE(6,'(a30,F9.3,F9.2,a)') nzone(zone),Cp(zone),abs(Cp(zone)/pDrag*100),' (%)'
!        end if
      end do
      WRITE(6,*) '=============='
      WRITE(6,'(a,F16.8,a)') ' Total Pressure Forcing',pForcing,' (N)'
      WRITE(6,'(a,F19.3)') ' Total Pressure Drag',pDrag
      


      ! ELEMENTS Shear Stress + Friction Force + Friction Coefficient CALCULATION:
      Cf(:)=0.d0 ; CfE(:,:)=0.d0 ; taux(:,:)=0.d0 ; tauy(:,:)=0.d0 ; tauz(:,:)=0.d0
      do zone=1,zmax ; do l=1,lzone(zone)
        uv(l,zone)=sqrt(u(l,zone)**2+v(l,zone)**2)
        uw(l,zone)=sqrt(u(l,zone)**2+w(l,zone)**2)
        vw(l,zone)=sqrt(v(l,zone)**2+w(l,zone)**2)

        
!        taux(l,zone) = 1.5!mu*u(l,zone)/dz
!        tauy(l,zone) = 1.5!mu*u(l,zone)/dy
!        tauz(l,zone) = 0.d0!mu*u(l,zone)/dx
        taux(l,zone) = mu*vw(l,zone)/dz
        tauy(l,zone) = mu*uw(l,zone)/dy
        tauz(l,zone) = mu*uv(l,zone)/dx

                                   ! XY                           ! XZ                           ! YZ
        FfE(l,zone) = taux(l,zone)*Pz_area(l,zone) + tauy(l,zone)*Py_area(l,zone) + tauz(l,zone)*Px_area(l,zone)  
        CfE(l,zone) = FfE(l,zone) / (0.5 * (Ubulk**2) * AreaFront)
      end do ; end do
      
      ! FACE Friction Force + Friction Coefficient CALCULATION:
      do zone=1,zmax ; do l=1,lzone(zone)
        Ff(zone) = Ff(zone) + FfE(l,zone)
        Cf(zone) = Cf(zone) + CfE(l,zone)
      end do ; end do
      
      ! TOTAL Friction CALCULATION:
      fDrag=0.d0 ; fForcing=0.d0
      do zone=1,zmax
        fForcing=fForcing+Ff(zone)
        fDrag=fDrag+Cf(zone)
      end do

      WRITE(6,*) 
      WRITE(6,*) 'FRICTION DRAG:'
      WRITE(6,*) '=============='
      do zone=1,zmax
        WRITE(6,'(a30,E14.4,F9.3,a)') nzone(zone),Cf(zone),Cf(zone)/fDrag*100,' (%)'
      end do
      WRITE(6,*) '=============='
      WRITE(6,'(a,E14.4,a)') ' Total Friction Forcing',fForcing,' (N)'
      WRITE(6,'(a,E14.4)') ' Total Friction Drag',fDrag
      WRITE(6,*) 
      WRITE(6,'(a,F29.3)') ' Total Drag Coefficient',pDrag+fDrag
      WRITE(6,'(a,F9.2,a)') ' Total Pressure Drag/Total Drag',pDrag/(pDrag+fDrag)*100,' (%)'
      WRITE(6,'(a,F9.2,a)') ' Total Friction Drag/Total Drag',fDrag/(pDrag+fDrag)*100,' (%)'
      WRITE(6,*) 
 

      ! ELEMENTS Pressure Force + Pressure Coefficient CALCULATION:
      Cl(:)=0.d0 ; ClE(:,:)=0.d0 
      do zone=1,zmax ; do l=1,lzone(zone)
        FlE(l,zone) = (p(l,zone) * Pz_area(l,zone))
        ClE(l,zone) = FlE(l,zone) / (0.5 * (Ubulk**2) * AreaTop)
      end do; end do

      
      ! FACE Pressure Force + Pressure Coefficient CALCULATION:
      Fl=0.d0 ; Cl=0.d0
      do zone=1,zmax ; do l=1,lzone(zone)
        Fl(zone) = Fl(zone) + FlE(l,zone)
        Cl(zone) = Cl(zone) + ClE(l,zone)
      end do ; end do

      ! TOTAL Friction CALCULATION:
      pDrag=0.d0 ; pForcing=0.d0
      do zone=1,zmax
        pForcing=pForcing+abs(Fl(zone))
        pLift=pLift+abs(Cl(zone))
      end do

      WRITE(6,*) 
      WRITE(6,'(a,E14.4,F20.1)') ' AreaTop, Pref:',AreaTop,Pref
      WRITE(6,*) 
      WRITE(6,*) 'PRESSURE LIFT:'
      WRITE(6,*) '=============='
      do zone=1,zmax
!        if(zone.le.3) then
          WRITE(6,'(a30,F9.3,F9.2,a)') nzone(zone),Cl(zone),abs(Cl(zone)/pLift*100),' (%)'
!        end if
      end do
      WRITE(6,*) '=============='
      WRITE(6,'(a,F16.8,a)') ' Total Pressure Forcing',pForcing,' (N)'
      WRITE(6,'(a,F19.3)') ' Total Lift Coefficient',pLift
      WRITE(6,*) 
      



      ipts=0 ; ipts2=0
      WRITE(6,'(a)') ' Reading '//trim(adjustl(fin2))//' ...'
      open(200,file=fin2)
        read(200,*) 
        do zone=1,zmax ; do l=1,lzone(zone)
          ipts=ipts+1
          read(200,'(a)') zoneT(ipts)
          do ll=1,3
            read(200,*) xE(ipts,ll),yE(ipts,ll),zE(ipts,ll)
          end do
          read(200,*) ! 
          read(200,*) ! Element: 1,2,3 
        end do ; end do
      close(200)

      ipts=0 ; ipts2=0
      WRITE(6,'(a)') ' Writing '//trim(adjustl(fout1))//' ...'
      open(200,file=fout1)
        write(200,'(a)') 'Variables=x,y,z,umag,p,taux,tauy,Body'
        do zone=1,zmax ; do l=1,lzone(zone)
          ipts=ipts+1
          write(200,'(a)') trim(adjustl(zoneT(ipts)))
          do ll=1,3
            write(200,'(7E14.4,I9)') xE(ipts,ll),yE(ipts,ll),zE(ipts,ll),umag(l,zone),p(l,zone),taux(l,zone),tauy(l,zone),zone
          end do
          write(200,*)
          write(200,'(a)') '      1       2       3' 
        end do ; end do
      
      ipts=0 ; ipts2=0
      WRITE(6,'(a)') ' Writing '//trim(adjustl(fout2))//' ...'
      open(200,file=fout2)
        write(200,'(a)') 'Variables=x,y,z,umag,p,tauYZ,tauXZ,tauXY,PYZ_area,PXZ_area,PXY_area,CpE,CfE,l,Body'
        do zone=1,zmax 
          write(200,'(a)') 'ZONE T="'//trim(adjustl(nzone(zone)))//'"'
          do l=1,lzone(zone)
            write(200,'(13E14.4,2I9)') x(l,zone),y(l,zone),z(l,zone),umag(l,zone),p(l,zone)  & !5
     &                                ,taux(l,zone),tauy(l,zone),tauz(l,zone)  &               !3
     &                                ,Px_area(l,zone),Py_area(l,zone),Pz_area(l,zone) &       !3
     &                                ,CpE(l,zone),CfE(l,zone),l,zone             !2
          end do 
        end do
      close(200)


      end program 
