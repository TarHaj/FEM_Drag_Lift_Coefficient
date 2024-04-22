      program Msh2dat
      implicit none

      INTEGER :: ref1,ref2,ref3,sn,ll,nbody,body,fNmax,fEmax
      INTEGER :: l,n,entity1,entity2,entity3,entity4,entities,magnitude,ipts2
      INTEGER :: nodepackages,node,ipts,dum,dumlines
      DOUBLE PRECISION :: Cxx,Cyy,Czz,norm,dx,dy,dz,AreaProjected
      DOUBLE PRECISION :: p1,p2,p3,a1,a2,a3,s,v,vx1,vy1,vz1,vx2,vy2,vz2,vx3,vy3,vz3,wtot
      INTEGER,pointer,dimension(:) :: filenodes,fileelements
      DOUBLE PRECISION,pointer,dimension(:,:) :: nx,ny,nz,tx,ty,tz,xc,yc,zc,area,Px_area,Py_area,Pz_area
      DOUBLE PRECISION,pointer,dimension(:,:) :: w,cos_thetaX,cos_thetaY,cos_thetaZ
      DOUBLE PRECISION,pointer,dimension(:,:) :: xf,yf,zf
      DOUBLE PRECISION,pointer,dimension(:) :: order
      INTEGER,pointer,dimension(:,:,:) :: c
      CHARACTER(300),pointer,dimension(:) :: filename
      CHARACTER(200) :: fn,chN,chE,chl,fileall,fileend
      
      nbody=6
      allocate(filenodes(nbody),fileelements(nbody),filename(nbody))

      filename(1)='Ahmed_FEM_Head_0.020m.msh'
!      filename(2)='Ahmed_FEM_Slope_0.004m.msh'
      filename(2)='Ahmed_FEM_Back_0.020m.msh'
      filename(3)='Ahmed_FEM_Top_0.020m.msh'
      filename(4)='Ahmed_FEM_Bot_0.020m.msh'
      filename(5)='Ahmed_FEM_Right_0.020m.msh'
      filename(6)='Ahmed_FEM_Left_0.020m.msh'
      fileall='Ahmed_FEM_All_0.020m.dat'
      fileend='Ahmed_FEM_Elements_0.020m.dat'



      WRITE(6,*) '##################################' 
      WRITE(6,*) 'PRE-PROCESSING FOR FEM CALCULATION' 
      WRITE(6,*) '##################################' 
      WRITE(6,*) '1.Export GMsh Different Surface' 
      WRITE(6,*) '  -> Delete Surface' 
      WRITE(6,*) '  -> Export (Export all element)' 
      WRITE(6,*) '2.Run this software' 
      WRITE(6,*) '3.Export files for Tecplot and FEM_DragCoefficient.f90' 
      WRITE(6,*) '  -> *_Elements_*m.dat' 
      WRITE(6,*) '  -> CloudPoints_*.dat' 


      WRITE(6,*)
      WRITE(6,*) 'Importing FEM Files'
      WRITE(6,*) '-------------------'
     
      fNmax=-1000 ; fEmax=-1000
      DO body=1,nbody
      
      fn=trim(adjustl(filename(body)))
      WRITE(6,*) 'Reading '//trim(adjustl(filename(body)))//'...'

      open(20,file=fn)
      do l=1,4
	read(20,*) 
      end do

      ! READ THE ENTITIES SECTION:
      read(20,*) entity1,entity2,entity3,entity4
      entities=entity1+entity2+entity3+entity4
      do l=1,entities
	read(20,*) 
      end do
      read(20,*) !#EndEntities'
      
      ! READ THE NODES SECTION:
      read(20,*) !#Nodes'
      read(20,*) nodepackages,filenodes(body)
      ipts=0
      do n=1,nodepackages
	read(20,*) dum,dum,dum,node
	do dumlines=1,node     ! read all node ID
	  read(20,*) 
	end do
	do l=1,node	      ! read all node cordinates
	  ipts=ipts+1	  
         read(20,*) dum
	end do
      end do
      read(20,*) !#EndNodes'

      magnitude=3
      ! READ THE CONNECTIVITY SECTION:
       read(20,'(a)') chN
       read(20,*) nodepackages,fileelements(body)
      
      fNmax=max(fNmax,filenodes(body))
      fEmax=max(fEmax,fileelements(body))
      close(20)
      END DO
      WRITE(6,*) '-------------------'
      
      WRITE(6,*)
      WRITE(6,*) 'Enter the axis reference x=1 ; y=2 ; z=3'
      read(*,*) ref1,ref2,ref3

      WRITE(6,*)
      WRITE(6,*) 'Enter Cx,Cy,Cz'
      read(*,*) Cxx,Cyy,Czz

      WRITE(6,*)
      WRITE(6,*) 'Enter dx,dy,dz'
      read(*,*) dx,dy,dz

      
!      WRITE(6,'(a,2I9)') 'fNmax, fEmax:',fNmax,fEmax 

      allocate(xf(fNmax,nbody),yf(fNmax,nbody),zf(fNmax,nbody))
      allocate(c(fEmax,5,nbody),order(fEmax))!,cdum(fileelements,3))



      DO body=1,nbody
      
      ! READ .MSH FILE
      fn=trim(adjustl(filename(body)))


      open(20,file=fn)
      do l=1,4
	read(20,*) 
      end do

      ! READ THE ENTITIES SECTION:
      read(20,*) entity1,entity2,entity3,entity4
      entities=entity1+entity2+entity3+entity4
      do l=1,entities
	read(20,*) 
      end do
      read(20,*) !#EndEntities'
      
      
      ! READ THE NODES SECTION:
      read(20,*) !#Nodes'
      read(20,*) nodepackages,filenodes(body)
      ipts=0
      do n=1,nodepackages
	read(20,*) dum,dum,dum,node
	do dumlines=1,node     ! read all node ID
	  read(20,*) 
	end do
	do l=1,node	      ! read all node cordinates
	  ipts=ipts+1	  
	  if(ref1.eq.1 .and. ref2.eq.2 .and. ref3.eq.3)then  !XYZ
	    read(20,*) xf(ipts,body),yf(ipts,body),zf(ipts,body)
	  
	  else if(ref1.eq.2 .and. ref2.eq.1 .and. ref3.eq.3)then  !YXZ
	    read(20,*) yf(ipts,body),xf(ipts,body),zf(ipts,body)
	  
	  else if(ref1.eq.1 .and. ref2.eq.3 .and. ref3.eq.2)then  !XZY
	    read(20,*) xf(ipts,body),zf(ipts,body),yf(ipts,body)
	  
	  else if(ref1.eq.3 .and. ref2.eq.1 .and. ref3.eq.2)then  !ZXY
	    read(20,*) zf(ipts,body),xf(ipts,body),yf(ipts,body)

	  else if(ref1.eq.2 .and. ref2.eq.3 .and. ref3.eq.1)then  !YZX
	    read(20,*) yf(ipts,body),zf(ipts,body),xf(ipts,body)
	  
	  else if(ref1.eq.3 .and. ref2.eq.2 .and. ref3.eq.1)then  !ZYX
	    read(20,*) zf(ipts,body),yf(ipts,body),xf(ipts,body)
	  end if
!	  read(20,*) xf(ipts),yf(ipts),zf(ipts)
	  xf(ipts,body)=xf(ipts,body)/1000+Cxx
	  yf(ipts,body)=yf(ipts,body)/1000+Cyy
	  zf(ipts,body)=zf(ipts,body)/1000+Czz
	end do
      end do
      read(20,*) !#EndNodes'

      magnitude=3
      ! READ THE CONNECTIVITY SECTION:
       read(20,'(a)') chN
       read(20,*) nodepackages,fileelements(body)
       ipts=0 ; ipts2=0
       do n=1,nodepackages
 	read(20,*) magnitude,dum,dum,node
 	ipts2=ipts2+node
 	do l=1,node	      ! read all node cordinates
 	  if(magnitude.ge.2) then 
 	    ipts=ipts+1
 	    order(ipts)=1+magnitude
 	    read(20,*) dum,c(ipts,1:order(ipts),body)
 	  else 
 	    read(20,*) 
 	  end if
 	end do
       end do
       read(20,*) !#EndNodes'
       close(20) 
      
      WRITE(chN,'(I20)') filenodes(body)
      WRITE(chE,'(I20)') ipts!fileelements

      fileelements(body)=ipts
      END DO !nbody
      
      
      WRITE(6,*) 
      WRITE(6,*) 'FEM NODES AND ELEMENTS'
      WRITE(6,*) '----------------------'
      ! Export the filename
      fn=trim(adjustl(fileall))
!      WRITE(6,*) 'Writing '//trim(adjustl(fn))//'... | Total Nodes:'//trim(adjustl(chN))
      
      open(20,file=fn)
      write(20,'(a)') 'Variables=x,y,z,body'
      DO body=1,nbody
      WRITE(chN,'(I20)') filenodes(body)
      WRITE(chE,'(I20)') fileelements(body)
      WRITE(6,'(a)') ' '//trim(adjustl(filename(body)))//'   NODES: '//trim(adjustl(chN))//'   ELEMENTS: '//trim(adjustl(chE))
      write(20,'(a)') 'ZONE T="'//trim(adjustl(filename(body)))//'", NODES='//trim(adjustl(chN))  &
     &//', ELEMENTS='//trim(adjustl(chE))//', ZONETYPE=FETRIANGLE, DATAPACKING=POINT'
      do l=1,filenodes(body)
	write(20,'(3E16.4,I12)') xf(l,body),yf(l,body),zf(l,body),body
      end do
      write(20,*)
      do l=1,fileelements(body)
	write(20,'(10I10)') c(l,1:order(l),body)
      end do
      END DO ! nbody
      close(20)
      WRITE(6,*) '----------------------'
      WRITE(6,*) 
      WRITE(6,'(a)') ' '//trim(adjustl(fn))//'    exported'




      allocate(nx(fEmax,nbody),ny(fEmax,nbody),nz(fEmax,nbody))
      allocate(tx(fEmax,nbody),ty(fEmax,nbody),tz(fEmax,nbody))
      allocate(xc(fEmax,nbody),yc(fEmax,nbody),zc(fEmax,nbody))
      allocate(area(fEmax,nbody))!,w(fEmax,3,nbody))
      allocate(Px_area(fEmax,nbody),Py_area(fEmax,nbody),Pz_area(fEmax,nbody))
      allocate(cos_thetaX(fEmax,nbody),cos_thetaY(fEmax,nbody),cos_thetaZ(fEmax,nbody))
      
      DO body=1,nbody
      do l=1,fileelements(body)
        
        p1=c(l,1,body) ; p2=c(l,2,body) ; p3=c(l,3,body)
        
        ! Calculate the centroid of each triangle
        xc(l,body) = (xf(p1,body) + xf(p2,body) + xf(p3,body)) / 3.d0
        yc(l,body) = (yf(p1,body) + yf(p2,body) + yf(p3,body)) / 3.d0
        zc(l,body) = (zf(p1,body) + zf(p2,body) + zf(p3,body)) / 3.d0


!        WRITE(6,'(a,I12,3F12.3)') 'l,xc,yc,zc:',l,xc(l),yc(l),zc(l)

        ! Length of the size of the triangle:
        a1 = sqrt((xf(p2,body) - xf(p1,body))**2 + (yf(p2,body) - yf(p1,body))**2 + (zf(p2,body) - zf(p1,body))**2)
        a2 = sqrt((xf(p3,body) - xf(p2,body))**2 + (yf(p3,body) - yf(p2,body))**2 + (zf(p3,body) - zf(p2,body))**2)
        a3 = sqrt((xf(p1,body) - xf(p3,body))**2 + (yf(p1,body) - yf(p3,body))**2 + (zf(p1,body) - zf(p3,body))**2)

        ! Calculate half the perimeter of the triangle:
        s = (a1 + a2 + a3) / 2.0
    
        ! Calculate the area of the triangle using Heron's formula:
        area(l,body) = sqrt(s * (s - a1) * (s - a2) * (s - a3))

        ! Calculate the vectors from the centroid to the vertices:
        vx1 = xf(p2,body) - xf(p1,body)
        vy1 = yf(p2,body) - yf(p1,body)
        vz1 = zf(p2,body) - zf(p1,body)
      
        vx2 = xf(p3,body) - xf(p1,body)
        vy2 = yf(p3,body) - yf(p1,body)
        vz2 = zf(p3,body) - zf(p1,body)

        ! Calculate the components of the normal vector:
        nx(l,body) = vy1*vz2 - vy2*vz1
        ny(l,body) = vz1*vx2 - vz2*vx1
        nz(l,body) = vx1*vy2 - vx2*vy1
        norm = sqrt(nx(l,body)**2 + ny(l,body)**2 + nz(l,body)**2)
        
        
        ! Calculate projected area
        
        ! X-axis:
        cos_thetaX(l,body)=abs(nx(l,body))/sqrt(nx(l,body)**2+ny(l,body)**2+nz(l,body)**2)
        Px_area(l,body)=0.5*sqrt(nx(l,body)**2+ny(l,body)**2+nz(l,body)**2)*cos_thetaX(l,body)

        ! Y-axis:
        cos_thetaY(l,body)=abs(ny(l,body))/sqrt(nx(l,body)**2+ny(l,body)**2+nz(l,body)**2)
        Py_area(l,body)=0.5*sqrt(nx(l,body)**2+ny(l,body)**2+nz(l,body)**2)*cos_thetaY(l,body)

        ! Z-axis:
        cos_thetaZ(l,body)=abs(nz(l,body))/sqrt(nx(l,body)**2+ny(l,body)**2+nz(l,body)**2)
        Pz_area(l,body)=0.5*sqrt(nx(l,body)**2+ny(l,body)**2+nz(l,body)**2)*cos_thetaZ(l,body)

        nx(l,body) = xc(l,body) + nx(l,body)/norm*dx
        ny(l,body) = yc(l,body) + ny(l,body)/norm*dy
        nz(l,body) = zc(l,body) + nz(l,body)/norm*dz

      
        ! Calculate the tangential vector of the centroid:
        tx(l,body) = (xf(p2,body) - xf(p1,body)) + (xf(p3,body) - xf(p1,body))
        ty(l,body) = (yf(p2,body) - yf(p1,body)) + (yf(p3,body) - yf(p1,body))
        tz(l,body) = (zf(p2,body) - zf(p1,body)) + (zf(p3,body) - zf(p1,body))

        norm = sqrt(tx(l,body)**2 + ty(l,body)**2 + tz(l,body)**2)

        tx(l,body) = xc(l,body) + tx(l,body)/norm*dx!*xc(l)
        ty(l,body) = yc(l,body) + ty(l,body)/norm*dy!*yc(l)
        tz(l,body) = zc(l,body) + tz(l,body)/norm*dz!*zc(l)
        
        ! Calculate the normalised weight of each triangle point:
!        w(l,1,body) = sqrt((xc(l,body) - xf(p1,body))**2 + (yc(l,body) - yf(p1,body))**2 + (zc(l,body) - zf(p1,body))**2)
!        w(l,2,body) = sqrt((xc(l,body) - xf(p2,body))**2 + (yc(l,body) - yf(p2,body))**2 + (zc(l,body) - zf(p2,body))**2)
!        w(l,3,body) = sqrt((xc(l,body) - xf(p3,body))**2 + (yc(l,body) - yf(p3,body))**2 + (zc(l,body) - zf(p3,body))**2)
        
!        wtot=sum(w)
!        w(l,:,body) = (wtot-w(l,:,body))/(2*wtot)

!        WRITE(6,'(a,I12,4E14.4,3F9.3)') 'l,nx,ny,nz,area,weight:',l,tx(l),ty(l),tz(l),area(l),w(l,1:3)
      end do

      END DO ! nbody
    
      

!      fn='TFE_Normal_'//trim(adjustl(filename(:sn-4)))//'.dat'
!      open(200,file=fn)
!        write(200,'(a)') 'Variables=x,y,z'
!        do l=1,fileelements(body)
!          write(chl,'(I12)') l 
!          write(200,'(a)') 'ZONE T=N'//trim(adjustl(chl))
!          write(200,'(3E14.4)') xc(l,body),yc(l,body),zc(l,body)
!          write(200,'(3E14.4)') nx(l,body),ny(l,body),nz(l,body)
!        end do
!      close(200)
      
!      fn='TFE_Tangential_'//trim(adjustl(filename(:sn-4)))//'.dat'
!      open(200,file=fn)
!        write(200,'(a)') 'Variables=x,y,z'
!        do l=1,fileelements(body)
!          write(chl,'(I12)') l 
!          write(200,'(a)') 'ZONE T=T'//trim(adjustl(chl))
!          write(200,'(3E14.4)') xc(l,body),yc(l,body),zc(l,body)
!          write(200,'(3E14.4)') tx(l,body),ty(l,body),tz(l,body)
!        end do
!      close(200)

      
      fn=trim(adjustl(fileend))
      open(200,file=fn)
        write(200,'(a)') 'Variables=x,y,z,area,Px_area,Py_area,Pz_area,Px_cos,Py_cos,Pz_cos'
        DO body=1,nbody
        do l=1,fileelements(body)
          write(chl,'(I12)') l 
          write(200,'(a)') 'ZONE T=Element'//trim(adjustl(chl))//', NODES=3, ELEMENTS=1, ZONETYPE=FETRIANGLE, DATAPACKING=POINT'     
          
          p1=c(l,1,body) ; p2=c(l,2,body) ; p3=c(l,3,body)
          write(200,'(20E14.4)') xf(p1,body),yf(p1,body),zf(p1,body),area(l,body),   &
     &     Px_area(l,body),Py_area(l,body),Pz_area(l,body),cos_thetaX(l,body),cos_thetaY(l,body),cos_thetaZ(l,body)
          write(200,'(20E14.4)') xf(p2,body),yf(p2,body),zf(p2,body),area(l,body),   &
     &     Px_area(l,body),Py_area(l,body),Pz_area(l,body),cos_thetaX(l,body),cos_thetaY(l,body),cos_thetaZ(l,body)
          write(200,'(20E14.4)') xf(p3,body),yf(p3,body),zf(p3,body),area(l,body),   &
     &     Px_area(l,body),Py_area(l,body),Pz_area(l,body),cos_thetaX(l,body),cos_thetaY(l,body),cos_thetaZ(l,body)

          write(200,*) 

          ll=1
          write(200,'(3I14)') ll,ll+1,ll+2
          
        end do
        END DO
      close(200)
      WRITE(6,'(a)') ' '//trim(adjustl(fn))//'    exported'

      fn='CloudPoints_'//trim(adjustl(fileall))
      open(200,file=fn)
        write(200,'(a)') 'Variables=x,y,z,area,Px_area,Py_area,Pz_area,element'
        DO body=1,nbody
        write(200,'(a)') 'ZONE T='//trim(adjustl(filename(body)))
        do l=1,fileelements(body)
          write(200,'(7E14.4,I9)') nx(l,body),ny(l,body),nz(l,body),area(l,body),Px_area(l,body),Py_area(l,body),Pz_area(l,body),l
        end do
        END DO
      close(200)
      WRITE(6,'(a)') ' '//trim(adjustl(fn))//'    exported'
      WRITE(6,*) 
      
      AreaProjected=0
      do l=1,fileelements(1)
        AreaProjected=AreaProjected+Px_area(l,1)
      end do
      WRITE(6,'(a,F9.6,F9.3)') ' Projected Area Head:',AreaProjected,AreaProjected/0.112421*100

      AreaProjected=0
      do l=1,fileelements(2)
        AreaProjected=AreaProjected+Px_area(l,2)
      end do
      WRITE(6,'(a,F9.6,F9.3)') ' Projected Area Back:',AreaProjected,AreaProjected/0.112421*100


      end 
