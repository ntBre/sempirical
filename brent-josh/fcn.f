      SUBROUTINE FCN(M,N,X,FVEC)
      USE PARMS
      INTEGER nstruc ,M,N,I, istat, ngeom, eRnI, eRnF, nfc,
     + JJ
      DOUBLE PRECISION FVEC(M), RMSD, geomPr, fcpr,x(n)
      DOUBLE PRECISION E1(M), kcal2cm, NumX
      double precision enemin(nmolec)
      double precision, allocatable :: etemp(:)
      double precision, allocatable:: xyzc(:,:) 
      CHARACTER xx, adummy, readenergy*35,fileName1*40, optfile*40
      character readgeometry*50
      allocate(xyzc(maxatoms,3),etemp(maxval(nstruct)))    
      CALL FLUSH(6)
      OPEN(4,FILE='params.dat')
      DO I=1,N
       WRITE(4,'(A15,1X,A3,2x, F20.15)') APARM(I), ASYM(I), X(I)
      END DO
      REWIND(4)
      CLOSE(4)
          CALL SYSTEM('./run-1step.sh >2 /dev/null')
          CALL SYSTEM('python read-precise.py 60000 > ENERGIES.DAT')
          OPEN(654,file='ENERGIES.DAT')
       RMSD=0.0
         jj=0
      do i=1, nmolec
        do j=1, nstruct(i)
           jj=jj+1
               read(654,*) e1(jj)
               etemp(j)=e1(jj)
         enddo
         do j=jj-nstruct(i)+1, jj
           e1(j)=(e1(j)-minval(etemp))*8065.5401069
           fvec(j)=(EAB(j)-E1(j))
           RMSD=RMSD+fvec(j)**2
         enddo
      enddo
      rewind(654)
      close(654)
      CALL SYSTEM('rm ENERGIES.DAT')
       do i=1,nmolec
          write (optfile,'(A28, i0.2, A4)')"structureFolder/StructureOpt
     +", i, ".out"
          open(22,file=optfile)
        do 
         read (22,'(A50)', iostat = istat) readgeometry
          if (istat .ne. 0) exit
          if (readgeometry(28:50).eq.'  CARTESIAN COORDINATES') then
          do j=1, natoms(i)
            read (22,*) adummy, adummy, (xyzc(j,k),k=1,3)
           enddo
          endif
        enddo
       close(22)
         do j = 1, ngeoms(i)
          jj=jj+1
          if(igeoms(i,j,1).eq.2) then
            dist=0.
            jone=igeoms(i,j,2)
            jtwo=igeoms(i,j,3)
            do k=1,3
            dist=dist+(xyzc(jone,k)-xyzc(jtwo,k))**2
            enddo
            e1(jj)=sqrt(dist)
          endif
          if(igeoms(i,j,1).eq.3) then
            dist1=0.
            dist2=0.
            dist3=0.
            jone=igeoms(i,j,2)
            jtwo=igeoms(i,j,3)
            jthr=igeoms(i,j,4)
            do k=1,3
          dist1=dist1+(xyzc(jone,k)-xyzc(jtwo,k))**2
          dist2=dist2+(xyzc(jtwo,k)-xyzc(jthr,k))**2
          dist3=dist3+(xyzc(jone,k)-xyzc(jthr,k))**2
            enddo
          ang=(dist3-dist2-dist1)/(-2.*sqrt(dist1)*sqrt(dist2)) 
            e1(jj)=acos(ang)*180.0/acos(-1.0)
          endif
         enddo
       enddo
       rewind(76)
       jj=0
       do i=1,nmolec
          do j=1,nstruct(i)
            jj=jj+1
           write(76,'(I6,3(F9.5,2x))') jj,eab(jj), e1(jj),eab(jj)-e1(jj)
          enddo
       enddo
         jj=sum(nstruct)
      write(6,*)"Differences in geometries"
      write(6,'(A20,A15,A12,A13)')"        Type", "Ab. Ini.", 
     + "Calc.", "Difference"
      do i=1,nmolec
         if(ngeoms(i).gt.0)    write(6,'(A10,I2)') 'Molecule #',i
         do j=1,ngeoms(i)
           jj=jj+1
           if(igeoms(i,j,1).eq.2) then
           write(6,'(A4,10x,2(I2,1x),3x,3(f11.3,1x))')'Bond',
     .      igeoms(i,j,2),igeoms(i,j,3),eab(jj),e1(jj),eab(jj)-e1(jj)
                fvec(jj)=eab(jj)-e1(jj)
           elseif (igeoms(i,j,1).eq.3) then
           write(6,'(A5,9x,3(I2,1x),3(f11.3,1x))')'Angle',igeoms(i,j,2),
     .          igeoms(i,j,3),igeoms(i,j,4),
     .          eab(jj), e1(jj), eab(jj)-e1(jj)
                fvec(jj)=eab(jj)-e1(jj)
           endif
         enddo
         write(6,*) 
       enddo
      WRITE(6,'(A4,1x,F17.12)') 'RMSD',SQRT(RMSD/REAL(m-sum(ngeoms)))
      WRITE(6,*)
      CALL FLUSH(6)
      do i=1, nweights
         jj=0
         if (iweight(i,1).eq.1) then 
            do j=1,nmolec
               do k=1,nstruct(j) 
               jj=jj+1
                 if(j.eq.iweight(i,2)) then
                   fvec(jj) = fvec(jj)*weight(i)
                 endif
               enddo
            enddo
         else if (iweight(i,1).eq.2) then
            do j=1,nmolec
              do k=1,nstruct(j)
                jj=jj+1
                if(jj.ge.iweight(i,2).and.jj.le.iweight(i,3)) then
                  fvec(jj)=fvec(jj)*weight(i)
                endif
              enddo
            enddo     
         else if (iweight(i,1).eq.3) then
           do j=1,nmolec
             do k=1,nstruct(j)
               jj=jj+1
               do kk=3,iweight(i,2)+2
                 if(jj.eq.iweight(i,kk)) then
                   fvec(jj) = fvec(jj)*weight(i)
                 endif 
               enddo
             enddo
           enddo
         endif
       enddo
       rewind(88)
       do jj=1,sum(nstruct)+sum(ngeoms) 
       write(88,*) fvec(jj)
       enddo
       rmsdstep = SQRT(RMSD/REAL(m-sum(ngeoms)))
      if(rmsdstep.lt.RMSDMIN) then
      RMSDMIN = rmsdstep
      open (25, file="bestparams.out")
        write(25,*) 'RMSD= ', rmsdstep
        write(25,*) 
        do i=1, N
         write (25,*) APARM(i), ASYM(I), X(i)
        enddo
       write (25,*) 
      close(25)
      endif
      deallocate (xyzc,etemp)     
      write(6,*) 'i get out of fcn'
      RETURN
      END
