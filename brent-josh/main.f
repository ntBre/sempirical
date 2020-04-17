      MODULE PARMS
      character, allocatable ::  APARM(:)*7, ASYM(:)*2, atsym(:,:)*2
      double precision, allocatable :: eab(:), eabmin(:), weight(:)
      integer, allocatable :: natoms(:), ngeoms(:), nstruct(:)
      integer, allocatable :: igeoms(:,:,:), iweight(:,:)
      integer nmolec, maxatoms, nweights
      double precision rmsdmin
      END MODULE

      PROGRAM main
      use parms
      implicit double precision(a-h,o-z)
      integer info, n, mims, lda, lwmax, inform, lwork
      integer jj, fvecsize, e1size, ia0
      parameter (mims=3)
      parameter (lda=mims)
      parameter (lwmax=100000)
      integer, allocatable ::  iwa(:), charges(:)
      double precision coord(3) 
      double precision, allocatable :: etemp(:)
      double precision, allocatable :: X(:)    
      double precision, allocatable :: fvec(:),wa(:)
      double precision, allocatable :: Atm(:), xyz(:,:,:)
      double precision W(MIMS), WORK(LWMAX)
      character, allocatable ::  filegeom(:)*30 
      character, allocatable :: filename(:)*20 
      character :: filenames*50, method*5, optfile*40,string*60,adummy
      EXTERNAL FCN
      ! fcn.f

      IFLAGREAD=0 ! Flag to read abinitioenergies
      TOL = 0.000001 ! Tolerance for stopping (probably never get there).
      LWA = 1153000 ! 
      au2ang = 0.529177
      rmsdmin=10000.
      !!!  reading from main input file
      read(5,*) n   ! The number of parameters to be optimized
      read(5,*) method
      read(5,*) nmolec ! The number of molecules to be studied.
      allocate(X(N), IWA(N), WA(LWA),asym(n),aparm(n))
      allocate (natoms(nmolec))
      read(5,*) (natoms(i),i=1,nmolec)
      maxatoms = maxval(natoms)
      allocate(xyz(nmolec,maxatoms,3),filename(nmolec),charges(nmolec))
      allocate (atsym(nmolec,maxatoms),ngeoms(nmolec),filegeom(nmolec))
      allocate (igeoms(nmolec,10,5),nstruct(nmolec))
      do i = 1,nmolec
         read(5,*) charges(i)
         read(5,*) nstruct(i)
         read(5,*) filename(i)
         read(5,*) filegeom(i), ia0
         do j = 1, natoms(i)
            read(5,*) atsym(i,j),(xyz(i,j,k),k=1,3)
         enddo
         if(ia0.eq.0) then
           xyz=xyz/au2ang
         endif 
         read(5,*) ngeoms(i)
         do j = 1, ngeoms(i)
            read (5,*) idimen
            backspace(5)
            read(5,*) (igeoms(i,j,k),k=1,idimen+1)
         enddo
      enddo
        read(5,*) nweights
        allocate(iweight(nweights,sum(nstruct)), weight(nweights))
        iweight=0
        do i=1, nweights
          read(5,*) idum
          if(idum.eq.2) then
            backspace (5)
            read(5,*) (iweight(i,j),j=1,3), weight(i)
          elseif (idum.eq.1) then
            backspace(5)
            read(5,*) (iweight(i,j),j=1,2), weight(i)
         elseif (idum.eq.3) then 
            backspace(5)
            read(5,*) (iweight(i,j),j=1,2) 
            backspace(5)
            nats=iweight(i,2)
            write(6,*) 'nats',nats
            read(5,*) (iweight(i,j),j=1,nats+2), weight(i)
            write(6,*) 
         else
            write(6,*) 'Inappropriate weighting factor'
          endif
        enddo
      OPEN(72,FILE='params.initial')
        do i=1,n 
        READ(72,*) aparm(i), asym(i), x(i)
        END DO
      CLOSE(72)
      open(55,file='run-1step.sh')
      write(55,*)'#!/bin/tcsh'
      write (55,'(A43,i6,a19)')"./run-mopac.sh structureFolder/Structure
     . 1 ",sum(nstruct), " 16"
      nstructnum=1
      do i=1,nmolec
        open (20,file=filegeom(i))
        do j=1, nstruct(i)
         write (fileNames,'(A25, i0.5, A4)')"structureFolder/Structure",
     + nstructnum, ".mop"
          nstructnum=nstructnum+1
          open (21, file = fileNames)
          read (20,*) adummy
          write (21,'(a82,i5,a1,a4)')'threads=1 XYZ A0 scfcrt=1.D-21 
     .aux(precision=9) external=params.dat 1SCF charge=',
     .charges(i),' ',method
          write (21,*) "MOLECULE #", i
          write (21,*)
          if(j.eq.1) then
          write (optfile,'(A28, i0.2, A4)')"structureFolder/StructureOpt
     +", i, ".mop"
          write (55,*) "mopac ", optfile, " 2> /dev/null"
          open(22,file=optfile)
          write (22,'(a80,i5,a1,a4)')'threads=1 XYZ A0 scfcrt=1.D-21 
     .aux(precision=9) external=params.dat charge=',
     .charges(i),' ',method
          write (22,*) "OPT MOLECULE #", i
          write (22,*)
         endif 
          do jj=1, natoms(i)
            read (20,*) (Coord(k), k=1, 3)
            if(ia0.eq.0) then 
            coord= coord/au2ang
            endif
            write (21,*) atsym(i,jj), (coord(k), k=1,3) 
            if(j.eq.1) write (22,*) atsym(i,jj), (coord(k), k=1,3) 
          enddo 
          close (21)
        enddo 
        close (20)
      enddo
      close(55)
      E1size=Sum(ngeoms)+sum(nstruct)
      m=e1size
      allocate(eab(m),eabmin(nmolec),fvec(m),etemp(maxval(nstruct)))
      jj=0
      do i = 1, nmolec
         open(1,file=filename(i))
         do j = 1, nstruct(i)
         jj=jj+1
         read(1,*)   eab(jj)
         etemp(j)= eab(jj)
         enddo
         eabmin(i)=minval(etemp)
         close(1)
         do j = jj-nstruct(i)+1 , jj
            eab(j)=(eab(j)-eabmin(i))*219474.63
            write(89,*) eab(j)
         enddo
      enddo
      do i = 1, nmolec
         do j = 1, ngeoms(i)
          jj=jj+1
          if(igeoms(i,j,1).eq.2) then
            dist=0.
            do k=1,3
            dist=dist+(xyz(i,igeoms(i,j,2),k)-xyz(i,igeoms(i,j,3),k))**2
            enddo
            eab(jj)=sqrt(dist)*au2ang
          endif
          if(igeoms(i,j,1).eq.3) then
            dist1=0.
            dist2=0.
            dist3=0.
            do k=1,3
          dist1=dist1+(xyz(i,igeoms(i,j,2),k)-xyz(i,igeoms(i,j,3),k))**2
          dist2=dist2+(xyz(i,igeoms(i,j,3),k)-xyz(i,igeoms(i,j,4),k))**2
          dist3=dist3+(xyz(i,igeoms(i,j,2),k)-xyz(i,igeoms(i,j,4),k))**2
            enddo
          ang=(dist3-dist2-dist1)/(-2.*sqrt(dist1)*sqrt(dist2)) 
            eab(jj)=acos(ang)*180.0/acos(-1.0)
          endif
         enddo
      enddo      
      JJ = 1
       CALL LMDIF1(FCN,M,N,X,FVEC,TOL,INFO,IWA,WA,LWA)
       WRITE(6,*) 'Final parameters'
       do i=1,n
       write(6,*) aparm(i),asym(i),x(i)
       enddo
       WRITE(6,*) 'INFO',info
      deallocate (X,IWA,WA,asym,aparm)
      deallocate (natoms,etemp)
      deallocate (xyz,filename,charges)
      deallocate (atsym,ngeoms,filegeom)
      deallocate (igeoms,nstruct)
      deallocate(eab,eabmin,fvec)
      END
      SUBROUTINE PRINT_MATRIX(DESC, M, N, A, LDA)
      Character*(*) DESC
      Integer M, N, LDA
      Double Precision A(LDA, *)
      Integer I, J
      write(*,*)
      write(*,*) DESC
      do I = 1, M
       write(*,9998) (A(I, J), J=1, N)
      enddo
 9998 format(11(:,1X,F27.12))
      return
      end
