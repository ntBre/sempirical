      MODULE PARMS
      character, allocatable ::  APARM(:)*7, ASYM(:)*2, atsym(:,:)*2
      double precision, allocatable :: eab(:), eabmin(:), weight(:)
      integer, allocatable :: natoms(:), ngeoms(:), nstruct(:)
      integer, allocatable :: igeoms(:,:,:), iweight(:,:)
      integer nmolec, maxatoms, nweights
      double precision rmsdmin
      END MODULE

!
!   This code optimizes parameters
!   It calls MOPAC SRP-RM1 with some parameters,
!   examines the RMSD of the frequencies calculated with respect to abinitio
!   methods and minimizes it calling minpack routines.
! 
!
!   Diego Troya 27 July 2004
!   Diego Troya 11 January 2007. Adapted to Ar+C2F6
!   Diego Troya & Josh Layfield January 31 2007
!   Josh Layfield May 20 2012 Update for Freqs
!   Joe Arend and Josh Layfield September 2016
      program main
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


!     COMMON /IRE/ IFLAGREAD
!     COMMON /ABINI/ EAB 
!     COMMON /SYMS/ ASYM, APARM, NATOMS, JJ
!
!     
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

!
!    Add weighting factors to different regions of the fit
!    The first line tells you how many different weights (nweights)
!    The next (nweights) lines are structured as follows.  
!     First column: 1: The weight is applied to a specific molecule
!     First column: 2: The weight is applied to a range of points  
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
 
!
!   Read the initial parameters, which go to the array X
!
      OPEN(72,FILE='params.initial')
        do i=1,n 
        READ(72,*) aparm(i), asym(i), x(i)
        END DO
      CLOSE(72)

!
!     Build run-1step.sh
!
      open(55,file='run-1step.sh')
      write(55,*)'#!/bin/tcsh'
      write (55,'(A43,i6,a19)')"./run-mopac.sh structureFolder/Structure
     . 1 ",sum(nstruct), " 16"




!
!   Build the structure files 
!
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

!
!     Build the EAB array - Energy AB initio
!     This array contains all the "true" values that we are using to
!     train the method
!    
      E1size=Sum(ngeoms)+sum(nstruct) ! add energy values
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
         do j = jj-nstruct(i)+1 , jj ! make energy relative 
            eab(j)=(eab(j)-eabmin(i))*219474.63 ! convert to wavenumbers ?!
            write(89,*) eab(j)
         enddo
      enddo
!
      do i = 1, nmolec ! add the geometry values
         do j = 1, ngeoms(i)
          jj=jj+1
          if(igeoms(i,j,1).eq.2) then
            dist=0.
            do k=1,3
            dist=dist+(xyz(i,igeoms(i,j,2),k)-xyz(i,igeoms(i,j,3),k))**2 ! handling weights?
            enddo
            eab(jj)=sqrt(dist)*au2ang
          endif
!
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
!
! FCN = ??
! m = # geoms + # structs, n = # parameters
! X = parameter array
! FVEC = ??
! TOL = tolerance, defined above
! INFO =  ??
! IWA = array same size as X
! WA = array size LWA
! LWA = 1153000 from above but idk why
       CALL LMDIF1(FCN,M,N,X,FVEC,TOL,INFO,IWA,WA,LWA)
!      
       WRITE(6,*) 'Final parameters'
       do i=1,n
       write(6,*) aparm(i),asym(i),x(i)
       enddo
       WRITE(6,*) 'INFO',info

!       open(91, file = 'structureFolder/Structureopt.out')
!        do
!         read (91,'(A51)', iostat = istat) readgeometry
!            if (istat .ne. 0) exit
!          if (readgeometry(30:50) .eq. 'CARTESIAN COORDINATES') then
!           do i=1, NATOMS
!            read (91,*) adummy, ASYM(i), xc(i), yc(i), zc(i)
!           enddo
!             exit
!          endif
!        enddo
!       close(91)
 
c Finding the moments of inertia and rotational constants
!     do i=1, natoms
!      if (ASYM(i) .eq. "C" .or. ASYM(i) .eq. "c") Atm(i) = 12.01
!      if (ASYM(i) .eq. "H" .or. ASYM(i) .eq. "h") Atm(i) = 1.008
c         write(*,*) Atm(i)
!     enddo

c Finding the total mass of the molecule
!     tm = 0.
!      do i=1, natoms
!       tm = tm + atm(i)
!      enddo
!     write (*,*)
!     write (*,*) 'Total mass of the molecule = ', tm

c finding the moments of inertia Diagonal terms
!       ixx = 0.
!       ixy = 0. 
!       ixz = 0.
!       iyx = 0.
!       iyy = 0.
!       iyz = 0.
!       izx = 0.
!       izy = 0.
!       izz = 0.
!     do i=1, natoms
!        ixx = (atm(i)*(yc(i)**2+zc(i)**2))+ixx
!        iyy = (atm(i)*(xc(i)**2+zc(i)**2))+iyy
!        izz = (atm(i)*(xc(i)**2+yc(i)**2))+izz
!     enddo

c finding the moments of inertia off-diagonal terms
!     do i=1, natoms
!        ixy = (-atm(i)*xc(i)*yc(i))+ixy
!        ixz = (-atm(i)*xc(i)*zc(i))+ixz
!        iyx = (-atm(i)*xc(i)*yc(i))+iyx
!        iyz = (-atm(i)*yc(i)*zc(i))+iyz
!        izx = (-atm(i)*xc(i)*zc(i))+izx
!        izy = (-atm(i)*yc(i)*zc(i))+izy
!     enddo

c Putting moment of inertia into i matrix
!     im(1,1) = ixx
!     im(1,2) = iyx
!     im(1,3) = izx
!     im(2,1) = ixy
!     im(2,2) = iyy
!     im(2,3) = izy
!     im(3,1) = ixz
!     im(3,2) = iyz
!     im(3,3) = izz

!     CALL PRINT_MATRIX( 'Moment of Inertia Tensor',MIMS,MIMS,im, LDA)

!     LWORK = -1
!     CALL DSYEV( 'Vectors', 'Lower', MIMS, im,
!    .       LDA, W, WORK, LWORK, INFORM )
!     LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!     write(*,*) 'LWORK=',LWORK

*
*     Solve eigenproblem.
*
!     CALL DSYEV( 'Vectors', 'Lower', MIMS, im,
!    .   LDA, W, WORK, LWORK, INFORM )
*
*     Check for convergence.

!     write(*,*) 'INFO =', INFORM
!     IF( INFORM.GT.0 ) THEN
!        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
!        STOP
!     END IF
*
*     Print eigenvalues.
!     CALL PRINT_MATRIX( 'Eigenvalues', 1, MIMS, W, 1 )
*     Print eigenvectors.
c      CALL PRINT_MATRIX( 'Eigenvectors(stored columnwise)',MIMS,MIMS,im,
c     $                   LDA )
c      STOP
c      write (*,*)
c      write (*,*) 'Energy in joules for each rotation'


!       do i=1, 3
c        Ej(i) = ((6.626070040E-34)**2)/
c     a(8*acos(-1.0)*w(i)*1.60539040E-27*1E-20)

cc This is where the rotational wavenumber calculation needs to be fixed
!         Ej(i) = 3.46371101E-22/W(i)
c        write (*,*) Ej(i)
!       enddo

c     write (*,*)
c     write (*,*) 'Energy in wave numbers for each rotation (cm^-1)'
!      do i=1, 3
!       Ewn(i) = Ej(i)*(5.03445E22)
!        write (*,*) Ewn(i)
!       enddo
!      write (*,*)


c This takes the frequencies from inteder-new.out and puts it into 
c the readable form of allfreqs.out

c     CALL system('./grabfreq.sh')
      
c This will read the freqs that intder puts out
c Eventually, it will also write the freqs that spectro puts out

c     open(35, file="allfreqs.out")
c      read (35,*) HarmonicFrequencies, adummy, adummy
c      write (*,*) HarmonicFrequencies
c       do i=1,6
c        read (35,*) hfreq(i)
c        write (*,*) hfreq(i)
c       enddo
c      read (35,*) FundamentalFrequencies
c      write (*,*) FundamentalFrequencies
c     close (35)

      deallocate (X,IWA,WA,asym,aparm)
      deallocate (natoms,etemp)
      deallocate (xyz,filename,charges)
      deallocate (atsym,ngeoms,filegeom)
      deallocate (igeoms,nstruct)
      deallocate(eab,eabmin,fvec)
c     deallocate(iweight,weights)
            

      END



c     Auxiliary routine: printing a matrix.



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

