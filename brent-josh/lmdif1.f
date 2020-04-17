      subroutine lmdif1(fcn,m,n,x,fvec,tol,info
     .,iwa,wa,lwa)
      integer m,n,info,lwa
      integer iwa(n)
      double precision tol
      double precision x(n),fvec(m),wa(lwa)
      external fcn
      integer maxfev,mode,mp5n,nfev,nprint
      double precision epsfcn,factor,ftol,gtol,xtol,zero
      data factor,zero /1.0d-1,0.0d0/
      info = 0
      WRITE(6,*) N,M,TOL,LWA, m*n + 5*n + m
      if (n .le. 0 .or. m .lt. n .or. tol .lt. zero
     *    .or. lwa .lt. m*n + 5*n + m) go to 10
      maxfev = 20000*(n + 1)
      ftol = tol
      xtol = tol
      gtol = zero
      epsfcn = zero 
      mode = 1
      nprint = 0
      mp5n = m + 5*n
      call lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,wa(1),
     1           mode,factor,nprint,info,nfev,wa(mp5n+1),m,iwa,
     2           wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
      if (info .eq. 8) info = 4
   10 continue
      return
      end
