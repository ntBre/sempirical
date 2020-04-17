      subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
      integer m,n,ldfjac,iflag
      double precision epsfcn
      double precision x(n),fvec(m),fjac(ldfjac,n),wa(m)
      integer i,j
      double precision eps,epsmch,h,temp,zero
      double precision dpmpar
      data zero /0.0d0/
      epsmch = dpmpar(1)
      eps = dsqrt(dmax1(epsfcn,epsmch))
      do 20 j = 1, n
         temp = x(j)
         h = eps*dabs(temp)
         if (h .eq. zero) h = eps
         x(j) = temp + h
         call fcn(m,n,x,wa,j)
         if (iflag .lt. 0) go to 30
         x(j) = temp
         do 10 i = 1, m
            fjac(i,j) = (wa(i) - fvec(i))/h
   10       continue
   20    continue
   30 continue
      return
      end
