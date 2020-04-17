      subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
      integer n,ldr
      integer ipvt(n)
      double precision r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa(n)
      integer i,j,jp1,k,kp1,l,nsing
      double precision cos,cotan,p5,p25,qtbpj,sin,sum,tan,temp,zero
      data p5,p25,zero /5.0d-1,2.5d-1,0.0d0/
      do 20 j = 1, n
         do 10 i = j, n
            r(i,j) = r(j,i)
   10       continue
         x(j) = r(j,j)
         wa(j) = qtb(j)
   20    continue
      do 100 j = 1, n
         l = ipvt(j)
         if (diag(l) .eq. zero) go to 90
         do 30 k = j, n
            sdiag(k) = zero
   30       continue
         sdiag(j) = diag(l)
         qtbpj = zero
         do 80 k = j, n
            if (sdiag(k) .eq. zero) go to 70
            if (dabs(r(k,k)) .ge. dabs(sdiag(k))) go to 40
               cotan = r(k,k)/sdiag(k)
               sin = p5/dsqrt(p25+p25*cotan**2)
               cos = sin*cotan
               go to 50
   40       continue
               tan = sdiag(k)/r(k,k)
               cos = p5/dsqrt(p25+p25*tan**2)
               sin = cos*tan
   50       continue
            r(k,k) = cos*r(k,k) + sin*sdiag(k)
            temp = cos*wa(k) + sin*qtbpj
            qtbpj = -sin*wa(k) + cos*qtbpj
            wa(k) = temp
            kp1 = k + 1
            if (n .lt. kp1) go to 70
            do 60 i = kp1, n
               temp = cos*r(i,k) + sin*sdiag(i)
               sdiag(i) = -sin*r(i,k) + cos*sdiag(i)
               r(i,k) = temp
   60          continue
   70       continue
   80       continue
   90    continue
         sdiag(j) = r(j,j)
         r(j,j) = x(j)
  100    continue
      nsing = n
      do 110 j = 1, n
         if (sdiag(j) .eq. zero .and. nsing .eq. n) nsing = j - 1
         if (nsing .lt. n) wa(j) = zero
  110    continue
      if (nsing .lt. 1) go to 150
      do 140 k = 1, nsing
         j = nsing - k + 1
         sum = zero
         jp1 = j + 1
         if (nsing .lt. jp1) go to 130
         do 120 i = jp1, nsing
            sum = sum + r(i,j)*wa(i)
  120       continue
  130    continue
         wa(j) = (wa(j) - sum)/sdiag(j)
  140    continue
  150 continue
      do 160 j = 1, n
         l = ipvt(j)
         x(l) = wa(j)
  160    continue
      return
      end
