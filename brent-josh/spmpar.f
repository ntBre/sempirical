      real function spmpar(i)
      integer i
      integer mcheps(2)
      integer minmag(2)
      integer maxmag(2)
      real rmach(3)
      equivalence (rmach(1),mcheps(1))
      equivalence (rmach(2),minmag(1))
      equivalence (rmach(3),maxmag(1))
      data rmach(1) /1.192091E-07/
      data rmach(2) /1.175495E-38/
      data rmach(3) /3.402822E+38/
      spmpar = rmach(i)
      return
      end
