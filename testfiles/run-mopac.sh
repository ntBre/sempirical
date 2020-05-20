#!/bin/tcsh
#module load mopac
if ($#argv != 4 ) then 
  echo 'Usage: ./run-structures.sh <stem> <first #> <final #> <batch size>'
else
  set stems = $1
  set i     = $2
  set j     = $3
  set kmax  = $4
  set pwd   = `pwd`
@ j = $j + 1
  while ($i < $j)
   set k = $i
   set kcount = 0
   while ($kcount < $kmax)
  set  x = `printf "%05d" $k`
  sed -e "s@example@$stems$x@" ex.sh | sbatch 
#   mopac $stems$x.mop >& /dev/null & 
@ k = $k + 1
@ kcount = $kcount + 1
    end
@ i = $i + $kmax
   wait
  end
endif
