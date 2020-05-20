package main

import (
	"fmt"
	"gonum.org/v1/gonum/mat"
	"io/ioutil"
	"os"
	"strings"
)

const (
	energyLine      = "TOTAL_ENERGY:EV="
	mopacTerminated = "not sure what this is yet"
	EVperHt         = 27.211385 // from http://www.ilpi.com/msds/ref/energyunits.html
)

var (
	Input [Nkeys]string
)

func MakeInp() {
	if _, err := os.Stat("inp/"); os.IsNotExist(err) {
		os.Mkdir("inp", 0755)
	} else {
		os.RemoveAll("inp/")
		os.Mkdir("inp", 0755)
	}
}

func WriteParams() {
	ioutil.WriteFile("params.dat", []byte(Input[Params]), 0755)
}

func main() {
	// read input
	ParseInfile("new.inp")
	geoms := ReadGfile(Input[Gfile])
	names := GetAtomNames()
	mop := Mopac{}
	queue := Slurm{}
	toRead := make([]string, 0)
	// write mopac input files
	WriteParams() // need to update this on each iteration
	MakeInp()
	for i, geom := range geoms {
		filename := fmt.Sprintf("inp/Structure%05d", i)
		toRead = append(toRead, filename)
		mop.WriteIn(filename+".mop", names, strings.Fields(geom))
		// run mopac
		queue.Write(filename+".pbs", filename+".mop")
		queue.Submit(filename + ".pbs")
	}
	// gather energies
	// expect -153.486630272141 from CCSD(T)-F12/cc-pVDZ-F12 calculation
	energies := make([]float64, 0)
	for i, file := range toRead {
		energy, err := mop.ReadOut(file + ".aux")
		if err != nil {
			fmt.Fprintln(os.Stderr, i)
			panic(err)
		}
		energies = append(energies, energy)
	}
	eVec := mat.NewVecDense(len(energies), energies)
	eVec.ScaleVec(1/EVperHt, eVec)
	fmt.Println(eVec)
/* Jacobian
let DEi = | (ab initio Ei) - (semi empirical Ei) |
let Pj = parameter j
J = [dDEi/dPj]
J = [dDE1/dP1 ... dDE1/dPn
      ...      .   ...
     dDEm/dP1 ... dDEm/dPn]
so to calculate the jacobian matrix, step each parameter Pj,
calculate the energy change in each DE, and that gives a column of J
calculate Jij by central differences:
1. calcule DEi(pj + delta), DEi(pj - delta)
2. dDEi = DEi(pj+delta) - DEi(pj-delta) / (2*delta)
delta P is given by some other delta * (1+|pj|), scale delta by the size of pj
want to use Broyden's method to not evaluate this jacobian all the time
*/
	// adjust parameters
}
