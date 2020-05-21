package main

import (
	"fmt"
	"io/ioutil"
	"os"
	"strings"

	"io"
	"math"

	"strconv"

	"time"

	"gonum.org/v1/gonum/mat"
)

const (
	energyLine      = "TOTAL_ENERGY:EV="
	mopacTerminated = "not sure what this is yet"
	EVtoHt          = 1 / 27.211385 // from http://www.ilpi.com/msds/ref/energyunits.html
	delta           = 1e-8          // roughly the same as from fortran
)

var (
	Input        [Nkeys]string
	mop          = Mopac{}
	queue        = Slurm{}
	jobs         []string
	abInit       []float64
	initParams   []float64
	paramHeaders []string
)

// Create the directory inp if it does not exist
// Delete it and recreate it otherwise
func MakeInp() {
	if _, err := os.Stat("inp/"); os.IsNotExist(err) {
		os.Mkdir("inp", 0755)
	} else {
		os.RemoveAll("inp/")
		os.Mkdir("inp", 0755)
	}
}

// Build the parameters file from the separate parameters and their headers
func MakeParams(params []float64, headers []string) (lines string) {
	for i := range params {
		lines += fmt.Sprintf("%s%20.12f\n", headers[i], params[i])
	}
	return
}

// Write the parameters to the file read by MOPAC
func WriteParams(params []float64, headers []string) {
	lines := MakeParams(params, headers)
	ioutil.WriteFile("params.dat", []byte(lines), 0755)
}

// Write MOPAC input files
func WriteInfiles(geoms, names []string) []string {
	toRead := make([]string, 0)
	for i, geom := range geoms {
		filename := fmt.Sprintf("inp/Structure%05d", i)
		toRead = append(toRead, filename)
		mop.WriteIn(filename+".mop", names, strings.Fields(geom))
		// run mopac
		queue.Write(filename+".pbs", filename+".mop")
	}
	return toRead
}

// Submit the jobs identified in files
func Submit(files []string) {
	// TODO try goroutine on submission
	for _, file := range files {
		queue.Submit(file + ".pbs")
	}
}

// Read energies from MOPAC .aux files
func ReadEnergies(files []string) []float64 {
	energies := make([]float64, 0)
	for _, file := range files {
		// TODO try goroutine on reading
		// closure around
		energy, err := mop.ReadOut(file + ".aux")
		for err != nil {
			// assume problem from reading before energy present
			fmt.Fprintln(os.Stderr, file+".aux", err, "sleeping")
			time.Sleep(time.Second)
			// convert to hartrees here
			energy, err = mop.ReadOut(file + ".aux")
		}
		energies = append(energies, energy*EVtoHt)
		// this, with passing in energies and file
		// would have to use a fully made slice for this and keep track of the
		// file numbers, since they would no longer come out in order
	}
	return energies
}

// Read ab initio energies from filename
func AbEnergies(filename string) (abInit []float64) {
	lines := ReadFile(filename)
	for _, line := range lines {
		if line != "" {
			f, _ := strconv.ParseFloat(line, 64)
			abInit = append(abInit, f)
		}
	}
	return
}

func init() {
	ParseInfile("new.inp")
	MakeInp()
	jobs = WriteInfiles(ReadGfile(Input[Gfile]), GetAtomNames())
	abInit = AbEnergies(Input[Efile])
	initParams, paramHeaders = FloatParams()
}

func SaveJac(jac *mat.Dense) {
	f, err := os.Create("jac.dat")
	if err != nil {
		panic(err)
	}
	PlotData(len(jobs), len(Input[Params]), jac, f)
	f.Close()
}

func main() {
	nparam := len(initParams)
	nstruct := len(jobs)
	// "convergence criteria"
	maxIter := 1

	// J^t J h = J^T[y-f(B)]
	// where J is the Jacobian, h is the step,
	// y is the exact solution, and f(B) is
	// the evaluation of f at the parameter vector B
	params := initParams
	for i := 0; i < maxIter; i++ {
		// semi-empirical energy column vector
		fB := mat.NewDense(nstruct, 1, SEnergy(params))
		// ab initio column vector
		y := mat.NewDense(nstruct, 1, abInit)
		jac := Jacobian(params)
		A := mat.NewDense(nparam, nparam, nil)
		b := mat.NewDense(1, nstruct, nil)
		// A = jac^T jac + lambda I <- neglecting this for now
		A.Mul(jac.T(), jac)
		// b = jac^T*[y - fB]
		b.Sub(y, fB)
		b.Mul(jac.T(), b)
		step := mat.NewDense(nparam, 1, nil)
		// Solve Ax = b for x, which is our step in the parameters
		step.Solve(A, b)
		// add step to params
		// calculate new energy
	}
}

// Calculate the semi-empirical energy as a function of params
func SEnergy(params []float64) []float64 {
	WriteParams(params, paramHeaders)
	Submit(jobs)
	return ReadEnergies(jobs) // in hartrees
}

// Calculate the Jacobian matrix by central differences
func Jacobian(params []float64) (jac *mat.Dense) {
	// no backing slice => nil
	nstruct := len(jobs)
	nparam := len(params)
	jac = mat.NewDense(nstruct, nparam, nil)
	// for each parameter <=> column in jac
	for col := 0; col < nparam; col++ {
		reset := params[col]
		params[col] += delta + delta*math.Abs(params[col])
		forward := SEnergy(params)
		params[col] = reset
		params[col] -= delta + delta*math.Abs(params[col])
		back := SEnergy(params)
		diff := make([]float64, len(forward))
		for i := range forward {
			// convert mopac results to hartree
			diff[i] = (forward[i] - back[i]) / (2 * delta)
		}
		// reset params
		params[col] = reset
		jac.SetCol(col, diff)
	}
	return
}

func PlotData(nstruct, nparam int, jac *mat.Dense, w io.Writer) {
	var lines string
	for i := 0; i < nstruct; i++ {
		for j := 0; j < nparam; j++ {
			lines += fmt.Sprintf("%5d%5d%20.12f\n", i, j, jac.At(i, j))
		}
	}
	w.Write([]byte(lines))
}
