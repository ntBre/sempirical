package main

import (
	"fmt"
	"io/ioutil"
	"os"
	"strings"

	"math"

	"strconv"

	"time"

	"bytes"

	"gonum.org/v1/gonum/mat"
)

const (
	energyLine      = "TOTAL_ENERGY:EV="
	mopacTerminated = "not sure what this is yet"
	EVtoHt          = 1 / 27.211385 // from http://www.ilpi.com/msds/ref/energyunits.html
	delta           = 1e-8          // roughly the same as from fortran
	maxIter         = 1
	lambda0         = 1.0 // initial levenberg-marquardt damping parameter
)

var (
	Input        [Nkeys]string
	mop          = Mopac{}
	queue        = Slurm{}
	jobs         []string
	abInit       []float64
	initParams   []float64
	paramHeaders []string
	lambda       = lambda0
	nu           = 1.1
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
	initParams, paramHeaders = FloatParams(Input[Params])
}

// generate the nxn identity matrix
func Eye(n int) *mat.Dense {
	m := mat.NewDense(n, n, nil)
	for i := 0; i < n; i++ {
		m.Set(i, i, 1.0)
	}
	return m
}

// return a diagonal matrix corresponding to the diagonal
// of the input
func Diag(m *mat.Dense) *mat.Dense {
	r, c := m.Caps()
	if r != c {
		err := "diag: Input matrix is not square, aborting"
		panic(err)
	}
	n := mat.NewDense(r, c, nil)
	for i := 0; i < r; i++ {
		n.Set(i, i, m.At(i, i))
	}
	return n
}

// Pretty print a mat.Dense matrix
func PrettyPrint(m *mat.Dense) string {
	r, c := m.Caps()
	var buf bytes.Buffer
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			fmt.Fprintf(&buf, "%g", m.At(i, j))
			if j != c-1 {
				fmt.Fprint(&buf, " ")
			}
		}
		fmt.Fprint(&buf, "\n")
	}
	return buf.String()
}

// calculate the root-mean-square for a mat.Dense column vector
func ColRMS(m *mat.Dense) float64 {
	r, _ := m.Caps()
	var total float64
	for i := 0; i < r; i++ {
		val := m.At(i, 0)
		total += val * val
	}
	mean := total / float64(r)
	return math.Sqrt(mean)
}

// returning the backing slice for a mat.Dense
func DenseSlice(m *mat.Dense) []float64 {
	r, c := m.Caps()
	params := make([]float64, 0)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			params = append(params, m.At(i, j))
		}
	}
	return params
}

// Add two float64 slices
func AddSlice(a, b []float64) (sum []float64) {
	la, lb := len(a), len(b)
	if la != lb {
		err := fmt.Sprintf("AddSlice: slices have different lengths,"+
			" aborting. len(a): %d, len(b): %d\n", la, lb)
		panic(err)
	}
	for i := range a {
		sum = append(sum, a[i]+b[i])
	}
	return
}

// Minimize the difference between the semi-empirical
// and ab initio energies by varying the empirical
// parameters with the Levenberg-Marquardt algorithm
func LevMar(nparam, nstruct int, params []float64) {
	var newparams []float64
	for i := 0; i < maxIter; i++ {
		// semi-empirical energy column vector
		fB := mat.NewDense(nstruct, 1, SEnergy(params))
		// ab initio column vector
		y := mat.NewDense(nstruct, 1, abInit)
		b := mat.NewDense(nstruct, 1, nil)
		b.Sub(y, fB)
		rmsd := ColRMS(b)
		fmt.Println("initial rmsd:", rmsd)
		jac := Jacobian(params)
		jacT := mat.DenseCopyOf(jac.T())
		A := mat.NewDense(nparam, nparam, nil)
		// A = jac^T jac + lambda diag(jac^t jac)
		JTJ := mat.NewDense(nparam, nparam, nil)
		// lambda * diag(jac^t jac)
		ld := Diag(JTJ)
		ld.Scale(lambda, ld)
		JTJ.Mul(jacT, jac)
		A.Add(JTJ, ld)
		// b = jac^T*[y - fB]
		RHS := mat.NewDense(nparam, 1, nil)
		RHS.Mul(jacT, b)
		fmt.Println(PrettyPrint(RHS))
		step := mat.NewDense(nparam, 1, nil)
		fmt.Println(PrettyPrint(step))
		// Solve Ax = b for x, which is our step in the parameters
		step.Solve(A, RHS)
		fmt.Println("params, step:", params, step)
		newparams = AddSlice(params, DenseSlice(step))
		fmt.Println("new params:", newparams)
		fBnew := mat.NewDense(nstruct, 1, SEnergy(newparams))
		bnew := mat.NewDense(nstruct, 1, nil)
		bnew.Sub(y, fBnew)
		fmt.Println("final rmsd:", ColRMS(bnew))
		// update lambda and nu:
		// try step with lambda, and lambda/nu
		// if neither is better, try lambda*nu^k until improves
		// with k increasing on each attempt
		// lambda becomes whatever is successful
		// broyden update (eq 19 from lm.pdf):
		// J = J + ((fBnew - fB - J*step)*step^T)/(step^T*step)
		// lm.pdf says only calculate the jacobian directly on:
		// - the first iteration
		// -- if i == 0
		// - every 2n iterations, where n is # parameters
		// -- if i % 2nparam == 0
		// - iterations where rmsd increases
		// -- if rmsdNew > rmsd
		// if fBnew "better" than before keep the step
		// fmt.Printf("%v, %v, %v, %v, %v, %v\n", fB, y, jac, A, b, step)
	}
}

func main() {
	// J^t J h = J^T[y-f(B)]
	// where J is the Jacobian, h is the step,
	// y is the exact solution, and f(B) is
	// the evaluation of f at the parameter vector B
	LevMar(len(initParams), len(jobs), initParams)
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
