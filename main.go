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
	maxIter         = 5
	lambda0         = 1e-8 // initial levenberg-marquardt damping parameter
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
		lines += fmt.Sprintf("%s%25.12f\n", headers[i], params[i])
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
			time.Sleep(time.Second)
			// convert to hartrees here
			energy, err = mop.ReadOut(file + ".aux")
		}
		energies = append(energies, energy*EVtoHt)
		// assure we are looking at the new files each time
		os.Remove(file + ".aux")
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

// Update the Jacobian using Broyden's method
// jac = jac + ((seColnew - seCol - jac*step)*step^T)/(step^T*step)
func Broyden(nstruct, nparam int, jac, step, seColnew, seCol *mat.Dense) *mat.Dense {
	diff := mat.NewDense(nstruct, 1, nil)
	tempJac := mat.DenseCopyOf(jac)
	// j*step
	prod := mat.NewDense(nstruct, 1, nil)
	prod.Mul(tempJac, step)
	// fBnew - fb
	diff.Sub(seColnew, seCol)
	// fBnew - fB - J*step
	diff.Sub(diff, prod)
	// (fBnew - fB - j*step)*step^T
	prod2 := mat.NewDense(nstruct, nparam, nil)
	prod2.Mul(diff, step.T())
	num := mat.NewDense(1, 1, nil)
	num.Mul(step.T(), step)
	// num should be 1 by 1, so take first element in slice version
	prod2.Scale(1/DenseSlice(num)[0], prod2)
	jac.Add(jac, prod2)
	return jac
}

// Calculate new parameter values based on the Levenberg-Marquardt algorithm.
// b is the righthand side of Ax=b, where x are the parameters.
// b = J^T * [y - f(B)], where y is the exact solution and f(B) is
// the function of the parameters
func NewParams(nstruct, nparam int, lscale float64, params []float64,
	b, jacTjac *mat.Dense) (newparams []float64, step *mat.Dense) {

	ld := Eye(nparam)
	ld.Scale(lscale, ld)
	A := mat.NewDense(nparam, nparam, nil)
	A.Add(jacTjac, ld)
	step = mat.NewDense(nparam, 1, nil)
	err := step.Solve(A, b)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
	}
	return AddSlice(params, DenseSlice(step)), step
}

// Calculate the RMSD associated with a set of parameters
func RMSD(nstruct int, params []float64,
	seCol, aiCol *mat.Dense) (rmsd float64, diff *mat.Dense) {

	diff = mat.NewDense(nstruct, 1, nil)
	diff.Sub(aiCol, seCol)
	return ColRMS(diff), diff
}

// Minimize the difference between the semi-empirical
// and ab initio energies by varying the empirical
// parameters with the Levenberg-Marquardt algorithm
func LevMar(nparam, nstruct int, params []float64) {

	var (
		jac      *mat.Dense
		seCol    *mat.Dense
		seColnew *mat.Dense
		step     *mat.Dense
	)

	// ab initio energy column vector
	aiCol := mat.NewDense(nstruct, 1, abInit)

	for i := 0; i < maxIter; i++ {

		// calculate and print the RMSD
		seCol = mat.NewDense(nstruct, 1, SEnergy(params))
		rmsd, diff := RMSD(nstruct, params, seCol, aiCol)
		fmt.Println("initial rmsd:", rmsd)

		// every 2n iterations, calculate the Jacobian by finite differences
		if i%(2*nparam) == 0 {
			jac = Jacobian(params)
		} else {
			// otherwise update by Broyden's method
			jac = Broyden(nstruct, nparam, jac, step, seColnew, seCol)
		}

		// jac transpose
		jacT := mat.DenseCopyOf(jac.T())

		// jac transpose jac
		jacTjac := mat.NewDense(nparam, nparam, nil)
		jacTjac.Mul(jacT, jac)

		// b = jac^T*[aiCol - seCol]
		b := mat.NewDense(nparam, 1, nil)
		b.Mul(jacT, diff)

		// update the parameters
		newparams, newStep := NewParams(nstruct, nparam, lambda, params, b, jacTjac)
		newparamsnu, newStepnu := NewParams(nstruct, nparam,
			lambda/nu, params, b, jacTjac)

		// check rmsd with lambda and lambda/nu dampening
		newRMSD, _ := RMSD(nstruct, newparams, aiCol,
			mat.NewDense(nstruct, 1, SEnergy(newparams)))
		newRMSDnu, _ := RMSD(nstruct, newparamsnu, aiCol,
			mat.NewDense(nstruct, 1, SEnergy(newparamsnu)))

		// if both are greater, try lambda*nu^k until a k that improves it
		if newRMSD > rmsd && newRMSDnu > rmsd {
			tryRMSD := newRMSD
			tryParams := newparams
			tryStep := newStep
			k := 1
			for tryRMSD > rmsd {
				prod := 1.0
				for i := 0; i < k; i++ {
					prod *= nu
				}
				tryParams, tryStep = NewParams(nstruct, nparam, lambda*prod,
					params, b, jacTjac)
				tryRMSD, _ = RMSD(nstruct, tryParams, aiCol,
					mat.NewDense(nstruct, 1, SEnergy(tryParams)))
				k++
			}
			newRMSD = tryRMSD
			newparams = tryParams
			newStep = tryStep
		} else if newRMSDnu < newRMSD {
			newRMSD = newRMSDnu
			newparams = newparamsnu
			newStep = newStepnu
			lambda = lambda / nu
		}
		// else lambda stays the same and use newRMSD as it was
		fmt.Println("final rmsd:", newRMSD)
		params = newparams
		step = newStep
		seColnew = mat.NewDense(nstruct, 1, SEnergy(params))
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
