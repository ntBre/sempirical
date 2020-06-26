/*
sempirical optimizes semi-empirical parameters in MOPAC to minimize
the difference between the semi-empirical energies and input ab initio
energies
Requires:
  - new.inp
  - ab initio energy file, named in new.inp
  - geometry file, named in new.inp
*/
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
	mopacTerminated = "END OF MOPAC FILE"
	EVtoHt          = 1 / 27.211385 // from http://www.ilpi.com/msds/ref/energyunits.html
	HtToCm          = 219474.5459784
	EvToCm          = 8065.541154
	delta           = 1e-8 // roughly the same as from fortran
	epsilon         = 1e-5 // required convergence in Ht for ~1 cm-1
	lambda0         = 1e-2 // from Marquardt63
	nuDown          = 5
	nuUp            = 1.5
	broyden         = true // toggle broyden update
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
	nu           = 10.0 // from Marquardt63
	maxIter      = 100
	logfile, _   = os.Create("semp.log")
	parlog, _    = os.Create("params.log")
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
	parlog.Write([]byte(lines + "\n"))
}

// Write MOPAC input files with the given geometries
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
		energy, err := mop.ReadAux(file + ".aux")
		for err != nil {
			// this indicates a bad step, so should propagate this error
			// and revise the step higher up
			if e := mop.CheckOut(file + ".out"); e != nil {
				// assume error caused by bad step so give a very large energy
				// and log the error
				fmt.Fprintln(os.Stderr, e)
				energy, err = 0.0, nil
				break
			}
			if err == ErrFinishedButNoEnergy {
				// also assume this is caused by bad step
				fmt.Fprintln(os.Stderr, err)
				energy, err = 0.0, nil
				break
			}
			// else assume problem from reading before energy present
			time.Sleep(time.Second)
			energy, err = mop.ReadAux(file + ".aux")
		}
		// convert to hartrees here
		energies = append(energies, energy*EVtoHt)
		// convert to cm-1
		// energies = append(energies, energy*EvToCm)
		// ensure we are looking at the new files each time
		// should move this outside the closure if I make this a goroutine
		// to avoid too many sys calls
		os.Remove(file + ".aux")
		os.Remove(file + ".out")
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
			// in hartrees
			abInit = append(abInit, f)
			// in cm-1
			// abInit = append(abInit, f*HtToCm)
		}
	}
	return
}

func init() {
	ParseInfile("new.inp")
	if Input[MaxIter] != "" {
		maxIter, _ = strconv.Atoi(Input[MaxIter])
	}
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
	// j*step
	prod := mat.NewDense(nstruct, 1, nil)
	prod.Mul(jac, step)
	diff := mat.NewDense(nstruct, 1, nil)
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
	fmt.Fprintln(logfile, "num: ", DenseSlice(num))
	fmt.Fprintln(logfile, "jac: ", DenseSlice(jac))
	// think I just want negative scaled prod2 as jac, that looks reasonable
	// prod2.Scale(-1/DenseSlice(num)[0], prod2)
	// fmt.Fprintln(logfile, "-prod2: ", DenseSlice(prod2))
	// return prod2
	prod2.Scale(1/DenseSlice(num)[0], prod2)
	jac.Add(jac, prod2)
	return jac
}

// Wrapper for calculating new parameters and the RMSD.
// Returns the new params, the step that generated them,
// the rmsd, the difference matrix that it was calculated
// from, and the matrix of new semi-empirical energies
func TryStep(nstruct, nparam int, lam float64, params []float64,
	b, jacTjac, aiCol *mat.Dense) ([]float64, *mat.Dense,
	float64, *mat.Dense, *mat.Dense) {
	newparams, step := NewParams(nstruct, nparam, lam, params, b, jacTjac)
	seColNew := mat.NewDense(nstruct, 1, SEnergy(newparams))
	rmsd, diff := RMSD(nstruct, newparams, seColNew, aiCol)
	return newparams, step, rmsd, diff, seColNew
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
		numJacs   int
		newparams []float64
		rmsd      float64
		newRMSD   float64
		jac       *mat.Dense
		diff      *mat.Dense
		step      *mat.Dense
		seCol     *mat.Dense
		seColNew  *mat.Dense
	)

	// ab initio energy column vector
	aiCol := mat.NewDense(nstruct, 1, abInit)

	// initial calculations
	jac = Jacobian(params)
	seCol = mat.NewDense(nstruct, 1, SEnergy(params))
	rmsd, diff = RMSD(nstruct, params, seCol, aiCol)
	i := 1
	for rmsd > epsilon && i < maxIter {
		// calculate and print the RMSD
		if i > 1 {
			rmsd = newRMSD
		}
		fmt.Println("initial rmsd:", rmsd)

		// jac transpose
		jacT := mat.DenseCopyOf(jac.T())

		// jac transpose jac
		jacTjac := mat.NewDense(nparam, nparam, nil)
		jacTjac.Mul(jacT, jac)

		// b = jac^T*[aiCol - seCol]
		b := mat.NewDense(nparam, 1, nil)
		b.Mul(jacT, diff)

		// Try with lambda/nu
		newparams, step, newRMSD, diff, seColNew = TryStep(nstruct, nparam, lambda/(nuDown*nu),
			params, b, jacTjac, aiCol)

		// if both are greater, try lambda*nu^k until a k that improves it
		// Marquardt63 case i pg 8
		if newRMSD < rmsd {
			lambda /= nuDown * nu
		} else {
			// Try with lambda
			// case ii, successful step, no change
			newparams, step, newRMSD, diff, seColNew = TryStep(nstruct, nparam, lambda,
				params, b, jacTjac, aiCol)
			if newRMSD >= rmsd {
				try := 1
				// 6/23
				// need to be better than the previous, not just equal to it
				for newRMSD >= rmsd {
					if newRMSD-rmsd == 0 {
						// 6/23
						// give up and run a numerical jacobian
						// might want to hold the starting lambda value and reset if this happens
						break
					}
					// if we make it here, lambda has to be updated this way
					// so just go ahead and modify lambda
					lambda *= nuUp * nu
					newparams, step, newRMSD, diff, seColNew = TryStep(nstruct,
						nparam, lambda, params, b, jacTjac, aiCol)
					b.Mul(jacT, diff)
					fmt.Printf("\ttry = %d, rmsd: %g, lambda: %g\n", try, newRMSD,
						lambda)
					try++
				}
			}
		}
		fmt.Printf("-> final rmsd: %g, lambda: %g\n", newRMSD, lambda)
		fmt.Printf("     Delta(%d): %g\n", i, newRMSD-rmsd)
		params = newparams
		// every 2n iterations, calculate the Jacobian by finite differences
		// or if broyden is disabled or if delta too small
		// 6/23
		// dont run numjac every 2n times if it's going fine, other criterion better
		// i%(2*nparam) == 0 ||
		if !broyden || newRMSD-rmsd == 0 {
			numJacs++
			fmt.Println("running numerical jacobian")
			jac = Jacobian(params)
			// if change is too small, try taking a bigger step
			// could also try altering delta? for bigger jacobian step
			// maybe after too many num jacs in a row
			// two might be too many
			if newRMSD-rmsd == 0 {
				// shrink lambda faster if multiple num jacs in a row
				// 6/23 - add nuDown to product
				lambda /= float64(numJacs) * nu * nuDown
				// 6/23
				// try resetting lambda, not working otherwise
				// lambda = lambda0
				// full reset is bad because then it loops right away.
				// at least if you vary it by a different amount from nuUp, it
				// hits different values
			}
			// I think lambda should be adjusted in here
			// or count times through here with no change so we can break
			// see last run for infinite num jac loop
		} else {
			// otherwise update by Broyden's method
			jac = Broyden(nstruct, nparam, jac, step, seColNew, seCol)
			// and reset numjacs
			numJacs = 0
		}
		seCol = seColNew
		i++
	}
}

func main() {
	// J^t J h = J^T[y-f(B)]
	// where J is the Jacobian, h is the step,
	// y is the exact solution, and f(B) is
	// the evaluation of f at the parameter vector B
	LevMar(len(initParams), len(jobs), initParams)
	fmt.Println("exiting")
	logfile.Close()
	parlog.Close()
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
