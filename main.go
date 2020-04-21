package main

import (
	_ "gonum.org/v1/gonum/mat"
	"log"
	"github.com/maorshutman/lm"
	"fmt"
	"math"
)

const (
	nparam            = "NPARAM"
	method            = "METHOD"
	nmolec            = "NMOLEC"
	charge            = "CHARGE"
	nstruct           = "NSTRUCT"
	efile             = "EFILE"
	gfile             = "GFILE"
	units             = "UNITS"
	geometry          = "GEOMETRY"
	parameters        = "PARAMS"
	threads           = 1            // number of threads in mopac input
	scfcrt            = "1.D-21"     // scf convergence criteria for mopac inp
	auxprcsn          = 9            // aux precision in mopac inp
	paramsfile        = "params.dat" // file for current mopac params
	BaseMopFilename   = "inp/Structure"
	cm1               = 219474.5459784 // wavenumbers per hartree
	cm1ev             = 8065.541154    //wavenumbers per ev
	maxretries        = 3              // max number of times to try running mopac
	errTooManyRetries = MopacErr("Too many retries")
	errFilesNotFound  = MopacErr(".aux and .out files not found")
	eps               = 1.0e-6 // tolerance
)

func main() {
	inp := ReadInp("new.inp")
	geom := ReadGeoms("file07")
	energies := ReadEnergies("energy.dat")
	// why am i using relative energies when the result is not relative
	// energyVector := mat.NewVecDense(len(energies), energies)
	// emin := mat.Min(energyVector)
	// use relative energies
	// for i, _ := range energies {
	// 	energies[i] = (energies[i] - emin) * cm1
	// }
	Fcn := func(dst, x []float64)  {
		inp.Param.Values = x
		inp.Param.Write(paramsfile)
		var jobs []Job
		var sqDev float64
		for i, _ := range energies {
			filename := WriteMopacIn(inp, geom[i], i+1)
			jobnum := SlurmSubmit(filename)
			jobs = append(jobs, Job{filename, jobnum, "queued", 0})
		}
		njobs := len(jobs)
		for njobs > 0 {
			for i, job := range jobs {
				if job.Status != "done" {
					f, err := ReadMopacOut(&job)
					switch err {
					case nil:
						dst[i] = energies[i] - f.(float64)*cm1ev
						sqDev += math.Pow(dst[i], 2)
						job.Status = "done"
						njobs -= 1
					case errFilesNotFound:
						continue
					case errTooManyRetries:
						log.Fatal(err)
					}
				}
			}
		}
		fmt.Println(math.Sqrt(sqDev / float64(inp.Nstruct))) // print rmsd?
	}
	fcnJac := lm.NumJac{Func: Fcn}
	Prob := lm.LMProblem{
		Dim: inp.Nparam,
		Size: inp.Nstruct,
		Func: Fcn,
		Jac: fcnJac.Jac,
		InitParams: inp.Param.Values,
		Tau: 1e-6,
		Eps1: 1e-8,
		Eps2: 1e-8}
	biggsResults, biggsErr := lm.LM(Prob, &lm.Settings{Iterations: 1, ObjectiveTol: 1e-16})
	fmt.Println(biggsResults)
	fmt.Println(biggsErr)
}
