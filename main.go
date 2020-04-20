package main

import (
	"fmt"
	"github.com/maorshutman/lm"
	"gonum.org/v1/gonum/mat"
	"log"
	"time"
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
)

func main() {
	inp := ReadInp("new.inp")
	geom := ReadGeoms("file07")
	energies := ReadEnergies("energy.dat")
	energyVector := mat.NewVecDense(len(energies), energies)
	emin := mat.Min(energyVector)

	LossFunction := func(dst, x []float64) {
		inp.Param.Values = x
		inp.Param.Write("inp/" + paramsfile)
		var jobs []Job
		dst = make([]float64, len(geom))
		for i, _ := range energies {
			// for i, _ := range []int{1, 2, 3, 4, 5} {
			energies[i] = (energies[i] - emin) * cm1 // ab initio energies now in wavenumbers
			filename := WriteMopacIn(inp, geom[i], i+1)
			jobnum := SlurmSubmit(filename)
			jobs = append(jobs, Job{filename, jobnum, "queued", 0})
			time.Sleep(100 * time.Millisecond)
		}
		njobs := len(jobs)
		for njobs > 0 {
			for i, job := range jobs {
				if job.Status != "done" {
					f, err := ReadMopacOut(job)
					switch err {
					case nil:
						dst[i] = energies[i] - f.(float64)*cm1ev
						job.Status = "done"
						njobs -= 1
					case errFilesNotFound:
						continue
					case errTooManyRetries:
						log.Fatal(err)
					}
				}
			}
			time.Sleep(1000 * time.Millisecond)
		}
		inp.Param.Write("inp/" + paramsfile)
	}

	biggsNumJac := lm.NumJac{Func: LossFunction}
	biggsProb := lm.LMProblem{
		Dim:        inp.Nparam,
		Size:       inp.Nstruct,
		Func:       LossFunction,
		Jac:        biggsNumJac.Jac,
		InitParams: inp.Param.Values,
		Tau:        1e-3, // 6 -> 3
		Eps1:       1e-8,
		Eps2:       1e-8}
	biggsResults, biggsErr := lm.LM(biggsProb, &lm.Settings{Iterations: 100, ObjectiveTol: 1e-16})
	fmt.Println(biggsResults)
	fmt.Println(biggsErr)
}
