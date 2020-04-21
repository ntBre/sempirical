package main

import (
	_ "gonum.org/v1/gonum/mat"
	"log"
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

func Fcn(inp *Input, energies []float64, geom [][]float64) (fvec []float64) {
	inp.Param.Write(paramsfile)
	var jobs []Job
	fvec = make([]float64, len(geom))
	for i, _ := range energies {
		filename := WriteMopacIn(*inp, geom[i], i+1)
		jobnum := SlurmSubmit(filename)
		jobs = append(jobs, Job{filename, jobnum, "queued", 0})
	}
	njobs := len(jobs)
	for njobs > 0 {
		for i, job := range jobs {
			if job.Status != "done" {
				f, err := ReadMopacOut(job)
				switch err {
				case nil:
					fvec[i] = energies[i] - f.(float64)*cm1ev
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
	return fvec
}

func Enorm(fvec []float64) float64 {
	sum := 0.0
	for _, v := range fvec {
		sum += math.Pow(v, 2)
	}
	return math.Sqrt(sum)
}

func Lmdif(inp Input, energies []float64, geom [][]float64) {
	// fvec := Fcn(inp, energies, geom)
	// fnorm := Enorm(fvec)
}

func Fdjac2(inp *Input, energies []float64, geom [][]float64, fvec []float64) [][]float64 {
	fjac := make([][]float64, len(fvec))
	for k, _ := range fjac {
		fjac[k] = make([]float64, len(inp.Param.Values))
	}
	for j, _ := range inp.Param.Values {
		temp := inp.Param.Values[j]
		h := eps * math.Abs(temp)
		inp.Param.Values[j] = inp.Param.Values[j] + h
		wa := Fcn(inp, energies, geom)
		inp.Param.Values[j] = temp
		for i, _ := range fvec {
			fjac[i][j] = (wa[i] - fvec[i]) / h
		}
	}
	return fjac

}

func Qrfac(fjac [][]float64) {
	//return work arrays? maybe up to 3 arrays
}

// func main() {
// 	inp := ReadInp("new.inp")
// 	geom := ReadGeoms("file07")
// 	energies := ReadEnergies("energy.dat")
// 	energyVector := mat.NewVecDense(len(energies), energies)
// 	emin := mat.Min(energyVector)
// 	for i, _ := range energies {
// ab initio energies now in wavenumbers
// 		energies[i] = (energies[i] - emin) * cm1
// 	}

// fcn -> Lossfunction
// fortran fcn is my lossfunction
// }
