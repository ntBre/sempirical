package main

import (
	"bufio"
	"errors"
	"fmt"
	"gonum.org/v1/gonum/mat"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path"
	"regexp"
	"strconv"
	"strings"
	"time"
        "github.com/maorshutman/lm"
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

type MopacErr string

func (e MopacErr) Error() string {
	return string(e)
}

type Geometry struct {
	Atoms  []string
	Coords [][]float64
}

type Parameters struct {
	Names   []string
	Targets []string
	Values  []float64
}

func (p Parameters) Write(filename string) error {
	var n int
	var e error
	infile, err := os.OpenFile(filename,
		os.O_CREATE|os.O_WRONLY, 0755)
	defer infile.Close()
	if err != nil {
		log.Fatal(err)
	}
	for i, _ := range p.Names {
		line := fmt.Sprintf("%20s%10s%20.12f\n", p.Names[i],
			p.Targets[i], p.Values[i])
		n, e = infile.WriteString(line)
	}
	if e != nil {
		return errors.New(fmt.Sprintf("Error after writing %d bytes to %s", n, filename))
	}
	return nil
}

type Input struct {
	Nparam  int // len of param
	Method  string
	Nmolec  int
	Charge  int
	Nstruct int // len of efile/gfile read
	Efile   string
	Gfile   string
	Units   string
	Geom    Geometry
	Param   Parameters
}

func String2Array(s string) []string {
	re := regexp.MustCompile(`\s+`)
	fix := re.ReplaceAllString(s, " ")
	re2 := regexp.MustCompile(`^\s*`)
	fix2 := re2.ReplaceAllString(fix, "")
	split := strings.Split(fix2, " ")
	return split
}

func NotBlank(s string) bool {
	b, _ := regexp.Match(`^\s*$`, []byte(s))
	if !b {
		return true
	}
	return false
}

func NewParams(lines []string) Parameters {
	params := Parameters{}
	for _, val := range lines {
		split := String2Array(val)
		params.Names = append(params.Names, string(split[0]))
		params.Targets = append(params.Targets, string(split[1]))
		f, _ := strconv.ParseFloat(split[2], 64)
		params.Values = append(params.Values, f)
	}
	return params
}

func NewGeometry(lines []string) Geometry {
	geom := Geometry{}
	for _, val := range lines {
		split := String2Array(val)
		geom.Atoms = append(geom.Atoms, string(split[0]))
		tmp := make([]float64, 0)
		for _, coord := range split[1:4] {
			f, _ := strconv.ParseFloat(coord, 64)
			tmp = append(tmp, f)
		}
		geom.Coords = append(geom.Coords, tmp)
	}
	return geom
}

func NewInput(lines []string) Input {
	var start, end int
	var gflag, endflag bool
	params := Input{}
	for i := 0; i < len(lines); i++ {
		line := lines[i]
		switch {
		case strings.Contains(line, "="):
			line = strings.ReplaceAll(line, " ", "")
			split := strings.Split(line, "=")
			val, err := strconv.Atoi(split[1])
			if err != nil {
				switch split[0] {
				case method:
					params.Method = split[1]
				case efile:
					params.Efile = split[1]
				case gfile:
					params.Gfile = split[1]
				case units:
					params.Units = split[1]
				case geometry:
					start = i
					gflag = true
				case parameters:
					start = i
					gflag = false
				}
			} else {
				switch split[0] {
				case nparam:
					params.Nparam = val
				case nmolec:
					params.Nmolec = val
				case charge:
					params.Charge = val
				case nstruct:
					params.Nstruct = val
				}
			}
		case strings.Contains(line, "{"):
			start = i
		case strings.Contains(line, "}"):
			end = i
			endflag = true
		}
		if endflag {
			if gflag {
				params.Geom = NewGeometry(lines[start+1 : end])
			} else {
				params.Param = NewParams(lines[start+1 : end])
			}
			endflag = false
		}
	}
	return params
}

func ReadFile(filename string) (lines []string) {
	file, err := os.Open(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if len(line) > 0 {
			lines = append(lines, line)
		}
	}
	return lines
}

func ReadInp(filename string) Input {
	return NewInput(ReadFile(filename))
}

func ReadGeoms(filename string) [][]float64 {
	geoms := make([][]float64, 0)
	lines := ReadFile(filename)
	tmp := make([]float64, 0)
	for _, line := range lines[1:] {
		if strings.Contains(line, "#") {
			geoms = append(geoms, tmp)
			tmp = make([]float64, 0)
		} else {
			split := String2Array(line)
			for _, v := range split {
				if NotBlank(v) {
					f, _ := strconv.ParseFloat(v, 64)
					tmp = append(tmp, f)
				}
			}
		}
	}
	geoms = append(geoms, tmp)
	return geoms
}

func ReadEnergies(filename string) []float64 {
	lines := ReadFile(filename)
	energies := make([]float64, 0)
	for _, line := range lines {
		if NotBlank(line) {
			f, _ := strconv.ParseFloat(line, 64)
			energies = append(energies, f)
		}
	}
	return energies
}

func MakeMopacIn(inp Input, geom []float64) (mopfile []string) {
	mopfile = append(mopfile,
		fmt.Sprintf("threads=%d XYZ A0 scfcrt=%s"+
			" aux(precision=%d) external=%s 1SCF charge= %d %s",
			threads, scfcrt, auxprcsn, paramsfile, inp.Charge,
			inp.Method),
		"MOLECULE # 1", "")
	for i, v := range inp.Geom.Atoms {
		mopfile = append(mopfile, fmt.Sprintf("%s %.10f %.10f %.10f",
			v, geom[3*i], geom[3*i+1], geom[3*i+2]))
	}
	return mopfile
}

func WriteMopacIn(inp Input, geom []float64, n int) string {
	err := os.Mkdir("inp", 0755)
	// if err != nil {
	// 	log.Println("Warning: directory 'inp' exists, overwriting")
	// }
	lines := MakeMopacIn(inp, geom)
	writelines := strings.Join(lines, "\n")
	filename := BaseMopFilename + fmt.Sprintf("%05d.mop", n)
	err = ioutil.WriteFile(filename, []byte(writelines), 0755)
	if err != nil {
		log.Fatal(err)
	}
	return filename
}

func SlurmSubmit(filename string) int {
	cmd := exec.Command("sbatch", "--job-name=qff", "--ntasks=1",
		"--cpus-per-task=1", "--mem=1gb ",
		"--output=test.out", "--error=test.err",
		"/home/qc/mopac2016/mopac.sh", filename)
	out, err := cmd.Output()
	if err != nil {
		log.Fatal(err)
	}
	val := strings.Split(string(out), " ")
	i, _ := strconv.Atoi(strings.TrimSpace(val[len(val)-1]))
	return i
}

func ReadMopacOut(job Job) (interface{}, error) {
	var eline []string
	re := regexp.MustCompile(`mop`)
	auxfile := re.ReplaceAllString(job.Name, "aux")
	re = regexp.MustCompile(`aux`)
	outfile := re.ReplaceAllString(auxfile, "out")
	infile := re.ReplaceAllString(auxfile, "mop")
	if _, err := os.Stat(auxfile); !os.IsNotExist(err) {
		lines := ReadFile(auxfile)
		for _, line := range lines {
			if strings.Contains(line, "TOTAL_ENERGY") {
				re := regexp.MustCompile(`D`)
				fix := re.ReplaceAllString(line, "E")
				eline = strings.Split(fix, "=")
				f, _ := strconv.ParseFloat(eline[len(eline)-1], 64)
				return f, nil // f in ev
			}
		}
	} else if _, err := os.Stat(outfile); !os.IsNotExist(err) {
		lines := ReadFile(outfile)
		for _, line := range lines {
			if strings.Contains(line, "MOPAC DONE") &&
				job.Retries < maxretries {
				return SlurmSubmit(infile),
					errors.New("Resubmitting job")
			} else if job.Retries >= maxretries {
				return nil, errTooManyRetries
			}
		}
	}
	return nil, errFilesNotFound
}

func WriteOutfile(filename string) error {
	// should take the basename from main method and write
	// the iteration number, rmsd, maybe something else if needed
	return nil
}

func Basename(filename string) string {
	file := path.Base(filename)
	re := regexp.MustCompile(path.Ext(file))
	basename := re.ReplaceAllString(file, "")
	return basename
}

func LossFunction(dst, x []float64) {
	inp := ReadInp("new.inp")
	geom := ReadGeoms("file07")
	energies := ReadEnergies("energy.dat")
	energyVector := mat.NewVecDense(len(energies), energies)
	inp.Param.Values = x
	inp.Param.Write("inp/" + paramsfile)
	emin := mat.Min(energyVector)
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

type Job struct {
	Name    string
	Number  int
	Status  string
	Retries int
}

func main() {
	// this inp is different from the inp inside
	// hence, the inside one is not updating and the params are not
	// either, rendering it useless
	// the initparams below should be being updated but are not written
	// and the params in LossFunction are reread from the input
	// on each iteration
	inp := ReadInp("new.inp")
	biggsNumJac := lm.NumJac{Func: LossFunction}
	biggsProb := lm.LMProblem{
	Dim:        inp.Nparam,
 	Size:       inp.Nstruct,
 	Func:       LossFunction,
 	Jac:        biggsNumJac.Jac,
 	InitParams: inp.Param.Values,
	Tau:        1e-3, // 6 -> 3
 	Eps1:       1e-8,
 	Eps2:       1e-8,}
biggsResults, biggsErr := lm.LM(biggsProb, &lm.Settings{Iterations: 100, ObjectiveTol: 1e-16})
fmt.Println(biggsResults)
fmt.Println(biggsErr)
}
