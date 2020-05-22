package main

import (
	"math"
	"reflect"
	"testing"

	"gonum.org/v1/gonum/mat"
)

func TestAbEnergies(t *testing.T) {
	got := AbEnergies("testfiles/energy.min")
	want := []float64{-153.506766246758,
		-153.506744023064,
		-153.506744023097,
		-153.506677292842,
		-153.506677292897,
		-153.506736347881}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func TestMakeParams(t *testing.T) {
	// rely on FloatParams working
	params, headers := FloatParams(Input[Params])
	got := MakeParams(params[:2], headers[:2])
	want := "FN11 H      0.024184000000\n" +
		"FN11 C      0.046302000000\n"
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func TestSEnergy(t *testing.T) {
	holdJobs := jobs
	jobs = []string{"inp/Structure00000", "inp/Structure00001"}
	params, _ := FloatParams(Input[Params])
	got := SEnergy(params)
	want := []float64{-19.018233032391734, -19.018233032391734}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
	jobs = holdJobs
}

func TestJacobian(t *testing.T) {
	holdJobs := jobs
	jobs = []string{"inp/Structure00000", "inp/Structure00001"}
	params := []float64{1.0, 2.0}
	got := Jacobian(params)
	// forward and reverse differences are the same
	// for the dummy sbatch
	back := []float64{0, 0, 0, 0}
	want := mat.NewDense(len(jobs), len(params), back)
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
	jobs = holdJobs
}

func TestLevMar(t *testing.T) {
	holdJobs := jobs
	holdParams := initParams
	holdAb := abInit
	abInit = AbEnergies("testfiles/energy.2")
	jobs = []string{"inp/Structure00000", "inp/Structure00001", "inp/Structure00002"}
	params := []float64{1.0, 2.0}
	LevMar(len(params), len(jobs), params)
	// if got != want {
	// 	t.Errorf("got %v, wanted %v\n", got, want)
	// }
	jobs = holdJobs
	initParams = holdParams
	abInit = holdAb
}

func TestPrettyPrint(t *testing.T) {
	back := []float64{1, 2, 3, 4, 5, 6, 7, 8, 9}
	m := mat.NewDense(3, 3, back)
	got := PrettyPrint(m)
	want := "1 2 3\n4 5 6\n7 8 9\n"
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func TestDiag(t *testing.T) {
	back := []float64{1, 2, 3, 4, 5, 6, 7, 8, 9}
	m := mat.NewDense(3, 3, back)
	got := Diag(m)
	want := mat.NewDense(3, 3,
		[]float64{1, 0, 0, 0, 5, 0, 0, 0, 9})
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func TestEye(t *testing.T) {
	got := Eye(2)
	want := mat.NewDense(2, 2,
		[]float64{1, 0, 0, 1})
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func TestColRMS(t *testing.T) {
	m := mat.NewDense(2, 1,
		[]float64{3, 4})
	got := ColRMS(m)
	want := math.Sqrt(12.5)
	if got != want {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}
