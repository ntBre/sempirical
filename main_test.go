package main

import (
	"reflect"
	"regexp"
	"strconv"
	"testing"
	"fmt"
	"gonum.org/v1/gonum/mat"
)

var testinp = Input{
	Nparam:  27,
	Method:  "PM6",
	Nmolec:  1,
	Charge:  1,
	Nstruct: 10353,
	Efile:   "energy.dat",
	Gfile:   "file07",
	Units:   "ANG",
	Geom: Geometry{[]string{"H", "C", "C", "C", "C", "H", "H"},
		[][]float64{{0.0000000000, 0.0000000000, -5.0705374640},
			{0.0000000000, 0.0000000000, -3.0439358050},
			{0.0000000000, 0.0000000000, -0.7125467130},
			{0.0000000000, -1.2030187010, 1.8800877760},
			{0.0000000000, 1.2030187010, 1.8800877760},
			{0.0000000000, -3.1388743520, 2.5132649210},
			{0.0000000000, 3.1388743520, 2.5132649210}}},
	Param: Parameters{[]string{"FN11", "FN11", "XFAC_C", "HSP",
		"XFAC_C", "ALPB_C", "GPP", "USS", "GSP", "ZS", "GSS",
		"FN31", "GSS", "BETAS", "ZP", "FN31", "ZS", "FN21",
		"XFAC_H", "ALPB_C", "FN21", "ALPB_H", "UPP", "USS",
		"BETAP", "BETAS", "GP2"},
		[]string{"H", "C", "H", "C", "C", "H", "C", "H", "C", "H",
			"C", "C", "H", "C", "C", "H", "C", "C", "H", "C",
			"H", "H", "C", "C", "C", "H", "C"},
		[]float64{0.02418400, 0.04630200, 0.21650600, 0.71732200,
			0.81351000, 1.02780600, 10.77832600, -11.24695800,
			11.52813400, 1.26864100, 13.33551900, 1.33395900,
			14.44868600, -15.38523600, 1.70284100, 1.78601100,
			2.04755800, 2.10020600, 2.24358700, 2.61371300,
			3.05595300, 3.54094200, -39.93792000, -51.08965300,
			-7.47192900, -8.35298400, 9.48621200}}}

var testgeom = [][]float64{{0.0000000000, 0.0000000000, -5.0705374640,
	0.0000000000, 0.0000000000, -3.0439358050,
	0.0000000000, 0.0000000000, -0.7125467130,
	0.0000000000, -1.2030187010, 1.8800877760,
	0.0000000000, 1.2030187010, 1.8800877760,
	0.0000000000, -3.1388743520, 2.5132649210,
	0.0000000000, 3.1388743520, 2.5132649210},
	{-0.0000000000, 0.0000000000, -5.0706091104,
		-0.0000000000, 0.0000000000, -3.0440074514,
		0.0000000000, 0.0000000000, -0.7126183594,
		0.0140780719, -1.2029363254, 1.8800161296,
		-0.0140780719, 1.2029363254, 1.8800161296,
		-0.0053955418, -3.1386119859, 2.5134440370,
		0.0053955418, 3.1386119859, 2.5134440370}}

func TestReadInp(t *testing.T) {

	t.Run("no spaces in param names", func(t *testing.T) {
		got := ReadInp("new.inp")
		if !reflect.DeepEqual(got, testinp) {
			t.Errorf("got %v, wanted %v", got, testinp)
		}
	})

	t.Run("spaces in param names", func(t *testing.T) {
		got := ReadInp("spaces.inp")
		if !reflect.DeepEqual(got, testinp) {
			t.Errorf("got %v, wanted %v", got, testinp)
		}
	})

	t.Run("blank lines between params", func(t *testing.T) {
		got := ReadInp("returns.inp")
		if !reflect.DeepEqual(got, testinp) {
			t.Errorf("got %v, wanted %v", got, testinp)
		}
	})
}

func TestString2Array(t *testing.T) {
	got := String2Array("this is a string")
	want := []string{"this", "is", "a", "string"}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, want %v", got, want)
	}
}

func TestReadGeoms(t *testing.T) {
	got := ReadGeoms("test.geom")
	if !reflect.DeepEqual(got, testgeom) {
		t.Errorf("got %v, want %v", got, testgeom)
	}
}

func TestReadEnergies(t *testing.T) {
	got := ReadEnergies("test.nrg")
	want := []float64{-153.506766246758, -153.506744023064,
		-153.506744023097, -153.506677292842,
		-153.506677292897, -153.506736347881,
		-153.506736347898, -153.506714117998,
		-153.506714118182, -153.506714118077}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, want %v", got, want)
	}
}

func TestMakeMopacIn(t *testing.T) {
	got := MakeMopacIn(testinp, testgeom[0])
	want := []string{"threads=1 XYZ A0 scfcrt=1.D-21 aux(precision=9) external=params.dat 1SCF charge= 1 PM6", "MOLECULE # 1", "",
		"H 0.0000000000 0.0000000000 -5.0705374640",
		"C 0.0000000000 0.0000000000 -3.0439358050",
		"C 0.0000000000 0.0000000000 -0.7125467130",
		"C 0.0000000000 -1.2030187010 1.8800877760",
		"C 0.0000000000 1.2030187010 1.8800877760",
		"H 0.0000000000 -3.1388743520 2.5132649210",
		"H 0.0000000000 3.1388743520 2.5132649210"}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %#v, want %#v", got, want)
	}
}

func TestWriteMopacIn(t *testing.T) {
	got := WriteMopacIn(testinp, testgeom[0], 1)
	want := BaseMopFilename + "00001.mop"
	if got != want {
		t.Errorf("got %s, want %s", got, want)
	}
}

func TestSlurmSubmit(t *testing.T) {
	got := strconv.Itoa(SlurmSubmit(BaseMopFilename + "00001.mop"))
	re := regexp.MustCompile(`[0-9]{5,8}`)
	if !re.MatchString(got) {
		t.Errorf("got %q, wanted %v", got, re)
	}
}

// func TestReadMopacOut(t *testing.T) {
// 	job := Job{"test.aux", 11111, "queued", 0}
// 	t.Run("Job ran successfully", func(t *testing.T) {
// 		got, _ := ReadMopacOut(job)
// 		want := 0.277615868297292e+78
// 		if got != want {
// 			t.Errorf("got %f, wanted %f", got, want)
// 		}
// 	})

// 	t.Run("Job failed but output produced", func(t *testing.T) {
// 		job := Job{"fail.aux", 11111, "queued", 0}
// 		got, _ := ReadMopacOut(job)
// 		re := regexp.MustCompile(`[0-9]{5,8}`)
// 		strgot := strconv.Itoa(got.(int))
// 		if !re.MatchString(strgot) {
// 			t.Errorf("got %q, wanted %v", got, re)
// 		}
// 	})

// 	t.Run("Job failed, output produced, but out of retries", func(t *testing.T) {
// 		job := Job{"fail.aux", 11111, "queued", maxretries + 1}
// 		_, err := ReadMopacOut(job)
// 		want := errTooManyRetries
// 		if err != want {
// 			t.Errorf("got %s, wanted %s", err, want)
// 		}
// 	})

// 	t.Run("Job failed and no output", func(t *testing.T) {
// 		job := Job{"fail1.aux", 11111, "queued", 0}
// 		_, err := ReadMopacOut(job)
// 		if err == nil {
// 			t.Errorf("Wanted an error but didn't get one")
// 		}
// 	})
// }

// func TestWriteParams(t *testing.T) {
// 	t.Run("No error", func(t *testing.T) {
// 		err := testinp.Param.Write("test-params.dat")
// 		if err != nil {
// 			t.Errorf("Didn't want an error but got one")
// 		}
// 	})
// }

// func TestBasename(t *testing.T) {
// 	got := Basename("/home/brent/Projects/sempirical/main.inp")
// 	want := "main"
// 	if got != want {
// 		t.Errorf("got %s, wanted %s", got, want)
// 	}
// }

// func TestFcn(t *testing.T) {
// 	got := Fcn(&testinp, []float64{-153.506766246758, -153.506744023064},
// 		testgeom)
// 	want := []float64{4.1738645326326294e+06, 4.1738585662503745e+06}
// 	if !reflect.DeepEqual(got, want) {
// 		t.Errorf("got %v, wanted %v", got, want)
// 	}
// }
	
// func TestEnorm(t *testing.T) {
// 	got := Enorm([]float64{3, 4})
// 	want := 5.0
// 	if got != want {
// 		t.Errorf("got %f, wanted %f", got, want)
// 	}
// }

// func TestLmdif(t *testing.T) {
// 	got := Lmdif(testinp, []float64{-153.506766246758, -153.506744023064},
// 		testgeom)
// }

// func TestFdjac2(t *testing.T) {
// 	inp := testinp
// 	fvec := Fcn(&inp, []float64{-153.506766246758, -153.506744023064},
// 		testgeom)
// 	got := Fdjac2(&inp, []float64{-153.506766246758, -153.506744023064},
// 		testgeom, fvec)
// 	fmt.Println(got)
// }

func TestQrfac(t *testing.T) {
	// A := [][]float64{{12, -51, 4}, {6, 167, -68}, {-4, 24, -41}}
	// Ap := []float64{12, -51, 4, 6, 167, -68, -4, 24, -41}
	Ap2 := []float64{-1, -1, 1, 1, 3, 3, -1, -1, 5, 1, 3, 7}
	// R := [][]float64{{14, 21, -14}, {0, 175, -70}, {0, 0, -35}}
	var qr mat.QR
	var r mat.Dense
	matA := mat.NewDense(4, 3, Ap2)
	// got := Qrfac(A)
	// want := R
	fmt.Println(matA)
	qr.Factorize(matA)
	qr.RTo(&r)
	fmt.Println(r)
}
