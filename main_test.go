package main

import (
	"reflect"
	"testing"
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
	got := SlurmSubmit(BaseMopFilename + "00001.mop")
	want := "--job-name=qff --ntasks=1 --cpus-per-task=1 " +
		"--mem=1gb --output=/dev/null --error=/dev/null " +
		"/home/qc/mopac2016/mopac.sh inp/Structure00001.mop"
	if got != want {
		t.Errorf("got %q, want %q", got, want)
	}
}

func TestReadMopacOut(t *testing.T) {
	t.Run("Job ran successfully", func(t *testing.T) {
		got := ReadMopacOut("test.aux")
		want := 0.277615868297292E+78
		if got != want {
			t.Errorf("got %f, wanted %f", got, want)
		}
	})

	t.Run("Job failed but output produced", func(t *testing.T) {
		got := ReadMopacOut("fail.aux")
		want := "--job-name=qff --ntasks=1 --cpus-per-task=1 " +
		"--mem=1gb --output=/dev/null --error=/dev/null " +
		"/home/qc/mopac2016/mopac.sh fail.mop"
		if got != want {
			t.Errorf("got %s, wanted %s", got, want)
		}
	})

	t.Run("Job failed and no output", func(t *testing.T) {
		err := ReadMopacOut("fail1.aux")
		if err == nil {
			t.Errorf("Wanted an error but didn't get one")
		}
	})
}

func TestWriteParams(t *testing.T) {
	t.Run("No error", func(t *testing.T) {
		err := testinp.Param.Write("test-params.dat")
		if err != nil {
			t.Errorf("Didn't want an error but got one")
		}
	})
}
