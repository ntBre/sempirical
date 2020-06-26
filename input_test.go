package main

import (
	"reflect"
	"testing"
)

func TestParseInfile(t *testing.T) {
	ParseInfile("new.inp")
	filled := [Nkeys]string{
		"27", "PM6", "1", "1", "10353", "energy.dat", "file07", "ANGSTROMS", "",
		"H        0.0000000000        0.0000000000       -5.0705374640\n" +
			"C        0.0000000000        0.0000000000       -3.0439358050\n" +
			"C        0.0000000000        0.0000000000       -0.7125467130\n" +
			"C        0.0000000000       -1.2030187010        1.8800877760\n" +
			"C        0.0000000000        1.2030187010        1.8800877760\n" +
			"H        0.0000000000       -3.1388743520        2.5132649210\n" +
			"H        0.0000000000        3.1388743520        2.5132649210",
		"FN11    H    0.02418400    CORE-CORE    VDW    MULTIPLIER    1\n" +
			"FN11    C    0.04630200    CORE-CORE    VDW    MULTIPLIER    1\n" +
			"XFAC_C    H    0.21650600    XFAC    factor\n" +
			"HSP    C    0.71732200    EV    ONE-CENTER    INTEGRAL    (SP,SP)\n" +
			"XFAC_C    C    0.81351000    XFAC    factor\n" +
			"ALPB_C    H    1.02780600    ALPB    factor\n" +
			"GPP    C    10.77832600    EV    ONE-CENTER    INTEGRAL    (PP,PP)\n" +
			"USS    H    -11.24695800    EV    ONE-CENTER    ENERGY    FOR    S\n" +
			"GSP    C    11.52813400    EV    ONE-CENTER    INTEGRAL    (SS,PP)\n" +
			"ZS    H    1.26864100    AU    ORBITAL    EXPONENT    FOR    S\n" +
			"GSS    C    13.33551900    EV    ONE-CENTER    INTEGRAL    (SS,SS)\n" +
			"FN31    C    1.33395900    CORE-CORE    VDW    POSITION    1\n" +
			"GSS    H    14.44868600    EV    ONE-CENTER    INTEGRAL    (SS,SS)\n" +
			"BETAS    C    -15.38523600    EV    BETA    PARAMETER    FOR    S\n" +
			"ZP    C    1.70284100    AU    ORBITAL    EXPONENT    FOR    P\n" +
			"FN31    H    1.78601100    CORE-CORE    VDW    POSITION    1\n" +
			"ZS    C    2.04755800    AU    ORBITAL    EXPONENT    FOR    S\n" +
			"FN21    C    2.10020600    CORE-CORE    VDW    EXPONENT    1\n" +
			"XFAC_H    H    2.24358700    XFAC    factor\n" +
			"ALPB_C    C    2.61371300    ALPB    factor\n" +
			"FN21    H    3.05595300    CORE-CORE    VDW    EXPONENT    1\n" +
			"ALPB_H    H    3.54094200    ALPB    factor\n" +
			"UPP    C    -39.93792000    EV    ONE-CENTER    ENERGY    FOR    P\n" +
			"USS    C    -51.08965300    EV    ONE-CENTER    ENERGY    FOR    S\n" +
			"BETAP    C    -7.47192900    EV    BETA    PARAMETER    FOR    P\n" +
			"BETAS    H    -8.35298400    EV    BETA    PARAMETER    FOR    S\n" +
			"GP2    C    9.48621200    EV    ONE-CENTER    INTEGRAL    (PP*,PP*)"}
	if !reflect.DeepEqual(filled, Input) {
		t.Errorf("\ngot %#v\nwad %#v\n", Input, filled)
	}
}

func TestGetAtomNames(t *testing.T) {
	ParseInfile("new.inp")
	got := GetAtomNames()
	want := []string{"H", "C", "C", "C", "C", "H", "H"}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func TestReadGfile(t *testing.T) {
	got := ReadGfile("testfiles/min07")
	want := []string{
		"0.0000000000        0.0000000000       -5.0705374640\n" +
			"0.0000000000        0.0000000000       -3.0439358050\n" +
			"0.0000000000        0.0000000000       -0.7125467130\n" +
			"0.0000000000       -1.2030187010        1.8800877760\n" +
			"0.0000000000        1.2030187010        1.8800877760\n" +
			"0.0000000000       -3.1388743520        2.5132649210\n" +
			"0.0000000000        3.1388743520        2.5132649210\n",
		"-0.0000000000        0.0000000000       -5.0706091104\n" +
			"-0.0000000000        0.0000000000       -3.0440074514\n" +
			"0.0000000000        0.0000000000       -0.7126183594\n" +
			"0.0140780719       -1.2029363254        1.8800161296\n" +
			"-0.0140780719        1.2029363254        1.8800161296\n" +
			"-0.0053955418       -3.1386119859        2.5134440370\n" +
			"0.0053955418        3.1386119859        2.5134440370\n"}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %q, wanted %q\n", got, want)
	}
}

func TestFloatParams(t *testing.T) {
	params, headers := FloatParams(Input[Params])
	want := []float64{
		0.02418400, 0.04630200, 0.21650600,
		0.71732200, 0.81351000, 1.02780600,
		10.77832600, -11.24695800, 11.52813400,
		1.26864100, 13.33551900, 1.33395900,
		14.44868600, -15.38523600, 1.70284100,
		1.78601100, 2.04755800, 2.10020600,
		2.24358700, 2.61371300, 3.05595300,
		3.54094200, -39.93792000, -51.08965300,
		-7.47192900, -8.35298400, 9.48621200}
	if !reflect.DeepEqual(params, want) {
		t.Errorf("got %v, wanted %v\n", params, want)
	}
	want2 := []string{
		"FN11 H", "FN11 C", "XFAC_C H",
		"HSP C", "XFAC_C C", "ALPB_C H",
		"GPP C", "USS H", "GSP C",
		"ZS H", "GSS C", "FN31 C",
		"GSS H", "BETAS C", "ZP C",
		"FN31 H", "ZS C", "FN21 C",
		"XFAC_H H", "ALPB_C C", "FN21 H",
		"ALPB_H H", "UPP C", "USS C",
		"BETAP C", "BETAS H", "GP2 C"}
	if !reflect.DeepEqual(headers, want2) {
		t.Errorf("got %q, wanted %q\n", headers, want2)
	}
}
