package main

import (
	"reflect"
	"testing"
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
	params, headers := FloatParams()
	got := MakeParams(params[:2], headers[:2])
	want := "FN11 H      0.024184000000\n" +
		"FN11 C      0.046302000000\n"
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}
