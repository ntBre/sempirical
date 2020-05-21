package main

import (
	"reflect"
	"strings"
	"testing"
)

func TestMakeHead(t *testing.T) {
	got := mop.MakeHead()
	want := []string{"threads=1 XYZ ANGSTROMS scfcrt=1.D-21 aux(precision=9)" +
		" external=params.dat 1SCF charge=1 PM6",
		"MOLECULE # 1", ""}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %q, wanted %q\n", got, want)
	}
}

func TestMakeIn(t *testing.T) {
	names := GetAtomNames()
	coords := strings.Fields(ReadGfile("testfiles/min07")[0])
	got := mop.MakeIn(names, coords)
	want := []string{"threads=1 XYZ ANGSTROMS scfcrt=1.D-21 aux(precision=9)" +
		" external=params.dat 1SCF charge=1 PM6",
		"MOLECULE # 1", "",
		"H 0.0000000000 0.0000000000 -5.0705374640",
		"C 0.0000000000 0.0000000000 -3.0439358050",
		"C 0.0000000000 0.0000000000 -0.7125467130",
		"C 0.0000000000 -1.2030187010 1.8800877760",
		"C 0.0000000000 1.2030187010 1.8800877760",
		"H 0.0000000000 -3.1388743520 2.5132649210",
		"H 0.0000000000 3.1388743520 2.5132649210"}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func TestReadOut(t *testing.T) {
	t.Run("file not found", func(t *testing.T) {
		_, err := mop.ReadOut("testfiles/not_a_real_file.aux")
		if err != ErrFileNotFound {
			t.Errorf("got %s, wanted %s\n", err, ErrFileNotFound)
		}
	})
	t.Run("blank file", func(t *testing.T) {
		_, err := mop.ReadOut("testfiles/blank.aux")
		if err != ErrBlankOutput {
			t.Errorf("got %s, wanted %s\n", err, ErrBlankOutput)
		}
	})
	t.Run("one line error", func(t *testing.T) {
		_, err := mop.ReadOut("testfiles/e1line.aux")
		if err != ErrFileContainsError {
			t.Errorf("got %s, wanted %s\n", err, ErrFileContainsError)
		}
	})
	t.Run("energy found", func(t *testing.T) {
		got, err := mop.ReadOut("testfiles/test.aux")
		want := 0.277615868297292e+78
		if got != want {
			t.Errorf("got %v, wanted %v\n", got, want)
		}
		if err != nil {
			t.Errorf("got error %q, but didn't want one", err)
		}
	})
}
