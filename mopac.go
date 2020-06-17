package main

import (
	"errors"
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"runtime"
	"strconv"
	"strings"
)

var (
	brokenFloat            = math.NaN()
	ErrEnergyNotFound      = errors.New("Energy not found in MOPAC output")
	ErrFileNotFound        = errors.New("MOPAC output file not found")
	ErrEnergyNotParsed     = errors.New("Energy not parsed in MOPAC output")
	ErrFinishedButNoEnergy = errors.New("MOPAC output finished but no energy found")
	ErrFileContainsError   = errors.New("MOPAC output file contains an error")
	ErrBlankOutput         = errors.New("MOPAC output file exists but is blank")
	ErrInputGeomNotFound   = errors.New("Geometry not found in input file")
	ErrTimeout             = errors.New("Timeout waiting for signal")
)

type Mopac struct{}

func (m Mopac) MakeHead() (headLines []string) {
	return []string{"threads=1 XYZ " + Input[Units] + " scfcrt=1.D-21 aux(precision=9) SPARKLE " +
		"external=params.dat 1SCF charge=" + Input[Charge] + " " + Input[Method],
		"MOLECULE # 1", ""}
}

func (m Mopac) MakeIn(names []string, coords []string) (inputLines []string) {
	body := make([]string, 0)
	for i, _ := range names {
		tmp := make([]string, 0)
		tmp = append(tmp, names[i])
		for _, c := range coords[3*i : 3*i+3] {
			tmp = append(tmp, c)
		}
		body = append(body, strings.Join(tmp, " "))
	}
	return append(append(inputLines, m.MakeHead()...), body...)
}

func (m Mopac) WriteIn(filename string, names []string, coords []string) {
	lines := m.MakeIn(names, coords)
	writelines := strings.Join(lines, "\n")
	err := ioutil.WriteFile(filename, []byte(writelines), 0755)
	if err != nil {
		panic(err)
	}
}

func (m Mopac) ReadAux(filename string) (result float64, err error) {
	// copied from go-cart molpro, update accordingly
	runtime.LockOSThread()
	if _, err = os.Stat(filename); os.IsNotExist(err) {
		runtime.UnlockOSThread()
		return brokenFloat, ErrFileNotFound
	}
	err = ErrEnergyNotFound
	result = brokenFloat
	lines := ReadFile(filename)
	// ASSUME blank file is only created when PBS runs
	// blank file has a single newline
	if len(lines) == 1 {
		if strings.Contains(strings.ToUpper(lines[0]), "ERROR") {
			return result, ErrFileContainsError
		}
		return result, ErrBlankOutput
	}
	for _, line := range lines {
		if strings.Contains(strings.ToUpper(line), "ERROR") {
			return result, ErrFileContainsError
		}
		if strings.Contains(line, energyLine) {
			split := strings.Split(line, "=")
			if len(split) == 2 {
				err = nil
				result, err = strconv.ParseFloat(
					strings.Replace(split[1], "D", "E", -1), 64)
				if err != nil {
					err = ErrEnergyNotParsed
				}
			}
		}
		if strings.Contains(line, mopacTerminated) && err != nil {
			err = ErrFinishedButNoEnergy
		}
	}
	runtime.UnlockOSThread()
	return result, err
}

// Check a MOPAC output file for existence and errors
func (m Mopac) CheckOut(filename string) error {
	runtime.LockOSThread()
	if _, err := os.Stat(filename); os.IsNotExist(err) {
		runtime.UnlockOSThread()
		return nil
	}
	lines := ReadFile(filename)
	for _, line := range lines {
		if strings.Contains(strings.ToUpper(line), "ERROR") {
			return fmt.Errorf("CheckOut: error %q in file %s", line,
				filename)
		}
	}
	return nil
}
