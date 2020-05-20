package main

import (
	"regexp"
	"strings"
)

type Key int

const (
	Nparam Key = iota
	Method
	Nmolec
	Charge
	Nstruct
	Efile
	Gfile
	Units
	Geometry
	Params
	Nkeys
)

type Regexp struct {
	Expr *regexp.Regexp
	Name Key
}

func ParseInfile(filename string) {
	lines := ReadFile(filename)
	Keywords := []Regexp{
		Regexp{regexp.MustCompile(`(?i)nparam=`), Nparam},
		Regexp{regexp.MustCompile(`(?i)method=`), Method},
		Regexp{regexp.MustCompile(`(?i)nmolec=`), Nmolec},
		Regexp{regexp.MustCompile(`(?i)charge=`), Charge},
		Regexp{regexp.MustCompile(`(?i)nstruct=`), Nstruct},
		Regexp{regexp.MustCompile(`(?i)efile=`), Efile},
		Regexp{regexp.MustCompile(`(?i)gfile=`), Gfile},
		Regexp{regexp.MustCompile(`(?i)units=`), Units},
	}
	geom := regexp.MustCompile(`(?i)geometry={`)
	params := regexp.MustCompile(`(?i)params={`)
	for i := 0; i < len(lines); {
		if len(lines[i]) > 0 && lines[i][0] == '#' {
			i++
			continue
		}
		if geom.MatchString(lines[i]) {
			i++
			geomlines := make([]string, 0)
			for !strings.Contains(lines[i], "}") {
				geomlines = append(geomlines, lines[i])
				i++
			}
			Input[Geometry] = strings.Join(geomlines, "\n")
		} else if params.MatchString(lines[i]) {
			i++
			paramlines := make([]string, 0)
			for !strings.Contains(lines[i], "}") {
				paramlines = append(paramlines, lines[i])
				i++
			}
			Input[Params] = strings.Join(paramlines, "\n")
		} else {
			for _, kword := range Keywords {
				if kword.Expr.MatchString(lines[i]) {
					split := strings.Split(lines[i], "=")
					Input[kword.Name] = split[len(split)-1]
				}
			}
			i++
		}
	}
}

func GetAtomNames() (names []string) {
	lines := strings.Split(Input[Geometry], "\n")
	for _, line := range lines {
		if len(line) >= 1 {
			split := strings.Fields(line)
			names = append(names, split[0])
		}
	}
	return
}
			
func ReadGfile(filename string) []string {
	lines := strings.Join(ReadFile(filename), "\n")
	// skip empty string before first split
	return strings.Split(lines, "# GEOMUP #################\n")[1:]
}
