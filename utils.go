package main

import (
	"io/ioutil"
	"strings"
)

func ReadFile(filename string) []string {
	text, err := ioutil.ReadFile(filename)
	if err != nil {
		panic(err)
	}
	split := strings.Split(string(text), "\n")
	for i, line := range split {
		split[i] = strings.TrimSpace(line)
	}
	return split
}
