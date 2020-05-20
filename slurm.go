package main

import (
	"io/ioutil"
	"os/exec"
	"strings"
	"time"
)

type Slurm struct{}

func (s Slurm) MakeHead() []string {
	return []string{
		"#!/bin/bash",
		"#SBATCH --job-name=sempirical",
		"#SBATCH --ntasks=1",
		"#SBATCH --cpus-per-task=1",
		"#SBATCH -o /dev/null",
		"#SBATCH --no-requeue",
		"#SBATCH --mem=50mb"}
}

func (s Slurm) Make(filename string) (lines []string) {
	body := []string{"/home/qc/mopac2016/mopac.sh " + filename}
	return append(append(lines, s.MakeHead()...), body...)
}

func (s Slurm) Write(qfile, chemfile string) {
	lines := s.Make(chemfile)
	writelines := strings.Join(lines, "\n")
	err := ioutil.WriteFile(qfile, []byte(writelines), 0755)
	if err != nil {
		panic(err)
	}
}

func (s Slurm) Submit(filename string) {
	_, err := exec.Command("sbatch", filename).Output()
	// have to use sbatch because srun grabs a whole node
	// and runs interactively
	for err != nil {
		time.Sleep(time.Second)
		_, err = exec.Command("sbatch", filename).Output()
	}
}
