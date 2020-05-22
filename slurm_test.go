package main

import (
	"reflect"
	"testing"
)

var (
	slurm = Slurm{}
)

func TestMakeSlurmHead(t *testing.T) {
	got := slurm.MakeHead()
	want := []string{
		"#!/bin/bash",
		"#SBATCH --job-name=sempirical",
		"#SBATCH --ntasks=1",
		"#SBATCH --cpus-per-task=1",
		"#SBATCH -o /dev/null",
		"#SBATCH --no-requeue",
		"#SBATCH --mem=50mb"}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func TestMakeSlurm(t *testing.T) {
	got := slurm.Make("inp/Structure00005.mop")
	want := make([]string, 0)
	want = append(append(want, slurm.MakeHead()...),
		"/home/qc/mopac2016/mopac.sh inp/Structure00005.mop")
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func TestSubmit(t *testing.T) {
	got := slurm.Submit("test")
	want := "test"
	if got != want {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}
