## 

Running this image requires a local directory called `.julia/` in the working directory.    

## Build

sudo singularity build julia.sif julia.def 1> build.log 2>&1

## Usage

E.g.:
```./julia.sif myscript.jl options...```