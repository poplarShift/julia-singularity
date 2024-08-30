## Julia environment as singularity image file

Running this image will create a directory called `.julia/` in the working directory.    

## Build

sudo singularity build julia.sif julia.def 1> build.log 2>&1

## Usage

E.g.:
```./julia.sif myscript.jl options...```