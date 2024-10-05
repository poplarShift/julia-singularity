#!/usr/bin/env julia
using Pkg
Pkg.add([
    Pkg.PackageSpec(name="ComponentArrays", version="0.14.1"),
    Pkg.PackageSpec(name="CSV", version="0.10.11"),
    Pkg.PackageSpec(name="DataFrames", version="1.6.1"),
    Pkg.PackageSpec(name="DataFramesMeta", version="0.14.1"),
    Pkg.PackageSpec(name="DiffEqParamEstim", version="2.0.1"),
    Pkg.PackageSpec(name="DiffEqSensitivity", version="6.79.0"),
    Pkg.PackageSpec(name="DifferentialEquations", version="7.7.0"),
    Pkg.PackageSpec(name="Interpolations", version="0.14.7"),
    Pkg.PackageSpec(name="ForwardDiff", version="0.10.36"),
    Pkg.PackageSpec(name="ModelingToolkit", version="8.50.0"),
    Pkg.PackageSpec(name="Optimization", version="3.14.0"),
    Pkg.PackageSpec(name="OptimizationBBO", version="0.1.4"),
    Pkg.PackageSpec(name="OptimizationCMAEvolutionStrategy", version="0.1.3"),
    Pkg.PackageSpec(name="OptimizationOptimJL", version="0.1.8"),
    Pkg.PackageSpec(name="OptimizationPolyalgorithms", version="0.1.1"),
    Pkg.PackageSpec(name="Plots", version="1.39.0"),
    Pkg.PackageSpec(name="RecursiveArrayTools", version="2.38.10"),
    Pkg.PackageSpec(name="SciMLSensitivity", version="7.29.0"),
    Pkg.PackageSpec(name="StructuralIdentifiability", version="0.4.10"),
    Pkg.PackageSpec(name="Zygote", version="0.6.64"),
])