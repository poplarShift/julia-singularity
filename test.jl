#!/usr/bin/env julia
using Pkg
Pkg.installed()

using DifferentialEquations, RecursiveArrayTools, Plots, DiffEqParamEstim, Optimization, ForwardDiff, OptimizationOptimJL, OptimizationBBO, Plots, CSV, DataFrames, Random, Zygote, OptimizationPolyalgorithms, SciMLSensitivity, ModelingToolkit, LinearAlgebra, OptimizationCMAEvolutionStrategy

