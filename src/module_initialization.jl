# Pkg Module
using .Sys
using .Threads
using CSV
using DataFrames
using DataStructures
using Dates
using FilesystemDatastructures
using ForwardDiff               # Automatic differentiation package
using HDF5
using Interpolations
using LaTeXStrings
using LinearAlgebra
using Logging
using LoggingExtras
using NearestNeighbors
using Pkg
using Plots
using Printf
using PyCall
using QuadGK
using Statistics
using StatsPlots
using StatsBase

"""
    initialize_pyplot_backend()
Initialize the built-in pyplot backend.
"""
function initialize_pyplot_backend()
    push!(pyimport("sys")."path", dirname(pathof(PhantomRevealer)))
    prplt = pyimport("pyplot_backend")
    return prplt
end