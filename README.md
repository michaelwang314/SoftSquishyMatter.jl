# SoftSquishyMatter.jl
A simple package for simulating soft matter things.  At the moment, the package can simulate simple systems of interacting passive Brownian or systems of active particles (active Brownian and run-and-tumble).  There are plans for adding additional interactions, integrators, and more.  Also, hopefully the code can be further optimized in the future.  The goal is to create a package written in Julia that is easily customizable and readable (at the cost of some performance probably).

This is mostly just a fun/educational and (hopefully) long-term project!

<p align="center">
    <img src="https://github.com/michaelwang314/simulation-storage/blob/master/Example_gifs/Example_LennardJonesFluid.gif" width="300" height="300">
    <img src="https://github.com/michaelwang314/simulation-storage/blob/master/Example_gifs/Example_TwoTemperature.gif" width="300" height="300">
    <img src="https://github.com/michaelwang314/simulation-storage/blob/master/Example_gifs/Example_PolyDispersed.gif" width="300" height="300">
    <img src="https://github.com/michaelwang314/simulation-storage/blob/master/Example_gifs/Example_RunAndTumble.gif" width="300" height="300">
    <img src="https://github.com/michaelwang314/simulation-storage/blob/master/Example_gifs/Example_RunAndTumbleFunnels.gif" width="400" height="400">
</p>

More simulations can be found in [this repository](https://github.com/michaelwang314/simulation-storage).

## Installation
To download, open Julia command line, press `]`, and use `add`
```
(v1.4) pkg> add https://github.com/michaelwang314/SoftSquishyMatter.jl
```
or do
```
julia> using Pkg
julia> Pkg.add("https://github.com/michaelwang314/SoftSquishyMatter.jl")
```
You may also need to add [Plots.jl](http://docs.juliaplots.org/latest/).  I guess I should register this at some point so it's easier to add the package.

## Usage
Just add the line
```
using SoftSquishyMatter
```
There are several example simulations in the [Example folder](https://github.com/michaelwang314/SoftSquishyMatter.jl/tree/master/Examples).
