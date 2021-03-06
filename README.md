# Metaplex Networks

A collection of programs demonstrating the metaplex model.

## Installation
All the files are in Julia, so you'll need to install it from [here](https://julialang.org/).

The original code in the `old` folder was reworked into a [Julia Package](https://github.com/csimal/NetworkEpidemics.jl). You'll need it to run the notebooks. To install it, simply enter the following in the Julia REPL:
```julia
] add https://github.com/csimal/NetworkEpidemics.jl
```
Several community packages are also used. You can install them by running the `install.jl` file, either by entering `include("install.jl")` in the REPL, or by running `julia install.jl` in the command line.

## Usage
The notebooks are self contained, and generate the figures of the paper (to be published). The `old` folder contains the original code of the project. Ignore it unless you want to see simple implementations of the various models. The package made from it is considerably more complex and might not be as easy to understand.