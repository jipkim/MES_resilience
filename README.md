# MES_resilience

This code was developed under Julia v0.6 by Jip Kim in 2018.
The following packages must be installed:

  - Gurobi
  - Graphs
  - JuMP
  - JLD  

This code includes two execution files:

PH.jl: Solving the proposed optimization problem using the progressive hedging algorithm
BF.jl: Solving the proposed optimization problem without any algorithmic enhancement

To run the code, execute PH.jl/BF.jl, or include() it from a Julia prompt.

The input data is given as follows:

  - load24data-PJM.csv: 24 hour load profile
  - DataMiner-Export_2018-01-10-121427.csv: 24 hour LMP data
  - Node.csv: Distribution network data
  - Line.csv: Distribution line data
  - Generator.csv: Distributed generation data
  - Storage1MWhx3.csv: ES unit data
  - Station.csv: ES unit charging station data  
  - We also specify additional input data in the beginning of PH.jl/BF.jl

The results of the solve are saved as result.jld file using JLD(Julia Data format) package.
