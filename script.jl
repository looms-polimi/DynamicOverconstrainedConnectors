#=
  Test script for DOCC
=#
using Revise
using Plots
using Test

import OM

global file = "DynamicOverconstrainedConnectors.mo"

function plotPartial(sols,
                     vars,
                     tspan,
                     file)
  #= Separate plot=#
  local subPlots = []
  for sol in sols
    p = plot(sol; vars = vars)
    push!(subPlots, p)
  end
  fp = plot(subPlots...)
  Plots.pdf(fp, file)
  #= Joint plot=#
  p = plot(first(sols); vars = vars, xlim = first(tspan), ylim=last(tspan))
  for i in 2:length(sols)
    p = plot!(sols[i]; vars = vars, xlim = first(tspan), ylim=last(tspan))
  end    
  Plots.pdf(p, file * "Combined")
end

include("powerGrids.jl")
include("incompressibleFluidNetworks.jl")

#solver = :(Rodas5())
#simulateAndPlotSystem("System2"; solver = solver)
#= Simulate System 3=#
#simulateAndPlotSystem("System3"; solver = solver)

#println("Testing System 2")
#res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System2", file; MSL = true)
#println(res)

# println("Testing System 3")
# res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System3", file; MSL = true)
# println(res)

# println("Testing System 4")
# res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System4", file; MSL = true)
# println(res)


