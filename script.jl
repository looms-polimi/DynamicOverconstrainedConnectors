#=
  Test script for DOCC
=#
using Revise
using Plots
using Test

import OM


file = "DynamicOverconstrainedConnectors.mo"

function simulateAndPlotSystem(modelName; solver)
  println("Testing $(modelName)")
  sol = OM.simulate("DynamicOverconstrainedConnectors.PowerGridsReal.$(modelName)",
                    file;
                    startTime = 0.0,
                    stopTime = 50.0,
                    MSL = true,
                    solver = solver)
  #= Ω =#
  p = plot(sol; vars = [OM.OMBackend.G1_omega])
  Plots.pdf(p, "./Plots/$(modelName)/$(modelName)G1_omega")
  p = plot(sol; vars = [OM.OMBackend.G2_omega])
  Plots.pdf(p, "./Plots/$(modelName)/$(modelName)G2_omega")
  #= θ =#
  p = plot(sol; vars = [OM.OMBackend.G1_theta])
  Plots.pdf(p, "./Plots/$(modelName)/$(modelName)G1_theta")
  p = plot(sol; vars = [OM.OMBackend.G2_theta])
  Plots.pdf(p, "./Plots/$(modelName)/$(modelName)G2_theta")
  #= Power output =#
  p = plot(sol; vars = [OM.OMBackend.G1_Pe])
  Plots.pdf(p, "./Plots/$(modelName)/$(modelName)G1_Pe")
  p = plot(sol; vars = [OM.OMBackend.G2_Pe])
  Plots.pdf(p, "./Plots/$(modelName)/$(modelName)G2_Pe")
  #= Primary frequency control power in per unit =#
  p = plot(sol; vars = [OM.OMBackend.G1_Pc])
  Plots.pdf(p, "./Plots/$(modelName)/$(modelName)G1_Pc")
  p = plot(sol; vars = [OM.OMBackend.G2_Pc])
  Plots.pdf(p, "./Plots/$(modelName)/$(modelName)G2_Pc")
  #= Active power output set point =#
  p = plot(sol; vars = [OM.OMBackend.G1_Ps])
  Plots.pdf(p, "./Plots/$(modelName)/$(modelName)G1_Ps")
  p = plot(sol; vars = [OM.OMBackend.G2_Ps])
  Plots.pdf(p, "./Plots/$(modelName)/$(modelName)G2_Ps")
  #= P variable  =#
    p = plot(sol; vars = [OM.OMBackend.L1_P])
  Plots.pdf(p, "./Plots/$(modelName)/$(modelName)L1_P")
  p = plot(sol; vars = [OM.OMBackend.L2_P])
  Plots.pdf(p, "./Plots/$(modelName)/$(modelName)L2_P")
end


#OM.LogBackend()
#println("Testing IF_EQ")
#sol = OM.simulate("IfEquationDer", "./IfEquationDer.mo"; MSL = true, startTime = 0.0, stopTime = 12.0)
#@assert last(last(sol.u)) == 43.99999999999987
#println("END Testing IF_EQ")
res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System1", file; MSL = true)
println(res)
solver = :(Rosenbrock23())
simulateAndPlotSystem("System1"; solver = solver)
#solver = :(Rodas5())
simulateAndPlotSystem("System2"; solver = solver)
#= Simulate System 3=#
simulateAndPlotSystem("System3"; solver = solver)

#println("Testing System 2")
#res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System2", file; MSL = true)
#println(res)

# println("Testing System 3")
# res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System3", file; MSL = true)
# println(res)

# println("Testing System 4")
# res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System4", file; MSL = true)
# println(res)


