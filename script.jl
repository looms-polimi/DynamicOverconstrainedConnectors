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
  p = plot(first(sols); vars = vars, limits = tspan)
  for i in 2:length(sols)
    p = plot!(sols[i]; vars = vars, limits = tspan)
  end    
  Plots.pdf(p, file * "Combined")
end

include("powerGrids.jl")
include("incompressibleFluidNetworks.jl")

function plotSystem4Special1()
  solver = :(Rodas5())  
  sol = OM.simulate("DynamicOverconstrainedConnectors.PowerGridsReal.System4",
                    file;
                    startTime = 0.0,
                    stopTime = 50,
                    MSL = true,
                    solver = solver);
  p1 = plot(sol[1]; vars = [OM.OMBackend.G2_port_v_re], tspan = (0.0, 10.0))
  p2 = plot(sol[2]; vars = [OM.OMBackend.G2_port_v_re], tspan = (10.0, 50.0))
  p = plot(p1, p2)
  Plots.pdf(p, "System4_G2_port_v_re_paper")
end



function plotSystem8Special1()
  solver = :(Rodas5())  
  sol = OM.simulate("DynamicOverconstrainedConnectors.PowerGridsReal.System8",
                    file;
                    startTime = 0.0,
                    stopTime = 50,
                    MSL = true,
                    solver = solver);
  p1 = plot(sol[1]; vars = [OM.OMBackend.G2_port_v_re], tspan = (0.0, 10.0))
  p2 = plot(sol[2]; vars = [OM.OMBackend.G2_port_v_re], tspan = (10.0, 50.0))
  p = plot(p1, p2)
  Plots.pdf(p, "System8_G2_port_v_re_paper")
end


function plotSystem8SpecialT3_Pab()
  solver = :(Rodas5())  
  sol = OM.simulate("DynamicOverconstrainedConnectors.PowerGridsReal.System8",
                    file;
                    startTime = 0.0,
                    stopTime = 50,
                    MSL = true,
                    solver = solver);
  p1 = plot(sol[1]; vars = [OM.OMBackend.T3_Pab], tspan = (0.0, 10.0))
  p2 = plot(sol[2]; vars = [OM.OMBackend.T3_Pab], tspan = (10.0, 50.0))
  p = plot(p1, p2)
  Plots.pdf(p, "System8_T3_Pab_v_re_paper")
end



function plotSystem8SpecialG3_port_v_re()
  solver = :(Rodas5())  
  sol = OM.simulate("DynamicOverconstrainedConnectors.PowerGridsReal.System8",
                    file;
                    startTime = 0.0,
                    stopTime = 50,
                    MSL = true,
                    solver = solver);
  p1 = plot(sol[1]; vars = [OM.OMBackend.G3_port_v_re], tspan = (0.0, 10.0))
  p2 = plot(sol[2]; vars = [OM.OMBackend.G3_port_v_re], tspan = (10.0, 50.0))
  p = plot(p1, p2)
  Plots.pdf(p, "System8_G3_port_v_re_paper")
end



function plotSystem7G3_v_re()
  solver = :(Rodas5())  
  sol = OM.simulate("DynamicOverconstrainedConnectors.PowerGridsReal.System7",
                    file;
                    startTime = 0.0,
                    stopTime = 50,
                    MSL = true,
                    solver = solver);
  p1 = plot(sol; vars = [OM.OMBackend.G3_port_v_re], tspan = (0.0, 10.0))
  p2 = plot(sol; vars = [OM.OMBackend.G3_port_v_re], tspan = (10.0, 50.0))
  p = plot(p1, p2)
  Plots.pdf(p, "System7_G3_port_v_re_paper")
end


function plotSystem7G2_v_re()
  solver = :(Rodas5())  
  sol = OM.simulate("DynamicOverconstrainedConnectors.PowerGridsReal.System7",
                    file;
                    startTime = 0.0,
                    stopTime = 50,
                    MSL = true,
                    solver = solver);
  p1 = plot(sol; vars = [OM.OMBackend.G2_port_v_re], tspan = (0.0, 10.0))
  p2 = plot(sol; vars = [OM.OMBackend.G2_port_v_re], tspan = (10.0, 50.0))
  p = plot(p1, p2)
  Plots.pdf(p, "System7_G2_port_v_re_paper")
end


function plotSystem7Special2()
  solver = :(Rodas5())  
  sol = OM.simulate("DynamicOverconstrainedConnectors.PowerGridsReal.System7",
                    file;
                    startTime = 0.0,
                    stopTime = 50,
                    MSL = true,
                    solver = solver);
  p1 = plot(sol; vars = [OM.OMBackend.G1_Pe, OM.OMBackend.G2_Pe, OM.OMBackend.G3_Pe], tspan = (0.0, 10.0))
  p2 = plot(sol; vars = [OM.OMBackend.G1_Pe, OM.OMBackend.G2_Pe, OM.OMBackend.G3_Pe], tspan = (10.0, 50.0))
  p = plot(p1, p2)
  Plots.pdf(p, "System7_G1_G2_G3_port_v_re_paper")
end


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


