"""
  Generates Plots for the Electrical System
"""
function simulateAndPlotPowerSystem(modelName; solver)
  println("Testing $(modelName)")
  @time sol = OM.simulate("DynamicOverconstrainedConnectors.PowerGridsReal.$(modelName)",
                          file;
                          startTime = 0.0,
                          stopTime = 50,
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
  #=OMEGA ref=#
  p = plot(sol; vars = [OM.OMBackend.G1_port_omegaRef])
  Plots.pdf(p, "./Plots/$(modelName)/$(modelName)G1_port_omegaRef")
  p = plot(sol; vars = [OM.OMBackend.G2_port_omegaRef])
  Plots.pdf(p, "./Plots/$(modelName)/$(modelName)G2_port_omegaRef")
  #=Omega and Omegaref=#
  p = plot(sol; vars = [OM.OMBackend.G2_omega, OM.OMBackend.G2_port_omegaRef])
  Plots.pdf(p, "./Plots/$(modelName)/$(modelName)G2_omega_and_port_omegaRef")
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
  if modelName == "System3"
    p = plot(sol; vars = [OM.OMBackend.T2_port_a_i_re])
    Plots.pdf(p, "./Plots/$(modelName)/$(modelName)T2_port_a_i_re")
    p = plot(sol; vars = [OM.OMBackend.T2_port_a_i_im])
    Plots.pdf(p, "./Plots/$(modelName)/$(modelName)T2_port_a_i_im")
  end
end

function simulateAndPlotPowerSystemSpecial(modelName; solver)
  println("Testing $(modelName)")
  local startTime = 0.0
  local stopTime = 50.0
  @time sols = OM.simulate("DynamicOverconstrainedConnectors.PowerGridsReal.$(modelName)",
                          file;
                          startTime = startTime,
                          stopTime = stopTime,
                          MSL = true,
                          solver = solver)
  local tspan = (startTime, stopTime)
  #= Ω =#
  p = plotPartial(sols, [OM.OMBackend.G1_omega], tspan, "./Plots/$(modelName)/$(modelName)G1_omega")
  p = plotPartial(sols, [OM.OMBackend.G2_omega], tspan, "./Plots/$(modelName)/$(modelName)G2_omega")
  #= Θ =#
  p = plotPartial(sols, [OM.OMBackend.G1_theta], tspan, "./Plots/$(modelName)/$(modelName)G1_theta")
  p = plotPartial(sols, [OM.OMBackend.G2_theta], tspan, "./Plots/$(modelName)/$(modelName)G2_theta")
  #= Ω Ref =#
  p = plotPartial(sols, [OM.OMBackend.G1_port_omegaRef], tspan, "./Plots/$(modelName)/$(modelName)G1_port_omegaRef")
  p = plotPartial(sols, [OM.OMBackend.G2_port_omegaRef], tspan, "./Plots/$(modelName)/$(modelName)G2_port_omegaRef")
  #= Power Output =#
  p = plotPartial(sols, [OM.OMBackend.G1_Pe], tspan, "./Plots/$(modelName)/$(modelName)G1_Pe")
  p = plotPartial(sols, [OM.OMBackend.G2_Pe], tspan, "./Plots/$(modelName)/$(modelName)G2_Pe")
  #= Primary frequency control power in per unit =#
  p = plotPartial(sols, [OM.OMBackend.G1_Pc], tspan, "./Plots/$(modelName)/$(modelName)G1_Pc")
  p = plotPartial(sols, [OM.OMBackend.G2_Pc], tspan, "./Plots/$(modelName)/$(modelName)G2_Pc")
  #= Active Power output set point =#
  p = plotPartial(sols, [OM.OMBackend.G1_Ps], tspan, "./Plots/$(modelName)/$(modelName)G1_Ps")
  p = plotPartial(sols, [OM.OMBackend.G2_Ps], tspan, "./Plots/$(modelName)/$(modelName)G2_Ps")
  #= P Variable for the loads =#
  p = plotPartial(sols, [OM.OMBackend.L1_P], tspan, "./Plots/$(modelName)/$(modelName)L1_P")
  p = plotPartial(sols, [OM.OMBackend.L2_P], tspan, "./Plots/$(modelName)/$(modelName)L2_P")
  #= The IM and REAL variables  =#
  p = plotPartial(sols, [OM.OMBackend.T2_port_a_i_re], tspan, "./Plots/$(modelName)/$(modelName)T2_port_a_i_re")
  p = plotPartial(sols, [OM.OMBackend.T2_port_a_i_im], tspan, "./Plots/$(modelName)/$(modelName)T2_port_a_i_im")
  p = plotPartial(sols, [OM.OMBackend.G2_port_v_re], tspan, "./Plots/$(modelName)/$(modelName)G2_port_v_re")
  if (modelName == "System8")
    p = plotPartial(sols, [OM.OMBackend.G3_port_v_re], tspan, "./Plots/$(modelName)/$(modelName)G3_port_v_re")
    p = plotPartial(sols, [OM.OMBackend.G3_omega], tspan, "./Plots/$(modelName)/$(modelName)G3_omega")
    p = plotPartial(sols, [OM.OMBackend.G3_theta], tspan, "./Plots/$(modelName)/$(modelName)G3_theta")
    p = plotPartial(sols, [OM.OMBackend.G3_port_omegaRef], tspan, "./Plots/$(modelName)/$(modelName)G3_port_omegaRef")
    p = plotPartial(sols, [OM.OMBackend.T3_Pab_v_re], tspan, "./Plots/$(modelName)/$(modelName)T3_Pab_v_re_paper")
  end
end

function PowerGridsSystem1()
  res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System1", file; MSL = true)
  write("./FlatModels/PowerGridsSystem1.mo", res)
  solver = :(Rodas5())
  simulateAndPlotPowerSystem("System1"; solver = solver)
end

function PowerGridsSystem2()
  res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System2", file; MSL = true)
  write("./FlatModels/PowerGridsSystem2.mo", res)
  solver = :(Rodas5())
  simulateAndPlotPowerSystem("System2"; solver = solver)
end

function PowerGridsSystem3()
  res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System3", file; MSL = true)
  write("./FlatModels/PowerGridsSystem3.mo", res)
  solver = :(Rodas5())
  simulateAndPlotPowerSystem("System3"; solver = solver)
end

function PowerGridsSystem4()
  res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System4", file; MSL = true)
  write("./FlatModels/PowerGridsSystem4.mo", res)
  solver = :(Rodas5())
  simulateAndPlotPowerSystemSpecial("System4"; solver = solver)
end

function PowerGridsSystem5()
  res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System5", file; MSL = true)
  write("./FlatModels/PowerGridsSystem5.mo", res)
  solver = :(Rodas5())
  simulateAndPlotPowerSystem("System5"; solver = solver)
end

function PowerGridsSystem6()
  res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System6", file; MSL = true)
  write("./FlatModels/PowerGridsSystem6.mo", res)
  solver = :(Rodas5())
  simulateAndPlotPowerSystemSpecial("System6"; solver = solver)
end

function PowerGridsSystem7()
  res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System7", file; MSL = true)
  write("./FlatModels/PowerGridsSystem7.mo", res)
  solver = :(Rodas5())
  simulateAndPlotPowerSystem("System7"; solver = solver)
end

function PowerGridsSystem8()
  res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System8", file; MSL = true)
  write("./FlatModels/PowerGridsSystem8.mo", res)
  solver = :(Rodas5())
  simulateAndPlotPowerSystemSpecial("System8"; solver = solver)
end

function PowerGridsSystem9()
  res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System9", file; MSL = true)
  write("./FlatModels/PowerGridsSystem9.mo", res)
  solver = :(Rodas5())
  @time sol = OM.simulate("DynamicOverconstrainedConnectors.PowerGridsReal.System9",
                          file;
                          startTime = 0.0,
                          stopTime = 50,
                          MSL = true,
                          solver = solver)
  if sol.retcode == :MaxIters
    println("System was singular at $(last(sol.t)) as expected")
  else
    println("System was not singular as expected at $(last(sol.t))")
  end
end


function PowerGridsSystem99()
  res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System99", file; MSL = true)
  write("./FlatModels/PowerGridsSystem99.mo", res)
  solver = :(Rodas5())
  @time sol = OM.simulate("DynamicOverconstrainedConnectors.PowerGridsReal.System99",
                          file;
                          startTime = 0.0,
                          stopTime = 50,
                          MSL = true,
                          solver = solver)
  if sol.retcode == :MaxIters
    println("System was singular at $(last(sol.t)) as expected")
  else
    println("System was not singular as expected at $(last(sol.t))")
  end
end

function PowerGridsSystem10()
  res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.PowerGridsReal.System10", file; MSL = true)
  write("./FlatModels/PowerGridsSystem10.mo", res)
  solver = :(Rodas5())
  simulateAndPlotPowerSystemSpecial("System10"; solver = solver)
end

function runPowergrids()
  PowerGridsSystem1()
  PowerGridsSystem2()
  PowerGridsSystem3()
  PowerGridsSystem4()
  PowerGridsSystem5()
  PowerGridsSystem6()
  PowerGridsSystem7()
  PowerGridsSystem8()
  PowerGridsSystem9()
end
