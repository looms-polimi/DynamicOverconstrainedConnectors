global file = "DynamicOverconstrainedConnectors.mo"
println("Running fluid Networks")

function simulateAndPlotFluid(modelName; solver, tspan = (0.0, 2.0))
  println("Testing $(modelName)")
  local startTime = first(tspan)
  local stopTime = last(tspan)
  @time sol = OM.simulate("DynamicOverconstrainedConnectors.IncompressibleFluid.$(modelName)",
                          file;
                          startTime = startTime,
                          stopTime = stopTime,
                          MSL = true,
                          solver = solver)
  @info "Simulation done"
  if modelName == "System1"
    p = plot(sol; vars = [OM.OMBackend.pump_dp])
    Plots.pdf(p, "./Plots/Fluid/$(modelName)/$(modelName)pump_dp")
    p = plot(sol; vars = [OM.OMBackend.pump_outlet_w])
    Plots.pdf(p, "./Plots/Fluid/$(modelName)/$(modelName)pump_outlet_w")
    p = plot(sol; vars = [OM.OMBackend.pump_n])
    #= The normalized rotational speed =#
    Plots.pdf(p, "./Plots/Fluid/$(modelName)/$(modelName)pump__n")
  elseif modelName == "System2" || modelName == "System3" || modelName == "System4"
    p = plot(sol; vars = [OM.OMBackend.pumpA_dp])
    Plots.pdf(p, "./Plots/Fluid/$(modelName)/$(modelName)pumpA_dp")
    p = plot(sol; vars = [OM.OMBackend.pumpA_outlet_w])
    Plots.pdf(p, "./Plots/Fluid/$(modelName)/$(modelName)pumpA_outlet_w")
    p = plot(sol; vars = [OM.OMBackend.pumpA_n])
    Plots.pdf(p, "./Plots/Fluid/$(modelName)/$(modelName)pumpA__n")
  else
    @info "Nothing to plot"
  end
end

function simulateAndPlotFluidSpecial(modelName; solver, tspan = (0.0, 2.0))
  println("Testing $(modelName)")
  local startTime = first(tspan)
  local stopTime = last(tspan)
  @time sols = OM.simulate("DynamicOverconstrainedConnectors.IncompressibleFluid.$(modelName)",
                           file;
                           startTime = startTime,
                           stopTime = stopTime,
                           MSL = true,
                           solver = solver)
  @info "Length:" length(sols)
  @info "Simulation done"
  p = plotPartial(sols, [OM.OMBackend.pumpA_dp], tspan, "./Plots/Fluid/$(modelName)/$(modelName)pumpA_dp")

  p = plotPartial(sols, [OM.OMBackend.pumpA_outlet_w], tspan, "./Plots/Fluid/$(modelName)/$(modelName)pumpA_outlet_w")
  p = plotPartial(sols, [OM.OMBackend.pumpA_outlet_w], tspan, "./Plots/Fluid/$(modelName)/$(modelName)pumpA__n")

  p = plotPartial(sols, [OM.OMBackend.valveAB_Kv], tspan, "./Plots/Fluid/$(modelName)/$(modelName)valveAB_Kv")
  p = plotPartial(sols, [OM.OMBackend.valveBA_Kv], tspan, "./Plots/Fluid/$(modelName)/$(modelName)valveBA_Kv")

  p = plotPartial(sols, [OM.OMBackend.valveAB_closed], tspan, "./Plots/Fluid/$(modelName)/$(modelName)valveAB_closed")
  p = plotPartial(sols, [OM.OMBackend.valveBA_closed], tspan, "./Plots/Fluid/$(modelName)/$(modelName)valveBA_closed")
end

function FluidNetworksSystem1()
  res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.IncompressibleFluid.System1", file; MSL = true)
  write("./FlatModels/FluidNetworksSystem1.mo", res)
  solver = :(Rodas5())
  simulateAndPlotFluid("System1"; solver = solver)
end

function FluidNetworksSystem2()
  res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.IncompressibleFluid.System2", file; MSL = true)
  write("./FlatModels/FluidNetworksSystem2.mo", res)
  solver = :(Rodas5())
  simulateAndPlotFluid("System2"; solver = solver, tspan = (0.0, 3.0))
end

function FluidNetworksSystem3()
  res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.IncompressibleFluid.System3", file; MSL = true)
  write("./FlatModels/FluidNetworksSystem3.mo", res)
  solver = :(Rodas5())
  simulateAndPlotFluid("System3"; solver = solver, tspan = (0.0, 5.0))
end

function FluidNetworksSystem4()
  res = OM.generateFlatModelica("DynamicOverconstrainedConnectors.IncompressibleFluid.System4", file; MSL = true)
  write("./FlatModels/FluidNetworksSystem4.mo", res)
  solver = :(Rodas5())
  simulateAndPlotFluidSpecial("System4"; solver = solver, tspan = (0.0, 5.0))
end

function runFluidSystems()
  FluidNetworksSystem1()
  FluidNetworksSystem2()
  FluidNetworksSystem3()
  #FluidNetworksSystem4()
end
