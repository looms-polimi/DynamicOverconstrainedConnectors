# DynamicOverconstrainedConnectors
Test cases for the proposed extension of overconstrained connectors in Modelica

The DynamicOverconstrainedConnectors package contains three sets of test cases for the proposed extension of the overconstrained connector semantics in Modelica,
that allows to change the network topology at runtime. The proposal is described in OpenModelica ticket OpenModelica/OpenModelica#6240.

## Case 1: Electro-mechanical power grid simulation with complex phasors

The first set of test cases is found in the **PowerGridsComplex** sub-package. The models in there are minimimalistic versions of electro-mechanical power grids, with
radical simplifying assumptions that reduce the complexity of the models, while retaining the relevant aspects for this proposal. The main simplifying assumptions are:
- purely inductive transmission lines
- idealized synchonous generators that impose a voltage at their port with fixed magnitude and a phase equal to the rotor angle
- droop-based primary frequency control on the generators
- the reference frame for the phasors is rigidly connected to the rotor of the generator that is selected as the root node in the overconstrained graph

**System1** demonstrates the components of the library: two generators G1 and G2 with local loads L1 and L2 are connected by a transmission line T with a susceptance B = -5 p.u. The system starts in steady-state. When load L2 is changed at time = 1, the system undergoes a transient with some damped oscillations and settles to a new steady state at a slightly higher frequency. Since the reference frame for the whole system has the same frequency, the voltage and current phasors settle to a steady value, allowing stiff variable step-size solvers to increase the step, which is essential for long-term simulation studies (e.g. voltage stability).

**System 2** splits the transmission network into three components: two parallel lines T1a and T1b, each with B = -5 p.u., series connected with a line T2 having B = -10 p.u. The overall network has the same equivalent susceptance of the line in System 1, i.e., B = -5 p.u.; hence, the transients of the frequencies G1.omega and G2.omega, as well as the power flows G1.Pc, G1.Pe, G2.Pc, and G2.Pe are exactly the same as in System 1. On the other hand, the parallel connection of T1a and T1B creates a loop in the connection graph, so it can be used to test the ability of the tool to break that loop correctly.

**System 3** has the same structure of System2 and tests what happens when the susceptance of line T2 is additionally brought to zero at time = 10, representing the behaviour when line breakers are tripped open. The system is effectively split into two separate synchronous islands, each fed by one generator; as the loads were unbalanced at the time of the splitting, one island accelerates and one decelerates, so there are two different steady state values for the generator frequency omega at the end of the transient. With the standard fixed overconstrained connection semantics, only one reference node is still used, so the voltage and current phasors of the other synchronous island start rotating at a frequency of 0.5 Hz. This causes the corresponding parts of the system Jacobian to change all the time, preventing variable 
step-size solvers from increasing the step size once the new frequency steady-state has been reached.

**System4** is built as System3, except that now the line breaker implements the proposed extension, i.e. it dynamically removes the unbreakeable branch between its to connectors when the susceptance B is brought to zero. This splits the formerly connected connection graph into two disconnected islands, each with its root node providing the reference frequency. As a consequence, once the new frequency steady state is reached, phasors in both islands will remain constant, allowing a stiff solver to take much longer steps. Notice that the phasor representation of currents and voltages in the connectors will be different, because of the change in the frames of reference, but the physically meaningful variables, namely the generator frequencies G1.omega and G2.omega, as well as the generator powers G1.Pe, G1.Pc, G2.Pe, and G2.Pc will be exactly the same as in System3.

**System5** is similar to System3, except that now the line being tripped open is T1a. Since the other parallel line keeps the two generators connected, there is no splitting of the synchronous system into two islands and both generators eventually reach the same steady-state frequency.

**System6** is the same as System5, except that the extended line breaker is used. Since there is not synchronous system splitting when the unbreakable branch of T1a is removed, the connection graph remains fully connected after the event at time = 10 and there is no need to reconfigure the root nodes and the equality chains in the connected sub-networks.

**System7** adds to System3 one more generator G3 and one more load L3 that are connected to G2 via an additional line T3. When the line T2 is tripped open, two synchronous islands are formed, one including G1, and one including G2 and G3.

**System8** is like System7, except that the extended line breaker is used for T2. In this case, G1.port should be selected as root node for the first connected sub-graph, while G3.port should be selected as root node for the second sub-graph, since it has lower priority number than G2. Also in this case, once both islands have reached steady-state, the phasors become constant allowing longer time steps; the physical transents of the generator frequencies G1.omega, G2.omega and G3.omega, as well as the powers Pe and Ps in the three generators, should be exactly the same as in System7.

**System9** describes a possible occurence in power grid simulation, namely the formation of islands without any available power generation. A generator G1 has a connected local load L1, and feeds two other loads L2 and L3 through the series connection of L1 and L2. If L1 is tripped open, the system is effectively split into two islands: the left one, containing G1 and L1, undergoes a frequency transient because of the power imbalance caused by the disconnection of the two loads L2 and L3, but can continue operating; the right one, containg the other two loads, is no longer connected to any generator. This causes the network system to have no solution after T1 is tripped opened, because obviously loads L1 and L2 have cannot take their 0.2 p.u. active power from anywhere. Hence, the simulation aborts and cannot continue.

**System10** demonstrates how this situation can be handled by dynamic overconstrained connectors. 
The idea is that also load ports can become root nodes of the connection graph, albeit with very high priority number, hence with very low priority. 
As long as there is at least one generator in their connected island, they won't be selected, because of the low priority. 
However, if there are no generators at all in the island, then one load will be selected as root note, and will set the voltage to zero and omegaRef to zero. 
All other loads in the same island will sense that omegaRef is now zero, and absorb zero current, instead of the prescribed P and Q active and reactive power. 
As a consequence, that island will be effectively switched off (black-out mode), but the model will remain solvable and the simulation of the active islands can carry on unhindered.

## Case 2: Electro-mechanical power grid simulation with real number representation of complex phasors

The second set of test cases is found in the **PowerGridsReal** sub-package. The models found here are identical to the models in PowerGridsComplex, except they are written with Real numbers only and explicit equations for the real and imaginary part of complex phasor equations, without using Complex operator records. This can be handy in experimental Modelica environments that do not support them.

## Case 3: Incompressible fluid networks

The third set of test cases is found in the **IncompressibleFluid** sub-package. In this case, dynamic overconstrained connectors are used to handle incompressible fluid networks.

As it is well-known, any closed incompressible fluid system needs a component, usually an expansion tank, to uniquely set the pressure level in the system and to allow for thermal expansion of the fluid, even if not explicitly modelled. Without that component, the overall pressure level of the system is undetermined. As the system is closed and the fluid mass in it is constant, the flow rate into that component will be computed to be zero at all times. 

Consider now the case where the incompressible fluid system can be split into two (or more) sub-systems, by closing strategically placed valves that effectively separate the two (or more) sub-circuits. This is achieved by simply setting to zero the values of their flow coefficient. However, when this happens, only one part of the circuit is effectively connected to the expansion tank model, so that its level of pressure is well-determined. The remaining part(s) of the circuit are not effectively connected to any, so as a result their level of pressure will become undetermined, or, in other words, their hydraulic equations will become singular. In real life, these isolated sub-circuits would actually risk blowing up, since there is no place for the fluid to go when some (unmodelled) thermal expansion occurs. Hence, a proper design of the circuit is such that every *potentially isolated* sub-circuit has one expansion tank.

This situation can be easily modelled by means of dynamic overconstrained connectors. Expansion tanks fix their port pressure, but only if that port is selected as the root node of a connected graph; otherwise, they behave like plugs, setting the port flow rate to zero, leaving the job of determining the pressure level of the circuit to another expansion tank in the same connected circuit. In this case, the only part of the dynamic overconstrained connection semantics that is needed is the selection of root nodes; strictly speaking, there is no need of overconstrained connector variables. However, the semantics of overconstrained connections in Modelica 3.5 is defined on the overconstrained connector variables, rather than on the connectors themselves, requiring to define them somehow. One option is to use them to carry an identification number for the conneted sub-circuit, which is then set by the active expansion tank.

The flow equations for pumps and valves are simple linear relationships, which ensure that the implicit algbraic equations systems that need to be solved to simulate the test cases are linear and have no issues with initialization. Using more accurate quadratic relationships would make the system more realistic, but would not affect the need of dynamic overconstrained connectors at all, so it is preferrable to keep things as simple as possible in that respect.

**System1** is a simple demo of the library, which only contains one loop and does not require to handle the connection graphs dynamically. In this system, a pump drives the flow through a circuit made by a valve and two connected pipes. The pressure level is determined by an expansion tank connected at the pump inlet. At time t = 1 the rotational speed of the pump is reduced, so both the flow rates and the pressure differences across the valve and pipes are reduced.

**System2** contains two loops, formed by a pump and three series-connected pipes, which are joined togeher by two valves, one connecting the inlets of the two intermediate pipes, one connecting the outlets of the two intermediate pipes. Each loop has an expansion tank at the pump inlet. The system starts in a symmetric operating point; as it is fully connected, only on expansion tank is activated as a root node, and the entire circuit has id = 1. At time t = 1 the rotational speed of the left pump is reduced, bringing the system to a new asymmetric operating point. At time = 2, one of the two valves is closed; this prevents the flow between the two sub-circuits, changing the flow pattern again, but still leaves the circuit fully connected. The simulation stops at time = 3, before the other valve is closed, hence it can be handled by the standard static overconstrained connector semantics.

**System3** is the same as System2, except that the simulation continues further; at time = 4, also the second valve is closed. From that point onwards, the algebraic system of equations determining the circuit pressures and flows becomes singular.

**System4** is the same as System3, except that the valve models are redeclared to use the proposed dynamic connection graphs. 
When both valves are closed at time = 4, the formely fully connected circuit is split into two separate circuits, one for loop A, with id = 1, and one for loop B, with id = 2, thanks to the deactivation of the unbreakable branches of the two valves and to the activation of tankB's port as a root node for the second circuit. 
Note that the algebraic system of equations determining the flows and pressures in the system still encompasses the entire circuit, albeit it effectively describes two independent circuits, thanks to the specific choice of some strategic coefficient values. 
However, thanks to the overconstrained connector formulation, this system no longer becomes singular, as it happened in the previous test case.
