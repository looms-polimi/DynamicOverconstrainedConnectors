# DynamicOverconstrainedConnectors
Test cases for the proposed extension of overconstrained connectors in Modelica

The DynamicOverconstrainedConnectors package contains three sets of test cases for the proposed extension of the overconstrained connector semantics in Modelica,
that allows to change the network topology at runtime. The proposal is described in OpenModelica ticket OpenModelica/OpenModelica#6240.

The first set of test cases is found in the **PowerGridsComplex** sub-package. The models in there are minimimalistic versions of electro-mechanical power grids, with
radical simplifying assumptions that reduce the complexity of the models, while retaining the relevant aspects for this proposal. The main simplifying assumptions are:
- purely inductive transmission lines
- idealized synchonous generators that impose a voltage at their port with fixed magnitude and a phase equal to the rotor angle
- droop-based primary frequency control on the generators
- the reference frame for the phasors is rigidly connected to the rotor of the generator that is selected as the root node in the overconstrained graph

**System1** demonstrates the components of the library: two generators G1 and G2 with local loads L1 and L2 are connected by a transmission line T with a susceptance B = -5 p.u. The system starts in steady-state. When load L2 is changed at time = 1, the system undergoes a transient with some damped oscillations and settles to a new steady state at a slightly higher frequency. Since the reference frame for the whole system has the same frequency, the voltage and current phasors settle to a steady value, allowing stiff variable step-size solvers to increase the step, which is essential for long-term simulation studies (e.g. voltage stability).

**System 2** splits the transmission network into three components: two parallel lines T1a and T1b, each with B = -5 p.u., series connected with a line T2 having B = -10 p.u. The overall network has the same equivalent susceptance of the line in System 1, i.e., B = -5 p.u.; hence, the transients of the frequencies G1.omega and G2.omega, as well as the power flows G1.Pc, G1.Pe, G2.Pc, and G2.Pe are exactly the same as in System 1.

**System 3** has the same structure of System2 and tests what happens when the susceptance of line T2 is additionally brought to zero at time = 10, representing the behaviour when line breakers are tripped open. The system is effectively split into two separate synchronous islands, each fed by one generator; as the loads were unbalanced at the time of the splitting, one island accelerates and one decelerates, so there are two different steady state values for the generator frequency omega at the end of the transient. With the standard fixed overconstrained connection semantics, only one reference node is still used, so the voltage and current phasors of the other synchronous island start rotating at a frequency of 0.5 Hz. This causes the corresponding parts of the system Jacobian to change all the time, preventing variable 
step-size solvers from increasing the step size once the new frequency steady-state has been reached.

**System4** is built as System3, except that now the line breaker implements the proposed extension, i.e. it dynamically removes the unbreakeable branch between its to connectors when the susceptance B is brought to zero. This splits the formerly connected connection graph into two disconnected islands, each with its root node providing the reference frequency. As a consequence, once the new frequency steady state is reached, phasors in both islands will remain constant, allowing a stiff solver to take much longer steps. Notice that the phasor representation of currents and voltages in the connectors will be different, because of the change in the frames of reference, but the physically meaningful variables, namely the generator frequencies G1.omega and G2.omega, as well as the generator powers G1.Pe, G1.Pc, G2.Pe, and G2.Pc will be exactly the same as in System3.

**System5** is similar to System3, except that now the line being tripped open is T1a. Since the other parallel line keeps the two generators connected, there is no splittng of the synchronous system into two islands and both generators eventually reach the same steady-state frequency.

**System6** is the same as System5, except that the extended line breaker is used. Since there is not synchronous system splitting when the unbreakable branch of T1a is removed, the connection graph remains fully connected after the event at time = 10 and there is no need to reconfigure the root nodes and the equality chains in the connected sub-networks.

**System7** adds to System3 one more small generator that is connected to G2 via an additional line T3. When the line T2 is tripped open, two synchronous islands are formed, one including G1, and one including G2 and G3.

**System8** is like System7, except that the extended line breaker is used for T2. In this case, G1.port should be selected as root node for the first connected sub-graph, while G2.port should be selected as root node for the second sub-graph, since it has lower priority than G3. Also in this case, once both islands have reached steady-state, the phasors become constant allowing longer time steps; the physical transents of the generator frequencies G1.omega, G2.omega and G3.omega, as well as the powers Pe and Ps in the three generators, should be exactly the same as in System7.
