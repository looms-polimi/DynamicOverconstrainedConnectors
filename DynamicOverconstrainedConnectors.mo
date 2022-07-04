package DynamicOverconstrainedConnectors
  package PowerGridsComplex "Power grid models using complex records"
    import SI = Modelica.SIunits;
    import CM = Modelica.ComplexMath;
    constant SI.Frequency f_n = 50 "Nominal grid frequency";
    constant SI.PerUnit pi = Modelica.Constants.pi;
    final constant SI.AngularVelocity omega_n = 2*pi*f_n;
  
    type ReferenceAngularSpeed "Overconstrained type for per-unit reference angular speed"
      extends SI.PerUnit;
      function equalityConstraint
        input ReferenceAngularSpeed omega1;
        input ReferenceAngularSpeed omega2;
        output SI.PerUnit residue[0] "No constraints";
      end equalityConstraint;
      function reconnectable
        input ReferenceFrequency f1;
        input ReferenceFrequency f2;
        output Boolean canReconnect;
      protected
        constant Real eps = 1e-6;
      algorithm
        canReconnect := abs(f1 - f2) < eps;
      end reconnectable;
    end ReferenceAngularSpeed;
    
    connector ACPort "Port for per unit 3-phase AC systems"
      SI.ComplexPerUnit v "Per unit voltage phasor referred to phase";
      flow SI.ComplexPerUnit i "Per unit current phasor referred to phase";
      ReferenceAngularSpeed omegaRef "Reference per-unit angular speed";
    end ACPort;
  
    model Load "AC load model"
      ACPort port;
      SI.PerUnit P = 0 "Active per unit power";
      Real Q = 0 "Reactive per unit power";
    equation
      port.v*CM.conj(port.i) = Complex(P,Q);
    end Load;

    partial model TransmissionLineBase "Purely inductive transmission line - base model"
      parameter SI.PerUnit B = -5.0 "Line series per unit susceptance";
      discrete SI.PerUnit B_act "Actual value of per unit susceptance including breaker status";
      Boolean closed "State of line breaker";
      Boolean open = false "Command to open the line breaker";
      Boolean close = false "Command to close the line breaker";
      ACPort port_a;
      ACPort port_b;
    initial equation
      closed = true;
      B_act = B;
    equation
      port_a.i + port_b.i = Complex(0);
      port_a.i = Complex(0,B_act)*(port_a.v - port_b.v);
      when open then
        closed = false;
        B_act = 0;
      elsewhen close then
        closed = true;
        B_act = B;
      end when;
    end TransmissionLineBase;
    
    model TransmissionLine "Purely inductive transmission line model - static overconstrained connector version"
      extends TransmissionLineBase;
    equation
      port_a.omegaRef = port_b.omegaRef;
      Connections.branch(port_a.omegaRef, port_b.omegaRef);
    end TransmissionLine;

    model TransmissionLineVariableBranch "Purely inductive transmission line model with time-varying connection branch"
      extends TransmissionLineBase;
    equation
  // This if-equation is not valid according to Modelica 3.5
  // Section 8.3.4 but it would according to this extension proposal
      if closed then
        port_a.omegaRef = port_b.omegaRef;
        Connections.branch(port_a.omegaRef, port_b.omegaRef);
      end if;
    end TransmissionLineVariableBranch;
    
     model Generator "Idealized synchronous generator with ideal voltage control and basic primary frequency control"
       parameter SI.PerUnit V = 1 "Fixed rotor per unit voltage magnitude";
       parameter SI.Time Ta = 10 "Acceleration time of turbogenerator Ta = J*omega^2/Pnom";
       parameter SI.PerUnit droop = 0.05 "Droop coefficient of primary frequency control";
       parameter Integer p = 0 "Potential root node priority";
       ACPort port;
       SI.PerUnit Ps = 1 "Active power output set point";
       SI.PerUnit Pc "Primary frequency control power in per unit";
       SI.PerUnit Pe "Electrical power output";
       SI.Angle theta(start = 0, fixed = true) "Machine angle relative to port.theta";
       SI.PerUnit omega(start = 1, fixed = true) "Per unit angular speed";
    equation
       der(theta) = (omega - port.omegaRef)*omega_n;
       Ta*omega*der(omega) = Ps + Pc - Pe;
       port.v = CM.fromPolar(V, theta);
       Pe = -CM.real(port.v*CM.conj(port.i));
       Pc = -(omega-1)/droop;
       // Connections.potentialRoot(port.omegaRef,p);
       Connections.potentialRoot(port.omegaRef);
       if Connections.isRoot(port.omegaRef) then
         port.omegaRef = omega;
       end if;
    end Generator;
    
    model System1 "Two generators, one line, fixed branches"
      Generator G1;
      Generator G2;
      Load L1(P = 1);
      Load L2(P = if time < 1 then 1 else 0.8);
      TransmissionLine T;
    equation
      connect(G1.port, L1.port);
      connect(G2.port, L2.port);
      connect(G1.port, T.port_a);
      connect(G2.port, T.port_b);
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System1;
    
    model System2 "Two generators, two-parallel and one series lines, fixed branches"
      Generator G1;
      Generator G2;
      Load L1(P = 1);
      Load L2(P = if time < 1 then 1 else 0.8);
      TransmissionLine T1a(B = -5.0);
      TransmissionLine T1b(B = -5.0);
      replaceable TransmissionLine T2(B = -10.0)
        constrainedby TransmissionLineBase;
    equation
      connect(G1.port, L1.port);
      connect(G2.port, L2.port);
      connect(G1.port, T1a.port_a);
      connect(G1.port, T1b.port_a);
      connect(T1a.port_b, T2.port_a);
      connect(T1b.port_b, T2.port_a);
      connect(G2.port, T2.port_b);
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System2;
    
    model System3 "Two generators, two-parallel and one series line with breaker, fixed branches"
      extends System2(T2(open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System3;
    
    model System4 "Two generators, two-parallel and one series line with breaker, dynamic branches"
      extends System3(
        redeclare TransmissionLineVariableBranch T2(B = -10.0, open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System4;
    
    model System5 "Two generators, two parallel (one with breaker) and one series, fixed branches"
      extends System2(T1b(open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System5;
    
    model System6 "Two generators, two-parallel and one series line with breaker, dynamic branches"
      extends System5(
        redeclare TransmissionLineVariableBranch T1b(B = -5.0, open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System6;
    
    model System7 "Three generators, two-parallel and two series lines with breaker, fixed branches"
      extends System3(G2(p = 1));
      Generator G3(Ta = 15);
      Load L3(P = 1);
      TransmissionLine T3;
    equation
      connect(G3.port, L3.port);
      connect(G3.port, T3.port_b);
      connect(G2.port, T3.port_a);
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System7;
    
    model System8 "Three generators, two-parallel and two series lines with breaker dynamic branches"
      extends System7(
        redeclare TransmissionLineVariableBranch T2(B = -10.0, open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System8;
  end PowerGridsComplex;

  package PowerGridsReal "Power grid models using real variables"
    import SI = Modelica.Units.SI;
    import CM = Modelica.ComplexMath;
    constant SI.Frequency f_n = 50 "Nominal grid frequency";
    constant SI.PerUnit pi = Modelica.Constants.pi;
    final constant SI.AngularVelocity omega_n = 2*pi*f_n;
  
    type ReferenceAngularSpeed "Overconstrained type for per-unit reference angular speed"
      extends SI.PerUnit;
      function equalityConstraint
        input ReferenceAngularSpeed omega1;
        input ReferenceAngularSpeed omega2;
        output SI.PerUnit residue[0] "No constraints";
      end equalityConstraint;
      function reconnectable
        input ReferenceFrequency f1;
        input ReferenceFrequency f2;
        output Boolean canReconnect;
      protected
        constant Real eps = 1e-6;
      algorithm
        canReconnect := abs(f1 - f2) < eps;
      end reconnectable;
    end ReferenceAngularSpeed;
    
    connector ACPort "Port for per unit 3-phase AC systems"
      SI.PerUnit v_re "Per unit voltage phasor referred to phase, real part";
      SI.PerUnit v_im "Per unit voltage phasor referred to phase, imaginary part";
      flow SI.PerUnit i_re "Per unit current phasor referred to phase, real part";
      flow SI.PerUnit i_im "Per unit current phasor referred to phase, imaginary part";
      ReferenceAngularSpeed omegaRef "Reference per-unit angular speed";
    end ACPort;
  
    model Load "AC load model"
      ACPort port;
      SI.PerUnit P = 0 "Active per unit power";
      Real Q = 0 "Reactive per unit power";
    equation
       port.v_re*port.i_re + port.v_im*port.i_im = P;
      -port.v_re*port.i_im + port.v_im*port.i_re = Q;
    end Load;

    partial model TransmissionLineBase "Purely inductive transmission line - base model"
      parameter SI.PerUnit B = -5.0 "Line series per unit susceptance";
      discrete SI.PerUnit B_act "Actual value of per unit susceptance including breaker status";
      Boolean closed "State of line breaker";
      Boolean open = false "Command to open the line breaker";
      Boolean close = false "Command to close the line breaker";
      ACPort port_a;
      ACPort port_b;
    initial equation
      closed = true;
      B_act = B;
    equation
      port_a.i_re + port_b.i_re = 0;
      port_a.i_im + port_b.i_im = 0;
      port_a.i_re = -B_act*(port_a.v_im - port_b.v_im);
      port_a.i_im =  B_act*(port_a.v_re - port_b.v_re);
      when open then
        closed = false;
        B_act = 0;
      elsewhen close then
        closed = true;
        B_act = B;
      end when;
    end TransmissionLineBase;
    
    model TransmissionLine "Purely inductive transmission line model"
      extends TransmissionLineBase;
    equation
      port_a.omegaRef = port_b.omegaRef;
      Connections.branch(port_a.omegaRef, port_b.omegaRef);
    end TransmissionLine;
    
    model TransmissionLineVariableBranch "Purely inductive transmission line model with time-varying connection branch"
    extends TransmissionLineBase;
    equation
      // This if-equation is not valid according to Modelica 3.5
    // Section 8.3.4 but it would according to this extension proposal
      if closed then
        port_a.omegaRef = port_b.omegaRef;
        Connections.branch(port_a.omegaRef, port_b.omegaRef);
      end if;
    end TransmissionLineVariableBranch;
    
    model Generator "Idealized synchronous generator with ideal voltage control and basic primary frequency control"
      parameter SI.PerUnit V = 1 "Fixed rotor per unit voltage magnitude";
      parameter SI.Time Ta = 10 "Acceleration time of turbogenerator Ta = J*omega^2/Pnom";
      parameter SI.PerUnit droop = 0.05 "Droop coefficient of primary frequency control";
      parameter Integer p = 0 "Potential root node priority";
      ACPort port;
      SI.PerUnit Ps = 1 "Active power output set point";
      SI.PerUnit Pc "Primary frequency control power in per unit";
      SI.PerUnit Pe "Electrical power output";
      SI.Angle theta(start = 0, fixed = true) "Machine angle relative to port.theta";
      SI.PerUnit omega(start = 1, fixed = true) "Per unit angular speed";
    equation
      der(theta) = (omega - port.omegaRef)*omega_n;
      Ta*omega*der(omega) = Ps + Pc - Pe;
      port.v_re = V*cos(theta);
      port.v_im = V*sin(theta);
      Pe = -(port.v_re*port.i_re + port.v_im*port.i_im);
      Pc = -(omega-1)/droop;
      // Connections.potentialRoot(port.omegaRef,p);
      Connections.potentialRoot(port.omegaRef);
      if Connections.isRoot(port.omegaRef) then
        port.omegaRef = omega;
      end if;
    end Generator;
    
    model System1 "Two generators, one line, fixed branches"
      Generator G1;
      Generator G2;
      Load L1(P = 1);
      Load L2(P = if time < 1 then 1 else 0.8);
      TransmissionLine T;
    equation
      connect(G1.port, L1.port);
      connect(G2.port, L2.port);
      connect(G1.port, T.port_a);
      connect(G2.port, T.port_b);
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System1;
    
    model System2 "Two generators, two-parallel and one series lines, fixed branches"
      Generator G1;
      Generator G2;
      Load L1(P = 1);
      Load L2(P = if time < 1 then 1 else 0.8);
      TransmissionLine T1a(B = -5.0);
      TransmissionLine T1b(B = -5.0);
      replaceable TransmissionLine T2(B = -10.0)
        constrainedby TransmissionLineBase;
    equation
      connect(G1.port, L1.port);
      connect(G2.port, L2.port);
      connect(G1.port, T1a.port_a);
      connect(G1.port, T1b.port_a);
      connect(T1a.port_b, T2.port_a);
      connect(T1b.port_b, T2.port_a);
      connect(G2.port, T2.port_b);
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System2;
    
    model System3 "Two generators, two-parallel and one series line with breaker, fixed branches"
      extends System2(T2(open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System3;
    
    model System4 "Two generators, two-parallel and one series line with breaker, dynamic branches"
      extends System3(
        redeclare TransmissionLineVariableBranch T2(B = -10.0, open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System4;
    
    model System5 "Two generators, two parallel (one with breaker) and one series, fixed branches"
      extends System2(T1b(open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System5;
    
    model System6 "Two generators, two-parallel and one series line with breaker, dynamic branches"
      extends System5(
        redeclare TransmissionLineVariableBranch T1b(B = -5.0, open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System6;
    
    model System7 "Three generators, two-parallel and two series lines with breaker, fixed branches"
      extends System3(G2(p = 1));
      Generator G3(Ta = 15);
      Load L3(P = 1);
      TransmissionLine T3;
    equation
      connect(G3.port, L3.port);
      connect(G3.port, T3.port_b);
      connect(G2.port, T3.port_a);
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System7;
    
    model System8 "Three generators, two-parallel and two series lines with breaker dynamic branches"
      extends System7(
        redeclare TransmissionLineVariableBranch T2(B = -10.0, open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System8;
  end PowerGridsReal;

annotation(uses(Modelica(version="3.2.3")));
end DynamicOverconstrainedConnectors;
