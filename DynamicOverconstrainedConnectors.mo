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
      SI.ComplexPerUnit v(re(start = 1)) "Per unit voltage phasor referred to phase";
      flow SI.ComplexPerUnit i "Per unit current phasor referred to phase";
      ReferenceAngularSpeed omegaRef "Reference per-unit angular speed";
    annotation(
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 1, grid = {2, 2}), graphics = {Rectangle(origin = {92, 3}, fillColor = {85, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-192, 97}, {8, -103}})}));
    end ACPort;
    
    model LoadBase "Base class for AC load models"
      ACPort port annotation(
        Placement(visible = true, transformation(origin = {-1.42109e-14, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-1.42109e-14, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SI.PerUnit P = 0 "Active per unit power";
      Real Q = 0 "Reactive per unit power";
    annotation(
        Icon(graphics = {Line(origin = {0, -20}, points = {{0, 20}, {0, -20}, {0, -20}}), Polygon(origin = {0, -70}, fillPattern = FillPattern.Solid, points = {{-40, 30}, {40, 30}, {0, -30}, {-40, 30}}), Text(origin = {-1, -120}, lineColor = {0, 0, 255}, extent = {{-81, 20}, {81, -20}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
    end LoadBase;
  
    model Load "AC load model"
      extends LoadBase;
    equation
      port.v*CM.conj(port.i) = Complex(P,Q);
    end Load;

    model LoadVariableRoot "AC load model with provision for unsupplied sub-networks"
      extends LoadBase;
    equation
      if port.omegaRef > 0 then 
        // the load is connected to a synchronous island with a generator as a root node
        port.v*CM.conj(port.i) = Complex(P,Q);
      elseif Connections.isRoot(port.omegaRef) then
        // the load is the root node of an island without any generator root node
        port.v = Complex(0);
      else
        // the load is connected to an island that doesn't have any generator, but is not the root node
        port.i = Complex(0);
      end if;
      // Very high priority number, gets selected only if there are no generators in the sub-graph
      Connections.potentialRoot(port.omegaRef, 10000);
      // Set omegaRef = 0 if the load is selected as root node, hence there are no generators in the sub-graph
      if Connections.isRoot(port.omegaRef) then
        port.omegaRef = 0;
      end if;
    end LoadVariableRoot;

    partial model TransmissionLineBase "Purely inductive transmission line - base model"
      parameter SI.PerUnit B = -5.0 "Line series per unit susceptance";
      discrete SI.PerUnit B_act "Actual value of per unit susceptance including breaker status";
      Boolean closed "State of line breaker";
      Boolean open = false "Command to open the line breaker";
      Boolean close = false "Command to close the line breaker";
      ACPort port_a annotation(
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ACPort port_b annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
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
    annotation(
        Icon(graphics = {Rectangle(origin = {-1, -1}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-59, 11}, {61, -11}}), Line(origin = {-80, 0}, points = {{-20, 0}, {20, 0}}), Line(origin = {80, 0}, points = {{-20, 0}, {20, 0}}), Text(origin = {0, 40}, lineColor = {0, 0, 255}, extent = {{-80, 20}, {80, -20}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
    end TransmissionLineBase;
    
    model TransmissionLine "Purely inductive transmission line model - static overconstrained connector version"
      extends TransmissionLineBase;
    equation
      port_a.omegaRef = port_b.omegaRef;
      Connections.branch(port_a.omegaRef, port_b.omegaRef);
    end TransmissionLine;

    model TransmissionLineDynamicBranch "Purely inductive transmission line model with time-varying connection branch"
      extends TransmissionLineBase;
    equation
// This if-equation is not valid according to Modelica 3.5
// Section 8.3.4 but it would according to this extension proposal
      if closed then
        port_a.omegaRef = port_b.omegaRef;
        Connections.branch(port_a.omegaRef, port_b.omegaRef);
      end if;
    end TransmissionLineDynamicBranch;
    
     model Generator "Idealized synchronous generator with ideal voltage control and basic primary frequency control"
       parameter SI.PerUnit V = 1 "Fixed rotor per unit voltage magnitude";
       parameter SI.Time Ta = 10 "Acceleration time of turbogenerator Ta = J*omega^2/Pnom";
       parameter SI.PerUnit droop = 0.05 "Droop coefficient of primary frequency control";
       parameter Integer p = 0 "Potential root node priority";
       ACPort port annotation(
         Placement(visible = true, transformation(origin = {-1.42109e-14, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-1.42109e-14, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
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
  Connections.potentialRoot(port.omegaRef,p);
       if Connections.isRoot(port.omegaRef) then
         port.omegaRef = omega;
       end if;
    annotation(
         Icon(graphics = {Ellipse(origin = {0, -1}, extent = {{-60, 61}, {60, -59}}), Text(origin = {0, 80}, lineColor = {0, 0, 255}, extent = {{-80, 20}, {80, -20}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
     end Generator;
    
    model System1 "Two generators, one line, fixed branches"
      Generator G1 annotation(
        Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Generator G2 annotation(
        Placement(visible = true, transformation(origin = {20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Load L1(P = 1) annotation(
        Placement(visible = true, transformation(origin = {-20, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Load L2(P = if time < 1 then 1 else 0.8) annotation(
        Placement(visible = true, transformation(origin = {20, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      TransmissionLine T annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(G1.port, L1.port) annotation(
        Line(points = {{-20, 0}, {-20, -10}}));
      connect(G2.port, L2.port) annotation(
        Line(points = {{20, 0}, {20, -10}}));
      connect(G1.port, T.port_a) annotation(
        Line(points = {{-20, 0}, {-10, 0}}));
      connect(G2.port, T.port_b) annotation(
        Line(points = {{20, 0}, {10, 0}}));
    annotation(experiment(StopTime = 50, Interval = 0.02),
  Diagram(coordinateSystem(extent = {{-40, -25}, {40, 15}})),
  Icon(coordinateSystem(extent = {{-40, -25}, {40, 15}})));
    end System1;
    
    model System2 "Two generators, two-parallel and one series lines, fixed branches"
      Generator G1 annotation(
        Placement(visible = true, transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Generator G2 annotation(
        Placement(visible = true, transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Load L1(P = 1)  annotation(
        Placement(visible = true, transformation(origin = {-40, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Load L2(P = if time < 1 then 1 else 0.8) annotation(
        Placement(visible = true, transformation(origin = {40, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      TransmissionLine T1a(B = -5.0) annotation(
        Placement(visible = true, transformation(origin = {-14, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      TransmissionLine T1b(B = -5.0) annotation(
        Placement(visible = true, transformation(origin = {-14, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      replaceable TransmissionLine T2(B = -10.0) annotation(
        Placement(visible = true, transformation(origin = {16, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)))
        constrainedby TransmissionLineBase ;
    equation
      connect(G1.port, L1.port) annotation(
        Line(points = {{-40, 0}, {-40, -10}}));
      connect(G2.port, L2.port) annotation(
        Line(points = {{40, 0}, {40, -10}}));
      connect(G1.port, T1a.port_a) annotation(
        Line(points = {{-40, 0}, {-30, 0}, {-30, 6}, {-24, 6}}));
      connect(G1.port, T1b.port_a) annotation(
        Line(points = {{-40, 0}, {-30, 0}, {-30, -6}, {-24, -6}}));
      connect(T1a.port_b, T2.port_a) annotation(
        Line(points = {{-4, 6}, {2, 6}, {2, 0}, {6, 0}}));
      connect(T1b.port_b, T2.port_a) annotation(
        Line(points = {{-4, -6}, {2, -6}, {2, 0}, {6, 0}}));
      connect(G2.port, T2.port_b) annotation(
        Line(points = {{26, 0}, {40, 0}}));
    annotation(experiment(StopTime = 50, Interval = 0.02),
  Diagram(coordinateSystem(extent = {{-50, -25}, {50, 15}})),
  Icon(coordinateSystem(extent = {{-50, -25}, {50, 15}})));
    end System2;
    
    model System3 "Two generators, two-parallel and one series line with breaker, fixed branches"
      extends System2(T2(open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System3;
    
    model System4 "Two generators, two-parallel and one series line with breaker, dynamic branches"
      extends System3(
        redeclare TransmissionLineDynamicBranch T2(B = -10.0, open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System4;
    
    model System5 "Two generators, two parallel (one with breaker) and one series, fixed branches"
      extends System2(T1b(open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System5;
    
    model System6 "Two generators, two-parallel and one series line with breaker, dynamic branches"
      extends System5(
        redeclare TransmissionLineDynamicBranch T1b(B = -5.0, open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System6;
    
    model System7 "Three generators, two-parallel and two series lines with breaker, fixed branches"
      extends System3(G2(p = 1));
      Generator G3(Ta = 15) annotation(
        Placement(visible = true, transformation(origin = {84, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Load L3(P = 1) annotation(
        Placement(visible = true, transformation(origin = {84, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      TransmissionLine T3  annotation(
        Placement(visible = true, transformation(origin = {62, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(G2.port, T3.port_a) annotation(
        Line(points = {{40, 0}, {52, 0}}));
      connect(T3.port_b, G3.port) annotation(
        Line(points = {{72, 0}, {84, 0}}));
      connect(G3.port, L3.port) annotation(
        Line(points = {{84, 0}, {84, -10}}));
    annotation(experiment(StopTime = 50, Interval = 0.02),
  Diagram(coordinateSystem(extent = {{-60, -25}, {100, 15}})),
  Icon(coordinateSystem(extent = {{-60, -25}, {100, 15}})));
    end System7;
    
    model System8 "Three generators, two-parallel and two series lines with breaker dynamic branches"
      extends System7(
        redeclare TransmissionLineDynamicBranch T2(B = -10.0, open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System8;
    
    model System9 "System with loads that can get disconnected"
      Generator G1  annotation(
        Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Load L1(P = 0.8) annotation(
        Placement(visible = true, transformation(origin = {-20, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      replaceable Load L2(P = 0.1) annotation(
        Placement(visible = true, transformation(origin = {10, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)))
        constrainedby LoadBase;
      replaceable Load L3(P = 0.1) annotation(
        Placement(visible = true, transformation(origin = {36, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)))
        constrainedby LoadBase;
      replaceable TransmissionLine T1(open = if time < 10 then false else true)
        constrainedby TransmissionLineBase  annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      TransmissionLine T2 annotation(
        Placement(visible = true, transformation(origin = {26, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(G1.port, L1.port) annotation(
        Line(points = {{-20, 0}, {-20, -10}}));
      connect(G1.port, T1.port_a) annotation(
        Line(points = {{-20, 0}, {-10, 0}}));
      connect(T1.port_b, L2.port) annotation(
        Line(points = {{10, 0}, {10, -10}}));
      connect(T1.port_b, T2.port_a) annotation(
        Line(points = {{10, 0}, {16, 0}}));
      connect(T2.port_b, L3.port) annotation(
        Line(points = {{36, 0}, {36, -10}}));
    annotation(experiment(StopTime = 50, Interval = 0.02),
    Diagram(coordinateSystem(extent = {{-40, -25}, {60, 15}})),
    Icon(coordinateSystem(extent = {{-40, -25}, {60, 15}})));
    end System9;
    
    model System10 "System with loads that can get disconnected, using dynamic overconstrained connectors"
      extends System9(
        redeclare LoadVariableRoot L2(P = 0.1),
        redeclare LoadVariableRoot L3(P = 0.1),
        redeclare TransmissionLineDynamicBranch T1(open = if time < 0 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System10;
  end PowerGridsComplex;

  package PowerGridsReal "Power grid models using real variables"
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
      SI.PerUnit v_re(start = 1) "Per unit voltage phasor referred to phase, real part";
      SI.PerUnit v_im "Per unit voltage phasor referred to phase, imaginary part";
      flow SI.PerUnit i_re "Per unit current phasor referred to phase, real part";
      flow SI.PerUnit i_im "Per unit current phasor referred to phase, imaginary part";
      ReferenceAngularSpeed omegaRef "Reference per-unit angular speed";
    end ACPort;
  
  model LoadBase "Base class for AC load models"
      ACPort port;
      SI.PerUnit P = 0 "Active per unit power";
      SI.PerUnit Q = 0 "Reactive per unit power";
    end LoadBase;
  
    model Load "AC load model"
      extends LoadBase;
    equation
       port.v_re*port.i_re + port.v_im*port.i_im = P;
      -port.v_re*port.i_im + port.v_im*port.i_re = Q;
    end Load;
    
    model LoadVariableRoot "AC load model"
      extends LoadBase;
    equation
      if port.omegaRef > 0 then 
        // the load is connected to a synchronous island with a generator as a root node
         port.v_re*port.i_re + port.v_im*port.i_im = P;
        -port.v_re*port.i_im + port.v_im*port.i_re = Q;
      elseif Connections.isRoot(port.omegaRef) then
        // the load is the root node of an island without any generator root node
        port.v_re = 0;
        port.v_im = 0;
      else
        // the load is connected to an island that doesn't have any generator, but is not the root node
        port.i_re = 0;
        port.i_im = 0;
      end if;
      // Very high priority number, gets selected only if there are no generators in the sub-graph
      Connections.potentialRoot(port.omegaRef, 10000);
      // Set omegaRef = 0 if the load is selected as root node, hence there are no generators in the sub-graph
      if Connections.isRoot(port.omegaRef) then
        port.omegaRef = 0;
      end if;
    end LoadVariableRoot;
  
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
    
    model TransmissionLineDynamicBranch "Purely inductive transmission line model with time-varying connection branch"
    extends TransmissionLineBase;
    equation
// This if-equation is not valid according to Modelica 3.5
// Section 8.3.4 but it would according to this extension proposal
      if closed then
        port_a.omegaRef = port_b.omegaRef;
        Connections.branch(port_a.omegaRef, port_b.omegaRef);
      end if;
    end TransmissionLineDynamicBranch;
    
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
  Connections.potentialRoot(port.omegaRef,p);
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
        redeclare TransmissionLineDynamicBranch T2(B = -10.0, open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System4;
    
    model System5 "Two generators, two parallel (one with breaker) and one series, fixed branches"
      extends System2(T1b(open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System5;
    
    model System6 "Two generators, two-parallel and one series line with breaker, dynamic branches"
      extends System5(
        redeclare TransmissionLineDynamicBranch T1b(B = -5.0, open = if time < 10 then false else true));
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
        redeclare TransmissionLineDynamicBranch T2(B = -10.0, open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System8;
  
    model System9 "System with loads that can get disconnected"
      Generator G1;
      Load L1(P = 0.8);
      replaceable Load L2(P = 0.1)
        constrainedby LoadBase;
      replaceable Load L3(P = 0.1)
        constrainedby LoadBase;
      replaceable TransmissionLine T1(open = if time < 10 then false else true)
        constrainedby TransmissionLineBase;
      TransmissionLine T2;
    equation
      connect(G1.port, L1.port);
      connect(G1.port, T1.port_a);
      connect(T1.port_b, L2.port);
      connect(T1.port_b, T2.port_a);
      connect(T2.port_b, L3.port);
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System9;
    
    model System10 "System with loads that can get disconnected, using dynamic overconstrained connectors"
      extends System9(
        redeclare LoadVariableRoot L2(P = 0.1),
        redeclare LoadVariableRoot L3(P = 0.1),
        redeclare TransmissionLineDynamicBranch T1(open = if time < 0 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System10;
  end PowerGridsReal;

  package IncompressibleFluid
    import SI = Modelica.SIunits;
    
    type CircuitIdentifier "Overconstrained type for circuit identifier"
      extends SI.PerUnit;
      function equalityConstraint
        input CircuitIdentifier id1;
        input CircuitIdentifier id2;
        output SI.PerUnit residue[0] "No constraints";
      end equalityConstraint;
    end CircuitIdentifier;
  
  
    connector FluidPort "Port for incompressible fluid flow"
      SI.Pressure p "Fluid pressure";
      flow SI.MassFlowRate w "Mass flow rate";
      CircuitIdentifier id "Circuit id number";
  annotation(
        Icon(graphics = {Ellipse(lineColor = {0, 0, 255}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}));
    end FluidPort;
  
    partial model FlowModel "Base class for flow models"
      FluidPort inlet annotation(
        Placement(visible = true, transformation(origin = {-68, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FluidPort outlet annotation(
        Placement(visible = true, transformation(origin = {-68, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SI.MassFlowRate w;
      SI.Pressure dp;
    equation
      w = inlet.w;
      inlet.w + outlet.w = 0;
    end FlowModel;
  
    model Pump "Simple pump model with linear flow characteristic law"
      extends FlowModel;
      parameter SI.Pressure dp0 = 1e5;
      parameter SI.MassFlowRate w0 = 1;
      
      SI.PerUnit n = 1 "Normalized rotational speed";
    equation
      dp = dp0*n^2 - w/w0*dp0*n;
      dp = outlet.p - inlet.p;
      Connections.branch(inlet.id, outlet.id);
      inlet.id = outlet.id;
    annotation(
        Icon(graphics = {Ellipse(extent = {{-60, 60}, {60, -60}}), Line(origin = {-77.2071, 1.20711}, points = {{17.2071, -1.20711}, {-14.7929, -1.20711}, {-16.7929, 0.792893}}), Line(origin = {76, 0}, points = {{-16, 0}, {16, 0}}), Polygon(origin = {14, 0}, points = {{-46, 50}, {-46, -50}, {46, 0}, {-46, 50}, {46, 0}, {-46, 50}}), Text(origin = {0, -80}, extent = {{-100, 20}, {100, -20}}, textString = "%name")}));
    end Pump;
  
    partial model BaseValve "Base model for simple on-off valve model"
      extends FlowModel;
      
      parameter SI.Pressure dp0 = 1e4;
      parameter SI.MassFlowRate w0 = 1;
      final parameter Real Kv0 = w0/dp0;
      
      discrete Real Kv "Flow coefficient";
      Boolean closed "State of line breaker";
      Boolean open = false "Command to open the line breaker";
      Boolean close = false "Command to close the line breaker";
    initial equation
      closed = false;
      Kv = Kv0;
    equation
      w = Kv*dp;
      when open then
        closed = false;
        Kv = Kv0;
      elsewhen close then
        closed = true;
        Kv = 0;
      end when;
      dp = inlet.p - outlet.p;
    annotation(
        Icon(graphics = {Polygon(origin = {0, 10}, points = {{-40, 30}, {40, -52}, {40, 30}, {-40, -50}, {-40, -50}, {-40, 30}}), Line(origin = {-69, -1.84746}, points = {{-29, 1}, {29, 1}}), Line(origin = {68.7458, -1.84746}, points = {{-29, 1}, {29, 1}}), Text(origin = {0, -80}, extent = {{-100, 20}, {100, -20}}, textString = "%name")}));
    end BaseValve;
    
    model Valve "Valve model with fixed connection branch"
      extends BaseValve;
    equation
      Connections.branch(inlet.id, outlet.id);
      inlet.id = outlet.id;
    end Valve;
    
    model ValveDynamicBranch "Valve model with dynamic connection branch"
      extends BaseValve;
    equation
      if closed then
        Connections.branch(inlet.id, outlet.id);
        inlet.id = outlet.id;
      end if;
    end ValveDynamicBranch;
  
    model Pipe "Simple pipe model with linear pressure loss"
      extends FlowModel;
      parameter SI.Pressure dp0 = 1e4;
      parameter SI.MassFlowRate w0 = 1;
    equation
      w = w0/dp0*dp;
      Connections.branch(inlet.id, outlet.id);
      inlet.id = outlet.id;
      dp = inlet.p - outlet.p;
    annotation(
        Icon(graphics = {Rectangle(fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-60, 20}, {60, -20}}), Line(origin = {-80, 0}, points = {{-20, 0}, {20, 0}}), Line(origin = {78, 0}, points = {{-18, 0}, {18, 0}}), Text(origin = {0, -40}, extent = {{-100, 20}, {100, -20}}, textString = "%name")}));
    end Pipe;
  
    model ExpansionTank "Ideal expansion tank, fixes the pressure of the circuit it is attached to"
      parameter SI.Pressure p0 "Fixed pressure";
      parameter Integer id = 0 "Circuit id";
      parameter Integer priority = 0 "Priority for selection as root node";
  
      FluidPort inlet  annotation(
        Placement(visible = true, transformation(origin = {-26, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation    
      Connections.potentialRoot(inlet.id, priority);
      // If the tank is selected as root node, it determines the circuit pressure and id
      // otherwise it behaves like a plug
      if Connections.isRoot(inlet.id) then
        inlet.p = p0;
        inlet.id = id;
      else
        inlet.w = 0;
      end if;
    annotation(
        Icon(graphics = {Line(origin = {0, -60}, points = {{0, 20}, {0, -20}, {0, -20}}), Polygon(origin = {0, 10}, points = {{-40, 50}, {-40, -50}, {40, -50}, {40, 50}, {-40, 50}}), Rectangle(origin = {0, -7}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{40, 33}, {-40, -33}}), Text(origin = {0, 80}, extent = {{-100, 20}, {100, -20}}, textString = "%name")}));
    end ExpansionTank;
    
    model System1 "Simple loop for testing"
      Pump pump(n = if time < 1 then 1 else 0.8) annotation(
        Placement(visible = true, transformation(origin = {-8, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Valve valve(w0 = 1, dp0 = 1e4) annotation(
        Placement(visible = true, transformation(origin = {22, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Pipe pipe1(w0 = 1, dp0 = 5e4) annotation(
        Placement(visible = true, transformation(origin = {40, -20}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Pipe pipe2(w0 = 1, dp0 = 5e4) annotation(
        Placement(visible = true, transformation(origin = {8, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      ExpansionTank tank annotation(
        Placement(visible = true, transformation(origin = {-18, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(pump.outlet, valve.inlet) annotation(
        Line(points = {{2, -6}, {12, -6}}, color = {0, 0, 255}));
      connect(valve.outlet, pipe1.inlet) annotation(
        Line(points = {{32, -6}, {40, -6}, {40, -10}}, color = {0, 0, 255}));
      connect(pipe1.outlet, pipe2.inlet) annotation(
        Line(points = {{40, -30}, {40, -34}, {18, -34}}, color = {0, 0, 255}));
      connect(pipe2.outlet, pump.inlet) annotation(
        Line(points = {{-2, -34}, {-26, -34}, {-26, -6}, {-18, -6}}, color = {0, 0, 255}));
    connect(tank.inlet, pump.inlet) annotation(
        Line(points = {{-18, 6}, {-18, -6}}, color = {0, 0, 255}));
    annotation(experiment(StopTime = 2, Interval = 0.02),
        Diagram(coordinateSystem(extent = {{-40, 40}, {60, -60}})));
    end System1;
    
    model System2 "Network with two separable loops"
      Pump pumpA(n = if time < 1 then 1 else 0.8) annotation(
        Placement(visible = true, transformation(origin = {-66, 10}, extent = {{-10, 10}, {10, -10}}, rotation = 90)));
      Pipe pipe1A(w0 = 1, dp0 = 3e4) annotation(
        Placement(visible = true, transformation(origin = {-44, 40}, extent = {{10, 10}, {-10, -10}}, rotation = 180)));
      Pipe pipe2A(w0 = 1, dp0 = 3e4) annotation(
        Placement(visible = true, transformation(origin = {-20, 10}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Pipe pipe3A(w0 = 1, dp0 = 3e4) annotation(
        Placement(visible = true, transformation(origin = {-44, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      Pump pumpB annotation(
        Placement(visible = true, transformation(origin = {60, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Pipe pipe1B(w0 = 1, dp0 = 3e4) annotation(
        Placement(visible = true, transformation(origin = {40, 40}, extent = {{-10, 10}, {10, -10}}, rotation = 180)));
      Pipe pipe2B(w0 = 1, dp0 = 3e4) annotation(
        Placement(visible = true, transformation(origin = {20, 10}, extent = {{-10, 10}, {10, -10}}, rotation = 270)));
      Pipe pipe3B(w0 = 1, dp0 = 3e4) annotation(
        Placement(visible = true, transformation(origin = {42, -20}, extent = {{10, -10}, {-10, 10}}, rotation = 180)));
      replaceable Valve valveAB(close = if time < 2 then false else true)  annotation(
        Placement(visible = true, transformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)))
        constrainedby BaseValve;
      replaceable Valve valveBA(close = if time < 4 then false else true) annotation(
        Placement(visible = true, transformation(origin = {0, -20}, extent = {{10, -10}, {-10, 10}}, rotation = 0)))
        constrainedby BaseValve;
      ExpansionTank tankA(id = 1, p0 = 2e5) annotation(
        Placement(visible = true, transformation(origin = {-90, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ExpansionTank tankB(id = 2, p0 = 2e5) annotation(
        Placement(visible = true, transformation(origin = {90, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(pumpA.outlet, pipe1A.inlet) annotation(
        Line(points = {{-66, 20}, {-66, 40}, {-54, 40}}, color = {0, 0, 255}));
  connect(pipe1A.outlet, pipe2A.inlet) annotation(
        Line(points = {{-34, 40}, {-20, 40}, {-20, 20}}, color = {0, 0, 255}));
  connect(pipe2A.outlet, pipe3A.inlet) annotation(
        Line(points = {{-20, 0}, {-20, -20}, {-34, -20}}, color = {0, 0, 255}));
  connect(pipe3A.outlet, pumpA.inlet) annotation(
        Line(points = {{-54, -20}, {-66, -20}, {-66, 0}}, color = {0, 0, 255}));
  connect(pumpB.outlet, pipe1B.inlet) annotation(
        Line(points = {{60, 20}, {60, 40}, {50, 40}}, color = {0, 0, 255}));
  connect(pipe1B.outlet, pipe2B.inlet) annotation(
        Line(points = {{30, 40}, {20, 40}, {20, 20}}, color = {0, 0, 255}));
  connect(pipe2B.outlet, pipe3B.inlet) annotation(
        Line(points = {{20, 0}, {20, -20}, {32, -20}}, color = {0, 0, 255}));
  connect(pipe3B.outlet, pumpB.inlet) annotation(
        Line(points = {{52, -20}, {60, -20}, {60, 0}}, color = {0, 0, 255}));
  connect(pipe1A.outlet, valveAB.inlet) annotation(
        Line(points = {{-34, 40}, {-10, 40}}, color = {0, 0, 255}));
  connect(pipe1B.outlet, valveAB.outlet) annotation(
        Line(points = {{30, 40}, {10, 40}}, color = {0, 0, 255}));
  connect(pipe3A.inlet, valveBA.outlet) annotation(
        Line(points = {{-34, -20}, {-10, -20}}, color = {0, 0, 255}));
  connect(valveBA.inlet, pipe3B.inlet) annotation(
        Line(points = {{10, -20}, {32, -20}}, color = {0, 0, 255}));
  connect(tankA.inlet, pumpA.inlet) annotation(
        Line(points = {{-90, 8}, {-90, 0}, {-66, 0}}, color = {0, 0, 255}));
  connect(tankB.inlet, pumpB.inlet) annotation(
        Line(points = {{90, 8}, {90, 0}, {60, 0}}, color = {0, 0, 255}));
    annotation(experiment(StopTime = 3, Interval = 0.02),
        Diagram(coordinateSystem(extent = {{-100, 60}, {100, -40}})));
    end System2;

    

    model System3 "Same as System 2 but with longer simulation time, exposing singularity"
      extends System2;
annotation(experiment(StopTime = 5, Interval = 0.02));
    end System3;
    
    model System4 "Same as System 3 but using dynamic overconstrained connectors"
      extends System2(redeclare ValveDynamicBranch valveAB(close = if time < 2 then false else true),
                      redeclare ValveDynamicBranch valveBA(close = if time < 4 then false else true));
    annotation(experiment(StopTime = 5, Interval = 0.02));
    end System4;

  /*  
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
        redeclare TransmissionLineDynamicBranch T2(B = -10.0, open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System4;
    
    model System5 "Two generators, two parallel (one with breaker) and one series, fixed branches"
      extends System2(T1b(open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System5;
    
    model System6 "Two generators, two-parallel and one series line with breaker, dynamic branches"
      extends System5(
        redeclare TransmissionLineDynamicBranch T1b(B = -5.0, open = if time < 10 then false else true));
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
        redeclare TransmissionLineDynamicBranch T2(B = -10.0, open = if time < 10 then false else true));
    annotation(experiment(StopTime = 50, Interval = 0.02));
    end System8;
    */
  end IncompressibleFluid;

annotation(uses(Modelica(version="3.2.3")));
end DynamicOverconstrainedConnectors;
