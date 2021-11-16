%% Symbolic math toolbox test
syms R1 R2 L C v(t) vC(t) i1 i2 i3 vL
vars = [i1; i2; i3; v(t); vC(t)];

eqs = [i1 + i2 + i3 == 0;
       i1 == vL/R1;
       i2 == (vL - v(t) - vC(t))/R2;
       vL == L*diff(i3,t);
       i2 == C*diff(vC(t),t)];
   
%    [eqs, vars, newVars] = reduceDifferentialOrder(eqs, vars)