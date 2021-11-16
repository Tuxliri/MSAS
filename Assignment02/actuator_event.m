function [value,isterminal,direction] = actuator_event(t,y,data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
xmax = data.actuator.cc_max;
x = y(2);
value = xmax-x;
isterminal = 1;
direction = -1;
end

