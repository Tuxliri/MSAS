function u = command(t,data)
% Distributor control command

t1 = data.command.t1;
t2 = data.command.t2;

if t <= t1
    u = 0;
elseif t > t1 && t < t2
    u = (t - t1)/(t2 - t1);
elseif t >= t2
    u = 1;
end

end