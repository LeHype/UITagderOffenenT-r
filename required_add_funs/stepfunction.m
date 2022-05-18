function f = stepfunction(a,b,t)
sympref('HeavisideAtOrigin',1);
f = heaviside(t-a) - heaviside(t-b);
end