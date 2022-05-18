function y = smoothHeaviside(x,slope)
    if nargin<2
        slope = 10^4; %10^4
    end
    y = (atan(slope*x) + pi/2)/pi; 
end