function f = smoothStepFunction(a,b,t)
    f = smoothHeaviside(t-a) - smoothHeaviside(t-b);
end