function t=pmntanh(x)
% PMNTANH
% -------
% Fast hyperbolic tangent function to be used in
% neural networks instead of the tanh provided by MATLAB
t=1-2./(exp(2*x)+1);

