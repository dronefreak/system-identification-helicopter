function qvalue = Qvalue(state, action, policy)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2000-2002 
%
% Michail G. Lagoudakis (mgl@cs.duke.edu)
% Ronald Parr (parr@cs.duke.edu)
%
% Department of Computer Science
% Box 90129
% Duke University, NC 27708
% 
%
% qvalue = Qvalue(state, action, policy)
%   
% Returns Q^policy(state,action) = 
%             policy.basis(state,action)' * policy.weights
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  phi = feval(policy.basis, state, action);
  qvalue = phi' * policy.weights;

  
  return
