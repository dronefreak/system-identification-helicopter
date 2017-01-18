function [action, actionphi] = policy_function(policy, state)

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
% [action, actionphi] = policy_function(policy, state)
%
% Computes the "policy" at the given "state".
%
% Returns the "action" that the policy picks at that "state" and the
% evaluation ("actionphi") of the basis at the pair (state, action).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Exploration or not? 
  if (rand < policy.explore)
    
    %%% Pick one action in random
    action = randint(policy.actions);
    actionphi = feval(policy.basis, state, action);
    
  else
    
    %%% Pick the action with maximum Q-value
    bestq = -inf;
    besta = [];
    
    %%% Find first all actions with maximum Q-value
    for i=1:policy.actions
      
      phi = feval(policy.basis, state, i);
      q = phi' * policy.weights;
      
      if (q > bestq)
	bestq = q;
	besta = [i];
	actionphi = [phi];
      elseif (q == bestq)
	besta = [besta; i];
	actionphi = [actionphi, phi];
      end
      
    end
    
    %%% Now, pick one of them
    which = 1;                         % Pick the first (deterministic)
    %which = randint(length(besta));    % Pick randomly
    
    action = besta(which);
    actionphi = actionphi(:,which);
    
  end
  
  
  return
