function [steps, totdrew, toturew] = justexec(initial_state, simulator, ...
					      policy, maxsteps)

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
% [steps, totdrew, toturew] = justexec(initial_state, simulator,
%                                      policy, maxsteps)
%
% Executes one episode (of at most "maxsteps" steps) on the
% "simulator" starting at the "initial_state" and using the "policy"
% to select actions.
%
% It returns the total number of steps and the total discounted and
% undiscounted reward collected during the episode. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %%% Initialize the random number generator to a random state
  rand('state', sum(100*clock));

  %%% Initialize variables
  totdrew = 0;
  toturew = 0;
  steps = 0;
  endsim = 0;
  
  
  %%% Set initial state
  state = feval(simulator, initial_state);
  
  
  %%% Run the episode
  while ( (steps < maxsteps) & (~endsim) )
    
    steps = steps + 1;
    
    %%% Select action 
    action = policy_function(policy, state);
    
    %%% Simulate
    [nextstate, reward, endsim] = feval(simulator, state, action);
    
    %%% Update the total reward(s)
    totdrew = totdrew + (policy.discount)^(steps-1) * reward;
    toturew = toturew + reward;
    
    %%% Continue
    state = nextstate;
    
  end
  
  
  return
  
