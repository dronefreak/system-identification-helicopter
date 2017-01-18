function [new_results, totdrew, toturew] = execute(initial_state, ...
						  simulator, ...
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
% [new_results, totdrew, toturew] = execute(initial_state, simulator,
%                                           policy, maxsteps)
%
% Executes one episode (of at most "maxsteps" steps) on the
% "simulator" starting at the "initial_state" and using the "policy"
% to select actions.
%
% Returns all the samples collected during the episode and the total
% discounted and undiscounted reward.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Initialize the random number generator to a random state
  rand('state', sum(100*clock));

  %%% Initialize storage for new samples 
  empty_result.state = feval(simulator);
  empty_result.action = 0;
  empty_result.reward = 0.0;
  empty_result.nextstate = empty_result.state;
  empty_result.absorb = 0;
  
  results = repmat(empty_result, 1, maxsteps);
  
  
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
    
    %%% Record sample
    results(steps).state = state;
    results(steps).action = action;
    results(steps).reward = reward;
    results(steps).nextstate = nextstate;
    results(steps).absorb = endsim;
    
    %%% Update the total reward(s)
    totdrew = totdrew + (policy.discount)^(steps-1) * reward;
    toturew = toturew + reward;

    %%% Continue
    state = nextstate;
    
  end
  
  
  %%% Return the results
  new_results = results(1:steps);
  
  
  return
  
