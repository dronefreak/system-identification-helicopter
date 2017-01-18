function new_samples = collect_samples(domain, maxepisodes, maxsteps, policy) 	
  
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
% new_samples = collect_samples(domain, maxepisodes, maxsteps, policy)
%
% Collects samples from the "domain" using the "policy" by running at
% most "maxepisodes" episodes each of which is at most "maxsteps"
% steps long.
%
% Input:
% 
% domain - A string containing the name of the domain
%          This string should be the prefix for all related
%          functions, for example if domain is 'chain' functions
%          should be chain_initialize_policy, chain_simulator,
%          etc.
%
% maxepisodes - An integer indicating the maximum number of
%               episodes from which samples are collected.
% 
% maxsteps - An integer indicating the maximum number of steps of each
%            episode (an episode may finish earlier if an absorbing
%            state is encountered).
% 
% policy - (optional argument) A policy to be used for collecting the
%          samples. It should be given as a struct with the following
%          fields (at least):
%
%          explore  : Exploration rate (real number)
%          discount : Discount factor (real number)
%          actions  : Total numbers of actions, |A|
%          basis    : The function handle for the basis
%                     associated with this policy
%          weights  : A column array of weights 
%                     (one for each basis function)
%
%          If a policy is not provided, samples will be collected by a
%          purely random policy initialized with "explore"=1.0,
%          "discount" and "basis" some dummy values, and "actions" and
%          "weights" as suggested by the domain (in the
%          domain_initialize_policy function.
%
% Output:
% 
% new_samples - An array of the collected samples. Each entry
%               new_samples(i) has the following form:
%
%               new_samples(i).state     : Arbitrary description of state
%               new_samples(i).action    : An integer in [1,|A|]
%               new_samples(i).reward    : A real value
%               new_samples(i).nextstate : Arbitrary description
%               new_samples(i).absorb    : Absorbing nextstate? (0 or 1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%% Initialize some variables
  simulator = [domain '_simulator'];
  initialize_policy = [domain '_initialize_policy'];
  initialize_state = [domain '_initialize_state'];
  
  
  %%% Initialize a purely random policy (if one is not provided)
  if nargin<7
    policy = feval(initialize_policy, 1.0, 0.0, @fake_basis);
  end
  

  %%% Initialize storage for new samples 
  empty_result.state = feval(simulator);
  empty_result.action = 0;
  empty_result.reward = 0.0;
  empty_result.nextstate = empty_result.state;
  empty_result.absorb = 0;
  
  samples = repmat(empty_result, 1, maxepisodes*maxsteps);
  
  
  %%% Initialize variables
  nextslot = 1;
  episodes = 0;
  
  
  %%% Initialize simulator
  feval(simulator);
  
  
  %%% Main loop
  while (episodes < maxepisodes)
    
    %%% Select initial state
    initial_state = feval(initialize_state, simulator); 
    
    
    %%% Run one episode (up to the max number of steps)
    epi_samples = execute(initial_state, simulator, policy, maxsteps);
    
    
    %%% Store the new samples
    if length(epi_samples) > 0
      oldslot = nextslot; 
      nextslot = oldslot + length(epi_samples);
      samples(oldslot:nextslot-1) = epi_samples;
    end
    
    
    %%% Next episode
    episodes = episodes + 1;
    
  end
  
  
  %%% Return the new samples
  new_samples = samples(1:nextslot-1);
  clear samples
  
  
  return
