function [w, A, b] = lsq(samples, policy, new_policy)

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
% [w, A, b] = lsq(samples, policy, new_policy)
%
% Evaluates the "policy" using the set of "samples", that is, it
% learns a set of weights for the basis specified in new_policy to
% form the approximate Q-value of the "policy" and the improved
% "new_policy". The approximation is the fixed point of the Bellman
% equation.
%
% Returns the learned weights w and the matrices A and b of the
% linear system Aw=b. 
%
% See also lsqfast.m for a faster (batch) implementation. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Initialize variables
  howmany = length(samples);
  k = feval(new_policy.basis);
  A = zeros(k, k);
  b = zeros(k, 1);
  mytime = cputime;
  
  
  %%% Loop through the samples 
  for i=1:howmany
    
    %%% Compute the basis for the current state and action
    phi = feval(new_policy.basis, samples(i).state, samples(i).action);
    
    
    %%% Make sure the nextstate is not an absorbing state
    if ~samples(i).absorb
      
      %%% Compute the policy and the corresponding basis at the next state 
      nextaction = policy_function(policy, samples(i).nextstate);
      nextphi = feval(new_policy.basis, samples(i).nextstate, nextaction);
      
    else
      nextphi = zeros(k, 1);
    end
  
    
    %%% Update the matrices A and b
    A = A + phi * (phi - new_policy.discount * nextphi)';
    b = b + phi * samples(i).reward;
    
  end

  phi_time = cputime - mytime;
  disp(['CPU time to form A and b : ' num2str(phi_time)]);
  mytime = cputime;
  
  %%% Solve the system to find w
  rankA = rank(A);
  
  rank_time = cputime - mytime;
  disp(['CPU time to find the rank of A : ' num2str(phi_time)]);
  mytime = cputime;
  
  disp(['Rank of matrix A : ' num2str(rankA)]);
  if rankA==k
    disp('A is a full rank matrix!!!');
    w = A\b;
  else
    disp(['WARNING: A is lower rank!!! Should be ' num2str(k)]);
    w = pinv(A)*b;
  end
  
  solve_time = cputime - mytime;
  disp(['CPU time to solve Aw=b : ' num2str(solve_time)]);
  
  return
  
