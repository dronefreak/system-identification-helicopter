function [w, A, b] = lsqfast(samples, policy, new_policy, firsttime)
  
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
% [w, A, b] = lsqfast(samples, policy, new_policy, firsttime)
%
% Evaluates the "policy" using the set of "samples", that is, it
% learns a set of weights for the basis specified in new_policy to
% form the approximate Q-value of the "policy" and the improved
% "new_policy". The approximation is the fixed point of the Bellman
% equation.
%
% "firsttime" is a flag (0 or 1) to indicate whether this is the first
% time this set of samples is processed. Preprossesing of the set is
% triggered if "firstime"==1.
%
% Returns the learned weights w and the matrices A and b of the
% linear system Aw=b. 
%
% See also lsq.m for a slower (incremental) implementation. 
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  persistent Phihat;
  persistent Rhat;
  
  
  %%% Initialize variables
  howmany = length(samples);
  k = feval(new_policy.basis);
  A = zeros(k, k);
  b = zeros(k, 1);
  PiPhihat = zeros(howmany,k);
  mytime = cputime;
  
  
  %%% Precompute Phihat and Rhat for all subsequent iterations
  if firsttime == 1
    
    Phihat = zeros(howmany,k);
    Rhat = zeros(howmany,1);
    
    for i=1:howmany
      phi = feval(new_policy.basis, samples(i).state, samples(i).action);
      Phihat(i,:) = phi';
      Rhat(i) = samples(i).reward;
    end
    
  end
  
  
  %%% Loop through the samples 
  for i=1:howmany
    
    %%% Make sure the nextstate is not an absorbing state
    if ~samples(i).absorb
      
      %%% Compute the policy and the corresponding basis at the next state 
      nextaction = policy_function(policy, samples(i).nextstate);
      nextphi = feval(new_policy.basis, samples(i).nextstate, nextaction);
      PiPhihat(i,:) = nextphi';
      
    end
    
  end
  
  
  %%% Compute the matrices A and b
  A = Phihat' * (Phihat - new_policy.discount * PiPhihat);
  b = Phihat' * Rhat;
  
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
  
