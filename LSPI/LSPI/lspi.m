function [policy, all_policies] = lspi(domain, algorithm, maxiterations, ...
				       epsilon, samples, basis, ...
				       discount, initial_policy)

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
% [policy, all_policies] = lspi(domain, algorithm, maxiterations,
%                               epsilon, samples, basis, discount,
%                               initial_policy)
%
% LSPI : Least-Squares Policy Iteration
%
% Finds a good policy given a set of samples and a basis
%
% Input:
% 
% domain - A string containing the name of the domain
%          This string should be the prefix for all related
%          functions, for example if domain is 'chain' functions
%          should be chain_initialize_policy, chain_simulator,
%          etc.
%
% algorithm - This is a number that indicates which evaluation
%             algorithm should be used (see the paper):
%
%             1-lsq       : The regular LSQ (incremental)
%             2-lsqfast   : A fast version of LSQ (uses more space)
%             3-lsqbe     : LSQ with Bellman error minimization 
%             4-lsqbefast : A fast version of LSQBE (more space)
%
%             LSQ is the evaluation algorithm for regular
%             LSPI. Use lsqfast in general, unless you have
%             really big sample sets. LSQBE is provided for
%             comparison purposes and completeness.
%
% maxiterations - An integer indicating the maximum number of
%                 LSPI iterations 
% 
% epsilon - A small real number used as the termination
%           criterion. LSPI converges if the distance between
%           weights of consequtive iterations is less than
%           epsilon.
%
% samples - The sample set. This should be an array where each
%            entry samples(i) has the following form: 
%
%            samples(i).state     : Arbitrary description of state
%            samples(i).action    : An integer in [1,|A|]
%            samples(i).reward    : A real value
%            samples(i).nextstate : Arbitrary description
%            samples(i).absorb    : Absorbing nextstate? (0 or 1)
%
% basis - The function that computes the basis for a given pair
%         (state, action) given as a function handle
%         (e.g. @chain_phi) or as a string (e.g. 'chain_phi').
% 
% discount - A real number in (0,1] to be used as the discount factor
%
% initial_policy - (optional argument) An initial policy for
%                  LSPI. It should be given as a struct with the
%                  following fields (at least):
%
%                  explore  : Exploration rate (real number)
%                  discount : Discount factor (real number)
%                  actions  : Total numbers of actions, |A|
%                  basis    : The function handle for the basis
%                             associated with this policy
%                  weights  : A column array of weights 
%                             (one for each basis function)
%
%                  If initial_policy is not provided it is initialized
%                  to a policy with "explore"=0.0, "discount" and
%                  "basis" as provided above, and "actions" and
%                  "weights" as suggested by the domain (in the
%                  domain_initialize_policy function.
%
% Output:
%
% policy - The learned policy (same struct as above)
% 
% all_policies - A cell array of size (iterations+1) containing
%                all the intermediate policies at each LSPI
%                iteration, including the initial policy. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Initialize the random number generator to a random state
  rand('state', sum(100*clock));

  
  %%% Create a new policy
  initialize_policy = [domain '_initialize_policy'];
  policy = feval(initialize_policy, 0.0, discount, basis);
  %%% Consider this to be the initial_policy if one is not provided
  if nargin<8
    initial_policy = policy;
  end

  
  %%% Initialize policy iteration 
  iteration = 0;
  distance = inf;
  all_policies{1} = initial_policy;
  
  
  %%% If no samples, return
  if length(samples)==0
    disp('Warning: Empty sample set');
    return
  end
  
  
  %%% Main LSPI loop  
  while ( (iteration < maxiterations) & (distance > epsilon) )
       
    %%% Update and print the number of iterations
    iteration = iteration + 1;
    disp('*********************************************************');
    disp( ['LSPI iteration : ', num2str(iteration)] );
    if (iteration==1)
      firsttime = 1;
    else
      firsttime = 0;
    end
    
 
    %%% You can optionally make a call to collect_samples right here
    %%% to change/update the sample set. Make sure firsttime is set
    %%% to 1 if you do so.
    

    %%% Evaluate the current policy (and implicitly improve)
    %%% There are several options here - choose one
    if (algorithm == 1)
      policy.weights = lsq(samples, all_policies{iteration}, policy);
    elseif (algorithm == 2)
      policy.weights = lsqfast(samples, all_policies{iteration}, ...
			       policy, firsttime);
    elseif (algorithm == 3)
      policy.weights = lsqbe(samples, all_policies{iteration}, policy);
    elseif (algorithm == 4)
      policy.weights = lsqbefast(samples, all_policies{iteration}, ...
				 policy, firsttime);
    end

    
    %%% Compute the distance between the current and the previous policy
    l1 = length(policy.weights);
    l2 = length(all_policies{iteration}.weights);
    if (l1 == l2)
      difference = policy.weights - all_policies{iteration}.weights;
      LMAXnorm = norm(difference,inf);
      L2norm = norm(difference);
    else
      LMAXnorm = abs(norm(policy.weights,inf) - ...
		     norm(all_policies{iteration}.weights,inf));
      L2norm = abs(norm(policy.weights) - ...
		   norm(all_policies{iteration}.weights));
    end
    distance = L2norm;
      
      
    
    %%% Print some information 
    disp( ['   Norms -> Lmax : ', num2str(LMAXnorm), ...
	   '   L2 : ',            num2str(L2norm)] );
    
    
    %%% Store the current policy
    all_policies{iteration+1} = policy;
    
    
    %%% Depending on the domain, print additional info if needed
    feval([domain '_print_info'], all_policies);
    
  end
  
  
  %%% Display some info
  disp('*********************************************************');
  if (distance > epsilon) 
    disp(['LSPI finished in ' num2str(iteration) ...
	  ' iterations WITHOUT CONVERGENCE to a fixed point']);
  else
    disp(['LSPI converged in ' num2str(iteration) ' iterations']);
  end
  disp('********************************************************* ');
  
  
  return
