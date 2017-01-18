function tr = settrain(trparms,varargin)
%SETTRAIN set parameters for training algorithm.
%   It is only necessary to set parameters specific to the selected
%   training algorithm. In case parameters that are needed in the training
%   algorithm are not set, the called training function will automatically
%   set these parameters to the default values.
%
%   TRPARMS = SETTRAIN
%      Set all parameters to default values.
%
%   SETTRAIN(TRPARMS)
%      List all parameters.
%
%   TRPARMS = SETTRAIN(TRPARMS,'field1',value1,'field2',value2,...)
%      Set specific parameters
%         TRPARMS.field1 = value1;
%         TRPARMS.filed2 = value2;
%         etc.
%      If value = 'default', the parameter is set to the default value.
%
%   The following fields are valid:
%   Information displayed during training
%    infolevel - Display little information (0) or much (1)
%
%   Stopping criteria (all algorithms)
%    maxiter   - Maximum iterations.
%    critmin   - Stop if criterion is below this value.
%    critterm  - Stop if change in criterion is below this value.
%    gradterm  - Stop if largest element in gradient is below this value.
%    paramterm - Stop if largest parameter change is below this value.
%    NB: critterm, gradterm and paramterm must all be satisfied. 
%
%   Weight decay (all algorithms trained with the Levenberg-Marquardt alg.).
%    D         - Row vector containing the weight decay parameters. If D has
%                one element a scalar weight decay will be used. If D has two
%                elements, the first element will be used as weight decay for
%                the hidden-to-output layer while second will be used for the
%                input-to-hidden layer weights. For individual weight decays,
%                D must contain as many elements as there are weights in the
%                network.
%
%   Levenberg-Marquardt parameters
%    lambda    - Initial Levenberg-Marquardt parameter. 
%
%   Backprop parameters
%    eta       - Step size.
%    alph      - Momentum.
%
%   RPE parameters
%    method    - Training method ('ff', 'ct', 'efra').
%
%    Forgetting factor
%    fflambda  - Forgetting factor.
%    p0        - Covariance matrix is initialized to p0*I.
%
%    Constant trace
%    ctlambda  - Forgetting factor.
%    alpha_min - Min. eigenvalue of P matrix.
%    alpha_max - Max. eigenvalue of P matrix.
%
%    EFRA
%    eflambda  - Forgetting factor.
%    alpha     - EFRA parameter.
%    beta      - EFRA parameter.
%    delta     - EFRA parameter.
%
%   For recurrent nets
%    skip      - Do not use the first 'skip' samples for training.
%
%   For multi-output nets
%    repeat    - Number of times the IGLS procedure should be repeated.

%  Programmed by : Magnus Norgaard, IAU/IMM 
%  LastEditDate  : Dec. 29, 1999

% >>>>>>>>>>>>>>>>>>>>>   SET ALL PARAMETERS TO DEFAULT   <<<<<<<<<<<<<<<<<<<<<
% Information level
trd.infolevel = 0;

% Termination values
trd.maxiter = 500;
trd.critmin = 0;
trd.critterm = 1e-7;
trd.gradterm = 1e-4;
trd.paramterm = 1e-3;

% Weight decay
trd.D = 0;

% Levenberg-Marquardt parameters
trd.lambda = 1;

% Backprop parameters
trd.eta  = 1e-4;
trd.alph = 0;

% RPE parameters
trd.method = 'ff';
trd.fflambda = 0.995;
trd.p0 = 10;

trd.alpha_min = 1e-3;
trd.alpha_max = 1e1;

trd.eflambda = 0.995;
trd.alpha = 1;
trd.beta  = 0.001;
trd.delta = 0.001;

% For recurrent nets
trd.skip = 0;

% For multi-output nets
trd.repeat = 5;

% Default names
dnames = fieldnames(trd);

if nargin==0
  tr = trd;
  
% >>>>>>>>>>>>>>>>>>>>>>>>>>>   DISPLAY PROPERTIES   <<<<<<<<<<<<<<<<<<<<<<<<<<<  
elseif nargin==1,
  names = fieldnames(trparms);
  for idx=1:length(names),
    tmp = getfield(trparms,names{idx});
    if ischar(tmp),
      fprintf('%15s  = %s\n',names{idx},tmp);
    elseif (size(tmp,1)==1 | size(tmp,1)==1)
      if rem(tmp,1)==0,
        fprintf('%15s  = %d\n',names{idx},tmp);
      else
        fprintf('%15s  = %4.3e\n',names{idx},tmp);
      end
    else
       fprintf('%15s  = [%dx%d double]\n',names{idx},size(tmp,1),size(tmp,2));
    end
  end

% >>>>>>>>>>>>>>>>>>>>>>>>>   SET SPECIFIC PROPERTIES   <<<<<<<<<<<<<<<<<<<<<<<<
elseif nargin>=2,
  tr = trparms;
  if rem(length(varargin),2),
    error('You must specify an even number of properties.');
  end
  for idx=1:2:length(varargin)
  
    % Check if field is a string
    if ~isstr(varargin{idx})
      error('Property name must be a string.');
      
    % Check if field is illegal
    elseif(isempty(find(strcmp(lower(dnames),lower(varargin{idx})))))
       errstr = sprintf('%s ''%s''.','Unknown property name',varargin{idx});
       error(errstr);
       
    % Set field to default value if requested
    elseif(strcmp(lower(varargin{idx+1}),'default'))
      if strcmp(lower(varargin{idx}),'d'),
        tr = setfield(tr,'D',getfield(trd,'D'));
      else
        tr=setfield(tr,lower(varargin{idx}),getfield(trd,lower(varargin{idx})));
      end
    
    % Set field to specified value
    else
      if strcmp(lower(varargin{idx}),'d'),
        tr = setfield(tr,'D',varargin{idx+1});
      else
        tr = setfield(tr,lower(varargin{idx}),varargin{idx+1});
      end
    end
  end
end
