function []=drawnet(w1,w2,CancelVal,instr,outstr)
%  DRAWNET
%  -------
%          Draws a two layer feedforward neural network.
%
%          drawnet(W1,W2)
%          drawnet(W1,W2,CancelVal)
%          drawnet(W1,W2,CancelVal,instring,outstring)
%          DRAWNET draws the network specified by the weight matrices 
%          W1 and W2. Positive weights are represented by a solid line while 
%          a dashed line represents a negative weight. Only weights and biases
%          larger than 'CancelVal' are drawn. A bias is represented by a vertical
%          line through the neuron.
%
%  INPUT:
%  W1       : Input-to-hidden layer weights. The matrix dimension is
%             dim(W1) = [(# of hidden units) * (inputs + 1)] (the 1 is due to the bias)
%  W2       : hidden-to-output layer weights.
%             dim(W2) = [(outputs)  *  (# of hidden units + 1)]
%  CancelVal: (Optional) Draw only weights/biases exceeding this value
%  instring : (Optional) A cell structure containing in each cell a 
%             string to be assigned to the corresponding input. The number of
%             cells should thus match the number of inputs. If it is not present, or
%             it is empty {}, the inputs are simply numbered.
%  outstring: (Optional). A cell structure containing in each cell a string to be
%             assigned to the corresponding output. The number of cells should thus
%             match the number of outputs.
%
%  See also OBDPRUNE, OBSPRUNE, NNPRUNE

%  Original function programmed by Claus Svarer, EI/CONNECT. Current version is
%  modified by Magnus Norgaard, IAU/IMM, Technical University of Denmark.
%  LastEditDate: Jan. 15, 2000
[N1,N0]=size(w1);                % Dimension of input-to-hidden weight matrix
N0=N0-1;
[N2,dummy]=size(w2);             % Dimension of hidden-to-output weight matrix
if nargin<3 | isempty(CancelVal)
   CancelVal=0;
end
if nargin>3,
   if ~iscell(instr),
      error('Argument "instring" must be cell structure');
   elseif ~isempty(instr) & length(instr)~=N0
      error('The number of cells in "instring" must correspond to the number of inputs');
   end
   if isempty(instr)
      instr={};
   end
end
if nargin>4
   if ~iscell(outstr),
      error('Argument "outstring" must be cell structure');
   elseif ~isempty(outstr) & length(outstr)~=N2
      error('The number of cells in "outstring" must correspond to the number of outputs');
   end
   if isempty(outstr)
      outstr={};
   end
end
MaxNeu=max([N0 N1 N2]);
cla
LengthTres=0.025*MaxNeu;
axis([-0.1 2.1 0.5 MaxNeu+0.5]);
axis('off')

hold on
for i = 1:N0,
   plot(0,(MaxNeu/(N0+1))*i,'ko');
   if nargin<=3 | isempty(instr),
     text(-0.1,(MaxNeu/(N0+1))*i-0.0,sprintf('%g',i));
   else
     text(-0.4,(MaxNeu/(N0+1))*i-0.0,sprintf(instr{i}));
   end
end;
for i = 1:N1,
   plot(1,(MaxNeu/(N1+1))*i,'ro');
   if (w1(i,N0+1) ~= 0)   
      plot([1 1],[((MaxNeu/(N1+1))*i-LengthTres) ((MaxNeu/(N1+1))*i+LengthTres)],'r');
   end;
end;
for i = 1:N2,
   plot(2,(MaxNeu/(N2+1))*i,'ro');
   if (w2(i,N1+1) ~= 0)   
      plot([2 2],[((MaxNeu/(N2+1))*i-LengthTres) ((MaxNeu/(N2+1))*i+LengthTres)],'r');
   end;
   if nargin==5 & ~isempty(outstr),
      text(2.05,(MaxNeu/(N2+1))*i-0.0,sprintf(outstr{i}));
   end  
end;

MaxColorNo=7;
ColorMatrix = [zeros(MaxColorNo,2) [1:-1/(MaxColorNo-1):0]'];
for i=1:N0,
   for j=1:N1,
      colour_int = ceil(abs(w1(j,i))*2+eps);
      if (colour_int > MaxColorNo),
         colour_int = MaxColorNo;
      end;
      colour=ColorMatrix(colour_int,:);
      if (w1(j,i) > 0),
         connect='-';
      else
         connect='--';
      end;
      if (abs(w1(j,i)) >= CancelVal),
         plot([0 1],[(MaxNeu/(N0+1))*i (MaxNeu/(N1+1))*j],connect,'Color',colour);
      end;
   end;
end
for i=1:N1,
   for j=1:N2,
      colour_int = ceil(abs(w2(j,i))*2+eps);
      if (colour_int > MaxColorNo),
         colour_int = MaxColorNo;
      end;
      if (w2(j,i) > 0),
         connect='-';
      else
         connect='--';
      end;
      if (abs(w2(j,i)) >= CancelVal),
         plot([1 2],[(MaxNeu/(N1+1))*i (MaxNeu/(N2+1))*j],connect,'Color',colour);
      end;
   end;
end
hold off
