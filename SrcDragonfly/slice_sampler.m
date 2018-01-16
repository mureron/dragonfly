% ----------------------------------------------------------
% UNIVARIATE SLICE SAMPLER - stepping out (Neal, 2003)
% W: optimal value in the range (3,10)*std(x)
%    - see C.Planas and A.Rossi (2014)
% func(x,lambda): any unnormalized pdf  f(x) [-log]
% with lambda (optional) a vector of auxiliaty parameters
% to be passed to f( ).
% ----------------------------------------------------------
function [theta, fxsim, neval] = slice_sampler(objective_function,theta,thetaprior,sampler_options,varargin)

if sampler_options.rotated %&& ~isempty(sampler_options.V1),
    [theta, fxsim, neval] = rotated_slice_sampler(objective_function,theta,thetaprior,sampler_options,varargin{:});
    if isempty(sampler_options.mode), % jumping 
       return,
    end
end
    
theta=theta(:);
npar = length(theta);
neval = zeros(npar,1);
W1 = sampler_options.W1;

for it=1:npar,
    neval(it) = 0;
    W = W1(it); 
    xold  = theta(it);
   % XLB   = thetaprior(3);
   % XUB   = thetaprior(4);
    XLB   = thetaprior(it,1);
    XUB   = thetaprior(it,2);
   
    
    % -------------------------------------------------------
    % 1. DRAW Z = ln[f(X0)] - EXP(1) where EXP(1)=-ln(U(0,1))
    %    THIS DEFINES THE SLICE S={x: z < ln(f(x))}
    % -------------------------------------------------------
    fxold = -feval(objective_function,theta,varargin{:});
    neval(it) = neval(it) + 1;
    Z = fxold + log(rand(1,1));
    % -------------------------------------------------------------
    % 2. FIND I=(L,R) AROUND X0 THAT CONTAINS S AS MUCH AS POSSIBLE
    %    STEPPING-OUT PROCEDURE
    % -------------------------------------------------------------
    u = rand(1,1);
    L = max(XLB,xold-W*u);
    R = min(XUB,L+W);
    while(L > XLB)
        xsim = L;
        theta(it) = xsim;
         fxl = -feval(objective_function,theta,varargin{:});
       % fxl = log(feval(objective_function,theta,varargin{:}));
        neval(it) = neval(it) + 1;
        if (fxl <= Z)
            break;
        end
        L = max(XLB,L-W);
    end
    while(R < XUB)
        xsim = R;
        theta(it) = xsim;
        fxr = -feval(objective_function,theta,varargin{:});
        % fxr = log(feval(objective_function,theta,varargin{:}));
        neval(it) = neval(it) + 1;
        if (fxr <= Z)
            break;
        end
        R = min(XUB,R+W);
    end
    % ------------------------------------------------------
    % 3. SAMPLING FROM THE SET A = (I INTERSECT S) = (LA,RA)
    % ------------------------------------------------------
    fxsim = Z-1;
    while (fxsim < Z)
        u = rand(1,1);
        xsim = L + u*(R - L);
        theta(it) = xsim;
        fxsim = -feval(objective_function,theta,varargin{:});
        %fxsim = log(feval(objective_function,theta,varargin{:}));
        neval(it) = neval(it) + 1;
        if (xsim > xold)
            R = xsim;
        else
            L = xsim;
        end
    end
    
end

