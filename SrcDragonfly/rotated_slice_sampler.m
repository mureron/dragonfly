% ----------------------------------------------------------
% UNIVARIATE SLICE SAMPLER - stepping out (Neal, 2003)
% W: optimal value in the range (3,10)*std(x)
%    - see C.Planas and A.Rossi (2014)
% func(x,lambda): any unnormilized pdf  f(x)
% with lambda (optional) a vector of auxiliaty parameters
% to be passed to f( ).
% ----------------------------------------------------------
function [theta, fxsim, neval] = rotated_slice_sampler(objective_function,theta,thetaprior,sampler_options,varargin)

theta=theta(:);
npar = length(theta);
neval = zeros(npar,1);
W1=[];
if isfield(sampler_options,'WR'),
    W1 = sampler_options.WR;
end
V1 = sampler_options.V1;
if sampler_options.random_rotation,
    % % random rotation
    % V1=rand(npar,npar)*2-1;
    % V1(:,1)=V1(:,1)./norm(V1(:,1));
    % % Gram-Schmidt
    % for j=2:npar,
    %     for k=1:j-1,
    %         V1(:,j)=V1(:,j)-V1(:,k)'*V1(:,j)*V1(:,k);
    %     end
    %     V1(:,j)=V1(:,j)/norm(V1(:,j));
    % end
    angle = rand*pi/2-pi/4;
    V1 = ortnormbaserot(npar, angle);
end
if ~isempty(sampler_options.mode),
    mm = sampler_options.mode;
    n = length(mm);
    r = randperm(n);
    r=r(1);
    V1 = mm(r).m;
    V1=V1-theta;
    V1=V1/norm(V1);
%     %d = chol(mm(r).invhess);
%     %V1 = transpose(feval(sampler_options.proposal_distribution, transpose(mm(r).m), d, npar));
% 
%     V1=eye(npar);
%     V1=V1(:,randperm(npar));
%     for j=1:2,
%         V1(:,j)=mm(r(j)).m-theta;
%         V1(:,j)=V1(:,j)/norm(V1(:,j));
%     end
%     % Gram-Schmidt
%     for j=2:npar,
%         for k=1:j-1,
%             V1(:,j)=V1(:,j)-V1(:,k)'*V1(:,j)*V1(:,k);
%         end
%         V1(:,j)=V1(:,j)/norm(V1(:,j));
%     end    
    for j=1:n,
        distance(j)=sqrt(sum((theta-mm(j).m).^2));
    end
    [m, im] = min(distance);
    if im==r, 
        fxsim=[];
        return,
    else
        theta1=theta;
    end
end
npar=size(V1,2);
    
for it=1:npar,
    theta0 = theta;
    neval(it) = 0;
    xold  = 0;
   % XLB   = thetaprior(3);
   % XUB   = thetaprior(4);
    tb=sort([(thetaprior(:,1)-theta)./V1(:,it) (thetaprior(:,2)-theta)./V1(:,it)],2);
    XLB=max(tb(:,1));
    XUB=min(tb(:,2));  
    if isempty(W1),
        W = (XUB-XLB)*0.8; 
    else
        W = W1(it);
    end
        
    % -------------------------------------------------------
    % 1. DRAW Z = ln[f(X0)] - EXP(1) where EXP(1)=-ln(U(0,1))
    %    THIS DEFINES THE SLICE S={x: z < ln(f(x))}
    % -------------------------------------------------------
    
    fxold = -feval(objective_function,theta,varargin{:});
    %I have to be sure that the rotation is for L,R or for Fxold, theta(it)
    neval(it) = neval(it) + 1;
    Z = fxold + log(rand(1,1));
    % -------------------------------------------------------------
    % 2. FIND I=(L,R) AROUND X0 THAT CONTAINS S AS MUCH AS POSSIBLE
    %    STEPPING-OUT PROCEDURE
    % -------------------------------------------------------------
    u = rand(1,1);
    L = max(XLB,xold-W*u);
    R = min(XUB,L+W);
    
    %[L R]=slice_rotation(L, R, alpha);
    while(L > XLB)
        xsim = L;
        theta = theta0+xsim*V1(:,it);
        fxl = -feval(objective_function,theta,varargin{:});
        neval(it) = neval(it) + 1;
        if (fxl <= Z)
            break;
        end
        L = max(XLB,L-W);
    end
    while(R < XUB)
        xsim = R;
        theta = theta0+xsim*V1(:,it);
        fxr = -feval(objective_function,theta,varargin{:});
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
        theta = theta0+xsim*V1(:,it);
        fxsim = -feval(objective_function,theta,varargin{:});
        neval(it) = neval(it) + 1;
        if (xsim > xold)
            R = xsim;
        else
            L = xsim;
        end
    end
end

if ~isempty(sampler_options.mode),
    dist1=sqrt(sum((theta-mm(r).m).^2));
    if dist1>distance(r),
        theta=theta1;
        fxsim=[];
    end
end
end
