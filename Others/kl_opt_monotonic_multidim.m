% This function solves the optimization problem
% minimize max_{||z||<1} D(P(x)||P(x-z))
% subject to E [c(||X||)] <= C
%            P(x) = f(||x||), where f(r) is monotically non-increasing in r
% where x,y,z are vectors of some dimension.
% Solved using an interior point second-order method.
% 
%
% The inputs are:
%  1. dim is the number of dimensions of the vectors in the problem
%  2. n is the quantization level. Bins are of length 1/n
%  3. xmax determines the size of the interval considered for the function
%  (from 0 to xmax)
%  4. C is the cost value
%  5. c_exp is the cost exponent. That is, the cost function is
%        c(r)=r^c_exp (Default is 2)
%  6. tol is the tolerance for the solver. (Default is 1e-8)
%  7. verbose controls how much information is displayed. 0 does nothing. 
%     1 prints some data. 2 also plots the candidate distribution. (Default is 2)
%
% The outputs are:
%  1. primobj is the primal objective value.
%  2. x is the vector of values at the centers of the quantization bins.
%  3. p is the solution to the optimization problem. This is actually the
%     PDF of the radius, not the PDF of the vector distribution. (Note that
%     this returns the PDF, not a PMF of the discretized distribution. So
%     it will sum to n rather than 1.)
%  4. pscaled is scaled so that it gives the true PDF of the vector
%     distribution
%  5. px1 is the PDF of (the absolute value of) a single element of the 
%     vector distribution
%  6. samplefnc is a function that will produce one random sample from the
%     vector distribution when called with no arguments.
%  7. gausobj is the KL divergence of the Gaussian distribution with the
%     same cost

function [primobj,x,p,pscaled,px1,samplefnc,gausobj] = kl_opt_monotonic_multidim(dim,n,xmax,C,c_exp,tol,verbose)

if nargin < 5
    c_exp = 2;
end
if nargin < 6
    tol = 1e-8;
end
if nargin < 7
    verbose = 2;
end

N = ceil(n*xmax);
x = (1:N)'/n;
l = length(x);

% The two constraints (probability normalization, and the variance
% constraint) are handled by A*p=b as follows. Note that even though the
% variance constraint is an inequality, it will always be tight at optimal,
% so we treat it as an equality constraint.

a1 = x'.^c_exp/n;
a2 = ones(1,l)/n;
A = [a1;a2];
b = [C;1];

% We start with a guess for the distribution p.
lam = C^(-1/c_exp)/gamma(1+c_exp);
p = exp(-lam*x);

% make sure that p is normalized (it actually doesn't need to satisfy the
% cost constraint exactly initially -- the algorithm will ensure that)
p = p/(A(2,:)*p);

% construct matrix for monotonicity constraint given by D'*p >= 0
nd = l-1;
D = sparse(l,nd);
for i=1:l-1
    D(i,i) = 1;
    D(i+1,i) = -(x(i)/x(i+1))^(dim-1);
end

% the objective function will be sum_{i,j} Q(i,j)*p(i)*log(p(i)/p(j))+c'*p
Q = zeros(l); % matrix with weights for objective function
for i=1:l
    for j=1:l
        r = x(i);
        q = x(j);
        if r+q >= 1 && abs(r-q) <= 1 && i~=j
            Q(i,j) = (q/r)*(((r+q)^2-1)*(1-(r-q)^2)/4/r^2)^((dim-3)/2);
        end
    end
end
% scale Q so that it actually gives the correct value of the KL divergence
Q = Q*(exp(gammaln(dim/2)-gammaln((dim-1)/2))/sqrt(pi)/n^2);
c = -(dim-1)*sum(Q.*log(x./x'),2);
sumQ = sum(Q,2);

% calculate KL divergence for the Gaussian with the same cost
% this is the variance of the Gaussian
sg2 = (C*exp(gammaln(dim/2)-gammaln((dim+c_exp)/2)))^(2/c_exp)/2;
gausobj = 1/2/sg2;

iter = 0; % iteration count
t = 1; % this is the weighting parameter

% boolean to keep track of feasibility with respect to the linear equality constraints
isfeasible = false;

primobj = getobj(p,Q,c,sumQ);

while 1
    iter = iter+1;

    dvals = D'*p;
    fval = primobj-sum(log(dvals))/t;
    
    grad = getgrad(p,Q,c,sumQ); % gradient
    grad = grad-D*(1./dvals/t);
    
    if isfeasible
        rpri = zeros(2,1);
    else
        rpri = A*p-b; % primal feasibility residual
    end
    
    % The next section of code is typically the majority of the run-time.
    
    % Compute Hessian matrix
    H = getH(p,Q,sumQ);
    fullH = H+D*spdiags(1./dvals.^2/t,0,nd,nd)*D';
    
    % Now solve for the Newton search direction
    [R,psd_flag] = chol(fullH);
    if psd_flag == 0 % fullH really should be PSD, but sometimes Cholesky fails
        gtilde = R'\grad;
        AR = A/R;
        
        [U,S,V] = svd(AR,'econ');
        s = diag(S);
        vtilde = -gtilde-V*(1./s.*(U'*rpri))+V*(V'*gtilde);
        
        v = R\vtilde; % search direction
    end
    if psd_flag > 0 || (isfeasible && grad'*v > 0)
        if verbose >= 1
            disp('Cholesky failed');
        end
        % in this case, solve for the search direction the old fashioned way
        vw = [fullH A';A zeros(2)]\[-grad;-rpri];
        v = vw(1:l);
    end
    
    nwt_dec = -grad'*v/2; % Newton decrement: estimate of gap to optimality in the soft-max problem
    if isfeasible && nwt_dec < 0
        disp('Newton decrement was negative: doing simpler thing');
        v = -grad.*p.^2;
        v = v-pinv(A)*(A*v);
        nwt_dec = -grad'*v/2;
        if nwt_dec < 0
            % if we're still somehow going the wrong way, just flip the direction
            v = -v;
        end
        while any(p+v<0)
            v = v/2;
        end
    end
    
    % implementing a line search. We move to p+dst*v
    dst = 1;
    posdst = 0;
    while 1
        pnew = p+dst*v;
        newdvals = D'*pnew;
        if all(pnew > 0) && all(newdvals > 0) % the p vector better be all positive and monotonic
            posdst = max(posdst,dst);
            primobj = getobj(pnew,Q,c,sumQ);
            
            newfval = primobj-sum(log(newdvals))/t;
            
            if ~isfeasible % if we weren't feasible before, get as close as possible
                break;
            end
            if newfval < fval+.1*dst*(grad'*v)
                break;
            end
            if dst < 1e-16
                break;
            end
        end
        dst = dst/2;
    end

    if ~isfeasible && dst==1
        isfeasible = true;
        if verbose >= 1
            disp('Feasible!');
        end
    end

    % print some useful data
    if verbose >= 1
        fprintf('iter=%i  dst=%1.1e  primobj=%1.4f  nwt_dec=%1.2e  t=%1.1e\n',iter,dst,primobj,nwt_dec,t);
    end
    
    if verbose >= 2
        semilogy(x,p); % plot the candidate distribution
        grid on
        drawnow
    end
    
    p = pnew;

    if isfeasible
        % l/t is an estimate for the duality gap, if we are at exact
        % optimality for the soft-max problem. Thus nwt_dec+l/t is an
        % estimate for the total duality gap
        if (nwt_dec+l/t)/primobj < tol && nwt_dec > 0
            break;
        end
        % if we take a full Newton step or get extremely close to optimal, then increase t
        if dst == 1 || (nwt_dec/primobj < tol/2 && nwt_dec > 0) %|| dst < 1e-8
            % This is a fairly modest increase in t. This seems to make the
            % solver the most robust, if not necessarily the most efficient.
            t = t*1.25;
       end
    end

end

pscaled = exp(gammaln(dim/2)-log(2)-(dim/2)*log(pi)-(dim-1)*log(x)).*p;

% get PDF of single variable
px1 = zeros(l,1);
for i=1:l
    px1(i) = sum(p.*(max(0,x.^2-x(i)^2).^((dim-3)/2))./x.^(dim-1));
end
px1 = px1/sum(px1)*n;

% return a sampling function
CDF = cumsum(p/n);
samplefnc = @() sample(x,CDF,n,dim);

end


function obj = getobj(p,Q,c,sumQ)
% calculate the objective value

obj = sumQ'*(p.*log(p))+p'*(c-(Q*log(p)));

end

function grad = getgrad(p,Q,c,sumQ)
% Calculate the gradient for a shift of k

grad = (log(p)+1).*sumQ-Q*log(p)-(Q'*p)./p+c;


end

function H = getH(p,Q,sumQ)
% Calculate the Hessian matrix

l = length(p);
d = sumQ./p+(Q'*p)./(p.^2);

H = -Q./p';
H = (H+H')+diag(d);

end


function s = sample(x,CDF,n,dim)
% produce one sample from the distribution

i = find(rand<CDF,1);
r = x(i)-rand/n;
z = randn(dim,1);
s = z*(r/norm(z));

end