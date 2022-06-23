% trying to solve the optimization problem
% minimize max_{|a|<1} D(P(x)||P(x-a))
% subject to E X^2 <= D
% Solved using an interior point second-order method

% It actually solves the finite dimensional approximation:
% minimize max_{|a|\le n} \sum_{i=-xmax*n}^{(xmax-1)*n} p(i)*log(p(i)/p(i+a))
% subject to \sum_i p(i) = 1
%            \sum_i (i*n)^2*p(i) <= D
% In the above, n is the parameter that controls the quantization, and xmax
% determines the range of the distribution that is computed. 

% It actually actually solves a modification of the above optimization
% problem to make the maximum into a differentiable function. Specifically,
% the objective function is replaced with
% \exp(t*\sum_{i=-xmax*n}^{(xmax-1)*n} p(i)*log(p(i)/p(i+a)))
% This remains a convex function, and for large t, it gives the same
% solution as the above. The basic technique is an equality-constrained (to 
% deal with the two constraints) Newton minimization, and then when close 
% to optimal the t value is increased

std = .35; % standard deviations
D = std^2; % variance
n = round(50/std); % quantization level
xmax = round(20*std);

x = ((-xmax):1/n:xmax)';
l = length(x);

% The two constraints (probability normalization, and the variance
% constraint) are handled by A*p=b as follows. Note that even though the
% variance constraint is an inequality, it will always be tight at optimal,
% so we treat it as an equality constraint.
A = [x'.^2;ones(1,l)];
b = [D;1];

% we start with a guess, based on the Laplace distribution
p = exp(-abs(x)/sqrt(D/2));

% make sure that p satisfies the two consraints
p = p/sum(p);
A2 = A*A';
p = p-A'*(A2\(A*p))+A'*(A2\b);

% this is the value of the objective for a Gaussian mechanism -- it has a 
% nice closed form!
objg = 1/D; 

iter = 0;
t = 3; % this is the weighting thing

while 1
    iter = iter+1;
    allobj = getallobj(p,n);
    maxobj = max(allobj);
    evals = exp(t*(allobj-maxobj));
    fval = sum(evals)/t;
    
    allgrads = zeros(l,n);
    H = sparse(l,l); % Hessian matrix
    for a=1:n
        allgrads(:,a) = getgrad(p,a);
        if evals(a) > 1e-12
            H = H+evals(a)*getH(p,a);
        end
    end
    grad = allgrads*evals; % gradient
    
    fullH = full(H);
    fullH = fullH+allgrads*diag(t*evals)*allgrads';
    mymat = [fullH A';A zeros(2)];
    
    vw = mymat\[-grad;zeros(2,1)];
    v = vw(1:l); % search direction
    
    % implementing a simple line search. we move to p+r*v
    r = 1;
    while 1
        pnew = p+r*v;
        if min(pnew) < 0 % the p vector better be all positive
            r = r/2;
        else
            newallobj = getallobj(pnew,n);
            newfval = sum(exp(t*(newallobj-maxobj)))/t;
            if newfval < fval+.25*r*(grad'*v)
                break;
            else
                r = r/2;
                if r == 0
                    error('r=0 reached');
                end
            end
        end
    end
    
    p = pnew;
    [maxobj,amax] = max(newallobj);
    nm = abs(grad'*v); % measure of closeness to optimality
    % print some useful data
    fprintf('iter=%i  r=%1.1e  maxobj=%f  amax=%i  nm=%1.2e  t=%1.1e\n',iter,r,maxobj,amax,nm,t);
    semilogy(x,p); % plot the distribution
    title(['std=' num2str(std)]);
    drawnow
    if t < 1e8
        if nm < 1e-5 % if fairly close to optimality, increase t
            t = t*3;
        end
    else
        if nm < 1e-8
            save('optimal_p.mat','p','x','std','r','n');
            break;
        end
    end
end

return

function obj = getobj(p,n)
% Calculate objective value for a shift of n.
% Note that this really computes the sum of two KL divergences; i.e.
% D(p(x)||p(x+n))+D(p(x)||p(x-n)). This just makes things nice and
% symmetric (although it means there's a factor of two somewhere)
l = length(p);
obj = sum((p(1+n:l)-p(1:l-n)).*log(p(1+n:l)./p(1:l-n)));
end

function allobj = getallobj(p,n)
% calculate all objective values for different shifts
allobj = zeros(n,1);
for a=1:n
    allobj(a) = getobj(p,a);
end
end

function grad = getgrad(p,n)
% Calculate gradient for a shift of n
l = length(p);
rat = p(1+n:l)./p(1:l-n);
rati = p(1:l-n)./p(1+n:l);
grad = [zeros(n,1);log(rat)+1-rati]+[log(rati)+1-rat;zeros(n,1)];
end

function H = getH(p,n)
% Calculate Hessian matrix for a shift of n.
l = length(p);

d = [1./p(1:l-n)+p(1+n:l)./p(1:l-n).^2;zeros(n,1)]+[zeros(n,1);1./p(1+n:l)+p(1:l-n)./p(1+n:l).^2];
a = -1./p(1:l-n)-1./p(1+n:l);
H = spdiags([[a;zeros(n,1)] d [zeros(n,1);a]],[-n,0,n],l,l);

end