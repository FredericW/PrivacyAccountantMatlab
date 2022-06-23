c_exp = 1; % the cost function is abs(x)^c_exp
stdmax = 1.5;
stdmin = 0.1;
steps = 100;
stds = linspace(stdmin,stdmax,steps);

figure

primobj = zeros(size(stds));

for i=1:length(stds)

    mystd = stds(i);
    disp(mystd);
    
    n = min(500,round(50/mystd));
    xmax = max(2,mystd*12);
    C = mystd^c_exp;
    opt_type = 1;
    r = 0.9;
    tol = 1e-8;
    verbose = 1;
    
    primobj(i) = kl_opt(n,xmax,C,c_exp,opt_type,r,tol,verbose)/2;
end

filename = sprintf('comparison_L%d_cactus.csv',c_exp);
csvwrite(filename,primobj);
plot(stds,primobj,"-",'LineWidth',1.5,'DisplayName','Cactus')
hold on

if c_exp == 2
    gaus_obj = 1./(2*stds.^2);
    filename = sprintf('comparison_L2_gaus.csv');
    csvwrite(filename,gaus_obj);
    plot(stds,gaus_obj,"--",'LineWidth',1.5,'DisplayName','Gaussian')
    hold on
end

if c_exp == 1
    gaus_obj = 1./(pi*stds.^2);

    laplace_obj = -1+exp(-1./stds)+1./stds;
    filename = sprintf('comparison_L1_laplace.csv');
    csvwrite(filename,laplace_obj);
    plot(stds,laplace_obj,"--",'LineWidth',1.5,'DisplayName','Laplace')
    hold on

    x0 = -fzero(@(x)airy(1,x),-1);
    logpairy = chebfun(@(x)2*(log(airy(2*x0/(3)*abs(x)-x0)))-log(3*airy(-x0)^2),[-20,20]);    
    airy_obj = zeros(size(stds));
    for i=1:length(stds)
        if 1/stds(i) > 40
            airy_obj(i) = inf;
        else
            airy_obj(i) = sum(chebfun(@(x) exp(logpairy(x)).*(logpairy(x)-logpairy(x+1/stds(i))),[-20,20-1/stds(i)]));
        end
    end
    filename = sprintf('comparison_L1_airy.csv');
    csvwrite(filename,airy_obj);
    plot(stds,airy_obj,"-.",'LineWidth',1.5,'DisplayName','Airy')
    hold on
end

filename = sprintf('comparison_L%d_stds.csv',c_exp);
csvwrite(filename,stds);

axis([0 1.5 0 4]);
xlabel('Cost C');
ylabel('KL divergence');
grid on
legend
