function y = ma_ComputeRDP(ps,qs,orders)

logmgf = @(t) log(sum(qs.*((ps./qs).^t))); %log-mgf with a fudge value of 1e-30
rdp = arrayfun(logmgf,orders);
y = (rdp)./(orders-1);

end




