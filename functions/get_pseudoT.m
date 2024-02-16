function pseudoT = get_pseudoT(C,Ca,Cc,Lead_fields)
% Calculate pseudo T statistic
mu = 0.05;
Cr = C + mu*max(svd(C))*eye(size(C));
Cr_inv = inv(Cr);
pseudoT = zeros(1,size(Lead_fields,3));
infolength = 0;
for n = 1:size(Lead_fields,3)
    this_L = Lead_fields(:,:,n);
    W_v = inv((this_L'*Cr_inv*this_L))*(this_L'*Cr_inv);
    iPower_v = this_L'*Cr_inv*this_L;
    [v,d] = svd(iPower_v);
    [~,id] = min(diag(d));
    lopt = this_L*v(:,id); % turn to nAm amplitude
    w = (lopt'*Cr_inv/(lopt'*Cr_inv*lopt))';
    pseudoT(:,n) = (w'*Ca*w - w'*Cc*w)./(0.5.*(w'*Ca*w + w'*Cc*w));
    
    fprintf(repmat('\b',1,infolength));
    infolength = fprintf('Evaluating pseudo-T (%d/%d)\n',n,size(Lead_fields,3));
end
end