function [post] = vmp_xi_consensus_smfret(traceset,vbem,rates,k,kprime,maxIter,baseline,guess)
% say we're given some traces and some vbem outputs
% what's the best transition matrix?

    minIter = 1;
    intermediateIter = 50;
    tolerance = .001;
    converged = 0;
    evidence = 0;
    N = length(traceset);
    lengths = zeros(length(traceset),1);
    for rk = 1:length(traceset)
        lengths(rk,:) = length(traceset{rk});
    end
    % declarations

    logrates = log(rates);
    condition = ones(size(rates,1),1);
    for rk = 1:size(rates,2)
        condition = condition&rates(:,rk)>0&~isinf(rates(:,rk))&~isnan(rates(:,rk));
    end
    
    temp = logrates(condition,:);
    if k == 1
        p = exp(median(temp));
        coeff = ones(length(traceset),1);
        rule = 1;
        consensus_rates = p;
    else
        [coeff,p] = kmeans(temp,k);
        [~,rule] = sort(p(:,1));
        consensus_rates = exp(p(rule,:));
    end
    classes = length(unique(coeff));
    prior.mark.alpha = zeros(kprime,kprime,k);
    for rk = 1:classes
        temp = find(coeff==rule(rk));
        for rrk = 1:length(temp)
            try
                prior.mark.alpha(:,:,rk) = prior.mark.alpha(:,:,rk)+vbem{temp(rrk)}.mark.alpha;
            end
        end
        prior.mark.alpha(:,:,rk) = prior.mark.alpha(:,:,rk)/length(temp);
        prior.alpha(rk,:) = sum(length(temp));
    end
    % prior for the transition matrices - this part is absolutely crucial.
    % Uses the kmeans algorithm to make a guess for how many counts that
    % might be expected from the dataset.
    
    consensus_tm = rates2tm(consensus_rates,k,kprime);
    fulldata = cell2mat(traceset(:));
    % cool little function that moves between transition matrices and rates
    
    post.alpha = ones(kprime,kprime,k);
    post.mixture.alpha = ones(k,1)/k;
    post.mixture.mix = ones(k,1)/k;
    
    post.emissions.mix = ones(kprime,1)/kprime;
    post.emissions.alpha = ones(kprime,1)/kprime;
    post.b = ones(kprime,1)*var(fulldata)/kprime^2;
    post.a = ones(kprime,1)*10;
    post.m = guess(:,1);
    post.beta = 1e1*ones(kprime,1);
    post.ensemble.alpha = post.alpha(:,:,1);
    % model in the emissions and their priors
    
    temp = prior.alpha;
    temp2 = prior.mark;
    prior = post;
    prior.mixture.alpha = temp;
    prior.alpha = temp2.alpha/10;
    prior.emissions.alpha = length(fulldata)*post.emissions.mix;
    count = 1;
    
    vbem = standardize(post,vbem,traceset);
    if baseline
        for rk = 1:length(vbem)
            vbem{rk}.m = guess(:,rk);
        end
    end
    while ~converged
        disp(['Starting Iteration ' num2str(count)])
        post.record(count,:) = evidence;
        [phi,~,post.conprob,cgamma,~,evidence] = transmat_e_step(traceset,post,vbem,consensus_tm,k,kprime,evidence);

        % E step
        gamma = cell2mat(cgamma(:));
        for rk = 1:k
            temppost.alpha = post.alpha(:,:,rk);
            tempprior.alpha = prior.alpha(:,:,rk);
            [tempvar,tempevidence(rk,:)] = dirichlet_update(temppost,tempprior,phi(:,:,rk),0);
            post.alpha(:,:,rk) = tempvar.alpha;
            consensus_tm{rk} = tempvar.mix';
            % transition matrix mixture update
        end
        
        evidence = evidence+sum(tempevidence);
        [post.mixture,evidence] = dirichlet_update(post.mixture,prior.mixture,post.conprob,evidence);
        % mixture dirichlet update
        [post.emissions,evidence] = dirichlet_update(post.emissions,prior.emissions,gamma,evidence);
        
        if baseline
            for rrk = 1:length(traceset)
                for rk = 1:kprime
                    [vbem{rrk},evidence] = gauss_update(traceset{rrk},vbem{rrk},prior,cgamma{rrk}(:,rk)',rk,evidence);
                end
            end
        else
            for rk = 1:kprime
                [post,evidence] = gauss_update(fulldata,post,prior,gamma(:,rk)',rk,evidence);
            end
            % gaussian update

            vbem = standardize(post,vbem,traceset);
            % standardize the emissions
        end
        
        % M step
        if count==1
            post.record(count,:) = evidence;
        end
        if mod(count,4)==0
            if count < intermediateIter
                diagnostic = mean((diff(post.record(end-3:end)/abs(mean(post.record(end-3:end))))));
            else
                diagnostic = abs(mean((diff(post.record(end-3:end)/abs(mean(post.record(end-3:end)))))));
            end
            condition = diagnostic<tolerance&&diagnostic>-eps&&count>minIter;
            if condition||count>maxIter
                converged = 1;
            elseif isnan(evidence)
                converged = 1;
            end
        end
        count = count+1;
    end
    post.consensus_rates = tm2rates(consensus_tm);
    post.consensus_tm = consensus_tm;
    post.gamma = gamma;
    [~,epsc,~,~,post.vbem,~] = transmat_e_step(traceset,post,vbem,consensus_tm,k,kprime,evidence);       
    [post.ensemble,~] = dirichlet_update(post.ensemble,prior.ensemble,epsc,0);
    post.ensemble.tm{1} = post.ensemble.mix';
    post.ensemble.alpha = epsc'+1;
    post.ensemble.rates = tm2rates(post.ensemble.tm);
    for rk = 1:length(vbem) 
        temp = -log(1-post.vbem{rk}.mark.mix);
        post.vbem{rk}.rate = temp(~eye(kprime,kprime));
    end
    for rk = 1:length(vbem)
        vbem{rk} = post.vbem{rk};
    end
    
    parfor rk = 1:length(vbem)
        [~,viterbi] = vit(traceset{rk},vbem{rk},kprime);
        vbem{rk}.ideals = viterbi;
    end
    post.vbem = vbem;
    for rk = 1:k
        post.alpha(:,:,rk) = phi(:,:,rk)'+1;
        count = 1;
        for rrk = 1:size(post.alpha,1)
            temp = post.alpha(rrk,:,rk);
            temp = temp(:)';
            for rrrk = 1:length(temp)
                if rrk~=rrrk
                    try
                        post.rate_confidence_intervals{rk,count} = confidence_intervals(temp(rrrk),sum(temp)-temp(rrrk));
                        count = count+1;
                    catch
                        keyboard
                    end
                end
            end
        end
    end
    % transpose is for convention with the field
    % acquire confidence intervals
    
    post.evidence = evidence;
end

function [like, path] = vit(x,post,k)
    conprob = exp(bsxfun(@plus,-log(2*pi)+(psi(post.a')-...
        log(post.b'))-post.a'./(2*post.b').*(1./post.beta'+post.m'.^2),...
        bsxfun(@times,post.a'./(2*post.b'),(bsxfun(@plus,-x.^2,...
        2*bsxfun(@times,post.m',x))))));
    conprob = log(bsxfun(@times,1./sum(conprob,2),conprob));
% evaluate the conditionals
% second line introduces a normalization factor

    % everything that follows is jacked almost entirely from vbFRET
    % chmmviterbi.m by Jonathan Bronson
    
    % arbitrary value, since there is no predecessor to t=1
    bestPriorZ = zeros(length(x),k);
    z_hat = zeros(length(x),1);

    % Compute values for timesteps 2-end.
    % omega(zn)=ln(p(xn|zn))+max{ln(p(zn|zn-1))+omega(zn-)}
    % CB 13.68
    A = log(post.mark.mix);
    for t=2:length(x)
        [temp bestPriorZ(t,:)] = max(bsxfun(@plus,A,conprob(t-1,:)'));
        conprob(t,:) = conprob(t,:)+temp;
        % this is like, the only part I wrote lol
    end
    [logLikelihood z_hat(length(x),:)] = max(conprob(length(x),:));
    like = exp(logLikelihood);
    path(length(x),:) = post.m(z_hat(length(x),:));
    for t=(length(x)-1):-1:1
        z_hat(t,:) = bestPriorZ(t+1,z_hat(t+1),:);
        path(t,:) = post.m(z_hat(t,:));
    end
end

% Viterbi algorithm

function [vbem] = standardize(post,vbem,data)
    for rk = 1:length(vbem)
        vbem{rk}.a = post.a;
        vbem{rk}.b = post.b;
        vbem{rk}.m = post.m;
        vbem{rk}.mix = post.emissions.mix;
        vbem{rk}.beta = post.beta;
    end
end

function [pci] = confidence_intervals(a,b)
    top = @(x) betainc(max([0 min([1 x])]),a,b)-0.005;
    bottom = @(x) betainc(max([0 min([1 x])]),a,b)-0.995;
    pci(1,:) = -log(1-fzero(top,a/(a+b)));
    pci(2,:) = -log(1-fzero(bottom,a/(a+b)));
    % the marginal distribution of a dirichlet is a beta
end

function [update,lowerbound] = dirichlet_update(update,hyper,phi,lowerbound)
    if size(hyper.alpha,2) == 1
        message = sum(phi,1)';
        phi = bsxfun(@plus,hyper.alpha-1,message);
        update.alpha = phi+1;
        strength = 1;
    else
        update.alpha = phi+1+hyper.alpha;
%         hyper.strength = min(sum(xi.*(ones(size(xi,1),size(xi,2))-eye(size(xi,1),size(xi,2))),1));
% pre-programmed the message into the forward-backward algorithm
    end
% straight out of a Winn's thesis

    for rk = 1:size(update.alpha,2)
        update.mix(:,rk) = exp(psi(update.alpha(:,rk))-psi(sum(update.alpha(:,rk))));
    end
    update.mix = bsxfun(@times,update.mix,1./sum(update.mix,1));
    kl = @(alpha,beta) gammaln(sum(alpha))-gammaln(sum(beta))-sum(gammaln(alpha))+sum(gammaln(beta))+(alpha-beta)*(psi(alpha)-psi(sum(alpha)))';
    for rk = 1:size(phi,2)
    %         lowerbound = lowerbound-sum((update.alpha(:,rk)-sum(update.alpha(:,rk))).*(psi(update.alpha(:,rk))-psi(sum(update.alpha(:,rk)))));
        lowerbound = lowerbound-kl(update.alpha(:,rk)',hyper.alpha(:,rk)');
    end


end
%% Dirichlet node

function [conprob,xi,logscale,lowerbound] = obs_update(x,post,k,bound)
    conprob = bsxfun(@plus,-log(2*pi)+(psi(post.a')-...
        log(post.b'))-post.a'./(2*post.b').*(1./post.beta'+post.m'.^2),...
        bsxfun(@times,post.a'./(2*post.b'),(bsxfun(@plus,-x.^2,...
        2*bsxfun(@times,post.m',x)))));
    conprob = exp(bsxfun(@plus,conprob,-max(conprob,[],2)));
    conprob = bsxfun(@times,1./sum(conprob,2),conprob);
% evaluate the conditionals
% second line introduces a normalization factor

    forw = conprob;
    backw = conprob;
    scale = ones(length(x),1);
% initialize the probabilities
    
    forw(1,:) = forw(1,:).*post.mix';
    scale(1,:) = sum(forw(1,:));
    backw(end,:) = ones(1,k);    
% explicit boundary condition

    for rk = 2:length(x)
        forw(rk,:) = (forw(rk-1,:)*post.mark.mix').*conprob(rk,:);
        scale(rk,:) = sum(forw(rk,:));
        forw(rk,:) = forw(rk,:)./scale(rk,:);
% "parent" to "child" message
    end
% for a good discussion of the scale variable, see JWVDM
    
    for rrk = length(x)-1:1
        backw(rrk,:) = ((backw(rrk+1,:)).*conprob(rrk+1,:))*post.mark.mix./scale(rrk+1,:);
% "child" to "parent" message  
    end
    
%     backw = bsxfun(@times,backw,1./sum(backw,2));
% does this have to be normalized?
    conprob = forw.*backw;
    conprob = bsxfun(@times,conprob,1./sum(conprob,2));
% new conditionals are put together pointwise

% shouldn't actually have to renorm; do so to be careful about numerical
% precision
    
    xi = zeros(k,k);
    
    for rk = 1:k
        for rrk = 1:k
            xi(rk,rrk) = sum(1./scale(2:end,:).*backw(2:end,rk).*conprob(2:end,rk).*post.mark.mix(rk,rrk).*forw(1:end-1,rrk))/exp(prod(scale));
        end
    end
% transition pseudocounts, courtesy of Jan Willem van de Meent
% could consider ramming out one of the loops

    logscale = sum(log(scale));

% normalize
    lowerbound = zeros(k,1);
%     for rk = 1:k
%         lowerbound(rk,:) = sum(post.a(rk,:)./post.b(rk,:).*post.m(rk,:).*(x.*conprob(:,rk))-post.a(rk,:)./post.b(rk,:)*...
% 			(x.*conprob(:,rk)).^2/2+.5*(psi(post.a(rk,:))-log(post.b(rk,:))-post.a(rk,:)./...
% 			post.b(rk,:)*(post.m(rk,:)^2+1./post.beta(rk,:))+log(2*pi)));
%     end
    lowerbound = bound-sum(lowerbound);
% lowerbound calculation
end
%% observed node (white noise only)

function [phi,epsc,conprob2,conprob3,vbem,bound] = transmat_e_step(traceset,post,vbem,consensus_tm,k,N,bound)
    conprob = zeros(length(traceset),k);
    xi = zeros(N,N,length(vbem),k);
    gamma = cell(length(vbem),k);
    conprob3 = cell(length(vbem),1);
    for rk = 1:length(vbem)
        for rrk = 1:k
            gamma{rk,rrk} = zeros(length(traceset{rk}),N);
        end
        conprob3{rk} = zeros(length(traceset{rk}),N);
    end
    phi = zeros(N,N,k);
    parfor rk = 1:length(vbem)
        for rrk = 1:k
            vbem{rk}.mark.mix = consensus_tm{rrk}';
            [gamma{rk,rrk},xi(:,:,rk,rrk),conprob(rk,rrk),~] = obs_update(traceset{rk},vbem{rk},length(vbem{rk}.m),bound);
        end
    end
    temp = conprob;
    conprob = exp(bsxfun(@plus,conprob,-max(conprob,[],2)));
    conprob2 = bsxfun(@times,conprob,post.mixture.mix');
    conprob = bsxfun(@times,conprob,1./sum(conprob,2));
    conprob2 = bsxfun(@times,conprob2,1./sum(conprob2,2));
    bound = sum(sum(conprob2.*temp));
    % now we are in a position to calculate phi
    for rk = 1:k
        phi(:,:,rk) = post.alpha(:,:,rk)-1;
        for rrk = 1:length(vbem)
            phi(:,:,rk) = bsxfun(@plus,phi(:,:,rk),conprob2(rrk,rk)*xi(:,:,rrk,rk));
            conprob3{rrk} = conprob3{rrk}+conprob(rrk,rk)*gamma{rrk,rk};
        end
    end
    epsc = zeros(N,N);
    if k > 1
        instruction = dsp.ArrayVectorMultiplier;
        instruction.Dimension = 4;
        instruction2 = dsp.ArrayVectorMultiplier;
        instruction2.Dimension = 3;
        for rk = 1:length(vbem)
            psc = sum(step(instruction,xi(:,:,rk,:),conprob2(rk,:)),4);
            vbem{rk}.pseudocounts = xi(:,:,rk,:);
            epsc = epsc+psc;
            vbem{rk}.mark.alpha = sum(step(instruction2,post.alpha,conprob2(rk,:)./sum(conprob2,1)),3);
            [vbem{rk}.mark,~] = dirichlet_update(vbem{rk}.mark,vbem{rk}.mark,psc,0);
        end
    else
        for rk = 1:length(vbem)
            psc = xi(:,:,rk);
            vbem{rk}.pseudocounts = psc;
            epsc = epsc+psc;
            vbem{rk}.mark.alpha = post.alpha/length(traceset);
            [vbem{rk}.mark,~] = dirichlet_update(vbem{rk}.mark,vbem{rk}.mark,psc,0);
        end
    end
    % form the "expected" pseudocounts and the "expected" concentration
    % parameters for use as updates for individual vbem
    % technically this isn't a part of the e step
    
end
% e step

function [consensus_tm] = rates2tm(consensus_rates,k,N)
% relies heavily on the continuum limit
    trans = 1-exp(-consensus_rates');
    trans = trans(:);
    consensus_tm = cell(k,1);
    count = 1;
    for rk = 1:k
        for rrk = 1:N
            for rrrk = 1:N
                if rrk~=rrrk
                    consensus_tm{rk}(rrk,rrrk) = trans(count);
                    count = count+1;
                end
            end
            consensus_tm{rk}(rrk,rrk) = 1-sum(consensus_tm{rk}(rrk,:));
        end
    end
end
function [rates] = tm2rates(tm)
    j = size(tm{1},1);
    rates = zeros(length(tm),j^2-j);
    for rk = 1:length(tm)
        temp = -log(1-tm{rk}');
        rates(rk,:) = temp(~eye(j,j))';
    end
end

function [update, lowerbound] = gauss_update(x,update,hyper,conprob,index,bound)
    message_mu(1,:) = x(:)'.*hyper.a(index,:)./hyper.b(index,:);
    message_mu(2,:) = -ones(1,length(x)).*hyper.a(index,:)./(2*hyper.b(index,:));
    message_mu = bsxfun(@times,conprob,message_mu);
    phi_mu = [hyper.beta(index,:)*hyper.m(index,:);-hyper.beta(index,:)/2]+sum(message_mu,2);
    update.beta(index,:) = -2*phi_mu(2,:);
    update.m(index,:) = phi_mu(1,:)/update.beta(index,:);
%   mu update

    message_prec(1,:) = -1/2.*(x(:)'.^2-2.*x(:)'.*update.m(index,:)+update.m(index,:)^2+1/update.beta(index,:));
    message_prec(2,:) = ones(1,length(x))/2;
    message_prec = bsxfun(@times,conprob,message_prec);
    phi_prec = [-hyper.b(index,:);hyper.a(index,:)-1]+sum(message_prec,2);
    update.b(index,:) = -phi_prec(1,:);
    update.a(index,:) = phi_prec(2,:)+1;
%   prec update
%     kl_gamma = @(a1,b1,a2,b2) (b1(index)-1)*psi(b1(index))-log(a1(index))-...
%         b1(index)-gammaln(b1(index))+gammaln(b2(index))+b2(index)*log(a2(index))+...
%         (1-b2(index))*(psi(b1(index))+log(a1(index)))+b1(index)*a1(index)/a2(index);
%     kl_gauss = @(mu1,s1,mu2,s2) 1/2*log(s1(index)/s1(index))+((1./s1(index))+(mu1(index)-mu2(index)).^2)/(1./s2(index));
%     
%     lowerbound = bound-kl_gauss(update.m,update.beta,hyper.m,hyper.beta)-...
%         kl_gamma(update.a,update.b,hyper.a,hyper.b);
    lowerbound = -([update.beta(index,:).*update.m(index,:)-hyper.beta(index,:).*hyper.m(index,:),-update.beta(index,:)./...
        2+hyper.beta(index,:)./2]*[update.m(index,:);update.m(index,:).^2+1./update.beta(index,:)]+...
        1/2*(log(update.beta(index,:))-update.beta(index,:).*update.m(index,:).^2-...
        log(hyper.beta(index,:))+hyper.beta(index,:).*hyper.m(index,:).^2));
    lowerbound = bound-hyper.emissions.mix(index,:)*(lowerbound+[-update.b(index,:)+hyper.b(index,:),update.a(index,:)-hyper.a(index,:)]*...
        [update.a(index,:)./update.b(index,:);psi(update.a(index,:))-log(update.b(index,:))]+update.a(index,:).*log(update.b(index,:))-...
       gammaln(update.a(index,:))-hyper.a(index,:).*log(hyper.b(index,:))+gammaln(hyper.a(index,:)));
%     %   lowerbound calculation
end
%% Gaussian node
