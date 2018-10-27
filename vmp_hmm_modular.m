function [out,rec] = vmp_hmm_modular(data,k,initguess)
% vmp, compiled by Jason Hon; based on variational message passing by
% John Winn; vbFRET, by Jon Bronson; math pulled from ebFRET by Jan-Willem
% van de Meent/John Winn/Matthew Beal
% Inputs: 
% 1: data
% 2: # of expected emissions, which are bound to transitions
% 3: initial guess for the centers (if not provided, code will create)
    restarts = 2;
    out.bound = -Inf;
    out.likelihood = -Inf;
    sk=1;
    s = length(data);
    while sk<restarts+1
        converged = 0;
        count = 1;
        tolerance = 0.01;
        check = 6;
        maxIter = 1e1;
        % loop metrics

        rec.bound = zeros(maxIter+1,1);
        rec.trans = zeros(maxIter+1,1);

        [post,~,~] = initialize(data,k,initguess);
        prior = post;
        lifetime = -log(post.mark.mix);
        lifetime = lifetime(~~eye(k,k));
        lifetime = 1-lifetime./sum(lifetime);
        bound = length(data)/1e2;
        prior.mark.alpha = getprior(bound,bsxfun(@times,bound.*post.mark.mix,lifetime),post.mark.mix);
        % initialize

        while ~converged
            if count<5
                post.m = prior.m;
            end
            lowerbound = 0;

            [conprob,xi,lowerbound] = obs_update(data,post,k,lowerbound);
            % forward-backward algorithm


            for rk = 1:k
                [post,lowerbound] = gauss_update(data,post,prior,conprob(:,rk)',rk,lowerbound);
                % mixture parameter updates
            end
            [post,lowerbound] = dirichlet_update(post,prior,conprob,xi,lowerbound,s);
            % emissions update
            
            [post.mark,~] = dirichlet_update(post.mark,prior.mark,conprob,xi,lowerbound,s);
            % transitions update
            rec.bound(count,:) = lowerbound;
            rec.trans(count,:) = sum(sum(xi.*(ones(k,k)-eye(k,k))));
            % one diagnostic - how do the parameters change? other diagnostic -
            % how do the interclass transitions change?

            if mod(count,check)==0
               diagnostic = k*abs(mean((diff(rec.bound(count-5:count,:))))/length(data));
               diagnostic2 = abs(rec.trans(count-(check-1),:)-rec.trans(count,:))./rec.trans(count-(check-1),:);
               if (abs(diagnostic)+diagnostic2)/2<tolerance&&diagnostic>0||count>maxIter||diagnostic2<.001
                   converged = 1;
                   sk = sk+1;
               elseif sum(isnan(post.mark.mix(~eye(k,k))))>0
                   converged = 1;
                   sk = sk+1;
               end
            end
            % success if the last 5 iterations are within tolerance

            if isnan(rec.bound(count,:))||isinf(rec.bound(count,:))
                converged = 1;
                sk = sk+1;
            end
            % run screaming if NaNs or Infs appear

            count = count+1;
        end
        rec.bound = rec.bound(1:count-1,:);
        rec.trans = rec.trans(1:count-1,:);
        post.ideals = zeros(length(data),1);
        [~,post.ideals] = vit(data,post,k);
        [~,xi,~] = obs_update(data,post,k,0);
        post.mark.alpha = xi'+1;
        post.rate = -1./log(1-post.mark.mix');
        post.bound = rec.bound(end);
        if post.bound>out.bound&&~isinf(post.bound)
           out = post;
           out.bound = post.bound(end);
        end
    end
%     clear rec;
%     rec = 1;
end
% main body

function [hyper,conprob,xi] = initialize(x,numstat,initguess)
%     hyper.strength = length(x)/numstat;
    hyper.strength = max([1 length(x)/1e4]);
    hyper.mark.strength = hyper.strength;
% add a memory lag to the dirichlet posterior
    hyper.b = ones(numstat,1)*var(x)/numstat^2;
    hyper.a = ones(numstat,1)*10;
    hyper.beta = 1*ones(numstat,1);
% emissions are a gmm
    hyper.alpha = eps*ones(numstat,1);
    hyper.mark.alpha = eps*ones(numstat,numstat);
    hyper.mix = 1/numstat*ones(numstat,1);
    hyper.mark.mix = 1/numstat*ones(numstat,numstat);
    conprob = zeros(length(x),numstat);
% hmm/mixture parameters
    if isempty(initguess)
        x = x-tds(hyper.rsqr,hyper.rsqr*(x-mean(x)));
        thresh = multithresh(x,numstat-1);
        thresh = thresh(:);
        thresh = [min(x)-1;thresh;max(x)+1];
        for rk = 2:length(thresh)
            hyper.m(rk-1,:) = mean(x(x>thresh(rk-1)&x<thresh(rk)));
            hyper.a(rk-1,:) = 10./std(x(x>thresh(rk-1)&x<thresh(rk)))^2;
        end
    else
        initguess = initguess(:);
        eta = bsxfun(@times, hyper.a'./(2*hyper.b'), ((bsxfun(@plus, -initguess', x))).^2);
        for rk = 1:numstat
            conprob(:,rk) = hyper.mix(rk,:)*sqrt(hyper.a(rk,:)./(2*hyper.b(rk,:))).*...
                exp(-eta(:,rk))./sum(bsxfun(@times,sqrt(hyper.a'./(2*hyper.b')).*hyper.mix',exp(-eta)),2);
        end
        hyper.m = initguess(1:end);
    end
    hyper.m = sort(hyper.m+2*sqrt(hyper.b./hyper.a).*randn(numstat,1));
    eta = bsxfun(@times, hyper.a'./(2*hyper.b'), ((bsxfun(@plus, -hyper.m', x))).^2);
    for rk = 1:numstat
        conprob(:,rk) = hyper.mix(rk,:)*sqrt(hyper.a(rk,:)./(2*hyper.b(rk,:))).*...
            exp(-eta(:,rk))./sum(bsxfun(@times,sqrt(hyper.a'./(2*hyper.b')).*hyper.mix',exp(-eta)),2);
    end
    xi = zeros(numstat,numstat);
    for rk = 1:numstat
        for rrk = 1:numstat
           xi(rk,rrk) = sum(conprob(1:end-1,rk).*conprob(2:end,rrk));
        end
    end
    lowerbound = 0;
    [hyper,~] = dirichlet_update(hyper,hyper,conprob,xi,lowerbound,length(x));
    [hyper.mark,~] = dirichlet_update(hyper.mark,hyper.mark,conprob,xi,lowerbound,length(x));
    
    generator = 10+1e1*rand(numstat,numstat);
    generator = 1./generator;
    for rk = 1:numstat
        t = zeros(numstat,1);
        t(rk) = 1;
        generator(rk,rk) = -sum(generator(rk,~t));
    end
    hyper.mark.mix = expm(generator)';
    hyper.mark.mix = bsxfun(@times,hyper.mark.mix,1./sum(hyper.mark.mix,2));
end

function [conprob,xi,logscale] = obs_update(x,post,k,bound)
    conprob = bsxfun(@plus,-log(2*pi)+(psi(post.a')-...
        log(post.b'))-post.a'./(2*post.b').*(1./post.beta'+post.m'.^2),...
        bsxfun(@times,post.a'./(2*post.b'),(bsxfun(@plus,-x.^2,...
        2*bsxfun(@times,post.m',x)))));
    conprob = exp(bsxfun(@plus,conprob,-max(conprob,[],2)));
    conprob = bsxfun(@times,conprob,1./sum(conprob,2));


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
    
    for rrk = length(x)-1:-1:1
        backw(rrk+1,:) = backw(rrk+1,:)/sum(backw(rrk+1,:));
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
    
    xi = zeros(size(post.mark.mix,1),size(post.mark.mix,2));
    for rk = 1:size(post.mark.mix,1)
        for rrk = 1:size(post.mark.mix,2)
            xi(rk,rrk) = sum(1./scale(2:end,:).*backw(2:end,rk).*conprob(2:end,rk).*post.mark.mix(rk,rrk).*forw(1:end-1,rrk))/exp(prod(scale));
        end
    end
% transition pseudocounts, courtesy of Jan Willem van de Meent
% could consider ramming out one of the loops
    logscale = sum(log(scale));
% normalize
%     lowerbound = zeros(k,1);
%     for rk = 1:k
%         lowerbound(rk,:) = sum(post.a(rk,:)./post.b(rk,:).*post.m(rk,:).*(x.*conprob(:,rk))-post.a(rk,:)./post.b(rk,:)*...
% 			(x.*conprob(:,rk)).^2/2+.5*(psi(post.a(rk,:))-log(post.b(rk,:))-post.a(rk,:)./...
% 			post.b(rk,:)*(post.m(rk,:)^2+1./post.beta(rk,:))+log(2*pi)));
%     end
%     lowerbound = bound-sum(lowerbound);
% lowerbound calculation
end
%% observed node (white noise only)

function [update, lowerbound] = gauss_update(x,update,hyper,conprob,index,bound)
%     update = hyper;
    message_mu(1,:) = x(:)'.*update.a(index,:)./update.b(index,:);
    message_mu(2,:) = -ones(1,length(x)).*update.a(index,:)./(2*update.b(index,:));
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
    kl_gamma = @(a1,b1,a2,b2) (b1(index)-1)*psi(b1(index))-log(a1(index))-...
        b1(index)-gammaln(b1(index))+gammaln(b2(index))+b2(index)*log(a2(index))+...
        (1-b2(index))*(psi(b1(index))+log(a1(index)))+b1(index)*a1(index)/a2(index);
    kl_gauss = @(mu1,s1,mu2,s2) 1/2*log(s1(index)/s1(index))+((1./s1(index))+(mu1(index)-mu2(index)).^2)/(1./s2(index));
    
    lowerbound = bound-kl_gauss(update.m,update.beta,hyper.m,hyper.beta)-...
        kl_gamma(update.a,update.b,hyper.a,hyper.b);
   %   lowerbound calculation
end
%% Gaussian node

function [update, lowerbound] = dirichlet_update(update,hyper,conprob,xi,bound,s)
    if size(hyper.alpha,2) == 1
        message = sum(conprob,1)';
        hyper.strength = sum(message);
    else
        message = xi;
%         hyper.strength = min(sum(xi.*(ones(size(xi,1),size(xi,2))-eye(size(xi,1),size(xi,2))),1));
% pre-programmed the message into the forward-backward algorithm
    end
    phi = bsxfun(@plus,hyper.alpha-1,message);
    update.alpha = phi+1;
% straight out of a Winn's thesis

    for rk = 1:size(update.alpha,2)
        update.mix(:,rk) = exp(psi(update.alpha(:,rk))-psi(sum(update.alpha(:,rk))));
    end
    update.mix = bsxfun(@times,update.mix,1./sum(update.mix,1));

% the transition matrix is put on its side to acquiesce to the dimensions
% of the mixture components
    lowerbound = bound;
    
    kl = @(alpha,beta) gammaln(sum(alpha))-gammaln(sum(beta))-sum(gammaln(alpha))+sum(gammaln(beta))+(alpha-beta)*(psi(alpha)-psi(sum(alpha)))';
    for rk = 1:size(phi,2)
%         lowerbound = lowerbound-sum((update.alpha(:,rk)-sum(update.alpha(:,rk))).*(psi(update.alpha(:,rk))-psi(sum(update.alpha(:,rk)))));
        lowerbound = lowerbound-kl(update.alpha(:,rk)',hyper.alpha(:,rk)');
    end
   

end
%% Dirichlet node

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


function out = getprior(bound,guess,tm)
% Based on Tom Minka's technical paper
    tolerance = 1;
    converged = 0;
    out = guess;
    newout = guess;
    while ~converged
        for rk = 1:size(tm,2)
            newout(:,rk) = invpsi(psi(sum(out(:,rk)))+log(tm(:,rk)));
        end
        newout = newout.*bound./sum(sum(newout));
        if sum(sum(abs(newout-out)))<tolerance
            converged = 1;
        end
        out = newout;
    end
end

function Y=invpsi(X)
% Y = INVPSI(X)
%
% Inverse digamma (psi) function.  The digamma function is the
% derivative of the log gamma function.  This calculates the value
% Y > 0 for a value X such that digamma(Y) = X.
%
% This algorithm is from Paul Fackler:
% http://www4.ncsu.edu/~pfackler/
%
    L = 1;
    Y = exp(X);
    while L > 10e-8
        Y = Y + L*sign(X-psi(Y));
        L = L / 2;
    end
end  

