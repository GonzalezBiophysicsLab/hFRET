function best = vmp_hhmm(data,options)
    maxIter = options.maxIter;
    data = data(:);
    N = length(data);
    fulldata = cell2mat(data);
    tolerance = 1e-3;
    best.record.evidence = -Inf;
    if ~isfield(options,'restarts')
        options.restarts = 1;
    end
    if ~isfield(options,'baseline')
        options.baseline = 1==0;
    end
    for nk = 1:options.restarts
        [post] = initialize(data,N,options);
        prior = post;
        converged = 0;
        count = 1;
        clear record
        acc.gamma_all = zeros(length(fulldata),options.production);
        while ~converged
            disp(['Starting Iteration ' num2str(count)]);
            tempev = zeros(length(data),1);
            fastpost = post;
            try
                fastpost = rmfield(fastpost,'summaries');
            end
            parfor index = 1:length(data)
                temp = data{index};
                if options.baseline
                    conprob = gauss_prob_evaluate(temp,fastpost.emissions(index));
                else
                    conprob = gauss_prob_evaluate(temp,fastpost.emissions);
                end
                [Ts(index),tempev(index)] = fb_activation(fastpost,conprob,length(temp));
            end
            post.summaries = Ts;
            evidence = sum(tempev);
            for index = 1:length(data)
                for rk = 1:post.depth
                    for rrk = 1:post.path(rk,end)
                        if index == 1                        
                            acc.param_depth(rk).counts(rrk).a = Ts(index).param_depth(rk).g(rrk).a;
                            acc.param_depth(rk).counts(rrk).exit = Ts(index).param_depth(rk).g(rrk).exit;
                            acc.param_depth(rk).counts(rrk).pi = Ts(index).param_depth(rk).g(rrk).pi;
                        else
                            acc.param_depth(rk).counts(rrk).a = acc.param_depth(rk).counts(rrk).a+Ts(index).param_depth(rk).g(rrk).a;
                            acc.param_depth(rk).counts(rrk).exit = acc.param_depth(rk).counts(rrk).exit+Ts(index).param_depth(rk).g(rrk).exit;
                            acc.param_depth(rk).counts(rrk).pi = acc.param_depth(rk).counts(rrk).pi+Ts(index).param_depth(rk).g(rrk).pi;
                        end
                    end
                end
                % collect counts
                if index == 1
                    slice = 1:length(data{index});
                else
                    slice = (1:length(data{index}))+slice(end);
                end
                acc.g(index).gamma = zeros(length(data{index}),options.production);
                for rrk = 1:post.path(post.depth,end)
                    acc.g(index).gamma = acc.g(index).gamma+Ts(index).g(rrk).gamma;
                end
                if ~options.baseline
                    acc.gamma_all(slice,:) = acc.g(index).gamma;
                end
                % collect production occupancies
            end
            % run the forward-backward-activation algorithm on all time series
            for rk = 1:post.depth
                for rrk = 1:post.path(rk,end)
                    [post.param_depth(rk).hhmm_posterior(rrk).tm_exit,evidence] = dirichlet_update(post.param_depth(rk).hhmm_posterior(rrk).tm_exit,prior.param_depth(rk).hhmm_posterior(rrk).tm_exit,[acc.param_depth(rk).counts(rrk).a;acc.param_depth(rk).counts(rrk).exit],evidence);
                    [post.param_depth(rk).hhmm_posterior(rrk).pi,evidence] = dirichlet_update(post.param_depth(rk).hhmm_posterior(rrk).pi,prior.param_depth(rk).hhmm_posterior(rrk).pi,acc.param_depth(rk).counts(rrk).pi(:),evidence);
                    post.param_depth(rk).hhmm_params(rrk).tm = post.param_depth(rk).hhmm_posterior(rrk).tm_exit.mix(1:size(prior.param_depth(rk).hhmm_params(rrk).tm,1),:);
                    post.param_depth(rk).hhmm_params(rrk).exit = post.param_depth(rk).hhmm_posterior(rrk).tm_exit.mix(end,:);
                    post.param_depth(rk).hhmm_params(rrk).pi = post.param_depth(rk).hhmm_posterior(rrk).pi.mix(:)';
                end
            end
            % update markov chain parameters
            if options.baseline
                for index = 1:length(data)
                    for rk = 1:options.production
                        [post.emissions(index),evidence] = gauss_update(data{index},post.emissions(index),prior.emissions(index),acc.g(index).gamma(:,rk)',rk,evidence);
                    end
                    [post.emissions(index).mixture,~] = dirichlet_update(post.emissions(index).mixture,prior.emissions(index).mixture,acc.g(index).gamma,evidence);
                end
            else
                for rk = 1:options.production
                    [post.emissions,evidence] = gauss_update(fulldata,post.emissions,prior.emissions,acc.gamma_all(:,rk)',rk,evidence);
                end
                [post.emissions.mixture,~] = dirichlet_update(post.emissions.mixture,prior.emissions.mixture,acc.gamma_all,evidence);
            end
            % finite difference inversion; gradient gaussian used for
            % lowerbound calculation
            
            % update gaussian emissions

            for rk = 1:post.depth
                for rrk = 1:post.path(rk,end)
                    record.param_depth(rk).hhmm_params(rrk).exit(count,:) = post.param_depth(rk).hhmm_params(rrk).exit;
                    record.param_depth(rk).hhmm_params(rrk).pi(count,:) = post.param_depth(rk).hhmm_params(rrk).pi;
                    record.param_depth(rk).hhmm_params(rrk).tm(:,:,count) = post.param_depth(rk).hhmm_params(rrk).tm;
                end
            end
            record.evidence(:,count) = evidence;

            count = count+1;

            if mod(count,5)==0
                diagnostic = abs(diff(record.evidence));
                diagnostic = diagnostic(end)./abs(record.evidence(end));
                if diagnostic<tolerance||count>maxIter
                    post.record = record;
                    post.prior = prior;
                    if best.record.evidence(end)<post.record.evidence(end)||nk==1
                        best = post;
                    end
                    converged = 1;
                end
            end
        end
    end
    clear post
    fastpost = rmfield(best,'summaries');
    Ts = best.summaries;
    pathstate = cell(length(data),1);
    ideals = cell(length(data),1);
    parfor index = 1:length(data)
        temp = data{index};
        if options.baseline
            out = viterbi_activation(fastpost,Ts(index).g,length(temp),index);
        else
            out = viterbi_activation(fastpost,Ts(index).g,length(temp),1);
        end
        pathstate{index} = out.pathstate;
        ideals{index} = out.ideals;
    end
    best.pathstate = pathstate;
    best.ideals = ideals;
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
    lowerbound = -([update.beta(index,:).*update.m(index,:)-hyper.beta(index,:).*hyper.m(index,:),-update.beta(index,:)./...
        2+hyper.beta(index,:)./2]*[update.m(index,:);update.m(index,:).^2+1./update.beta(index,:)]+...
        1/2*(log(update.beta(index,:))-update.beta(index,:).*update.m(index,:).^2-...
        log(hyper.beta(index,:))+hyper.beta(index,:).*hyper.m(index,:).^2));
    lowerbound = bound-hyper.mixture.mix(index,:)*(lowerbound+[-update.b(index,:)+hyper.b(index,:),update.a(index,:)-hyper.a(index,:)]*...
        [update.a(index,:)./update.b(index,:);psi(update.a(index,:))-log(update.b(index,:))]+update.a(index,:).*log(update.b(index,:))-...
       gammaln(update.a(index,:))-hyper.a(index,:).*log(hyper.b(index,:))+gammaln(hyper.a(index,:)));
   %   lowerbound calculation
end
%% Gaussian node


function [update,lowerbound] = dirichlet_update(update,hyper,phi,bound)
    if size(hyper.alpha,2) == 1
        message = sum(phi,1)';
        phi = bsxfun(@plus,hyper.alpha-1,message);
        update.alpha = phi+1;
    else
        update.alpha = phi+1+hyper.alpha;
% pre-programmed the message into the forward-backward algorithm
    end
% straight out of a Winn's thesis

    for rk = 1:size(update.alpha,2)
        update.mix(:,rk) = exp(psi(update.alpha(:,rk))-psi(sum(update.alpha(:,rk))));
    end
    update.mix = bsxfun(@times,update.mix,1./sum(update.mix,1));
    lowerbound = bound;
  
    kl = @(alpha,beta) gammaln(sum(alpha))-gammaln(sum(beta))-sum(gammaln(alpha))+sum(gammaln(beta))+(alpha-beta)*(psi(alpha)-psi(sum(alpha)))';
    for rk = 1:size(phi,2)
        lowerbound = lowerbound-kl(update.alpha(:,rk)',hyper.alpha(:,rk)');
    end
   
% update the lowerbound
% just for the dimensions really

end
%% Dirichlet node

function [conprob] = gauss_prob_evaluate(x,post)
    conprob = exp(bsxfun(@plus,-log(2*pi)+(psi(post.a')-...
        log(post.b'))-post.a'./(2*post.b').*(1./post.beta'+post.m'.^2),...
        bsxfun(@times,post.a'./(2*post.b'),(bsxfun(@plus,-x.^2,...
        2*bsxfun(@times,post.m',x))))));
    conprob = bsxfun(@times,1./sum(conprob,2),conprob);
end

function [out,like] = fb_activation(post,conprob,T)
    scale = zeros(T+1,1);
    for rk = 1:post.depth
        for rrk = 1:post.path(rk,end)
            fb.param_depth(rk).currentnode(rrk).alpha_entry = ones(T,size(post.param_depth(rk).hhmm_params(rrk).pi,2));
            fb.param_depth(rk).currentnode(rrk).alpha_exit = fb.param_depth(rk).currentnode(rrk).alpha_entry;
            fb.param_depth(rk).currentnode(rrk).beta_entry = fb.param_depth(rk).currentnode(rrk).alpha_entry;
            fb.param_depth(rk).currentnode(rrk).beta_exit = fb.param_depth(rk).currentnode(rrk).alpha_entry;
            if rk == 1
                fb.param_depth(rk).currentnode(rrk).alpha_entry(1,:) = post.param_depth(rk).hhmm_params(rrk).pi;
            else
                parent = unique(post.path(rk-1,post.path(rk,:)==rrk));
                ca = mod(rrk,post.depth_num(rk))+post.depth_num(rk)*(mod(rrk,post.depth_num(rk))==0);
                fb.param_depth(rk).currentnode(rrk).alpha_entry(1,:) = fb.param_depth(rk-1).currentnode(parent).alpha_entry(1,ca).*post.param_depth(rk).hhmm_params(rrk).pi;
            end
        end
    end
    % alpha boundary conditions
    
    for time = 2:T
        for rk = post.depth:-1:1
            for rrk = 1:post.path(rk,end)
                if rk == post.depth
                    fb.param_depth(rk).currentnode(rrk).alpha_exit(time-1,:) = fb.param_depth(rk).currentnode(rrk).alpha_entry(time-1,:).*conprob(time-1,:);
                    scale(time-1,:) = scale(time-1,:)+sum(fb.param_depth(rk).currentnode(rrk).alpha_exit(time-1,:));
                else
                    if rk ~= 1
                        child = unique(post.path(rk+1,post.path(rk,:)==rrk));
                        for zk = 1:length(child)
                            fb.param_depth(rk).currentnode(rrk).alpha_exit(time-1,zk) = sum(fb.param_depth(rk+1).currentnode(child(zk)).alpha_exit(time-1,:).*post.param_depth(rk+1).hhmm_params(child(zk)).exit);
                        end
                    else
                        for zk = 1:post.path(rk+1,end)
                            fb.param_depth(rk).currentnode(rrk).alpha_exit(time-1,zk) = sum(fb.param_depth(rk+1).currentnode(zk).alpha_exit(time-1,:).*post.param_depth(rk+1).hhmm_params(zk).exit);
                        end
                    end
                end
            end
            if rk == post.depth 
                for rrk = 1:post.path(rk,end)
                    fb.param_depth(rk).currentnode(rrk).alpha_exit(time-1,:) = fb.param_depth(rk).currentnode(rrk).alpha_exit(time-1,:)./scale(time-1,:);
                end
            end
        end
        for rk = 1:post.depth
            if rk == 1
                for rrk = 1:post.path(rk,end)
                    fb.param_depth(rk).currentnode(rrk).alpha_entry(time,:) = fb.param_depth(rk).currentnode(rrk).alpha_exit(time-1,:)*post.param_depth(rk).hhmm_params(rrk).tm';
                end
            else
                for rrk = 1:post.path(rk,end)
                    parent = unique(post.path(rk-1,post.path(rk,:)==rrk));
                    ca = mod(rrk,post.depth_num(rk))+post.depth_num(rk)*(mod(rrk,post.depth_num(rk))==0);
                    fb.param_depth(rk).currentnode(rrk).alpha_entry(time,:) = fb.param_depth(rk-1).currentnode(parent).alpha_entry(time,ca).*post.param_depth(rk).hhmm_params(rrk).pi+fb.param_depth(rk).currentnode(rrk).alpha_exit(time-1,:)*post.param_depth(rk).hhmm_params(rrk).tm';
                end
            end
        end
    end
    % forward pass
    
    time = T+1;
    for rk = post.depth:-1:1
        for rrk = 1:post.path(rk,end)
            if rk == post.depth
                fb.param_depth(rk).currentnode(rrk).alpha_exit(time-1,:) = fb.param_depth(rk).currentnode(rrk).alpha_entry(time-1,:).*conprob(time-1,:);
                scale(time-1,:) = scale(time-1,:)+sum(fb.param_depth(rk).currentnode(rrk).alpha_exit(time-1,:));
            else
                child = unique(post.path(rk+1,post.path(rk,:)==rk));
                for zk = 1:length(child)
                    fb.param_depth(rk).currentnode(rrk).alpha_exit(time-1,zk) = sum(fb.param_depth(rk+1).currentnode(child(zk)).alpha_exit(time-1,:).*post.param_depth(rk+1).hhmm_params(child(zk)).exit);
                end
            end
        end        
        if rk == post.depth 
            for rrk = 1:post.path(rk,end)
                fb.param_depth(rk).currentnode(rrk).alpha_exit(time-1,:) = fb.param_depth(rk).currentnode(rrk).alpha_exit(time-1,:)./scale(time-1,:);
            end
        end
    end
    % forward pass cleanup!
    
    scale(T+1,:) = sum(fb.param_depth(1).currentnode.alpha_exit(T,:).*post.param_depth(1).hhmm_params(1).exit);
    % probability of the FDT state
    
    fb.param_depth(1).currentnode(1).beta_exit(T,:) = post.param_depth(1).hhmm_params(1).exit./scale(T+1,:);
    for rk = 2:post.depth
        for rrk = 1:post.path(rk,end)
            parent = unique(post.path(rk-1,post.path(rk,:)==rrk));
            ca = mod(rrk,post.depth_num(rk))+post.depth_num(rk)*(mod(rrk,post.depth_num(rk))==0);
            fb.param_depth(rk).currentnode(rrk).beta_exit(T,:) = fb.param_depth(rk-1).currentnode(parent).beta_exit(T,ca).*post.param_depth(rk).hhmm_params(rrk).exit;
        end
    end
    % beta variable boundary conditions
    
    for time = T-1:-1:1
        for rk = post.depth:-1:1
            for rrk = 1:post.path(rk,end)
                if rk == post.depth
                    fb.param_depth(rk).currentnode(rrk).beta_entry(time+1,:) = fb.param_depth(rk).currentnode(rrk).beta_exit(time+1,:).*conprob(time+1,:)./scale(time+1,:);
                else
                    child = unique(post.path(rk+1,post.path(rk,:)==rk));
                    for zk = 1:length(child)
                        fb.param_depth(rk).currentnode(rrk).beta_entry(time+1,zk) = sum(fb.param_depth(rk+1).currentnode(child(zk)).beta_entry(time+1,:).*post.param_depth(rk+1).hhmm_params(child(zk)).pi);
                    end
                end
            end
        end
        for rk = 1:post.depth
            for rrk = 1:post.path(rk,end)
                if rk == 1
                    fb.param_depth(rk).currentnode(rrk).beta_exit(time,:) = fb.param_depth(rk).currentnode(rrk).beta_entry(time+1,:)*post.param_depth(rk).hhmm_params(rrk).tm;
                else
                    parent = unique(post.path(rk-1,post.path(rk,:)==rrk));
                    ca = mod(rrk,post.depth_num(rk))+post.depth_num(rk)*(mod(rrk,post.depth_num(rk))==0);
                    fb.param_depth(rk).currentnode(rrk).beta_exit(time,:) = fb.param_depth(rk-1).currentnode(parent).beta_exit(time,ca).*post.param_depth(rk).hhmm_params(rrk).exit+fb.param_depth(rk).currentnode(rrk).beta_entry(time+1,:)*post.param_depth(rk).hhmm_params(rrk).tm;
                end
            end
        end
    end
    % backward pass
    
    time = 0;
    for rk = post.depth:-1:1
        for rrk = 1:post.path(rk,end)
            if rk == post.depth
                fb.param_depth(rk).currentnode(rrk).beta_entry(time+1,:) = fb.param_depth(rk).currentnode(rrk).beta_exit(time+1,:).*conprob(time+1,:)./scale(time+1,:);
            else
                child = unique(post.path(rk+1,post.path(rk,:)==rk));
                for zk = 1:length(child)
                    fb.param_depth(rk).currentnode(rrk).beta_entry(time+1,zk) = sum(fb.param_depth(rk+1).currentnode(child(zk)).beta_entry(time+1,:).*post.param_depth(rk+1).hhmm_params(child(zk)).pi);
                end
            end
        end        
    end
    % backward pass cleanup!
    
    for rk = 1:post.depth
        for rrk = 1:post.path(rk,end)
            for zk = 1:post.depth_num(rk+1)
                for zzk = 1:post.depth_num(rk+1)
                    out.param_depth(rk).g(rrk).a(zzk,zk) = sum(fb.param_depth(rk).currentnode(rrk).alpha_exit(1:end-1,zk).*post.param_depth(rk).hhmm_params(rrk).tm(zzk,zk).*fb.param_depth(rk).currentnode(rrk).beta_entry(2:end,zzk));
                    a(:,zzk,zk) = fb.param_depth(rk).currentnode(rrk).alpha_exit(1:end-1,zk).*post.param_depth(rk).hhmm_params(rrk).tm(zk,zzk).*fb.param_depth(rk).currentnode(rrk).beta_entry(2:end,zzk);
                end
            end
            if rk ~= 1
                parent = unique(post.path(rk-1,post.path(rk,:)==rrk));
                ca = mod(rrk,post.depth_num(rk))+post.depth_num(rk)*(mod(rrk,post.depth_num(rk))==0);
                out.param_depth(rk).g(rrk).pi = fb.param_depth(rk).currentnode(rrk).alpha_entry(1,:).*fb.param_depth(rk).currentnode(rrk).beta_entry(1,:)+sum(bsxfun(@times,fb.param_depth(rk-1).currentnode(parent).alpha_entry(2:end,ca),bsxfun(@times,fb.param_depth(rk).currentnode(rrk).beta_entry(2:end,:),post.param_depth(rk).hhmm_params(rrk).pi)),1);
                out.param_depth(rk).g(rrk).exit = fb.param_depth(rk).currentnode(rrk).alpha_exit(T,:).*fb.param_depth(rk).currentnode(rrk).beta_exit(T,:)+sum(bsxfun(@times,fb.param_depth(rk-1).currentnode(parent).beta_exit(2:end,ca),bsxfun(@times,fb.param_depth(rk).currentnode(rrk).alpha_exit(2:end,:),post.param_depth(rk).hhmm_params(rrk).exit)));
                b = bsxfun(@times,fb.param_depth(rk-1).currentnode(parent).beta_exit(2:end,ca),bsxfun(@times,fb.param_depth(rk).currentnode(rrk).alpha_exit(2:end,:),post.param_depth(rk).hhmm_params(rrk).exit));
            else
                out.param_depth(rk).g(rrk).pi = fb.param_depth(rk).currentnode(rrk).alpha_entry(1,:).*fb.param_depth(rk).currentnode(rrk).beta_entry(1,:);
                out.param_depth(rk).g(rrk).exit = fb.param_depth(rk).currentnode(rrk).alpha_exit(T,:).*fb.param_depth(rk).currentnode(rrk).beta_exit(T,:);
            end
        end
    end
    for rrk = 1:post.path(end,end)
        out.g(rrk).gamma = fb.param_depth(end).currentnode(rrk).alpha_exit.*fb.param_depth(end).currentnode(rrk).beta_exit;
    end
    % the precursor sufficient statistics. From here, the forward and
    % backward variables are invisible. These will be plugged directly into
    % the dirichlets, etc.
    
    like = sum(log(scale));
    % likelihood function. Done.
end

function [post] = initialize(x,N,options)
    dat = x;
    x = cell2mat(x(:));
    depth_vector = fliplr(options.depth_vector(:)');
    depth = length(options.depth_vector)+1;
    prod_states = options.production;
    guess = options.guess;
    count = ones(depth,1);
    active = zeros(depth,1);
    active(1) = 1;
    path = zeros(depth,prod(depth_vector));
    for rrk = 1:prod(depth_vector)
        path(:,rrk) = count;
        for rk = 1:length(count)-1
            if mod(count(rk),depth_vector(rk))==0&&count(rk)~=1&&active(rk)
                active(rk+1) = 1;
            else
                active(rk+1) = 0;
            end
        end
        count = count+active;
    end
    path = fliplr(path')';
    post.pathnum = size(path,2);
    post.depth = depth;
    post.path = path;
    depth_vector = [1 fliplr(depth_vector)];
    depth_vector(end+1) = prod_states;
    post.depth_num = depth_vector;
    for rk = depth:-1:1
        for rrk = 1:post.path(rk,end)
            post.param_depth(rk).hhmm_params(rrk).pi = ones(1,depth_vector(rk+1))./depth_vector(rk+1);
            % easy
            
            if rk ~= depth
                post.param_depth(rk).hhmm_params(rrk).tm = eye(depth_vector(rk+1),depth_vector(rk+1)).*(.1+rand(depth_vector(rk+1),depth_vector(rk+1)))/5+~eye(depth_vector(rk+1),depth_vector(rk+1));
            else
                post.param_depth(rk).hhmm_params(rrk).tm = ~eye(depth_vector(rk+1),depth_vector(rk+1)).*(.1+rand(depth_vector(rk+1),depth_vector(rk+1)))/5+eye(depth_vector(rk+1),depth_vector(rk+1));    
            end
            post.param_depth(rk).hhmm_params(rrk).exit = min(post.param_depth(rk).hhmm_params(rrk).tm,[],1)/10;
            temp = [post.param_depth(rk).hhmm_params(rrk).tm;post.param_depth(rk).hhmm_params(rrk).exit];
            post.param_depth(rk).hhmm_params(rrk).tm = bsxfun(@times,post.param_depth(rk).hhmm_params(rrk).tm,1./sum(temp,1));
            post.param_depth(rk).hhmm_params(rrk).exit = post.param_depth(rk).hhmm_params(rrk).exit./sum(temp,1);
            % can either exit up or within a family; hence such
            % normalization
            
            if rk == depth
                bound = length(x)/post.path(rk,end)/20;
            else
                bound = carryover.param_depth(rk).node(rrk)/20;
            end
            temp = [post.param_depth(rk).hhmm_params(rrk).tm;post.param_depth(rk).hhmm_params(rrk).exit];
            temp2 = -log(post.param_depth(rk).hhmm_params(rrk).tm);
            temp2 = temp2(~~eye(post.depth_num(rk+1),post.depth_num(rk+1)));
            temp2 = 1-temp2./sum(temp2);
            guess_tm_exit = bsxfun(@times,bound./post.path(rk,end).*temp,temp2');
            guess_pi = bound/post.path(rk,end)*post.param_depth(rk).hhmm_params(rrk).pi';
            post.param_depth(rk).hhmm_posterior(rrk).tm_exit.alpha = getprior(bound,guess_tm_exit,temp);
            post.param_depth(rk).hhmm_posterior(rrk).pi.alpha = getprior(bound,guess_pi,post.param_depth(rk).hhmm_params(rrk).pi');
            % using Minka's EM routine for the Dirichlet prior
        end        
        if rk ~= 1
            for zk = 1:post.path(rk-1,end)
                carryover.param_depth(rk-1).node(zk) = 0;
                child = unique(post.path(rk,post.path(rk-1,:)==zk));
                for zzk = 1:length(child)
                    carryover.param_depth(rk-1).node(zk) = carryover.param_depth(rk-1).node(zk)+sum(post.param_depth(rk).hhmm_posterior(child(zzk)).tm_exit.alpha(end,:));
                end
            end
        end
    end
    if options.baseline
        for rk = 1:length(dat)
            post.emissions(rk).mixture.mix = ones(depth_vector(end),1)/depth_vector(end);
            post.emissions(rk).mixture.alpha = post.emissions(rk).mixture.mix*length(dat{rk})/10;
            post.emissions(rk).a = ones(depth_vector(end),1);
            post.emissions(rk).b = var(dat{rk})/depth_vector(end)^2.*ones(depth_vector(end),1);
            post.emissions(rk).beta = 1e2*ones(depth_vector(end),1);
            post.emissions(rk).m = sort(guess(:,rk));
        end
    else
        post.emissions.mixture.mix = ones(depth_vector(end),1)/depth_vector(end);
        post.emissions.mixture.alpha = post.emissions.mixture.mix*length(x)/10;
        post.emissions.a = ones(depth_vector(end),1);
        post.emissions.b = var(x)/depth_vector(end)^2.*ones(depth_vector(end),1);
        post.emissions.beta = 1e2*ones(depth_vector(end),1);
        post.emissions.m = guess(:);
    end
end

function out = viterbi_activation(post,conprob,T,index)
    for rk = 1:post.depth
        for rrk = 1:post.path(rk,end)
            vit.param_depth(rk).currentnode(rrk).dE = ones(1,post.depth_num(rk+1));
            vit.param_depth(rk).currentnode(rrk).dB = ones(1,post.depth_num(rk+1));
            vit.param_depth(rk).currentnode(rrk).psiB = ones(T,post.depth_num(rk+1));
            if rk~=post.depth
                vit.param_depth(rk).currentnode(rrk).psiE = ones(T,post.depth_num(rk+1));
            end
        end
    end
    z = zeros(T,post.depth+1);
    d = ones(T,1);
    % initialize variables
    
    vit.param_depth(1).currentnode.dB = log(post.param_depth(1).hhmm_params.pi);
    for rk = 2:post.depth
        for rrk = 1:post.path(rk,end)
            parent = unique(post.path(rk-1,post.path(rk,:)==rrk));
            ca = mod(rrk,post.depth_num(rk))+post.depth_num(rk)*(mod(rrk,post.depth_num(rk))==0);
            vit.param_depth(rk).currentnode(rrk).dB = vit.param_depth(rk-1).currentnode(parent).dB(ca)+log(post.param_depth(rk).hhmm_params(rrk).pi);
        end
    end
    % Entry probability
    
    for time = 1:T
        for rrk = 1:post.path(end,end)
            vit.param_depth(post.depth).currentnode(rrk).dE = vit.param_depth(post.depth).currentnode(rrk).dB+log(conprob(rrk).gamma(time,:));
        end
        % pull in the state designations from FBA algorithm
        
        for rk = post.depth-1:-1:1
            for rrk = 1:post.path(rk,end)
                child = unique(post.path(rk+1,post.path(rk,:)==rrk));
                for zk = 1:length(child)
                    prob = vit.param_depth(rk+1).currentnode(child(zk)).dE+log(post.param_depth(rk+1).hhmm_params(child(zk)).exit);
                    [a,b] = nanmax(prob);
                    vit.param_depth(rk).currentnode(rrk).psiE(time,zk) = b;
                    vit.param_depth(rk).currentnode(rrk).dE(zk) = a;
                end
            end
        end
        % make a map of where the system would go were there a vertical
        % transition
        
        if time ~= T
            for zk = 1:post.path(2,end)
                prob = vit.param_depth(1).currentnode(1).dE+log(post.param_depth(1).hhmm_params(1).tm(:,zk))';
                [a,b] = nanmax(prob);
                vit.param_depth(1).currentnode(1).psiB(time+1,zk) = b;
                vit.param_depth(1).currentnode(1).dB(zk) = a;
            end
            for rk = 2:post.depth
                for rrk = 1:post.path(rk,end)
                    for zk = 1:size(post.param_depth(rk).hhmm_params(rrk).pi,2)
                        parent = unique(post.path(rk-1,post.path(rk,:)==rrk));
                        ca = mod(rrk,post.depth_num(rk))+post.depth_num(rk)*(mod(rrk,post.depth_num(rk))==0);
                        prob = vit.param_depth(rk).currentnode(rrk).dE+log(post.param_depth(rk).hhmm_params(rrk).tm(:,zk))';
                        prob(end+1) = vit.param_depth(rk-1).currentnode(parent).dB(ca)+log(post.param_depth(rk).hhmm_params(rrk).pi(zk));
                        [a,b] = nanmax(prob);
                        vit.param_depth(rk).currentnode(rrk).psiB(time+1,zk) = b;
                        vit.param_depth(rk).currentnode(rrk).dB(zk) = a;
                    end
                end
            end
        end
        % map for everything else
    end
    prob = vit.param_depth(1).currentnode(1).dE+log(post.param_depth(1).hhmm_params(1).exit);
    % exiting the F1T state
    [~,b] = max(prob);
    phi(1) = 1;
    phi(2) = b;
    for rk = 2:post.depth
        ca = mod(phi(rk),post.depth_num(rk))+post.depth_num(rk)*(mod(phi(rk),post.depth_num(rk))==0);
        phi(rk+1) = post.depth_num(rk).*(phi(rk)-1)+vit.param_depth(rk-1).currentnode(phi(rk-1)).psiE(T-1,ca);
    end
    % backtrack
    ca = mod(phi(post.depth+1),post.depth_num(post.depth+1))+post.depth_num(post.depth+1)*(mod(phi(post.depth+1),post.depth_num(post.depth+1))==0);
    z(T,:) = phi;
    meanpath = zeros(T,1);
    meanpath(T) = post.emissions(index).m(ca);
    d(T,:) = 1;
    for time = T-1:-1:1        
        for rk = post.depth:-1:1
            ca = mod(phi(rk+1),post.depth_num(rk+1))+post.depth_num(rk+1)*(mod(phi(rk+1),post.depth_num(rk+1))==0);
            if vit.param_depth(rk).currentnode(phi(rk)).psiB(time+1,ca) <= post.depth_num(rk+1)
                d(time) = rk+1;
                phi(rk+1) = post.depth_num(rk+1).*(phi(rk)-1)+vit.param_depth(rk).currentnode(phi(rk)).psiB(time+1,ca);
                break;
            end
        end
        for rk = d(time):post.depth
            ca = mod(phi(rk),post.depth_num(rk))+post.depth_num(rk)*(mod(phi(rk),post.depth_num(rk))==0);
            phi(rk+1) = post.depth_num(rk+1).*(phi(rk)-1)+vit.param_depth(rk-1).currentnode(phi(rk-1)).psiE(time,ca);
        end
        ca = mod(phi(post.depth+1),post.depth_num(post.depth+1))+post.depth_num(post.depth+1)*(mod(phi(post.depth+1),post.depth_num(post.depth+1))==0);
        meanpath(time) = post.emissions(index).m(ca);
        z(time,:) = phi;
    end
    out.pathstate = z;
    out.ideals = meanpath;
end

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

