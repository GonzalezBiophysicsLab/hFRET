function out = hhmm_simulate_trace(tms,em,tlen)

    depth = length(tms);
    % design: tms(depth).params(link).tm
    for rk = 1:depth
        if rk~=1
            depth_vector(rk,:) = length(tms(rk).params)/depth_vector(rk-1,:);
        else
            depth_vector(rk,:) = 1;
        end
        for rrk = 1:length(tms(rk).params)
            tms(rk).params(rrk).checkmat1 = [zeros(size(tms(rk).params(rrk).tm,1),1) cumsum(tms(rk).params(rrk).tm(:,1:end-1),2)];
            tms(rk).params(rrk).checkmat2 = cumsum(tms(rk).params(rrk).tm,2);
            tms(rk).params(rrk).ent_chk1 = [0 cumsum(tms(rk).params(rrk).entry(:,1:end-1),2)];
            tms(rk).params(rrk).ent_chk2 = cumsum(tms(rk).params(rrk).entry,2);
    % normalized test variable
        end
    end

    path = getpath(depth,depth_vector);
        
    signpost = 1;
    emstate = 1;
    d = depth;
    % start us off in a production state
    % keep track of it with a logical
    %
    % signpost keeps track of which node we're in
    out.trace = zeros(tlen,1);
    out.id = zeros(tlen,1);
    for rk = 1:tlen
        production = 1==0;
        while ~production
            roll = rand;
            position = find(roll>tms(d).params(signpost).checkmat1(emstate,:)&roll<tms(d).params(signpost).checkmat2(emstate,:));
            if position==size(tms(depth).params(signpost).tm,1)+1
                emstate = mod(signpost,depth_vector(d-1))+depth_vector(d-1)*(mod(signpost,depth_vector(d-1))==0);
                signpost = unique(path(d-1,path(d,:)==signpost));
                d = d-1;
                production = 1==0;
                % the "position" in the upper level transition matrix is an
                % "emission" in the lower level so in the case of an exit,
                % the state just passes up with no comment and prepares for
                % the next roll
            elseif position~=size(tms(depth).params(signpost).tm,1)+1&&d~=depth
                roll = rand;
                emstate = find(roll>tms(d).params(signpost).ent_chk1&roll<tms(d).params(signpost).ent_chk2);
                signpost = (signpost-1)*depth_vector(d)+position;
                % downward is easy, its just a bit of modular arithmetic
                
                production = 1==0;
                d = d+1;
                
                
            elseif position~=size(tms(depth).params(signpost).tm,1)+1&&d==depth
                emstate = position;
                out.trace(rk,:) = em.m(emstate)+em.sigma(emstate)*randn;
                out.id(rk,:) = em.m(emstate);
                out.state(rk,:) = (signpost-1)*length(em.m)+position;
                % output the path index at the signpost
                
                production = 1==1;
            end
        end
    end
end

function path = getpath(depth,depth_vector)
    depth_vector = fliplr(depth_vector(:)');
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
end

