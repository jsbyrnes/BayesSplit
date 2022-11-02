function model = fill_gradient_waveforms(model, allWfs, Parameters)

    del = zeros(length(model.vector),1);%overwrite it

    if Parameters.parallel

        parfor k = 1:length(model.vector)
    
            del(k, 1) = internal_loop(model, allWfs, Parameters, k);
    
        end

    else

        for k = 1:length(model.vector)
    
            del(k, 1) = internal_loop(model, allWfs, Parameters, k);
    
        end

    end

    model.del = del;

end

function [ del, dN, dE ] = internal_loop(model, allWfs, Parameters, index)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get the linear index of what will be changed
    [nevt, nsta] = size(allWfs);
    nsplit       = numel(model.A);%number of splitting parameters. Could be clustered

    scale  = 1;

    index0 = index;

    if nargin < 4
        index  = randi([1 model.nparam]);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %map the index to the structure fields
    evtind = [];
    staind = [];
    type   = [];
    field  = [];

    if nsta > 1 && Parameters.wavelet

        nwavelet     = numel(model.wavelet);%n in wavelets
        %change to the wavelet
        if index <= nwavelet
    
            [evtind, ~] = ind2sub(size(model.wavelet), index);
            staind      = 1:nsta;
            
            field = 'wavelet';
            type  = 'all';
    
        else
    
            index = index - nwavelet;
    
        end
    
        if Parameters.use_polarization
    
            %change to polarization
            if index <= nevt && isempty(field)
        
                staind = 1:nsta; 
                evtind = index;%maps directly
                field  = 'polarization';
                type   = 'all';
                scale  = Parameters.polarization_std;
        
            elseif isempty(field)
        
                index = index - nevt;
        
            end
    
        end
    
        %%%%amp
        if index <= nsta*nevt && isempty(field)
    
            [evtind, staind] = ind2sub(size(allWfs), index);
            
            field = 'amp';
            scale = 1;
    
        elseif isempty(field)
    
            index = index - nsta*nevt;
    
        end
    
        %%%%delay
        if index <= nsta*nevt && isempty(field)
    
            [evtind, staind] = ind2sub(size(allWfs), index);
            
            field = 'shift';
            scale = Parameters.Delaystd;
    
        elseif isempty(field)
    
            index = index - nsta*nevt;
    
        end
    
        %%%%tS
        if index <= nsta*nevt && isempty(field)
    
            [evtind, staind] = ind2sub(size(allWfs), index);
            
            field = 'dtS';
            scale = Parameters.dtS(2);
    
        elseif isempty(field)
    
            index = index - nsta*nevt;
    
        end

    end

    %change to A/dt
    if index <= nsplit && isempty(field)

        if Parameters.cluster

            staind = 1:nsta;

        else

            [staind, ~] = ind2sub(size(model.A), index);%dont need layer ind

        end

        evtind      = 1:nevt;

        field = 'A';
        scale = Parameters.dtstd;%same search size for each
        type  = 'all';

    elseif isempty(field)

        index = index - nsplit;

    end

    %change to B/fastdir
    if index <= nsplit && isempty(field)

        %A works for size on all splitting parameters
        if Parameters.cluster

            staind = 1:nsta;

        else

            [staind, ~] = ind2sub(size(model.A), index);%dont need layer ind

        end

        evtind      = 1:nevt;
            
        field = 'B';
        scale = Parameters.dtstd;%same search size for each

        type  = 'all';

    elseif isempty(field)

        index = index - nsplit;

    end

    %change to rotation
    if index <= nsplit && isempty(field)

        %A works for size on all splitting parameters
        [~, lind] = ind2sub(size(model.A), index);
        if Parameters.cluster

            staind = 1:nsta;

        else

            [staind, ~] = ind2sub(size(model.A), index);%dont need layer ind

        end
        evtind = 1:nevt;
        field  = 'fast_dir_rotation';
        scale = Parameters.rotation_std;

        type  = 'all';

    elseif isempty(field)

        index = index - nsplit;

    end

    if Parameters.use_tSA && index <= nsplit && isempty(field)

        %A works for size on all splitting parameters
        [~, lind] = ind2sub(size(model.A), index);
        if Parameters.cluster

            staind = 1:nsta;

        else

            [staind, ~] = ind2sub(size(model.A), index);%dont need layer ind

        end
        
        evtind = 1:nevt;
        field  = 'tSA';
        scale = Parameters.tSAstd;

        type  = 'all';

    elseif Parameters.use_tSA && isempty(field)

        index = index - nsplit;

    end

    %change to orientations
    if Parameters.use_orientations && nsta > 1 && Parameters.wavelet

        if index <= nsta && isempty(field)
    
            staind = index;%applies directly
            evtind = 1:nevt;
            
            field = 'sta_or';
            scale = Parameters.sta_err(index);
            type  = 'all'; 
    
        elseif isempty(field)
    
            index = index - nsta;
    
        end

    end
    
    %change to sigma
    if index <= nsta*nevt && isempty(field)

        [evtind, staind] = ind2sub(size(allWfs), index);
        
        field = 'sig';
        scale = Parameters.sig_range(2);
        type  = 'error';%not really but this is where it is evaluated

    elseif isempty(field)

        index = index - nsta*nevt;

    end

    %change to r
    if index == 1 && isempty(field)

        index = 1;%scalar
        
        evtind = 1:nevt;
        staind = 1:nsta;
        scale = Parameters.r_range(2);
        field = 'r';
        type  = 'error';%not really but this is where it is evaluated

    elseif isempty(field)

        index = index - 1;

    end

    %change to f
    if index == 1 && isempty(field)

        evtind = 1:nevt;
        staind = 1:nsta;
        field  = 'f';
        scale  = Parameters.f_range(2);
        type   = 'error';%not really but this is where it is evaluated

    elseif isempty(field)

        index = index - 1;

    end
    
    modelp = model;
    modeln = model;
    
    h = Parameters.h*scale;

    modelp.(field)(index) = modelp.(field)(index) + h;
    modeln.(field)(index) = modeln.(field)(index) - h;
    modelp                = apply_update(modelp, field, Parameters.t);
    modeln                = apply_update(modeln, field, Parameters.t);

    modelp = evaluate(modelp, allWfs, Parameters, evtind, staind, type);
    modeln = evaluate(modeln, allWfs, Parameters, evtind, staind, type);

    del = (((modelp.lpst) - (modeln.lpst))/(2*h))*scale;

    dN = cell(numel(allWfs), 1);
    dE = cell(numel(allWfs), 1);

    if nargout == 3

        %give the corresponding row of the Jacobian
        
        for k = 1:numel(allWfs)
 
            dN{k} = (modelp.N{k} - modeln.N{k})/(2*h);
            dE{k} = (modelp.E{k} - modeln.E{k})/(2*h);            

        end

    end

end

%    if abs(del) > 100
% 
%        v = (-5:0.1:5)*scale;
% 
%         for k = 1:length(v)
%             mn(k) = model;
%             mn(k).(field)(index) = model.(field)(index) + v(k);
%             mn(k)                = apply_update(mn(k), field, Parameters.t);
%             mn(k)                = evaluate(mn(k), allWfs, Parameters, evtind, staind, type);
%             lpvec(k) = mn(k).lpst;
%         end
%         plot(v, lpvec)
%       
%    end

%     model.G = zeros(model.nparam, model.nparam);
% 
%     ndata = model.nparam - numel(allWfs) - 2;
%     ncvar = model.nparam - ndata;
% 
%     Gdata = zeros(ndata, ndata);
%     Gcvar = zeros(ncvar, ncvar);
% 
%     %Gdata matrix
%     for k = 1:((ndata^2 - ndata)/2)
% 
%         [k1,k2] = ind2sub(size(Gdata), k);
% 
%         if k2 > k1
% 
%             continue
% 
%         end
% 
%         %get derivatives
%         dn1 = dN{k1};
%         de1 = dE{k1};
%         dn2 = dN{k2};
%         de2 = dE{k2};
% 
%         for j = 1:length(dn1)
% 
%             [j1,j2] = ind2sub(size(allWfs), j);
% 
%             Gdata(k1,k2) = model.G(k1,k2) + (dn2{j}'*model.Cinv*dn1{j})/exp(model.sig(j1,j2))^2;
%             Gdata(k1,k2) = model.G(k1,k2) + (de2{j}'*model.Cinv*de1{j})/exp(model.sig(j1,j2))^2;
% 
%         end
% 
%     end
% 
%     %Gcvar matrix
%     %covarience parameters, do numerically
%     Cinv = model.Cinv;%take this out of the structure
%     
%     %r
%     modelp = model;
%     modeln = model;
%     
%     modelp.r = modelp.r + Parameters.h;
%     modeln.r = modeln.r + Parameters.h;
%     
%     [~, Cp] = build_C(modelp, Parameters.t);
%     [~, Cn] = build_C(modeln, Parameters.t);
%     
%     dCr = (Cp - Cn)/(2*Parameters.h);
%     
%     %f
%     modelp = model;
%     modeln = model;
%     
%     modelp.f = modelp.f + Parameters.h;
%     modeln.f = modeln.f + Parameters.h;
%     
%     [~, Cp] = build_C(modelp, Parameters.t);
%     [~, Cn] = build_C(modeln, Parameters.t);
%     
%     dCf = (Cp - Cn)/(2*Parameters.h);
% 
%     for k = 1:((ncvar^2 - ncvar)/2)
% 
%         [k1,k2] = ind2sub(size(Gcvar), k);
% 
%         if k2 > k1
% 
%             continue
% 
%         end
% 
%         if k1 == k2 && k1 <= (ncvar - 2)
% 
%             %sigma levels
%             %only non-zero on the diagonal
%         
%             sig = exp(model.sig(k1));
%             model.G(k1,k2) = nsta*(4/sig^2);
% 
%         end
% 
%         if k2 == (ncvar - 1) && k1 <= (ncvar - 2)
% 
%             sig = exp(model.sig(k1));
%             model.G(k1,k2) = nsta*(4/sig^2);
% 
%         end
% 
%     end

