function model = evaluate(model, allWfs, Parameters, evt_ind, sta_ind, type, prelogf)
    
    if isempty(type)

        type = 'all';%not used, a good thing

    end

    t = Parameters.t;
    [nevt, nsta] = size(allWfs);

    if nargin < 7

        fs = 1/(t(2) - t(1));
    
        dFreq = fs/length(t); %frequency spacing
        fNyq=fs/2;   %Nyquist frequency
        
        %the next two steps build the frequency vector (freqs at which spectrum
        %was calculated) postivie and negative
        f=(0:length(t)-1)'*dFreq;
        f(f>fNyq)=f(f>fNyq)-fNyq*2;
    
        f = -2*pi*f;
        prelogf = log([1; f(2:end)]/(400*2*pi));

    end

    model = apply_model(model, allWfs, Parameters, evt_ind, sta_ind, prelogf);

    model.llh = -0.5*log(2*pi)*length(t) - ...
        (2*Parameters.n*model.logdet + 4*length(t)*sum(model.sig(:)))/2 ...
       - model.phi/2;
   
    model.lp = sum( ((model.sig(:) - Parameters.sig_range(1)).^2)/(Parameters.sig_range(2)^2));
    
    if nsta > 1 && Parameters.wavelet
         
        model.lp = model.lp - (sum(model.wavelet(:).^2) +...
            sum( ((model.dtS(:) - Parameters.dtS(1)).^2)/Parameters.dtS(2)^2) + ...
            sum( ((model.amp(:)).^2)) + ...
            sum( ((model.shift(:)).^2)/Parameters.Delaystd^2) )/2;

    end

    if Parameters.use_covarience

        model.lp = model.lp - (((model.r - Parameters.r_range(1))^2)/Parameters.r_range(2)^2 + ...
            ((model.f - Parameters.f_range(1))^2)/Parameters.f_range(2)^2)/2;
    
    end

    if Parameters.prior_information.use

        model.lp = model.lp - (...
                sum((( model.dt(:, 1) - Parameters.prior_information.dt(1)).^2)/Parameters.prior_information.dt_std(1)^2 ) ...
            +   sum( -cos(2*(model.fast_dir(:, 1) - Parameters.prior_information.phi(1)))/Parameters.prior_information.phi_std(1)^2) ...%note the negative sign
            +   sum((( model.fast_dir_rotation(:, 1) - Parameters.prior_information.rot(1)).^2)/Parameters.prior_information.rot_std(1)^2))/2;

        if length(Parameters.prior_information.dt) == 2

            model.lp = model.lp - (sum((( model.dt(:, 2) - Parameters.prior_information.dt(2)).^2)/Parameters.prior_information.dt_std(2)^2 ) ...
                + sum( -cos(2*(model.fast_dir(:, 2) - Parameters.prior_information.phi(2)))/Parameters.prior_information.phi_std(2)^2) ...
                + sum((( model.fast_dir_rotation(:, 2) - Parameters.prior_information.rot(2)).^2)/Parameters.prior_information.rot_std(2)^2))/2;

        end

        if Parameters.prior_information.cluster

            %times
            model.lp = model.lp - sum( ((model.dt(:, 1) - mean(model.dt(:, 1))).^2)/...
                (Parameters.prior_information.cluster_weight*Parameters.prior_information.dt_std(1))^2);

            %fast directions
            %get the distance from the center first
            fd_dist = circ_dist(2*model.fast_dir(:, 1), circ_mean(2*model.fast_dir(:, 1)))/2;
            model.lp = model.lp + sum( cos(2*(fd_dist))...
                /(Parameters.prior_information.cluster_weight*Parameters.prior_information.phi_std(1))^2);
            %rotations
            model.lp = model.lp - sum( ((model.fast_dir_rotation(:, 1) - mean(model.fast_dir_rotation(:, 1))).^2)/...
                (Parameters.prior_information.cluster_weight*Parameters.prior_information.rot_std(1))^2);

            if length(Parameters.prior_information.dt) == 2
                
                %times
                model.lp = model.lp - sum( ((model.dt(:, 2) - mean(model.dt(:, 2))).^2)/...
                    (Parameters.prior_information.cluster_weight*Parameters.prior_information.dt_std(2))^2);
    
                %fast directions
                %get the distance from the center first
                fd_dist = circ_dist(2*model.fast_dir(:, 2), circ_mean(2*model.fast_dir(:, 2)))/2;
                model.lp = model.lp - sum( -cos(2*(fd_dist))...
                    /(Parameters.prior_information.cluster_weight*Parameters.prior_information.phi_std(2))^2);
                %rotations
                model.lp = model.lp - sum( ((model.fast_dir_rotation(:, 2) - mean(model.fast_dir_rotation(:, 2))).^2)/...
                    (Parameters.prior_information.cluster_weight*Parameters.prior_information.rot_std(2))^2);

            end

        end

    else

%         model.lp = model.lp - (sum( ((model.fast_dir_rotation(:)).^2)/Parameters.rotation_std^2) + ...
%         sum( ((model.A(:)).^2)/Parameters.ABstd^2) + sum( ((model.B(:)).^2)/Parameters.ABstd^2))/2;
    %seperate priors for A&B biases the fast directions
        model.lp = model.lp - ( sum(((model.fast_dir_rotation(:)).^2)/Parameters.rotation_std^2) + ...
            sum( ((model.dt(:)).^2)/Parameters.dtstd^2) )/2;

    end

    if Parameters.use_polarization && nsta > 1

        kappa = 1/Parameters.polarization_std^2;
        model.lp = model.lp + sum(kappa*(cos(model.polarization))); %von mises for polarizations

    end

    if Parameters.tSAstd

        model.lp = model.lp - 0.5*sum((model.tSA(:).^2)/Parameters.tSAstd^2); %von mises for polarizations

    end

    %need to loop to do for stations since kappa can vary
    %and then damp the difference in phi to avoid aligned or
    %perpendicular splitting, also by station

    if Parameters.use_orientations && nsta > 1

        for k = 1:nsta
        
            model.lp = model.lp + cos(model.sta_or(k))/(Parameters.sta_err(k)^2);%zero mean
    
        end

    end

    model.lpst = model.llh + model.lp;

end

%%%%%%%%
%old search for alignment
%                     lags = 1:length(dN);
%                     for j = 1:length(lags)
% 
%                         %comutationally intensive
%                         phi(j) = md(circshift(No, lags(j)), circshift(Eo, lags(j)));
%     
%                     end
%     
%                     [R(2), idx] = min(phi);
%                     %find neighbors, assuming periodic
%                     if idx==1
%     
%                         R(1)=phi(end);
%                         R(3)=phi(2);
%     
%                     elseif idx==numel(phi)
%     
%                         R(1)=phi(end-1);
%                         R(3)=phi(1);
%     
%                     else
%     
%                         R(1)=phi(idx-1);
%                         R(3)=phi(idx+1);
%     
%                     end
%     
%                     c     = (R(3)-R(1))/(2*(2*R(2)-R(1)-R(3)));
%                     %lag   = mod(idx-1+floor(length(lags)/2),length(lags))-floor(length(lags)/2);%integer part
%                     delay = idx+c;%delay estimate
%                     N     = delay_continuous(No, Parameters.sample_rate, delay/Parameters.sample_rate );
%                     E     = delay_continuous(Eo, Parameters.sample_rate, delay/Parameters.sample_rate );


%double grid search for amplitude
%                     phi = []
%                     amps = (-6:0.05:3)
%                     for j = 1:length(amps)
%     
%                         phi(j) = md(exp(amps(j))*N, exp(amps(j))*E);
%     
%                     end
%     
%                     [mphi, idx] = min(phi);
% 
%                     if idx > 1
% 
%                         %refined search - needs to be extremely accurate
%                         phi = [];
%                         amps = (-0.5:0.05:0.5) + amps(idx);
%                         for j = 1:length(amps)
%         
%                             %computationally intensive
%                             phi(j) = md(exp(amps(j))*N, exp(amps(j))*E);
%         
%                         end
%                         [mphi, idx] = min(phi);        
%                         c     = (phi(idx+1) - phi(idx-1))/(2*(2*mphi - phi(idx-1) - phi(idx+1)));
%                         %ampi  = mod(idx - 1 + floor(length(amps)/2), length(amps)) - floor(length(amps)/2);%integer part
%                         amp   = amps(idx) + c*0.05;%amp estimate
% 
%                     else
% 
%                         amp = -6;%basically zero, its a huge mismatch. 
% 
%                     end
% 
%                     N = exp(amp)*N;
%                     E = exp(amp)*E;

%         if sum(model.dt(k, :)) > model.dt(k, 1)%check for a second layer
% 
%             %account for rotation
% %             dphi     = abs( (model.fast_dir(k, 1) + model.fast_dir_rotation(k,1)) ...
% %                 - (model.fast_dir(k, 2) - model.fast_dir_rotation(k,2)));
%             dphi     = abs( model.fast_dir(k, 1) - model.fast_dir(k, 2));
% 
%             if dphi > pi/2
% 
%                 dphi = dphi - pi/2;
% 
%             end
% 
%             dphi = dphi/(pi/2);
% 
%             model.lp = model.lp + (Parameters.beta - 1)*(log(dphi) + log(1 - dphi));
% 
%         end
