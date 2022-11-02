function model = make_starting_models(TD_parameters, allWfs)

    [nevt, nsta] = size(allWfs);

    %per evt parameters
    for k = 1:nevt

        if TD_parameters.wavelet

            model.wavelet(k, :) = randn(length(TD_parameters.t),1);%model.wavelet(k, :)/rms(model.wavelet(k, :));% + randn();

        else

            model.wavelet = dct(TD_parameters.syn_wavelet);

        end
            
        model.polarization(k,1) = normrnd(TD_parameters.polarization(k), TD_parameters.polarization_std);%fixed for synthetic test
        model.polarization(k,1) = mod(model.polarization(k), 2*pi);

    end

   for k = 1:nsta

        x = [allWfs(:, k).orientation_k];

        if isempty(x(~isnan(x)))
            sta_err(k,1) = 2*pi;
        else
            sta_err(k,1) = unique(x(~isnan(x)));
            sta_err(k,1) = sqrt(1/sta_err(k,1));
        end

    end

    if TD_parameters.prior_information.use

        %per station parameters
        if TD_parameters.cluster
    
            model.dt(1, 1) = randn(1,1)*TD_parameters.prior_information.dt_std(1) + TD_parameters.prior_information.dt(1);
            model.dt(1, 2) = randn(1,1)*TD_parameters.prior_information.dt_std(2) + TD_parameters.prior_information.dt(2);
            model.fast_dir(1,1) = randn(1,1)*TD_parameters.prior_information.phi_std(1) + TD_parameters.prior_information.phi(1);%in radians
            model.fast_dir(1,2) = randn(1,1)*TD_parameters.prior_information.phi_std(2) + TD_parameters.prior_information.phi(2);%in radians
            model   = make_AB(model);
            model   = make_dtphi(model);%enforce a pi/2 wrap
            model.fast_dir_rotation(1,1) = randn(1,1)*TD_parameters.prior_information.rot_std(1) + TD_parameters.prior_information.rot(1);%in radians
            model.fast_dir_rotation(1,2) = randn(1,1)*TD_parameters.prior_information.rot_std(2) + TD_parameters.prior_information.rot(2);%in radians
        
        else
    
            model.dt(:, 1) = randn(nsta,1)*TD_parameters.prior_information.dt_std(1) + TD_parameters.prior_information.dt(1);
            model.dt(:, 2) = randn(nsta,1)*TD_parameters.prior_information.dt_std(2) + TD_parameters.prior_information.dt(2);
            model.fast_dir(:,1) = randn(nsta,1)*TD_parameters.prior_information.phi_std(1) + TD_parameters.prior_information.phi(1);%in radians
            model.fast_dir(:,1) = randn(nsta,1)*TD_parameters.prior_information.phi_std(2) + TD_parameters.prior_information.phi(2);%in radians
            model   = make_AB(model);
            model   = make_dtphi(model);%enforce a pi/2 wrap
            model.fast_dir_rotation(:,1) = randn(nsta,1)*TD_parameters.prior_information.rot_std(1) + TD_parameters.prior_information.rot(1);%in radians
            model.fast_dir_rotation(:,2) = randn(nsta,1)*TD_parameters.prior_information.rot_std(2) + TD_parameters.prior_information.rot(2);%in radians
    
        end

    else

        %per station parameters
        if TD_parameters.cluster
    
            model.A = randn(1,TD_parameters.max_layers)*TD_parameters.ABstd;
            model.B = randn(1,TD_parameters.max_layers)*TD_parameters.ABstd;
            model   = make_dtphi(model);
            model.fast_dir_rotation = randn(1,TD_parameters.max_layers)*TD_parameters.rotation_std;%in radians
        
        else
    
            model.A = randn(nsta,TD_parameters.max_layers)*TD_parameters.ABstd;
            model.B = randn(nsta,TD_parameters.max_layers)*TD_parameters.ABstd;
            model   = make_dtphi(model);
            model.fast_dir_rotation = randn(nsta,TD_parameters.max_layers)*TD_parameters.rotation_std;%in radians
    
        end

    end

    model.sta_or            = randn(nsta,1)*TD_parameters.orientation_std;%approximation to von mises

    %per station-evt pairs
    model.sig      = randn(nevt, nsta)*TD_parameters.sig_range(2) + TD_parameters.sig_range(1);%zeros(nevt, nsta)*diff(TD_parameters.sig_range) + TD_parameters.sig_range(1);

    %for bookkeeping
    model.amp      = zeros(nevt, nsta);%log
    model.shift    = zeros(nevt, nsta)*TD_parameters.Delaystd;
    model.dtS      = zeros(nevt, nsta)*TD_parameters.dtS(2);%*TD_parameters.max_dtS;%ones(20, 1);% abs(randn(n,1)*TD_parameters.dt_std);

    %global parameters
    model.r = randn()*TD_parameters.r_range(2) + TD_parameters.r_range(1);%rand(n,1)*TD_parameters.max_dt;%log(0.9);%
    model.f = randn()*TD_parameters.f_range(2) + TD_parameters.f_range(1);%rand(n,1)*TD_parameters.max_dt;%log(0.5);

    model   = build_C(model, TD_parameters.t);
    model.T = 1;%tempering not enabled
    model   = evaluate(model, allWfs, TD_parameters, 1:nevt, 1:nsta, []);

    %this quantity is used a lot
    nwavelet = numel(model.wavelet);%n in wavelets
    %nevt is npolarization
    %nevt*nsta is ndts
    nsplit = numel(model.A);%number of splitting parameters. Could be clustered
    %nsta is orientations
    %nevt*nsta + 2 is error

    %model.nparam = nwavelet + nevt + nevt*nsta + nsplit*3 + nsta + nevt*nsta + 2;
    model.nparam = nwavelet + nevt + nevt*nsta + nsta + nevt*nsta + 2;

end

% n = floor(length(TD_parameters.t)/3);
% 
% for k = 1:(nevt*nsta)
% 
%     %model.amp(k) = log(max([ allWfs(k).north; allWfs(k).east ]));%peak is real peak, scaled down to damp gradient
%     %model.sig(k) = log(std([ allWfs(k).north(1:n); allWfs(k).east(1:n) ])/4);%~quarter variation is noise
% 
% end
% 
% model.amp = reshape(model.amp, size(allWfs));
% %model.sig = reshape(model.sig, size(allWfs));
