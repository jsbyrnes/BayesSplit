function model = devectorize_model(model, nevt, nsta, n, Parameters)

    %you act on the model vector, this fills it back out for evaluation
    %source parameters

    nl  = Parameters.max_layers;
    vec = model.vector;%will clip for simplicity

    if Parameters.wavelet && nsta > 1

        %source parameters
        model.wavelet = reshape(vec(1:(nevt*n)), [nevt n]);
        vec = vec((nevt*n + 1):end);

        if Parameters.use_polarization

            model.polarization = vec(1:nevt)*Parameters.polarization_std;
            vec = vec((nevt + 1):end);

        end

        %evt-station pairs
        model.amp = reshape(vec(1:nevt*nsta), [ nevt nsta]);
        vec = vec((nevt*nsta + 1):end);
        model.shift = reshape(vec(1:nevt*nsta), [ nevt nsta])*Parameters.Delaystd;
        vec = vec((nevt*nsta + 1):end);
        model.dtS = reshape(vec(1:nevt*nsta), [ nevt nsta])*Parameters.dtS(2) + Parameters.dtS(1);
        vec = vec((nevt*nsta + 1):end);

    end

    if Parameters.cluster

        model.A                 = reshape(vec(1:nl), [ 1 nl ])*Parameters.dtstd;
        vec                     = vec((nl + 1):end);
        model.B                 = reshape(vec(1:nl), [ 1 nl ])*Parameters.dtstd;
        vec                     = vec((nl + 1):end);

        model.fast_dir_rotation = reshape(vec(1:nl), [ 1 nl ])*Parameters.rotation_std;
        vec                     = vec((nl + 1):end);

        if Parameters.use_tSA
    
            model.tSA           = reshape(vec(1:nl), [ 1 nl ])*Parameters.tSAstd;
            vec                 = vec((nl + 1):end);
            
        end

    else

        model.A                 = reshape(vec(1:nsta*nl), [ nsta nl ])*Parameters.dtstd;
        vec                     = vec((nsta*nl + 1):end);
        model.B                 = reshape(vec(1:nsta*nl), [ nsta nl ])*Parameters.dtstd;
        vec                     = vec((nsta*nl + 1):end);

        model.fast_dir_rotation = reshape(vec(1:nsta*nl), [ nsta nl ])*Parameters.rotation_std;
        vec                     = vec((nsta*nl + 1):end);

        if Parameters.use_tSA
    
            model.tSA           = reshape(vec(1:nsta*nl), [ nsta nl ])*Parameters.tSAstd;
            vec                 = vec((nsta*nl + 1):end);
            
        end

    end

    model = make_dtphi(model);

    if Parameters.use_orientations && nsta > 1

        model.sta_or  = vec(1:nsta).*Parameters.sta_err;
        vec = vec((nsta + 1):end);

    end

    %error terms
    model.sig = reshape((vec(1:nevt*nsta)*Parameters.sig_range(2)) + Parameters.sig_range(1), [ nevt nsta ]);
    vec = vec((nevt*nsta + 1):end);

    if Parameters.use_covarience

        %model.r = min([ vec(1)*Parameters.r_range(2) + Parameters.r_range(1), log(1)]);
        %model.f = min([ vec(2)*Parameters.f_range(2) + Parameters.f_range(1), log(Parameters.sample_rate)]);
        model.r = vec(1)*Parameters.r_range(2) + Parameters.r_range(1);
        model.f = vec(2)*Parameters.f_range(2) + Parameters.f_range(1);

        %now apply parameters to useful model
        model = build_C(model, Parameters.t);

    else

        model.Cinv   = eye(n);
        model.logdet = 1;
        
    end

end