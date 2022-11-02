function model = vectorize_model(model, Parameters)
%non-dimensionalized model vector 

    model.vector = [];

    [~, nsta] = size(model.sig);

    %source parameters
    if Parameters.wavelet && nsta > 1

        model.vector = [ model.vector; model.wavelet(:)   ];
    
        if Parameters.use_polarization

            model.vector = [ model.vector; model.polarization/Parameters.polarization_std ];

        end

        %evt-station pairs
        model.vector = [ model.vector; model.amp(:)   ];%std of 1
        model.vector = [ model.vector; model.shift(:)/Parameters.Delaystd ];
        model.vector = [ model.vector; (model.dtS(:) - Parameters.dtS(1))/Parameters.dtS(2)   ];

    end

    model.vector = [ model.vector; model.A(:)/Parameters.dtstd ];%ignore the sqrt(2) you need here. It's just nondimensionalization
    model.vector = [ model.vector; model.B(:)/Parameters.dtstd ];

    model.vector = [ model.vector; model.fast_dir_rotation(:)/Parameters.rotation_std ];

    if Parameters.use_tSA

        model.vector = [ model.vector; model.tSA(:)/Parameters.tSAstd ];

    end

    if Parameters.use_orientations && nsta > 1 && Parameters.wavelet

        model.vector = [ model.vector; model.sta_or./Parameters.sta_err ];

    end

    %error terms
    model.vector = [ model.vector; (model.sig(:) - Parameters.sig_range(1))/Parameters.sig_range(2) ];
     
    if Parameters.use_covarience

        model.vector = [ model.vector; (model.r - Parameters.r_range(1))/Parameters.r_range(2) ];
        model.vector = [ model.vector; (model.f - Parameters.f_range(1))/Parameters.f_range(2) ];

    end

    model.nparam = length(model.vector);

end