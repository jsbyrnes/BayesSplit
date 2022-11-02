function mout = fit_waveforms(allWfs, Parameters)
%synmodel only used in debugging

    [nevt, nsta] = size(allWfs);

    disp('-> Building a starting model')
    model = make_starting_models_GA(Parameters, allWfs);

    model = vectorize_model(model, Parameters);

    if Parameters.solver_printout

        options                  = optimoptions('fminunc','Display','iter','Algorithm','quasi-newton',...
            'MaxIterations', 1e9, 'OptimalityTolerance', 1e-4, ...
            'MaxFunctionEvaluations', 1e9, 'SpecifyObjectiveGradient', true);
     
    else

        options                  = optimoptions('fminunc','Display','none','Algorithm','quasi-newton',...
            'MaxIterations', 1e9, 'OptimalityTolerance', 1e-4, ...
            'MaxFunctionEvaluations', 1e9, 'SpecifyObjectiveGradient', true);
        
    end

    %first attempt at fitting
    disp('-> Fitting a first model')

    f                        = @(x) fitwaveform_wrapper(x, model, allWfs, Parameters, -1);
    [x, ~, ~, fmin_out, del] = fminunc(f, model.vector, options);

    model.vector   = x;
    model          = devectorize_model(model, nevt, nsta, length(Parameters.t), Parameters);
    model          = evaluate(model, allWfs, Parameters, 1:nevt, 1:nsta, []);
    model.del      = del;
    model.fmin_out = fmin_out;

    disp([ '--> First solution at posterior density of ' num2str(model.lpst) ]);

    iter         = 1;
    rounds       = 0;
    reset_thresh = 0.1;%not important to fine tune, threshold for when you don't count as an improvement

    modelr = model;%for the model in this round

    while rounds < Parameters.reset_rounds
        
        if iter == 0

            modelr = make_starting_models_GA(Parameters, allWfs);
            modelr = vectorize_model(modelr, Parameters);
        
            iter = 1;

        end
    
        sig           = 10^(diff(Parameters.reset_size)*(iter/Parameters.reset_iter)...
            + Parameters.reset_size(1));
        modeln        = modelr;
        modeln.vector = modelr.vector + sig*randn(size(modelr.vector));
        modeln        = devectorize_model(modeln, nevt, nsta, length(Parameters.t), Parameters);

        modeln        = evaluate(modeln, allWfs, Parameters, 1:nevt, 1:nsta, []);
        disp([' -> Attempting a trial starting at posterior density of ' num2str(modeln.lpst) ]);
        
        f                        = @(x) fitwaveform_wrapper(x, modeln, allWfs, Parameters, -1);
        [x, ~, ~, fmin_out, del] = fminunc(f, modeln.vector, options);

        modeln.vector   = x;
        modeln          = devectorize_model(modeln, nevt, nsta, length(Parameters.t), Parameters);
        modeln          = evaluate(modeln, allWfs, Parameters, 1:nevt, 1:nsta, []);

        if (modeln.lpst - modelr.lpst) > reset_thresh

            disp([' -> Better trial solution found with posterior density of ' num2str(modeln.lpst) ]);
            modelr          = modeln;
            modelr.del      = del;
            modelr.fmin_out = fmin_out;

            iter = max([1 (iter - 1)]);

        elseif reset_thresh > (modeln.lpst - modelr.lpst) && (modeln.lpst - modelr.lpst) >= 0

            disp( '    Solution found within search tolerence');
            iter = iter + 1;
            disp(['    Counter raised to ' num2str(iter) ' of ' num2str(Parameters.reset_iter) ])

        else

            disp(['    Test solution rejected at posterior density of ' num2str(modeln.lpst) '']);                
            iter = iter + 1;
            disp(['    Counter raised to ' num2str(iter) ' of ' num2str(Parameters.reset_iter) ])

        end

        if iter == Parameters.reset_iter

            if (modelr.lpst - model.lpst) > reset_thresh

                disp(['--> Better solution over a round found with posterior density of ' num2str(modelr.lpst) ]);
                disp( '--> **Round counter set to zero**');
                
                model = modelr;

                rounds = 0;

            else

                rounds = rounds + 1;

                disp('    ')
                disp([ '--> Round counter raised to ' num2str(rounds)])

            end

            modelr = model;%and then reset it
            iter   = 0;%start over

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%
    %HMC search for errors

    if Parameters.get_errors

        f   = @(x) fitwaveform_wrapper(x, model, allWfs, Parameters, 1);
        smp = hmcSampler(f, model.vector, 'NumSteps', 5, 'MassVectorTuningMethod', 'hessian', 'StepSize', 1e-1);
    
        c = [];
    
        disp('-> Starting HMC search for errors');
        if Parameters.solver_printout
    
            smp          = tuneSampler(smp, 'VerbosityLevel', 2, 'NumPrint',1, 'NumStepSizeTuningIterations', 50, 'NumStepsLimit', 5);
            smp.NumSteps = 5;
    
            for k = 1:Parameters.nchains
    
                disp([ '     Chain #' num2str(k) ' of ' num2str(Parameters.nchains) ]);
                
                tmp = drawSamples(smp, 'Burnin', Parameters.burnin, 'NumSamples', Parameters.batch_size, ...
                    'ThinSize', Parameters.batch_thin ...
                , 'Verbosity', 1, 'NumPrint', 1, 'StartPoint', model.vector);
    
                c = [ c; tmp ];
    
            end
    
        else
    
            smp          = tuneSampler(smp, 'NumStepSizeTuningIterations', 50, 'NumStepsLimit', 5);
            smp.NumSteps = 5;
    
            for k = 1:Parameters.nchains
    
                disp([ '     Chain #' num2str(k) ' of ' num2str(Parameters.nchains) ]);
    
                tmp = drawSamples(smp, 'Burnin', Parameters.burnin, 'NumSamples', Parameters.batch_size, ...
                    'ThinSize', Parameters.batch_thin, 'StartPoint', model.vector);
    
                c = [ c; tmp ];
    
            end
    
        end
    
        [n, ~] = size(c);
    
        disp('   Building models from HMC search');
        mout = model;%save the optimimum
        for k = 2:n+1
    
            m        = model;
            m.vector = c(k-1, :)';
            m        = devectorize_model(m, nevt, nsta, length(Parameters.t), Parameters);
            mout(k)  = evaluate(m, allWfs, Parameters, 1:nevt, 1:nsta, []);
    
        end

    else

        mout = model;

    end
       
end