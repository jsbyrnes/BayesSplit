function [lpst, grad] = fitwaveform_wrapper(vector, model, allWfs, Parameters, sign)

    [nevt, nsta] = size(allWfs);

    model.vector = vector;
    model        = devectorize_model(model, nevt, nsta, length(Parameters.t), Parameters);
    
    model = evaluate(model, allWfs, Parameters, 1:nevt, 1:nsta, []);
    lpst  = sign*model.lpst;%fminunc is minimization

    if nargout > 1 % gradient required

        model = fill_gradient_waveforms(model, allWfs, Parameters);
        grad = sign*model.del;

    end

end
