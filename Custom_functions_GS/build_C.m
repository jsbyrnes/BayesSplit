function [model, C] = build_C(model, t)

    C = zeros(length(t), length(t));
    
    f = exp(model.f);
    r = exp(model.r);

    for i = 1:length(t)
    
        for j = 1:length(t)
    
            dt     = abs(t(j) - t(i));
            C(i,j) = exp(-r*dt)*cos(2*pi*f*dt);
    
        end
    
    end

    try

        model.Cinv   = inv(C);
        model.logdet = logdet(C, 'chol');

    catch

        model.logdet = -1e9;%occasionally occurs, hard to predict. Reject the model 
    
    end

end