function [ N, E, model ] = apply_operator_multi(model, f, evt_ind, sta_ind, prelogf, cluster)

    if nargin < 6

        cluster = false;

    end

    type    = 'all';

    f = -2*pi*f; %needs to be negative for the attenuation operator

    usewavelet = isfield(model, 'wavelet');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % attenuation and delay operators
     if nargin==2       % set w1 AND offset to default values
         w1=400;
         delay = 0;
     end
     
     if nargin==3       % set w1 OR  offset to default values
       if w1>10 % case where w1 is w1
         delay=0;
       else % case where w1 is offset
         delay=w1;
         w1=400*2*pi;
       end
     end
        
    if nargin < 5
        w1      = 400*2*pi;
        prelogf = log([1; f(2:end)]/w1);%This line is slow, precalculated
    
    end
    
    if usewavelet

        w           = model.wavelet(evt_ind, :)';

    else

        w                     = zeros(length(f), 1);
        w(round(length(f)/2)) = 1;%impulse

    end

    wavelet     = fft((w));%could save but it's not too bad computationally
    
    if ~strcmp(type, 'all')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % splitting operator saved from befor

        if cluster

            irF        = model.irF{evt_ind, 1};
            irS        = model.irS{evt_ind, 1};
            ending_phi = model.endphi(evt_ind, 1);

        else
    
            irF        = model.irF{evt_ind, sta_ind};
            irS        = model.irS{evt_ind, sta_ind};
            ending_phi = model.endphi(evt_ind, sta_ind);

        end

    else

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % splitting operator also in F domain
        % first getting the splitting parameters and determine if this is one
        % or two layer splitting

        if usewavelet

            polarization = model.polarization(evt_ind) + model.polarization0(evt_ind);

        else

            polarization = model.polarization0(evt_ind);

        end

        if cluster

            phi          = model.fast_dir(1, :);
            dphi         = model.fast_dir_rotation(1, :);
            delay        = model.dt(1, :);
            tSA          = model.tSA(1, :);

        else

            phi          = model.fast_dir(sta_ind, :);
            dphi         = model.fast_dir_rotation(sta_ind, :);
            delay        = model.dt(sta_ind, :);
            tSA          = model.tSA(sta_ind, :);

        end
            
        f = -1*f;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %This calculates the impulse response for up to two rotating layers. Could be
        %generalized to N layers but idk sks phases aren't that great. 
    
        [irF1, irS1]   = ir_rotating_layer(delay(1), phi(1), dphi(1), tSA(1), f, polarization, prelogf);%bottom layer
        
        if sum(delay) > delay(1)
    
            %impulse response for fast polarized waves
            [irF_F, irS_F] = ir_rotating_layer(delay(2), phi(2), dphi(2), tSA(2), f, phi(1) + dphi(1), prelogf);%top layer
        
            %impulse response for slow polarized waves
            [irF_S, irS_S] = ir_rotating_layer(delay(2), phi(2), dphi(2), tSA(2), f, phi(1) + dphi(1) + pi/2, prelogf);%top layer
                
            ff = real(ifft(fft((irF_F)).*fft(irF1)));
            fs = real(ifft(fft((irF_S)).*fft(irS1)));
            sf = real(ifft(fft((irS_F)).*fft(irF1)));
            ss = real(ifft(fft((irS_S)).*fft(irS1)));
        
            irF = ff + fs;
            irS = ss + sf;
    
            ending_phi = phi(2) + dphi(2);
    
        else
        
            irF = irF1;
            irS = irS1;
            
            ending_phi = phi(1) + dphi(1);
    
        end

        if cluster

            model.irF{evt_ind,    1} = irF;
            model.irS{evt_ind,    1} = irS;
            model.endphi(evt_ind, 1) = ending_phi;

        else

            model.irF{evt_ind,    sta_ind} = irF;
            model.irS{evt_ind,    sta_ind} = irS;
            model.endphi(evt_ind, sta_ind) = ending_phi;

        end

    end

    F = real(ifft(fft(irF).*wavelet));
    S = real(ifft(fft(irS).*wavelet));

    N = F*cos(ending_phi) - S*sin(ending_phi);%to north
    E = F*sin(ending_phi) + S*cos(ending_phi);%to east

end


    %wavelet_amp = max(abs(w));%to undo the attenuation's amplitude effect
    %shift       = model.shift(evt_ind, sta_ind);
    %tS          = exp(model.dtS(evt_ind, sta_ind));
    %scalar_amp  = exp(model.amp(evt_ind, sta_ind));

    %-------------------------------------------------------------------
    %Amp = exp(-f*tS/2);                 % amplitude spectrum
    %%%%need to normalize amp
    %Dt  = -(tS/pi)*prelogf; % dummy F=0, OK on next line
    %Dt  = -(tStar/pi)*log(f/w1); %Faster, but this line is still about 0.25 of the entire program
    %Dt(1) = -(tStar/pi)*log(1/w1);
    %Dt  = Dt - min(Dt) + shift;
    %Dt  = Dt - min(Dt);
    %Dph = f.*(Dt);                   % delta phase shift spectrum
    % attn opperator in freq with a delay in time
    %attnOperator = Amp.*(cos(Dph) + 1i*sin(Dph));%This line is about 1/10th of the entire program
    %attnOperator = Amp.*exp(1i*Dph);%operator for both the fast and slow directions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %pull it out so that at some point I can save them
    %F = real(ifft(wavelet.*attnOperator));
    %S = real(ifft(wavelet.*attnOperator));
    %F = real(ifft(fft(irF).*wavelet.*attnOperator));
    %S = real(ifft(fft(irS).*wavelet.*attnOperator));
    %F = real(ifft(wavelet.*attnOperator));

