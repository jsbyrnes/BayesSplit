function [model, N, E, dN, dE] = apply_model(model, allWfs, Parameters, evt_ind, sta_ind, prelogf)
%the last four optional outputs are really only useful when evt_ind and
%sta_ind are scalars

    f = Parameters.f;

    if nargin < 6

        w1      = 400*2*pi;
        prelogf = log([1; -2*pi*f(2:end)]/w1);%This line is slow, precalculated

    end

    [~, nsta] = size(allWfs);

    for k = evt_ind

        for kk = sta_ind

            if isnan(allWfs(k, kk).north)

                model.phi_vec(k, kk) = 0;
                model.sig(k, kk)     = 0;

            else
                    
                %rotate the data
                dN = cos(model.sta_or(kk))*allWfs(k, kk).north ...
                    - sin(model.sta_or(kk))*allWfs(k, kk).east;
                dE = sin(model.sta_or(kk))*allWfs(k, kk).north ...
                    + cos(model.sta_or(kk))*allWfs(k, kk).east;

                [No, Eo, model] = apply_operator_multi(model, Parameters.f, k, kk, prelogf, Parameters.cluster);

                if any(isnan(No)) || any(isnan(Eo))

                    model.phi_vec(k, kk) = 1e9; %aphysical model triggers nans
                    continue

                end
                
                if nsta > 1 && Parameters.wavelet

                    N = new_waveform(No, model.shift(k,kk), model.amp(k,kk), model.dtS(k,kk), f, prelogf);
                    E = new_waveform(Eo, model.shift(k,kk), model.amp(k,kk), model.dtS(k,kk), f, prelogf);

                else

                    amp = rms([  No; Eo ]);
                    N   = No/amp;
                    E   = Eo/amp;

                end
                
                model.N{k,kk} = N;
                model.E{k,kk} = E;

                if nsta > 1

                    model.phi_vec(k, kk) = ((N - dN)'*(model.Cinv/...
                        (exp(model.sig(k, kk))^2))*(N - dN) ...
                        + (E - dE)'*(model.Cinv/(exp(model.sig(k, kk))^2))...
                        *(E - dE));

                else

                    %cross convolution with difference
                    r = ifft(fft(N).*fft(dE)) - ifft(fft(E).*fft(dN));
                    model.phi_vec(k, kk) = r'*(model.Cinv/(exp(model.sig(k, kk))^2))*r;

                end

            end

        end

    end

    model.phi     = sum(model.phi_vec(:));

end

%%%%debugging plot
% subplot(311), sv = -5:0.01:5;
% for q = 1:length(sv), lvec(q) = lpst(new_waveform(No, sv(q), amp, tS, f, prelogf), new_waveform(Eo, sv(q), amp, tS, f, prelogf), sv(q), amp, tS);end
% plot(sv, lvec)
% xlabel('shift')
% subplot(312)
% a = -5:0.01:1;lvec=[];
% for q = 1:length(a), lvec(q) = lpst(exp(a(q))*No, exp(a(q))*Eo, shift, a(q), tS);end
% plot(a, lvec)
% xlabel('amp')
% subplot(313)
% tsv = -3:0.01:3;lvec = [];
% for q = 1:length(tsv), lvec(q) = lpst(new_waveform(No, shift, amp, tsv(q), f, prelogf), new_waveform(Eo, shift, amp, tsv(q), f, prelogf), shift, amp, tsv(q));end
% plot(tsv, lvec)
% xlabel('tS')

% 
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %align
%                 %[xc, lags] = xcorr([No; Eo], [dN; dE]);
%                 %[~, idx] = max(xc);
%                 %N = circshift(No, -lags(idx));
%                 %E = circshift(Eo, -lags(idx));
% 
%                 %This is  more computationally intensive but can be worth it
%                 lags = 1:length(dN);
%                 for j = 1:length(lags)
% 
%                     phi(j) = md(circshift(No, lags(j)), circshift(Eo, lags(j)));
% 
%                 end
% 
%                 mphi = min(phi);
% 
%                 [~, peaks_ind] = sort(phi - mphi);
% 
%                 idx_vector     = peaks_ind(1:4);
%                 delay_vec      = [];
%                 phi0_vec       = [];
% 
%                 for idx = 1:length(idx_vector)
%     
%                     limited_step = false;
% 
%                     N = circshift(No, lags(idx_vector(idx)));
%                     E = circshift(Eo, lags(idx_vector(idx)));
%                     %test at close points
%                     delay = 0;%starting guess
% 
%                     %weak grid search to start. Probably a saddle
%                     %point
% %                     x   = linspace(-2,2,10);
% %                     phi = [];
% %                     for q = 1:length(x)
% %                         phi(q) = md(delay_continuous(N, Parameters.sample_rate, x(q)/Parameters.sample_rate)...
% %                             , delay_continuous(E, Parameters.sample_rate, x(q)/Parameters.sample_rate));
% %                     end
% %
% %                    [~, id] = min(phi);
% %                    delay = x(id);
% 
%                     %delay = 0;
%                     h     = 5e-5;
%                     diff  = Inf;
%                     niter = 0;
%     
%                     while diff > h
%     
%                         phi0 = md(delay_continuous(N, Parameters.f, delay/Parameters.sample_rate)...
%                             , delay_continuous(E, Parameters.f, delay/Parameters.sample_rate));
%                         phip = md(delay_continuous(N, Parameters.f, (delay+h)/Parameters.sample_rate)...
%                             , delay_continuous(E, Parameters.f, (delay+h)/Parameters.sample_rate));
%                         phin = md(delay_continuous(N, Parameters.f, (delay-h)/Parameters.sample_rate)...
%                             , delay_continuous(E, Parameters.f, (delay-h)/Parameters.sample_rate));
%     
%                         d1 = (phip - phin)/(2*h);
%                         d2 = (phip + phin - 2*phi0)/(h^2);
% 
%                         step  = d1/d2;
% 
%                         if d2 > 0 && abs(step) <= 1 %blows up when second derivative is small
%                             
%                             delayn = delay - step;
%                             diff   = abs(delayn - delay);
%                             delay  = delayn;
%     
%                         elseif niter ~= 10 %newton's method won't work! psuedo grid search to reset the delay
% 
%                             phi0 = md(delay_continuous(N, Parameters.f, delay/Parameters.sample_rate)...
%                                 , delay_continuous(E, Parameters.f, delay/Parameters.sample_rate));
% 
%                             alpha = logspace(-3, -1, 10);
% 
%                             minphi = Inf;
%                             while minphi > phi0
% 
%     %                             delay_v = rand(10,1)*6 - 3;
%     %                             for q = 1:length(delay_v)
%     %                                 phi_search(q) = md(delay_continuous(N, Parameters.sample_rate, delay_v(q)/Parameters.sample_rate)...
%     %                                 , delay_continuous(E, Parameters.sample_rate, delay_v(q)/Parameters.sample_rate));
%     %                             end
%     %     
%     %                             [~, idx] = min(phi_search);
%     %     
%     %                             delay = delay_v(idx);
%         
%                                 %line search gradient descent
%                                 step_v = sign(d1)*alpha;
% %                                 step_v(step_v > 1 ) = [];
% %                                 step_v(step_v < -1) = [];
%                                 delay_v             = delay - step_v;
%         
%                                 phi_search          = [];
%                                 for q = 1:length(delay_v)
%                                     phi_search(q) = md(delay_continuous(N, Parameters.f, delay_v(q)/Parameters.sample_rate)...
%                                     , delay_continuous(E, Parameters.f, delay_v(q)/Parameters.sample_rate));
%                                 end
%     
%                                 [minphi, ii] = min(phi_search);
% 
%                                 if minphi < phi0
% 
%                                     diff  = abs(delay - delay_v(ii));
%                                     delay = delay_v(ii);%otherwise, marches into shit territory
% 
%                                 end
%             
%                                 alpha = alpha/5;
% 
%                                 if any(alpha < 1e-12)
% 
%                                     %pathological, so give up and grid
%                                     %search
%                                     x = -3:h:3;
%                                     phi_search = [];
%                                     for q = 1:length(x)
%                                         phi_search(q) = md(delay_continuous(N, Parameters.f, x(q)/Parameters.sample_rate)...
%                                             , delay_continuous(E, Parameters.f, x(q)/Parameters.sample_rate)); 
%                                     end
% 
%                                     [~, id] = min(phi_search);
%                                     delay = x(id);
%                                     break
% 
%                                 end
% 
%                             end
%     
%                         else
% 
%                             x = -2:(h*5):2;
%                             phi_search = [];
%                             for q = 1:length(x)
%                                 phi_search(q) = md(delay_continuous(N, Parameters.f, x(q)/Parameters.sample_rate)...
%                                     , delay_continuous(E, Parameters.f, x(q)/Parameters.sample_rate)); 
%                             end
%     
%                             [~, id] = min(phi_search);
%                             delay = x(id);
% 
%                         end
% 
%                         niter  = niter + 1;
% 
%                         if niter > 50
% 
%                             %pathological, so give up and grid
%                             %search
%                             x = -3:h:3;
%                             phi_search = [];
%                             for q = 1:length(x)
%                                 phi_search(q) = md(delay_continuous(N, Parameters.f, x(q)/Parameters.sample_rate)...
%                                     , delay_continuous(E, Parameters.f, x(q)/Parameters.sample_rate)); 
%                             end
% 
%                             [~, id] = min(phi_search);
%                             delay = x(id);
%                             break
%     
%                         end
% 
%                     end
% 
%                     phi0 = md(delay_continuous(N, Parameters.f, delay/Parameters.sample_rate)...
%                         , delay_continuous(E, Parameters.f, delay/Parameters.sample_rate));
%                     delay_vec(idx) = delay;
%                     phi0_vec(idx)  = phi0;
% 
%                 end
% 
%                 %kept the best one
%                 [~, idx] = min(phi0_vec);
% 
%                 N = circshift(No, lags(idx_vector(idx)));
%                 E = circshift(Eo, lags(idx_vector(idx)));
%                 N = delay_continuous(N, Parameters.f, delay_vec(idx)/Parameters.sample_rate);
%                 E = delay_continuous(E, Parameters.f, delay_vec(idx)/Parameters.sample_rate);
%                 
%                 %%%%%%%%%%%%%%%%%%%%%%%%
%                 %scale
%                 amp = 0;%starting guess
% 
%                 h = 5e-5;
% 
%                 diff = Inf;
%                 
%                 niter = 0;
% 
%                 %FD turns out to be efficient relative to analytical
%                 while diff > h
% 
%                     phi0 = md(exp(amp)*N, exp(amp)*E);
%                     phip = md(exp(amp + h)*N, exp(amp + h)*E);
%                     phin = md(exp(amp - h)*N, exp(amp - h)*E);
% 
%                     d1 = (phip - phin)/(2*h);
%                     d2 = (phip + phin - 2*phi0)/(h^2);
% 
%                     ampn = amp - d1/d2;
% 
%                     diff = abs(ampn - amp);
% 
%                     amp = ampn;
% 
%                     niter = niter + 1;
% 
%                     if amp < -4 || niter > 50
% 
%                        break
% 
%                     end
% 
%                 end
% 
%                 N             = exp(amp)*N;
%                 E             = exp(amp)*E;

%                 lpst = @(N,E, shift, amp, tS) -((N - dN)'*(model.Cinv/...
%                     (exp(model.sig(k, kk))^2))*(N - dN) ...
%                     + (E - dE)'*(model.Cinv/(exp(model.sig(k, kk))^2))...
%                     *(E - dE))/2 - (shift^2)/(2*Parameters.delay_std^2) ...
%                      - (amp^2)/(2) - ((tS - Parameters.dtS(1))^2)/(2*Parameters.dtS(2)^2);%other terms persist in function handle. -2 gives fit
%                 shift_vec = round(model.shift(k, kk)*ones(3,1)*...
%                     -1*Parameters.sample_rate);
% 
%                 %find up one cycle
%                 lp0 = lpst(No, Eo, 0, 0, 0);
% 
%                 peak1 = false;
%                 peak2 = false;
%                 shift = 0;
%                 start = false;
%                 while ~(peak1 && peak2)
% 
%                     %%%%% two possibilities. One is that you need cross the
%                     %%%%% optimum, so don't flag peaks yey
% 
%                     shift = shift + 1;
%                     llhp = lpst(circshift(No, shift), circshift(Eo, shift), 0, 0, 0);
% 
%                     if llhp > lp0 && ~start
% 
%                         lp0 = llhp;
%                         continue
% 
%                     end
% 
%                     %will trigger once you have climbed up at least once
%                     start = true;
% 
%                     if llhp > lp0
% 
%                         peak1 = true;
%                         lp0 = llhp;
%                         continue
% 
%                     end
% 
%                     if llhp < lp0 && peak1
% 
%                         peak2 = true;
%                         lp0 = llhp;
%                         continue
% 
%                     end
% 
%                     lp0 = llhp;
% 
%                 end
% 
%                 shift_vec(2) = shift-1;
% 
%                 %find up down cycle
%                 lp0 = lpst(No, Eo, 0, 0, 0);
% 
%                 peak1 = false;
%                 peak2 = false;
%                 shift = 0;
%                 start = false;
%                 while ~(peak1 && peak2)
% 
%                     %%%%% two possibilities. One is that you need cross the
%                     %%%%% optimum, so don't flag peaks yey
% 
%                     shift = shift - 1;
%                     llhp = lpst(circshift(No, shift), circshift(Eo, shift), 0, 0, 0);
% 
%                     if llhp > lp0 && ~start
% 
%                         lp0 = llhp;
%                         continue
% 
%                     end
% 
%                     %will trigger once you have climbed up at least once
%                     start = true;
% 
%                     if llhp > lp0
% 
%                         peak1 = true;
%                         lp0 = llhp;
%                         continue
% 
%                     end
% 
%                     if llhp < lp0 && peak1
% 
%                         peak2 = true;
%                         lp0 = llhp;
%                         continue
% 
%                     end
% 
%                     lp0 = llhp;
% 
%                 end
% 
%                 shift_vec(3) = shift+1;
%                 shift_vec    = -1*shift_vec/Parameters.sample_rate;
%                 shift_vec(abs(shift_vec) > 6*Parameters.delay_std) = [];
% 
%                 for i = 1:length(shift_vec)%do search from all three staring points
%                 
%                     shift = shift_vec(i);
%                     amp   = model.amp(k, kk);
%                     tS    = model.dtS(k, kk);%Parameters.dtS(1);
%                     %h     = 1e-6;
%                     %diff  = Inf;
%                     %niter = 0;
% 
%                     func = @(x) -lpst(new_waveform(No, x(1), x(2), x(3), f, prelogf), ...
%                                     new_waveform(Eo, x(1), x(2), x(3), f, prelogf), x(1), x(2), x(3));
%                 
%                     [x,fval,exitflag,output] = fminunc(func, [shift amp tS], options);
% 
%                     shift = x(1);
%                     amp   = x(2);
%                     tS    = x(3);

%                     keyboard
% 
%                     while diff > 1e-6
%     
%                         niter = niter + 1;
%                         if niter == 10
% 
%                             lims = (3*Parameters.delay_std);
%                             newshift = -lims:0.5:lims;
% 
%                             lpsearch = zeros(size(newshift));
%                             for q = 1:length(newshift)
%                                 lpsearch(q) = lpst(new_waveform(No, newshift(q), amp, tS, f, prelogf), ...
%                                     new_waveform(Eo, newshift(q), amp, tS, f, prelogf), newshift(q), amp, tS);
%                             end
% 
%                             [~, ii] = max(lpsearch);
%                             shift = newshift(ii);
% 
%                             continue
%                                         
%                         end
% 
%                         %three searches. Zero delay, +1 and -1 cycle
%     
%                         %notation. p is positive, n is negative, s is
%                         %shift, a is amp, t is t*
%                         lp0 = lpst(new_waveform(No, shift, amp, tS, f, prelogf), new_waveform(Eo, shift, amp, tS, f, prelogf), shift, amp, tS);
% 
%                         ps = lpst(new_waveform(No, shift + h, amp, tS, f, prelogf), new_waveform(Eo, shift + h, amp, tS, f, prelogf), shift + h, amp, tS);
%                         ns = lpst(new_waveform(No, shift - h, amp, tS, f, prelogf), new_waveform(Eo, shift - h, amp, tS, f, prelogf), shift - h, amp, tS);
%                         pa = lpst(new_waveform(No, shift, amp + h, tS, f, prelogf), new_waveform(Eo, shift, amp + h, tS, f, prelogf), shift, amp + h, tS);
%                         na = lpst(new_waveform(No, shift, amp - h, tS, f, prelogf), new_waveform(Eo, shift, amp - h, tS, f, prelogf), shift, amp - h, tS);
%                         pt = lpst(new_waveform(No, shift, amp, tS + h, f, prelogf), new_waveform(Eo, shift, amp, tS + h, f, prelogf), shift, amp, tS + h);
%                         nt = lpst(new_waveform(No, shift, amp, tS - h, f, prelogf), new_waveform(Eo, shift, amp, tS - h, f, prelogf), shift, amp, tS - h);
% 
%                         pspa = lpst(new_waveform(No, shift + h, amp + h, tS, f, prelogf), new_waveform(Eo, shift + h, amp + h, tS, f, prelogf), shift + h, amp + h, tS);
%                         nspa = lpst(new_waveform(No, shift - h, amp + h, tS, f, prelogf), new_waveform(Eo, shift - h, amp + h, tS, f, prelogf), shift - h, amp + h, tS);
%                         psna = lpst(new_waveform(No, shift + h, amp - h, tS, f, prelogf), new_waveform(Eo, shift + h, amp - h, tS, f, prelogf), shift + h, amp - h, tS);
%                         nsna = lpst(new_waveform(No, shift - h, amp - h, tS, f, prelogf), new_waveform(Eo, shift - h, amp - h, tS, f, prelogf), shift - h, amp - h, tS);
% 
%                         pspt = lpst(new_waveform(No, shift + h, amp, tS + h, f, prelogf), new_waveform(Eo, shift + h, amp, tS + h, f, prelogf), shift + h, amp, tS + h);
%                         nspt = lpst(new_waveform(No, shift - h, amp, tS + h, f, prelogf), new_waveform(Eo, shift - h, amp, tS + h, f, prelogf), shift - h, amp, tS + h);
%                         psnt = lpst(new_waveform(No, shift + h, amp, tS - h, f, prelogf), new_waveform(Eo, shift + h, amp, tS - h, f, prelogf), shift + h, amp, tS - h);
%                         nsnt = lpst(new_waveform(No, shift - h, amp, tS - h, f, prelogf), new_waveform(Eo, shift - h, amp, tS - h, f, prelogf), shift - h, amp, tS - h);
% 
%                         papt = lpst(new_waveform(No, shift, amp + h, tS + h, f, prelogf), new_waveform(Eo, shift, amp + h, tS + h, f, prelogf), shift, amp + h, tS + h);
%                         napt = lpst(new_waveform(No, shift, amp - h, tS + h, f, prelogf), new_waveform(Eo, shift, amp - h, tS + h, f, prelogf), shift, amp - h, tS + h);
%                         pant = lpst(new_waveform(No, shift, amp + h, tS - h, f, prelogf), new_waveform(Eo, shift, amp + h, tS - h, f, prelogf), shift, amp + h, tS - h);
%                         nant = lpst(new_waveform(No, shift, amp - h, tS - h, f, prelogf), new_waveform(Eo, shift, amp - h, tS - h, f, prelogf), shift, amp - h, tS - h);
% 
%                         ds = (ps - ns)/(2*h);
%                         da = (pa - na)/(2*h);
%                         dt = (pt - nt)/(2*h);
%     
%                         dss = (ps + ns - 2*lp0)/h^2;
%                         daa = (pa + na - 2*lp0)/h^2;
%                         dtt = (pt + nt - 2*lp0)/h^2;
%                 
%                         dsa = (pspa + nsna - nspa - psna)/(4*h^2);
%                         dst = (pspt + nsnt - nspt - psnt)/(4*h^2);
%                         dat = (papt + nant - napt - pant)/(4*h^2);
% 
%                         G = [ ds; da; dt ];
%                         H = [ dss dsa dst; dsa daa dat; dst dat dtt ];
% 
%                         step = H\G;
% 
%                         if ~any(abs(step) > 1) && ~any([dss daa dtt] > 0)% && any(sign(G) ~= sign(-step))%checks for quality of newton's search
%                             %flipped sign can be ok, but is also a sign of
%                             %weirdness
%                             shift0 = shift;
%                             amp0   = amp;
%                             tS0    = tS;
% 
%                             shift = shift - step(1);
%                             amp   = amp   - step(2);
%                             tS    = tS    - step(3);
% 
%                             lp   = lpst(new_waveform(No, shift, amp, tS, f, prelogf), ...
%                                 new_waveform(Eo, shift, amp, tS, f, prelogf), shift, amp, tS);
% 
%                             diff_lp = lp - lp0;
%                         
%                             diff = max(abs(step));
% 
%                             if diff_lp > 0
% 
%                                 continue %skip the grid search
% 
%                             end
% 
%                             %you missed so reset and either grid search or
%                             %gradient ascent search
%                             shift = shift0;
%                             amp   = amp0;
%                             tS    = tS0;
% 
%                         end
% 
%                         alpha = min([ max(abs(G)) 1]);
% 
%                         lp = -Inf;
%                         while lp < lp0
%     
%                             %line search gradient descent
%                             step  = sign(G).*abs(G).*alpha*(rand()*0.5 + 0.5);%rand keeps it from oscilatting
%     
%                             step(abs(step)>1) = 1;
%                             step              = abs(step).*sign(G);
% 
%                             lp = lpst(new_waveform(No, shift + step(1), amp + step(2), tS + step(3), f, prelogf), ...
%                                 new_waveform(Eo, shift + step(1), amp + step(2), tS + step(3), f, prelogf), ...
%                                 shift + step(1), amp + step(2), tS + step(3));
% 
%                             if lp > lp0
% 
%                                 diff  = max(abs(step));
%                                 shift = shift + step(1);
%                                 amp   = amp   + step(2);
%                                 tS    = tS    + step(3);
% 
%                             end
%         
%                             alpha = alpha/2; %always ends up finding something. Could become infinite in theory
% 
%                         end
%                        
%                         if niter > 100
%                             %keyboard
%                         end
% 
%                         if niter == 100
%                            disp(niter)
%                         end
%     
%                         if niter == 1000
%                            keyboard
%                         end
% 
%                     end
% 
%                     lp_vec(i) = lpst(new_waveform(No, shift, amp, tS, f, prelogf), new_waveform(Eo, shift, amp, tS, f, prelogf), shift, amp, tS);
% 
%                     shift_solution(i) = shift;
%                     amp_solution(i)   = amp;
%                     tS_solution(i)    = tS;
% 
%                 end
%                 %now you need to update sigma, a totally seperatable but
%                 %easy to solve for parameter. totally unimodal
%                 fit = ((N - dN)'*(model.Cinv)*(N - dN) ...
%                     + (E - dE)'*(model.Cinv)...
%                     *(E - dE));
%                 lpst = @(sig) -(fit/(exp(sig)^2) + (4*length(f)*sig))/2;%other terms persist in function handle. -2 gives fit
% 
%                 step = 1e9;
%                 h    = 1e-6;
%                 sig  = model.sig(k, kk);
%                 while step > 1e-6
% 
%                     lpst0 = lpst(sig);
%                     lpstp = lpst(sig + h);
%                     lpstn = lpst(sig - h);
% 
%                     d  = (lpstp - lpstn)/(2*h);
%                     dd = (lpstp + lpstn - 2*lpst0)/h^2;
% 
%                     step = -d/dd;
%                     if abs(step) > 1
% 
%                         step = sign(step);
% 
%                     end
% 
%                     sig = sig + step;
% 
%                 end
% 
%                 model.sig(k, kk)     = sig;
                %kept the best one
%                [~, idx] = max(lp_vec);
%                N = new_waveform(No, model.shift(k,kk), model.amp(k,kk), model.dtS(k,kk), f, prelogf);
%                E = new_waveform(Eo, model.shift(k,kk), model.amp(k,kk), model.dtS(k,kk), f, prelogf);

%                 model.amp(k,kk)      = amp_solution(idx);%log
%                 model.shift(k,kk)    = shift_solution(idx);
%                 model.dtS(k,kk)      = tS_solution(idx);%*Parameters.max_dtS;%ones(20, 1);% abs(randn(n,1)*Parameters.dt_std);

