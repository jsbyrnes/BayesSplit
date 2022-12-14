%%%%%%%%%this script prints out all off the events so that you can select
%%%%%%%%%what you want to do measurements on

clear
close all

addpath('Custom_functions_GS');
addpath('wfTools');
addpath('CircStat2012a');
addpath('irisFetch')

dataName            = 'SSIP';
run_name            = 'SSIP_picks_small';
scale               = 0.1;
include_depthphases = true;
radial_or_traces    = 2;%1 is radial, 2 is traces
time_after = 50;

%picking notes for SSIP
%SKS
%4,13,38,84,90,97,99,112,126,128,142,152,165,171,195,217,248,268,271,282,284,299,325,332,378,419
%SKKS
%13,88,99,107,128,152,195,296,299,419

%sSKS
%248

%smaller set, all SKS
% 4    38   248   271   282

j                = 282;
target_phase_use = { 'SKS' };

%parameters are amp, delay, fast direction, splitting time, src time, src, amp
TD_parameters.sig                 = 0.1;%always non-dimensionalized
TD_parameters.maxDelay            = 5;%in s
TD_parameters.n_iter              = 2e5;
TD_parameters.burn_in             = 1.5e5;
TD_parameters.keep_each           = 2.5e2;
TD_parameters.n_chains            = 48;
TD_parameters.n_cold              = 8;  
TD_parameters.maxCells            = 20;
TD_parameters.interp_style        = 'pchip';
TD_parameters.sig_range           = [ log(0.01),  0    ];%log. Traces are normalized to rms, so 0 is max possible
TD_parameters.r_range             = [ log(0.05),  0.25 ];%mean, std, log
TD_parameters.f_range             = [ log(0.2),   0.25 ];%mean, std, log
TD_parameters.max_dt              = 4;
TD_parameters.max_dtS             = 2;
TD_parameters.polarization_std    = 5*(2*pi)/360;
TD_parameters.print_on            = 1;
TD_parameters.north_names         = { 'BH1', 'HH1', 'BHN', 'HHN' };
TD_parameters.east_names          = { 'BH2', 'HH2', 'BHE', 'HHE' };
TD_parameters.debug               = 0;
TD_parameters.maxT                = 50;
TD_parameters.orientation_std     = 1;

%parameters for loading in real data
TD_parameters.total_time  = 200;%in seconds
TD_parameters.low_pass    = 1/12;
TD_parameters.sample_rate = 2;%larger by at least 2
TD_parameters.high_pass   = 1/30;%3/(TD_parameters.total_time);%in seconds

TD_parameters.dataName = dataName;

t               = (0:1/TD_parameters.sample_rate:(TD_parameters.total_time))';
TD_parameters.t = t;

load([ './FetchData/' dataName '.mat' ]);

allWfs = load_data(TD_parameters, E(j), S, target_phase_use, TD_parameters.total_time/2);%([43,44,61])

if any(~isnan([allWfs(:).latitude]))

    [ del, az ] = distance([allWfs(:).latitude], [allWfs(:).longitude], E(j).PreferredLatitude,...
        E(j).PreferredLongitude);

    snrlist = [allWfs(:).snr];
    evt_snr = mean(snrlist(~isnan(snrlist)));

    if radial_or_traces == 1

        figure(1)
        clf
        hold on
    
        for k = 1:length(allWfs)
    
            if ~isnan(allWfs(1, k).north)
    
                plot(t, scale*sqrt(allWfs(1, k).north.^2 + allWfs(1, k).east.^2) + del(k), 'k')
        
                for kk = 1:length(allWfs(1, k).phase_list)
        
                    %check if its a depth phase
                    phasename = allWfs(1, k).phase_list{kk};
                    if phasename(1)=='s' || phasename(1)=='p'
    
                        if ~include_depthphases
    
                            continue
    
                        end
    
                    end
    
                    plot([ allWfs(1, k).phase_times(kk), allWfs(1, k).phase_times(kk) ], [ del(k), (del(k) + 1)], 'k--')
                    text(allWfs(1, k).phase_times(kk) + 0.5, del(k) + 1, allWfs(1, k).phase_list{kk})
        
                end
    
            end
    
        end
        xlim([0 max(t)])
        ylim([min(del)-1, max(del)+1])
        title('Radial Traces')

        figure(2)
        clf
        hold on
    
        for k = 1:length(allWfs)
    
            if ~isnan(allWfs(1, k).north)
    
                plot(t, scale*sqrt(allWfs(1, k).north.^2 + allWfs(1, k).east.^2) + k, 'k')
        
                for kk = 1:length(allWfs(1, k).phase_list)
        
                    %check if its a depth phase
                    phasename = allWfs(1, k).phase_list{kk};
                    if phasename(1)=='s' || phasename(1)=='p'
    
                        if ~include_depthphases
    
                            continue
    
                        end
    
                    end
    
                    plot([ allWfs(1, k).phase_times(kk), allWfs(1, k).phase_times(kk) ], [ k, (k+ 1)], 'k--')
                    text(allWfs(1, k).phase_times(kk) + 0.5, k + 1, allWfs(1, k).phase_list{kk})
        
                end
    
            end
    
        end
        xlim([0 max(t)])
        ylim([0, length(allWfs)+1])
    
        title('Radial traces')

    elseif radial_or_traces == 2

        figure(1)
        clf
        subplot(121)
        hold on
    
        phasenames_printed = false;

        for k = 1:length(allWfs)
    
            if ~isnan(allWfs(1, k).north)
    
                plot(t, scale*allWfs(1, k).north + del(k), 'k')
        
                for kk = 1:length(allWfs(1, k).phase_list)
        
                    %check if its a depth phase
                    phasename = allWfs(1, k).phase_list{kk};
                    if phasename(1)=='s' || phasename(1)=='p'
    
                        if ~include_depthphases
    
                            continue
    
                        end
    
                    end
    
                    plot([ allWfs(1, k).phase_times(kk), allWfs(1, k).phase_times(kk) ], [ del(k), (del(k) + 1)], 'k--')
        
                    if ~phasenames_printed
                    
                        text(allWfs(1, k).phase_times(kk) + 0.5, del(k) + 1, allWfs(1, k).phase_list{kk})

                    end

                end
    
                phasenames_printed = true;

            end
    
        end
        xlim([0 max(t)])
        ylim([min(del)-1, max(del)+1])
        title('North component')
    
        subplot(122)
        hold on
    
        phasenames_printed = false;
        
        for k = 1:length(allWfs)
    
            if ~isnan(allWfs(1, k).north)
    
                plot(t, scale*allWfs(1, k).east + del(k), 'k')
        
                for kk = 1:length(allWfs(1, k).phase_list)
        
                    %check if its a depth phase
                    phasename = allWfs(1, k).phase_list{kk};
                    if phasename(1)=='s' || phasename(1)=='p'
    
                        if ~include_depthphases
    
                            continue
    
                        end
    
                    end
    
                    plot([ allWfs(1, k).phase_times(kk), allWfs(1, k).phase_times(kk) ], [ del(k), (del(k) + 1)], 'k--')
    
                    if ~phasenames_printed
                    
                        text(allWfs(1, k).phase_times(kk) + 0.5, del(k) + 1, allWfs(1, k).phase_list{kk})

                    end

                end

                phasenames_printed = true;
    
            end
    
        end
        title('East component')
        xlim([0 max(t)])
        ylim([min(del)-1, max(del)+1])
        sgtitle([ 'Mean snr:' num2str(evt_snr) ]);
    
        figure(2)
        clf
        subplot(121)
        hold on
    
        phasenames_printed = false;

        for k = 1:length(allWfs)
    
            if ~isnan(allWfs(1, k).north)
    
                plot(t, scale*allWfs(1, k).north + k, 'k')
        
                for kk = 1:length(allWfs(1, k).phase_list)
        
                    %check if its a depth phase
                    phasename = allWfs(1, k).phase_list{kk};
                    if phasename(1)=='s' || phasename(1)=='p'
    
                        if ~include_depthphases
    
                            continue
    
                        end
    
                    end
    
                    plot([ allWfs(1, k).phase_times(kk), allWfs(1, k).phase_times(kk) ], [ k, (k+ 1)], 'k--')

                    if ~phasenames_printed
                    
                        text(allWfs(1, k).phase_times(kk) + 0.5, k + 1, allWfs(1, k).phase_list{kk})

                    end
                            
                end

                phasenames_printed = true;
    
            end
    
        end
        xlim([0 max(t)])
        ylim([0, length(allWfs)+1])
    
        title('North component')
    
        subplot(122)
        hold on
    
        phasenames_printed = false;

        for k = 1:length(allWfs)
    
            if ~isnan(allWfs(1, k).north)
    
                plot(t, scale*allWfs(1, k).east + k, 'k')
        
                for kk = 1:length(allWfs(1, k).phase_list)
        
                    %check if its a depth phase
                    phasename = allWfs(1, k).phase_list{kk};
                    if phasename(1)=='s' || phasename(1)=='p'
    
                        if ~include_depthphases
    
                            continue
    
                        end
    
                    end
    
                    plot([ allWfs(1, k).phase_times(kk), allWfs(1, k).phase_times(kk) ], [ k, (k + 1)], 'k--')
        
                    if ~phasenames_printed
                    
                        text(allWfs(1, k).phase_times(kk) + 0.5, k + 1, allWfs(1, k).phase_list{kk})

                    end

                end

                phasenames_printed = true;
                
            end
    
        end
        title('East component')
        xlim([0 max(t)])
        ylim([0, length(allWfs)+1])
        sgtitle([ 'Mean snr:' num2str(evt_snr) ]);

    end

end

drawnow

k = str2num(input(['pick on for Figure 2 North (1) or east ' ...
    '(2) forwards, or North (3) or East (4) backwards: '],'s'));

figure(2)
if k == 1 || k == 3

    subplot(121)

elseif k == 2 || k == 4

    subplot(122)

end

x          = ginput(1);

if k == 1

    subplot(121)
    plot( [ x(1) x(1) ], [ 0 length(allWfs)+1], 'k', 'LineWidth', 2);
    plot( [ x(1) x(1) ] + time_after, [ 0 length(allWfs)+1], 'k', 'LineWidth', 2);

    pre_picked = TD_parameters.total_time/2 - x(1);

elseif k == 2

    subplot(122)
    plot( [ x(1) x(1) ], [ 0 length(allWfs)+1], 'k', 'LineWidth', 2);
    plot( [ x(1) x(1) ] + time_after, [ 0 length(allWfs)+1], 'k', 'LineWidth', 2);

    pre_picked = TD_parameters.total_time/2 - x(1);

elseif k == 3

    subplot(121)
    plot( [ x(1) x(1) ], [ 0 length(allWfs)+1], 'k', 'LineWidth', 2);
    plot( [ x(1) x(1) ] - time_after, [ 0 length(allWfs)+1], 'k', 'LineWidth', 2);

    pre_picked = TD_parameters.total_time/2 - x(1) - time_after;

elseif k == 4

    subplot(122)
    plot( [ x(1) x(1) ], [ 0 length(allWfs)+1], 'k', 'LineWidth', 2);
    plot( [ x(1) x(1) ] - time_after, [ 0 length(allWfs)+1], 'k', 'LineWidth', 2);

    pre_picked = TD_parameters.total_time/2 - x(1) - time_after;

end

drawnow

k          = str2num(input('Keep (1)? ', 's'));
pre_picked = TD_parameters.total_time/2 - x(1);

if k == 1

    if exist([ run_name '.mat' ]) == 0
    
        evt_ind         = j;
        target_phase{1} = target_phase_use;
        pre             = pre_picked;
    
    elseif exist([ run_name '.mat' ]) == 2
    
        load([ run_name '.mat' ])
        evt_ind(length(evt_ind) + 1) = j;
        target_phase{length(target_phase) + 1} = target_phase_use;
        pre(length(pre) + 1)             = pre_picked;
    
    end
    
    save([ run_name '.mat' ], 'evt_ind', 'target_phase', 'pre')

end
