%%%%%%%%%this script prints out all off the events so that you can select
%%%%%%%%%what you want to do measurements on

clear
close all

addpath('Custom_functions_RF');
addpath('wfTools');
addpath('CircStat2012a');
addpath('irisFetch')

dataName            = 'TA_Shen2012';
scale               = 0.1;
radial_or_traces    = 2;%2 is radial, 2 is traces
time_after = 50;

j                = 15;
target_phase_use = { 'P' };

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
TD_parameters.vertical_names      = { 'BHZ', 'HHZ' };
TD_parameters.north_names         = { 'BH1', 'HH1', 'BHN', 'HHN' };
TD_parameters.east_names          = { 'BH2', 'HH2', 'BHE', 'HHE' };
TD_parameters.debug               = 0;
TD_parameters.maxT                = 50;
TD_parameters.orientation_std     = 1;

%parameters for loading in real data
TD_parameters.total_time  = 50;%in seconds
TD_parameters.low_pass    = 2;
TD_parameters.sample_rate = 10;%larger by at least 2
TD_parameters.high_pass   = 1/5;%3/(TD_parameters.total_time);%in seconds

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

    figure(1)
    clf
    subplot(131)
    hold on
    phasenames_printed = false;

    for k = 1:length(allWfs)

        if ~isnan(allWfs(1, k).Z)

            plot(t, scale*allWfs(1, k).Z + del(k), 'k')
    
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
    title('Vertical component')

    subplot(132)
    hold on

    phasenames_printed = false;

    for k = 1:length(allWfs)

        if ~isnan(allWfs(1, k).north)

            plot(t, scale*allWfs(1, k).R + del(k), 'k')
    
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
    title('Radial component')

    subplot(133)
    hold on

    phasenames_printed = false;
    
    for k = 1:length(allWfs)

        if ~isnan(allWfs(1, k).north)

            plot(t, scale*allWfs(1, k).T + del(k), 'k')
    
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
    title('Transverse component')
    xlim([0 max(t)])
    ylim([min(del)-1, max(del)+1])
    sgtitle([ 'Mean snr:' num2str(evt_snr) ]);

end
