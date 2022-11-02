%%%%%%%%%this script prints out all off the events so that you can select
%%%%%%%%%what you want to do measurements on

clear
close all
clc

addpath('Custom_functions_GS');
addpath('wfTools');
addpath('CircStat2012a');
addpath('irisFetch')

dataName            = 'SSIP';
scale               = 0.01;
radial_or_traces    = 2;
include_depthphases = 0;

%parameters are amp, delay, fast direction, splitting time, src time, src, amp
Parameters.sig                 = 0.1;%always non-dimensionalized
Parameters.maxDelay            = 5;%in s
Parameters.n_iter              = 3e4;
Parameters.burn_in             = 2e4;
Parameters.keep_each           = 25;
Parameters.n_chains            = 64;
Parameters.n_cold              = 12;  
Parameters.maxCells            = 20;
Parameters.interp_style        = 'pchip';
Parameters.sig_range           = [ log(0.01), 0.5 ];%log. Traces are normalized so 0 is ~max, goes a bit high.
Parameters.r_range             = [ log(0.01), 1 ];%mean, std, log
Parameters.f_range             = [ log(0.05), 0.5 ];%mean, std, log
Parameters.max_dt              = 4;
Parameters.max_dtS             = 4;
Parameters.polarization_std    = 5*pi/180;
Parameters.orientation_std     = 0.1*pi/180;%default value
Parameters.rotation_std        = 15*pi/180;%30 degree std for layer rotation
Parameters.print_on            = 1e2;
Parameters.north_names         = { 'BH1', 'HH1', 'BHN', 'HHN' };
Parameters.east_names          = { 'BH2', 'HH2', 'BHE', 'HHE' };
Parameters.debug               = 0;
Parameters.maxT                = 10;

%parameters for loading in real data
Parameters.phase       = 'SKKS';
Parameters.total_time  = 200;
Parameters.low_pass    = 0.15;
Parameters.sample_rate = 1;%larger by at least 2
Parameters.high_pass   = 1/50;%in seconds

Parameters.dataName = dataName;

t               = (0:1/Parameters.sample_rate:(Parameters.total_time))';
Parameters.t = t;

load([ './FetchData/' dataName '.mat' ]);

mkdir([ './FetchData/' dataName '.' Parameters.phase ]);

for j = 1:length(E)

    try

        allWfs = load_data(Parameters, E(j), S, { Parameters.phase }, Parameters.total_time/2);

    catch

        continue

    end

    if any(~isnan([allWfs(:).latitude]))

        [ del, az ] = distance([allWfs(:).latitude], [allWfs(:).longitude], E(j).PreferredLatitude,...
            E(j).PreferredLongitude);

        evt_del(j) = mean(del(~isnan(del)));
        evt_azi(j) = circ_mean(az(~isnan(del))'*pi/180)*180/pi;
        
        snrlist    = [allWfs(:).snr];
        evt_snr(j) = mean(snrlist(~isnan(snrlist)));

        h = figure('visible', 'off');

        if radial_or_traces == 1
    
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
                        text(allWfs(1, k).phase_times(kk) + 0.5*scale, del(k) + sacle, allWfs(1, k).phase_list{kk})
            
                    end
        
                end
        
            end
            xlim([0 max(t)])
            ylim([min(del) - scale, max(del) + scale])
            title('Radial Traces')
        
        elseif radial_or_traces == 2
    
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
                        
                            text(allWfs(1, k).phase_times(kk) + 0.5*scale, del(k) + scale, allWfs(1, k).phase_list{kk})

                        end
            
                    end
        
                    phasenames_printed = true;

                end
        
            end
            xlim([0 max(t)])
            ylim([min(del) - 2*scale, max(del) + 2*scale])
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

                            text(allWfs(1, k).phase_times(kk) + 0.5*scale, del(k) + scale, allWfs(1, k).phase_list{kk})
            
                        end

                    end
        
                phasenames_printed = true;

                end
        
            end
            title('East component')
            xlim([0 max(t)])
            ylim([min(del) - 2*scale, max(del) + 2*scale])
            %sgtitle([ 'Mean snr:' num2str(evt_snr) ]);
            
        end

        print(h, [ './FetchData/' dataName '.' Parameters.phase '/' num2str(j) '.pdf' ], '-dpdf', '-bestfit')
        close(h)
        
    end

end

save([ dataName '.' Parameters.phase '.' num2str(Parameters.low_pass) 'Hz.mat'], 'evt_del', 'evt_azi', 'evt_snr')
