%%%%%%%%
%Script is designed to use command line arguments to guide what data gets
%loaded
%Test if the test is synthetic by looking for the words "Syn" in the first
%three letters
%%%%%%%%
rng('shuffle')
Parameters.t  = (0:1/Parameters.sample_rate:(Parameters.total_time))';

%frequency vector for FFTs
fs = 1/(Parameters.t(2) - Parameters.t(1));
dFreq = fs/length(Parameters.t); %frequency spacing
fNyq=fs/2;   %Nyquist frequency

if ~exist('station')

    station = [];

end

%the next two steps build the frequency vector (freqs at which spectrum
%was calculated) postivie and negative
f=(0:length(Parameters.t)-1)'*dFreq;
f(f>fNyq)=f(f>fNyq)-fNyq*2;
Parameters.f = f;

if ~any(~(Parameters.dataName(1:3) == 'Syn'))

    %SynName and sig are a command line arguments
    
    if exist(['./' Parameters.dataName '/' Parameters.dataName '.mat' ])

        oldrun                  = load(['./' Parameters.dataName '/' Parameters.dataName '.mat' ], 'model', 'allWfs', 'synmodel', 'Parameters');
        allWfs                  = oldrun.allWfs;
        synmodel                = oldrun.synmodel;
        Parameters.polarization = oldrun.Parameters.polarization;%to avoid overwriting everything

    else

        %validate the forward problem here!!!
        [ allWfs, synmodel, Parameters ] = load_data_syn(Parameters);

    end

    %[ allWfs, t ] = load_data(Parameters, E(10), S);

else

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp([ 'Getting the data for ' Parameters.dataName ]);

    if exist(['./' Parameters.name '/' Parameters.name 'Data.mat' ])

        oldrun                     = load(['./' Parameters.name '/' Parameters.name 'Data.mat' ], 'allWfs', 'Parameters');
        allWfs                     = oldrun.allWfs;
        Parameters.polarization = oldrun.Parameters.polarization;%to avoid overwriting everything

    else

        load([ './FetchData/' Parameters.dataName '.mat' ]);

        if isempty(station)

            allWfs = load_data(Parameters, E(evt_ind), S, target_phase, pre);

        else

            allWfs = load_data(Parameters, E(evt_ind), S(station), target_phase, pre);

        end
    
        lat = [allWfs(:).latitude];
        lon = [allWfs(:).longitude];
    
        for k = 1:length(evt_ind)
    
            [~,Parameters.polarization(k,1)] = distance(nanmean(lat), nanmean(lon), E(evt_ind(k)).PreferredLatitude,...
                E(evt_ind(k)).PreferredLongitude);
    
        end
    
        Parameters.polarization = mod(Parameters.polarization*pi/180, 2*pi);        

        mkdir(['./' Parameters.name '/'])
        save(['./' Parameters.name '/' Parameters.name 'Data.mat' ], 'allWfs', 'Parameters');

    end

end

%pre-save kappas and orientations
[nevt, nsta] = size(allWfs);

for k = 1:nsta

    x   = [allWfs(:, k).orientation_k];

    if isempty(x(~isnan(x)))
        Parameters.sta_err(k,1) = 2*pi;
    else
        Parameters.sta_err(k,1) = unique(x(~isnan(x)));
        Parameters.sta_err(k,1) = sqrt(1/Parameters.sta_err(k,1));
    end

end

%count up how many traces you have (missing data are nans)
n = 0;
for k = 1:numel(allWfs)

    if ~any(isnan(allWfs(k).north))

        n = n + 1;

    end

end

Parameters.n = n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
model  = fit_waveforms(allWfs, Parameters);

disp('-> Saving results')

%remove large fields that can be easily generated
msave = [];
for k = 1:length(model)

    m = rmfield(model(k), 'Cinv');
    m = rmfield(m,        'irF');
    m = rmfield(m,        'irS');
    m = rmfield(m,        'N'  );
    m = rmfield(m,        'E'  );

    msave = [ msave; m ];

end
model = msave;

mkdir(['./' Parameters.name '/'])
if ~any(~(Parameters.name(1:3) == 'Syn'))

    save([ './' Parameters.name '/' Parameters.name num2str(rand(), 4) '.' Parameters.solver ...
        '.mat' ], 'allWfs', 'Parameters', 'model', 'synmodel')

else

    save([ './' Parameters.name '/' Parameters.name num2str(rand(), 4) '.' Parameters.solver ...
        '.mat' ], 'allWfs', 'Parameters', 'model')

end
