# CochlearProcessing
%% Enter file path here:
filepath = '10.Single Male W- Low Chatter.m4a';
newfilepath = '';

%% Call Phase 1 processing function call
[resampled_data, resampled_freq] = process_signal(filepath, newfilepath);

%% Setup
lcfs = [105,200,300,400,510,630,770,920,1075,1260,1480,1710,1990,2310,2675,3125,3650,4350,5250,6350];
ucfs = [200,300,400,510,630,770,920,1080,1265,1480,1720,1990,2310,2690,3125,3675,4350,5250,6350,7995];

base_size = size(resampled_data);
duration = base_size(1)/resampled_freq;
t = 0:(1/resampled_freq):duration;

%% Task 5: function call for each filter
filtered_signals = zeros(base_size(1), 20);
for i=1:20
    passband = [lcfs(i),ucfs(i)];
    filtered_signals(:,i) = filter_data(passband, resampled_data);
end

%% Task 6
%figure, plot(t(1,1:end-1), filtered_signals(:,1));
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered Signal: 105 to 200 Passband');
%figure, plot(t(1,1:end-1), filtered_signals(:,20));
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered Signal: 6350 to 7995 Passband');

%% Tasks 7&8: envelope extraction for each passband
enveloped_signals = zeros(base_size(1), 20);
for i=1:20
    filtered_signals(:,i) = abs(filtered_signals(:,i));
    enveloped_signals(:,i) = envelop_data(filtered_signals(:,i));
end

%% Task 9
%figure, plot(t(1,1:end-1), enveloped_signals(:,1));
xlabel('Time (s)');
ylabel('Amplitude');
title('Enveloped Signal: 105 to 200 Passband');
%figure, plot(t(1,1:end-1), enveloped_signals(:,20));
xlabel('Time (s)');
ylabel('Amplitude');
title('Enveloped Signal: 6350 to 7995 Passband');

%% Task 10
cos_signals = zeros(base_size(1), 20);
for i=1:20
    cf = sqrt(lcfs(i) * ucfs(i))

    %%soundsc(cos(cf*2*pi*t), cf);
    %%a = cos(cf*2*pi*time(1,1:end));
    cos_signals(:,i) = cos(cf*2*pi*t(1,1:end-1));
    %%figure, plot(t(1,1:end-1), cos_signals(:,i))
end

%% Task 11
amp_mods = zeros(base_size(1), 20);
for i=1:20
    amp_mods(:,i) = cos_signals(:,i) .* enveloped_signals(:,i);
    %figure, plot(t(1,1:end-1), amp_mods(:,i));
end

%% Task 12
summed_signal = zeros(base_size(1), 1);
for i=1:20
    summed_signal(:, 1) = summed_signal(:, 1) + amp_mods(:, i);
end
abs_signal = zeros(base_size(1), 1);
abs_signal(:,1) = abs(summed_signal(:,1));
norm_val = max(abs_signal(:,1));
summed_signal = summed_signal(:,1) ./ norm_val;
soundsc(summed_signal, resampled_freq)
%% Phase 1 processing function
function [resampled_data, resampled_freq] = process_signal(filepath, newfilepath)
    cos_freq = 1000;
    resampled_freq = 16000;
    
    %% 3.1: Read file stored at filepath into Matlab
    [audiodata, rate] = audioread(filepath);
    
    %% 3.2: Check if signal is stereo or mono; convert to mono if needed
    sizevector = size(audiodata);
    monodata = zeros(sizevector(1), 1);
    if sizevector(2) == 2
        for i=1:sizevector(1)
            monodata(i) = audiodata(i,1) + audiodata(i,2);
        end
    else
        monodata = audiodata;
    end
    
    %% 3.3: Play the sound in Matlab
    %sound(audiodata, rate);
    
    %% 3.4: Write the sound to a new file
    %audiowrite(newfilepath, audiodata, rate);
    
    %% 3.5: Plot the sound waveform as a fct. of the sample number
    %plot(monodata)

    %% 3.6: Resample the signal to 16 kHz if necessary
    if rate > resampled_freq
        [p,q] = rat(resampled_freq/rate);
        resampled_data = resample(monodata, p, q);
    else
        resampled_data = monodata;
    end
    %%Resampled data is plotted here
    %figure, plot(resampled_data);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    title('Original Signal');

    
    %% 3.7: Generate cosine signal, play sound, and plot 2 waveforms
%     size_resampled = size(resampled_data);
%     duration = size_resampled(1)/resampled_freq;
%     t = 0:(1/size_resampled(1)):duration;
%     
%     soundsc(cos(cos_freq*2*pi*t), cos_freq);
%     
%     t_waveform = t(1:2/cos_freq*size_resampled(1));
%     figure, plot(t_waveform, cos(cos_freq*2*pi*t_waveform))
end

%% Bandpass filter (Task 5)
function [filtered_data] = filter_data(passband, data)
    Fs = 16000;  % Sampling Frequency

    N   = 10;    % Order
    Fc1 = passband(1);   % First Cutoff Frequency
    Fc2 = passband(2) ; % Second Cutoff Frequency

    % Construct an FDESIGN object and call its BUTTER method.
    h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
    Hd = design(h, 'butter');

    filtered_data = filter(Hd, data);
end

%% Lowpass filter for enveloping (Task 8)
function [enveloped_data] = envelop_data(x)
    Fs = 16000;  % Sampling Frequency

    N  = 10;   % Order
    Fc = 400;  % Cutoff Frequency

    % Construct an FDESIGN object and call its BUTTER method.
    h  = fdesign.lowpass('N,F3dB', N, Fc, Fs);
    Hd = design(h, 'butter');
    enveloped_data = filter(Hd, x);
end
