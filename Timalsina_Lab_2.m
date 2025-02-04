%% 3.1. Part A : Analyzing audio signal and downsampling

% (2) Read the WAV file and get its original sampling frequency
[signal, Fs] = audioread('human_voice.wav'); 
disp(['Original Sampling Frequency: ', num2str(Fs), ' Hz']);

% (3) Plot the original signal
t = (0:length(signal)-1) / Fs; % Time vector
figure;
plot(t, signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Audio Signal');
grid on;

% (4) Downsample the audio to 8kHz manually (without inbuilt functions)
Fs_new = 8000; 
downsample_factor = round(Fs / Fs_new); 
signal_downsampled = signal(1:downsample_factor:end); 
Fs_ds = Fs / downsample_factor; 

% (5) Count the number of samples after downsampling
num_samples_downsampled = length(signal_downsampled);
disp(['Number of audio samples: ', num2str(num_samples_downsampled)]);

% (6) Plot the downsampled signal
t_ds = (0:length(signal_downsampled)-1) / Fs_new;
figure;
plot(t_ds, signal_downsampled);
xlabel('Time (s)');
ylabel('Amplitude');
title('Downsampled Audio Signal (8 kHz)');
grid on;

% (7) Compare a section of the original and downsampled signal
% It can be seen from the graph that the sound becomes a little dampend
% when it is downsampled by 8Khz. If you look closely at the amplitude 0.6
% and -0.4, you can see a little difference between the two graph. We lost
% some of high frequencies during the downsampling.


%% 3.2 Part B: RMS, Cross-correlation, Sound Localization


% (1) Calculate RMS values for each audio signal
[M1, Fs1] = audioread('M1.wav');
[M2, Fs2] = audioread('M2.wav');
[M3, Fs3] = audioread('M3.wav');

rms_M1 = sqrt(mean(M1.^2));
rms_M2 = sqrt(mean(M2.^2));
rms_M3 = sqrt(mean(M3.^2));

fprintf('RMS Values:\n');
fprintf('M1: %.4f\n', rms_M1);
fprintf('M2: %.4f\n', rms_M2);
fprintf('M3: %.4f\n', rms_M3);

% (2) Determine which microphone is closer to the sound source
if rms_M1 > rms_M2
    fprintf('Microphone M1 is closer to the sound source.\n');
else
    fprintf('Microphone M2 is closer to the sound source.\n');
end

% (3) Compute Time Delay Using Cross-Correlation (Without Inbuilt Functions)
M1 = M1 - mean(M1);
M2 = M2 - mean(M2);

len_x = length(M1);
len_y = length(M2);
corr_length = len_x + len_y - 1;
cross_corr = zeros(1, corr_length);

for m = 1:corr_length
    sum_corr = 0;
    for n = 1:len_x
        if (n + m - len_x) > 0 && (n + m - len_x) <= len_y
            sum_corr = sum_corr + M1(n) * M2(n + m - len_x);
        end
    end
    cross_corr(m) = sum_corr;
end

[max_corr, max_index] = max(cross_corr);
time_lag = (max_index - len_x) / Fs1; 

fprintf('Time delay between M1 and M2: %.6f seconds\n', time_lag);

% (4) Compute Angle Theta for Robot Heading Correction
r = 0.1; 
speed_of_sound = 343; 
d1 = time_lag * speed_of_sound; 
d2 = sqrt(d1^2 + (2*r)^2 - 2*d1*(2*r)*cos(pi/2));

theta = atan((d2 - d1) / (2*r));
theta_degrees = rad2deg(theta);

fprintf('Robot must turn by an angle of %.4f degrees.\n', theta_degrees);


%% 3.3 Part C: Frequency Spectrum

% (1) Load and Plot the Audio Signal
[audio, Fs] = audioread('Cafe_with_noise.wav');
t = (0:length(audio)-1) / Fs; 

figure;
plot(t, audio);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Audio Signal - Cafe with Noise');
grid on;

% (2) Frequency Domain Analysis (FFT)
N = length(audio); 
freqs = linspace(-Fs/2, Fs/2, N); 
audio_fft = fftshift(fft(audio)); 

figure;
plot(freqs, abs(audio_fft) / N);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of Cafe_with_noise.wav');
grid on;

% Human voice typically lies in the range of 300 Hz to 3 kHz.
% Noise in a cafe often contains high frequencies (> 3 kHz).

% (3) Implement Low-Pass Filter to Remove Noise
fc = 3000; 
order = 6; 
[b, a] = butter(order, fc / (Fs / 2), 'low');
filtered_audio = filtfilt(b, a, audio);
figure;
plot(t, filtered_audio);
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered Audio Signal (Voice Signal)');
grid on;

% Play Original and Filtered Audio
% sound(audio, Fs); % Play original
% pause(length(audio)/Fs + 1);
% sound(filtered_audio, Fs); % Play filtered

% (4) Save Filtered Audio File
audiowrite('Filtered_Voice.wav', filtered_audio, Fs);



