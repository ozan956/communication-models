%% RECORDING SOUND

clear all
close all
clc

Fs = 40000;
td = 1/Fs;
T=4;
t=td:td:T;

nBits = 16; % no of beats to represent each sample L 
nChannels = 1;% MONO CHANNEL

ID = -1 ;%default aduio input device
recObj = audiorecorder(Fs, nBits, nChannels, ID);
disp('Start Speaking first message for 4 seconds');
recordblocking(recObj,4);
disp('End of the recording of the first message');
myRecording = getaudiodata(recObj);
message = getaudiodata(recObj);
filename = 'q1_m1.wav'; % name of file 
audiowrite(filename,message,Fs);
play(recObj);

%% DOWNSAMPLING THE MESSAGE

fs=8000; % sample the (approximately) continious time signal at a rate Fs.

ts = 1/fs;
N = ts/td; % N should be an integer

%DOWN SAMPLING
s_message = downsample(message,N); % downsampling at a rate of 8kHz.
s_message = upsample(s_message,N); % to standarize original message and downsampled message. 

%PLOT IN TIME DOMAIN
figure(1)
subplot(211)
plot(t,message,'b')
title('(a)')
xlabel('Time - Seconds'); ylabel('Amplitude'); legend('Original Signal')
xlim([1.1 1.12])
hold on
subplot(212)
stem(t,s_message,'r');
xlim([1.1 1.12])
title('(b)')
xlabel('Time - Seconds')
ylabel('Amplitude')
legend('Sampled Signal');

%% SPECTRUM ANALYZER

sa = dsp.SpectrumAnalyzer('SampleRate', Fs, ...
    'PlotAsTwoSidedSpectrum',true,'NumInputPorts',2, ...
    'ChannelNames',{'Original Signal','Downsampled Signal'});
sa(message,s_message);
release(sa);

%% LOWPASS FILTER
n=length(message);
f = (-(n-1)/2:(n-1)/2)*(Fs/n);

%CREATING LOW PASS FILTER
fcutoff = 3000;
fstop = 4000;
lpfit = designfilt('lowpassfir', 'PassbandFrequency', fcutoff,'StopbandFrequency',fstop, ...
    'SampleRate', Fs);
fvtool(lpfit);

%FILTERING
message_rec = n*filter(lpfit,s_message);

%PLOTTING IN TIME DOMAIN
figure(4)
subplot(211)
plot(t,message); 
xlim([1.1 1.12])
legend('Original Signal');
xlabel('Time - Seconds')
ylabel('Amplitude')
title('(a)')
grid on
subplot(212)
plot(t,message_rec);
xlim([1.1 1.12])
grid on
xlabel('Time - Seconds')
ylabel('Amplitude')
legend('Recovered Signal');
title('(b)')

%FFT
fre_message_rec =fftshift(fft(message_rec,n));
fre_message =fftshift(fft(message,n));

%PLOTTING IN FREQUECNY DOMAIN
figure(5)
subplot(211)
stem(f,abs(fre_message)/n,'bo');
xlabel('Frequency')
ylabel('Amplitude')
legend('Original Signal Spectrum');
title('(a)')
grid on
hold on 
subplot(212)
stem(f,abs(fre_message_rec/n), 'r-o');
grid on
xlabel('Frequency')
ylabel('Amplitude')
legend('Recovered Signal Spectrum');
title('(b)')

%PLAYING SOUNDS
disp('playing');
soundsc(message,Fs);
pause(4);
soundsc(message_rec,Fs);


%% SAME OPERATION WITH SAMPLING RATE OF 4 KHZ

fs=4000; % sample the (approximately) continious time signal at a rate Fs.

ts = 1/fs;
N = ts/td;% N should be an integer

s_message = downsample(message,N);
s_message = upsample(s_message,N);

figure(1)
subplot(211)
plot(t,message,'b')
xlabel('Time - Seconds')
legend('Original Message')
hold on
subplot(212)
stem(t,s_message,'r');
xlabel('Time - Seconds')
legend('Sampled Message');

n=length(message);
f = (-(n-1)/2:(n-1)/2)*(Fs/n);

fcutoff = 1500;
fstop = 2000;
lpfit = designfilt('lowpassfir', 'PassbandFrequency', fcutoff,'StopbandFrequency',fstop, ...
    'SampleRate', Fs);
%fvtool(lpfit);


message_rec = n*filter(lpfit,s_message);

sa = dsp.SpectrumAnalyzer('SampleRate', Fs, ...
    'PlotAsTwoSidedSpectrum',true,'NumInputPorts',2, ...
    'ChannelNames',{'Original Signal','Downsampled Signal'});
sa(message,s_message);
release(sa)

figure(4)
subplot(211)
plot(t,message);
xlim([1.1 1.12])
legend('original message');
grid on
hold on
subplot(212)
plot(t,message_rec);
xlim([1.1 1.12])
grid on
legend('recovered message');

fre_message_rec =fftshift(fft(message_rec,n));
fre_message =fftshift(fft(message,n));

figure(5)
subplot(211)
stem(f,abs(fre_message)/n,'bo');
legend('Original Message Spectrum');
grid on
hold on 
subplot(212)
stem(f,abs(fre_message_rec/n), 'r-o');
grid on
legend('Recovered Message Spectrum');
ylim([0 0.6]);

disp('playing');
soundsc(message,Fs);
pause(4)
soundsc(message_rec,Fs);