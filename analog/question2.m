%% RECORDING MESSAGES 

clear all 
close all
clc 

% FIRST MESSAGE
Fs = 40000;
nBits = 16; % no of beats to represent each sample L 
nChannels = 1;% MONO CHANNEL

ID = -1 ;%default aduio input device
recObj = audiorecorder(Fs, nBits, nChannels, ID);
disp('Start Speaking first message for 4 seconds');
recordblocking(recObj,4);
disp('End of the recording of the first message');
myRecording = getaudiodata(recObj);
m1 = getaudiodata(recObj);
filename = 'q2_m1.wav'; % name of file 
audiowrite(filename,m1,Fs);
play(recObj);

%SECOND MESSAGE
disp('Hit enter to record the second message')
pause
disp('Start speaking second message for 4 seconds')
recordblocking(recObj,4);
disp('End of Recording the second message');
myRecording = getaudiodata(recObj);
m2 = getaudiodata(recObj);
filename = 'q2_m2.wav';
audiowrite(filename,m2,Fs);

%% IMPLEMENTING QAM

%VARIABLES
ts = 1/ Fs;
T = 4;
t =ts:ts:T;
t = t';
n = length(m1);
f = (-(n-1)/2:(n-1)/2)*(Fs/n);
fc = 8000;

%MODULATION
s_qam = m1.*cos(2*pi*fc*t) + m2.*sin(2*pi*fc*t);

%SPECTRUM ANALYZE
sa = dsp.SpectrumAnalyzer('SampleRate', Fs, ...
    'PlotAsTwoSidedSpectrum',true,'NumInputPorts',3, ...
    'ChannelNames',{'Message Signal 1','Message Signal 2','QAM Modulated Singal'});
sa(m1,m2,s_qam);
release(sa);




%% GAUSSIAN CHANNEL NOISE

N = 1e-6;  % 1 µW noise power
w = sqrt(N)*randn(size(s_qam));
s_qam = s_qam + w;

sa = dsp.SpectrumAnalyzer('SampleRate', Fs, ...
    'PlotAsTwoSidedSpectrum',true,'NumInputPorts',3, ...
    'ChannelNames',{'Message Signal 1','Message Signal 2','QAM Modulated Singal'});
sa(m1,m2,s_qam);
release(sa);

%% QAM DEMODULATOR

%FIRST STEP OF DEMODULATION
m1_rec = 2*cos(2*pi*fc*t).*s_qam;
m2_rec = 2*sin(2*pi*fc*t).*s_qam;

% SPECTRUM ANALYZE AFTER DEMODULATION
sa = dsp.SpectrumAnalyzer('SampleRate', Fs, ...
    'PlotAsTwoSidedSpectrum',true,'NumInputPorts',3, ...
    'ChannelNames',{'Message Signal 1','Message Signal 2','QAM Modulated Singal'});
sa(m1_rec,m2_rec,s_qam);
release(sa);


% LOW PASS FILTER
f_cutoff = 6000; %Hz
f_stop = 8000;
lpFit = designfilt('lowpassfir','PassbandFrequency',f_cutoff,'StopbandFrequency', f_stop, 'SampleRate', Fs);

% VISIUALIZATION OF LOW PASS FILTER
fvtool(lpFit) % visualize filter figure 4

%FILTERING
m_rec1 = filter(lpFit,m1_rec);
m_rec2 = filter(lpFit,m2_rec);

%PLOT THE MESSAGES ON TIME DOMAIN
figure(6)
subplot(411)
plot(t,m1)
xlim([1.1 1.12])
title('(a)')
xlabel('Time - Seconds')
ylabel('Amplitude')
legend('Original Signal 1')
grid on
subplot(412)
plot(t,m_rec1)
xlim([1.1 1.12])
title('(b)')
xlabel('Time - Seconds')
ylabel('Amplitude')
legend('Reconstruced Signal 1')
grid on
subplot(413)
plot(t,m2)
xlim([1.1 1.12])
title('(c)')
xlabel('Time - Seconds')
ylabel('Amplitude')
legend('Original Signal 2')
grid on
subplot(414)
plot(t,m_rec2)
xlim([1.1 1.12])
title('(d)')
xlabel('Time - Seconds')
ylabel('Amplitude')
legend('Reconstruced Signal 2')
grid on


%FFT
n= length(m1); % n is the input signal
fre_m1 = fftshift(fft(m1,n)); % computes the fourier transform (y-axis)
fre_m2 = fftshift(fft(m2,n)); % computes the fourier transform (y-axis)
fre_recm1 = fftshift(fft(m_rec1,n)); % computes the fourier transform (y-axis)
fre_recm2 = fftshift(fft(m_rec2,n)); % computes the fourier transform (y-axis)
f= (-(n-1)/2 : (n-1)/2)*(Fs/n); % generate the discrete frequency vector (x-axis)


%PLOT THE MESSAGES ON FREQUENCY DOMAIN 
figure(7)
subplot(211)
stem(f,abs(fre_m1)/n,'bo');
title('(a)')
legend('Original Message 1 Spectrum', 'Location', 'southwest');
xlabel('Frequency')
ylabel('Amplitude')
grid on 
subplot(212)
stem(f,abs(fre_recm1)/n,'ro');
title('(b)')
grid on 
legend('Reconstruced Message 1 Spectrum','Location', 'southwest');
xlabel('Frequency')
ylabel('Amplitude')
figure(8)
subplot(211)
stem(f,abs(fre_m2)/n,'bo');
title('(a)')
legend('Original Message 2 Spectrum', 'Location', 'southwest');
xlabel('Frequency')
ylabel('Amplitude')
grid on
subplot(212)
stem(f,abs(fre_recm2)/n,'ro');
title('(b)')
grid on
legend( 'Reconstruced Message 2 Spectrum','Location', 'southwest');
xlabel('Frequency')
ylabel('Amplitude')

%PLAYING
soundsc(m1,Fs);
pause(4)
soundsc(m_rec1,Fs);
pause(4)
soundsc(m2,Fs);
pause(4)
soundsc(m_rec2,Fs);

%% POWER IS 1mW VERSION

% QAM MODULATION
ts = 1/ Fs;
T = 4;
t =ts:ts:T;
t = t';
n = length(m1);
f = (-(n-1)/2:(n-1)/2)*(Fs/n);

fc = 8000;

s_qam = m1.*cos(2*pi*fc*t) + m2.*sin(2*pi*fc*t);


sa = dsp.SpectrumAnalyzer('SampleRate', Fs, ...
    'PlotAsTwoSidedSpectrum',true,'NumInputPorts',3, ...
    'ChannelNames',{'Message 1','Message 2','QAM Modulated Singal'});
sa(m1,m2,s_qam);
release(sa);

% GAUSSIAN NOISE
N = 1e-3;  % 1 mW noise power
w = sqrt(N)*randn(size(s_qam));
s_qam = s_qam + w;


% QAM DEMODULATION % FILTERING
m1_rec = 2*cos(2*pi*fc*t).*s_qam;

m2_rec = 2*sin(2*pi*fc*t).*s_qam;

% SPECTRUM ANALYZE AFTER DEMODULATION
sa = dsp.SpectrumAnalyzer('SampleRate', Fs, ...
    'PlotAsTwoSidedSpectrum',true,'NumInputPorts',3, ...
    'ChannelNames',{'Message 1','Message 2','qam'});
sa(m1_rec,m2_rec,s_qam);
release(sa);


% LOW PASS FILTER
f_cutoff = 6000; %Hz
f_stop = 8000;
lpFit = designfilt('lowpassfir','PassbandFrequency',f_cutoff,'StopbandFrequency', f_stop, 'SampleRate', Fs);

% VISIUALIZATION OF LOW PASS FILTER
%fvtool(lpFit) % visualize filter figure 4

m_rec1 = filter(lpFit,m1_rec);
m_rec2 = filter(lpFit,m2_rec);
% 

%PLOT THE MESSAGES ON TIME DOMAIN
figure(6)
subplot(411)
plot(t,m1)
xlim([1.1 1.12])
title('(a)')
xlabel('Time - Seconds')
ylabel('Amplitude')
legend('Original Signal 1')
grid on
subplot(412)
plot(t,m_rec1)
xlim([1.1 1.12])
title('(b)')
xlabel('Time - Seconds')
ylabel('Amplitude')
legend('Reconstruced Signal 1')
grid on
subplot(413)
plot(t,m2)
xlim([1.1 1.12])
title('(c)')
xlabel('Time - Seconds')
ylabel('Amplitude')
legend('Original Signal 2')
grid on
subplot(414)
plot(t,m_rec2)
xlim([1.1 1.12])
title('(d)')
xlabel('Time - Seconds')
ylabel('Amplitude')
legend('Reconstruced Signal 2')
grid on

%FFT
n= length(m1); % n is the input signal
fre_m1 = fftshift(fft(m1,n)); % computes the fourier transform (y-axis)
fre_m2 = fftshift(fft(m2,n)); % computes the fourier transform (y-axis)
fre_recm1 = fftshift(fft(m_rec1,n)); % computes the fourier transform (y-axis)
fre_recm2 = fftshift(fft(m_rec2,n)); % computes the fourier transform (y-axis)
f= (-(n-1)/2 : (n-1)/2)*(Fs/n); % generate the discrete frequency vector (x-axis)


% PLOT THE MESSAGES ON FREQUENCY DOMAIN 
figure(1)
subplot(211)
stem(f,abs(fre_m1)/n,'bo');
hold on 
grid on 
subplot(212)
stem(f,abs(fre_recm1)/n,'ro');
grid on 
legend('Original Message 1 Spectrum', 'Recovered Message 1 Spectrum','Location', 'southwest');
xlabel('Frequency')

figure(2)
subplot(211)
stem(f,abs(fre_m2)/n,'bo');
hold on 
grid on
subplot(212)
stem(f,abs(fre_recm2)/n,'ro');
grid on
legend('Original Message 2 Spectrum', 'Recovered Message 2 Spectrum','Location', 'southwest');
xlabel('Frequency')

%PLAYING
soundsc(m1,Fs);
pause(4)
soundsc(m_rec1,Fs);
pause(4)
soundsc(m2,Fs);
pause(4)
soundsc(m_rec2,Fs);