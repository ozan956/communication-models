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
filename = 'q3_m1.wav'; % name of file 
audiowrite(filename,message,Fs);
play(recObj);

%% USB MODULATOR 


Fc = 8000;
n = length(message);
m_h = imag(hilbert(message));
t = t';

% CONSTRUCTING SIGNAL
s_usb = message.*cos(2*pi*Fc*t) - m_h .* sin(2*pi*Fc*t); 
s_dsb = message.*cos(2*pi*Fc*t);
fre_susb = fftshift(fft(s_usb,n));

%ANALYZE OF SIGNALS
sa = dsp.SpectrumAnalyzer('SampleRate', Fs, ...
    'PlotAsTwoSidedSpectrum',true,'NumInputPorts',3, ...
    'ChannelNames',{'Original Message','USB Modulated','DSB Modulated'});
sa(message,s_usb,s_dsb);
release(sa);

%% GAUSSIAN CHANNEL NOISE

N = 1e-3;  % 1 µW noise power
w = sqrt(N)*randn(size(s_usb));
s_usb = s_usb + w;

%ANALYZE OF NOISE ADDED SIGNALS
sa = dsp.SpectrumAnalyzer('SampleRate', Fs, ...
    'PlotAsTwoSidedSpectrum',true,'NumInputPorts',2, ...
    'ChannelNames',{'Original Message','USB Modulated'});
sa(message,s_usb);
release(sa);
%% USB DEMODULATOR

%VARIABLES
f = (-(n-1)/2:(n-1)/2)*(Fs/n);
r = 2*cos(2*pi*Fc*t).*s_usb;
fre_r = fftshift(fft(r,n));
fre_o = fftshift(fft(message,n));
r = r + 1e-3;

%LOW PASS FILTER
f_cutoff = 6000; %Hz
f_stop = 8000;
lpFit = designfilt('lowpassfir','PassbandFrequency',f_cutoff,'StopbandFrequency', f_stop, 'SampleRate', Fs);
fvtool(lpFit) % visualize filter figure 4

%FFT
m_rec = filter(lpFit,r);
fre_m_rec = fftshift(fft(m_rec,n));

%PLOT THE MESSAGES ON TIME DOMAIN
figure(6)
subplot(211)
plot(t,message)
xlim([1.1 1.12])
title('(a)')
xlabel('Time - Seconds')
ylabel('Amplitude')
legend('Original Signal 1')
grid on
subplot(212)
plot(t,m_rec)
xlim([1.1 1.12])
title('(b)')
xlabel('Time - Seconds')
ylabel('Amplitude')
legend('Reconstruced Signal 1')
grid on

%PLOT THE MESSAGES ON FREQUENCY DOMAIN 
figure(7)
subplot(211)
stem(f,abs(fre_o)/n,'bo');
title('(a)')
legend('Original Message 1 Spectrum', 'Location', 'southwest');
xlabel('Frequency')
ylabel('Amplitude')
grid on 
subplot(212)
stem(f,abs(fre_m_rec)/n,'ro');
title('(b)')
grid on 
legend('Reconstruced Message 1 Spectrum','Location', 'southwest');
xlabel('Frequency')
ylabel('Amplitude')

%LISTENING
soundsc(message,Fs);
pause(4)
soundsc(m_rec,Fs);

%% CHANNEL NOISE IS 1mW VERSION 

% USB / LSB / SSB MODULATION
Fc = 8000;
n = length(message);
m_h = imag(hilbert(message));

s_usb = message.*cos(2*pi*Fc*t) - m_h .* sin(2*pi*Fc*t);
s_lsb = message.*cos(2*pi*Fc*t) + m_h .* sin(2*pi*Fc*t);
s_ssb = message.*cos(2*pi*Fc*t);

fre_susb = fftshift(fft(s_usb,n));

sa = dsp.SpectrumAnalyzer('SampleRate', Fs, ...
    'PlotAsTwoSidedSpectrum',true,'NumInputPorts',3, ...
    'ChannelNames',{'am','USB','SSB'});
sa(message,s_usb,s_ssb);
release(sa);


% GAUSSIAN NOISE
N = 1e-3;  % 1 mW noise power
w = sqrt(N)*randn(size(s_usb));
s_usb = s_usb + w;


% LSB DEMODULATION
f = (-(n-1)/2:(n-1)/2)*(Fs/n);
r = 2*cos(2*pi*Fc*t).*s_usb;
fre_r = fftshift(fft(r,n));
figure(3)
subplot(211)
stem(f,abs(fre_r)/n,'bo');

f_cutoff = 5000; %Hz
f_stop = 10000;

lpFit = designfilt('lowpassfir','PassbandFrequency',f_cutoff,'StopbandFrequency', f_stop, 'SampleRate', Fs);

%fvtool(lpFit) % visualize filter figure 4

m_rec = filter(lpFit,r);
fre_m_rec = fftshift(fft(m_rec,n));
figure(3)
subplot(212)
stem(f,abs(fre_m_rec)/n,'bo');

%PLOT THE MESSAGES ON TIME DOMAIN
figure(6)
subplot(211)
plot(t,message)
xlim([1.1 1.12])
title('(a)')
xlabel('Time - Seconds')
ylabel('Amplitude')
legend('Original Signal 1')
grid on
subplot(212)
plot(t,m_rec)
xlim([1.1 1.12])
title('(b)')
xlabel('Time - Seconds')
ylabel('Amplitude')
legend('Reconstruced Signal 1')
grid on

soundsc(message,Fs);
pause(4)
soundsc(m_rec,Fs);
