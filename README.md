# Analog & Digital Communication Models
Analog &amp; Digital implementations of communication techniques including AM, PM, FM, MPAM, MQAM, MPSK and MFSK.

## Analog Part 
Analog communication models have examined and tested with using case study method.
### Cases
#### Case 1
a) In MATLAB, read the following sentence and record your own voice for 4 seconds at 
a rate 40 kHz and generate a “message.wav” file.
“Use a pencil to write the first draft.”

b) Down sample your message signal at a sampling rate 8 kHz, show the original signal 
and its discrete time samples for a very short period of time.

c) Using Spectrum Analyzer, show the two sided spectrum of the original signal and the 
sampled signal.

d) Define a low pass filter with suitable cutoff frequencies to recover the original low 
pass filter from its discrete time samples. Listen the original message signal and 
recovered message signal at 40 kHz and comment on your results.

e) Repeat the above steps for the new sampling rate 4 kHz. What happens at 4 kHz? 
Comment on your results.

#### Case 2
a) In MATLAB, read the following sentences and record your own voice for 4 seconds at 
a rate 40 kHz and generate “message1.wav” and “message2.wav” files.
“The two met while playing on the sand.”
“This is a grand season for hikes on the road.”

b) In MATLAB, implement an analog QAM modulator with carrier frequency 8 kHz. 
Show the two-sided QAM spectrum using Spectrum Analyzer.

c) Add Gaussian channel noise whose power is 1 micro-Watts.

d) Implement QAM demodulator and recover message 1 and message 2. Listen the 
original message signals and recovered message signals and comment on your results.

e) Repeat the above steps with channel noise whose power is 1 milli Watts. Comment on 
your results.
#### Case 3
In this question, you are NOT allowed to use “ssbmod” and “ssbdemod” functions.
a) In Matlab, read the following sentence and record your own voice for 4 seconds at a 
rate 40 kHz and generate “message.wav” file,
“We find joy in the simplest things. “

b) In Matlab, implement an USB modulator with carrier frequency 8 kHz. Show the two 
sided USB spectrum and DSB spectrum using Spectrum Analyzer.

c) Add Gaussian channel noise whose power is 1 micro Watts.

d) Implement USB/LSB demodulators and recover the message. Listen the original 
message signals and recovered message signals and comment on your results.

e) Repeat the above steps with channel noise whose power is 1 milliWatts. Comment on 
your results.

## Digital Part
Simulink models for each model has implemented and their test &amp; compare codes have added.
