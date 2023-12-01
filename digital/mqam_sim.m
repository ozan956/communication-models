clear all
close all
clc

M_Vect = [4 16 64];
EbN0dB_Vec = 3:2:15;
c = 1;

while(c<=length(M_Vect))


    M = M_Vect(c);
    L = log2(M);
    W = 1e6;%channel bandwidth
    Rs = W; %bandpass system
    Ts = 1/Rs;
    Rb = L*Rs;
    Tb = 1/Rb;
    m = sqrt(M); %to be used for Pb of m-PAM.

    b = 1;
    while b<=length(EbN0dB_Vec)
        EbN0dB = EbN0dB_Vec(b);
        sim('mqam_model_inclass.slx')
        Pb_sim(c,b) = ErrorVec(1);
        EbN0 = 10^(EbN0dB/10);
        Ps = (2*(m-1)/m) * qfunc(sqrt((6*log2(m)/((m^2)-1))*EbN0));
        Pb_theo(c,b) = Ps/log2(m);
        b=b+1;

    end

    c = c+1;

end

figure(1)
semilogy(EbN0dB_Vec,Pb_theo(1,:),'b-o');
hold on;
semilogy(EbN0dB_Vec,Pb_theo(2,:),'k-s');
semilogy(EbN0dB_Vec,Pb_theo(3,:),'m-v');

semilogy(EbN0dB_Vec,Pb_sim(1,:),'r*');
semilogy(EbN0dB_Vec,Pb_sim(2,:),'gx');
semilogy(EbN0dB_Vec,Pb_sim(3,:),'c+');
grid on;
legend('4QAM - Theoretical', '16QAM -Theoretical', '64QAM -Theoretical', ...
    '4QAM - Simulated', '16QAM -Simulated', '64QAM -Simulated','Location','southwest');
xlabel('Eb/N_0 dB');
ylabel('BER');