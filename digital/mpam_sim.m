clear all
close all
clc

M_Vect = [16];
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
    b = 1;

    while b<=length(EbN0dB_Vec)
        EbN0dB = EbN0dB_Vec(b);
        sim('mpam_model_inclass.slx')
        Pb_sim(c,b) = ErrorVec(1);
        EbN0 = 10^(EbN0dB/10);
        Ps = (2*(M-1)/M) * qfunc(sqrt( (6*log2(M)) / ( (M^2)-1 ) * EbN0));
        Pb_theo(c,b) = Ps/log2(M);
        b=b+1;

    end

    c = c+1;

end

figure(1)
semilogy(EbN0dB_Vec,Pb_theo(1,:),'b-o');
hold on;
%semilogy(EbN0dB_Vec,Pb_theo(2,:),'k-s');
%semilogy(EbN0dB_Vec,Pb_theo(3,:),'m-v');

semilogy(EbN0dB_Vec,Pb_sim(1,:),'r*');
%semilogy(EbN0dB_Vec,Pb_sim(2,:),'gx');
%semilogy(EbN0dB_Vec,Pb_sim(3,:),'c+');
grid on;
legend('BPAM - Theoretical', '4PAM -Theoretical', '8PAM -Theoretical', ...
    'BPAM - Simulated', '4PAM -Simulated', '8PAM -Simulated','Location','southwest');
