clear all
close all
clc

M_Vect = [2 4 8];
EbN0dB_Vec = 2:15;
c = 1;

while(c<=length(M_Vect))


    M = M_Vect(c);
    k = log2(M);
    Rb = 10e3;
    Rs = Rb/k;
    W = M*Rs;

    Ts = 1/Rs;
    Tb = 1/Rb;

    b = 1;
    while b<=length(EbN0dB_Vec)
        EbN0dB = EbN0dB_Vec(b);
        sim('mfsk_model_inclass.slx')
        Pb_sim(c,b) = ErrorVec(1);
        EbN0 = 10^(EbN0dB/10);
        EsN0 = k*EbN0;

        Ps = ((M-1)/2)*exp(-EsN0/2);
        Pb_theo(c,b) = ((M/2)/(M-1))*(Ps);
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
title('E_b/N_0 vs BER');
legend('BFSK - Theoretical','4FSK - Theoretical', '8FSK -Theoretical', ...
    'BFSK - Simulated','4FSK - Simulated', '8FSK -Simulated','Location','southwest');
xlabel('E_b / N_0 dB');
ylabel('BER - P_b');