clear all
close all
clc

W = 1e6;%channel bandwidth
Rs = W; %bandpass system
Ts = 1/Rs;

M_Vect = [4 8];
EbN0dB_Vec = 3:2:15;
c = 1;

while(c<=length(M_Vect))


    M = M_Vect(c);
    k = log2(M);

    Rb = k*Rs;
    Tb = 1/Rb;

    b = 1;
    while b<=length(EbN0dB_Vec)
        EbN0dB = EbN0dB_Vec(b);
        sim('mpsk_model_inclass.slx')
        Pb_sim(c,b) = ErrorVec(1);
        EbN0 = 10^(EbN0dB/10);
        EsN0 = k*EbN0;

        Ps = 2* qfunc(sqrt(2*EsN0)*sin(pi/M));
        Pb_theo(c,b) = Ps/k;
        b=b+1;

    end

    c = c+1;

end

figure(1)
semilogy(EbN0dB_Vec,Pb_theo(1,:),'b-o');
hold on;
semilogy(EbN0dB_Vec,Pb_theo(2,:),'k-s');

semilogy(EbN0dB_Vec,Pb_sim(1,:),'r*');
semilogy(EbN0dB_Vec,Pb_sim(2,:),'gx');

grid on;
legend('4PSK - Theoretical', '8PSK -Theoretical', ...
    '4PSK - Simulated', '8PSK -Simulated','Location','southwest');
xlabel('E_b / N_0 dB');
ylabel('BER - P_b');