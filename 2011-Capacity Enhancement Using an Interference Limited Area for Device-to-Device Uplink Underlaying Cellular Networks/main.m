M = 4;
K = 100;
R_cell = 500;
R_ILA = 100;
SNR = linspace(0,20,11);
Pcb_average = 10.^(SNR/10);

C_cn_conventional = M*log2(1+Pcb_average*log2(K));
C_cn_ILA = M*log2(1+Pcb_average*log2(K*(1-R_ILA/R_cell)));

C_cn_conventional_RRS = M*log2(1+Pcb_average/M*log2(K));
C_cn_ILA_RRS = M*log2(1+Pcb_average/M*log2(K*(1-R_ILA/R_cell)));

plot(SNR,C_cn_conventional,'-^k');
hold on;
plot(SNR,C_cn_ILA,'-vk');
hold on;
plot(SNR,C_cn_conventional_RRS,':^k');
hold on;
plot(SNR,C_cn_ILA_RRS,':^k');
axis([0 20 0 70]);
grid on;