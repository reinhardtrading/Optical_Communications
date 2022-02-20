%% Analyzing non-linear effects on a dispersion managed N*100 span fiber
clear
clc
%% initializzing gstate
Nsymb       = 1024;            % Number of symbols
Nt          = 16;              % Number of discrete points per symbol
Symbolrate  = 10;              % Symbol rate in Gb/s
Nsamp       = Nsymb*Nt;
fs          = Symbolrate*Nt;
inigstate(Nsamp,fs);
%% %% Transmission Fiber Parameters
Nspan         = 10;
RDPS          = 30;            % Residual dispersion per span, alternate with DM30, DM70
ft.length     = 100000;       % Length of trasmitting fiber in metres
ft.disp       = 17;           % Dispersion
ft.lambda     = 1550;
ft.ismanakov  = true;
ft.pmdpar     = 0; 
ft.alphadB    = 0.2;          % Attenuation          
ft.slope      = 0.057;               
ft.aeff       = 80;           % Effective Area
ft.n2         = 2.5e-20;      % nonlinear index [m^2/W]
ft.dzmax      = 2E4;          % maximum SSFM step size [m]
%ft.trace      = true;         % show information on screen
%% %% Compensating Fiber Parameters
fc = ft;
fc.alphadB    = 0.6;
fc.length     = 1e3;
fc.aeff       = 20; 
fc.dzmax      = 2E4;
fc.slope      = 0.057;
fc.disp       = (RDPS*1e3 - ft.disp*ft.length)/fc.length;
%% %% Amplifier Parameter
amp.gain = (ft.length*ft.alphadB + fc.length*fc.alphadB)*1e-3; % gain [dB]
amp.f    = 6;      % Amplifier noise figure
%% %% Transmitter side
modfor      = '16qam';          % Modulation format
Pdbm        = 1;               % Power in dbm 
tx.rolloff  = 0.01; 
tx.emph     = 'asin';
lam         = 1550;             % Wavelength in nm
%% %% Receiver side
rx.modformat = modfor;
rx.eftype = 'rootrc';
rx.ebw = 0.5;
rx.epar = tx.rolloff;
rx.type = 'bin';
rx.obw = Inf;
rx.oftype = 'gauss';
rx.sync.type = 'da';
%% %% Channel Part
Plin = 10.^(Pdbm/10);           % [mW]
E = lasersource(Plin,lam);
[Ex,Ey] = pbs(E);               % spliting into two orthogonal polarizations

rng(1);
patx = pattern1(Nsymb,'rand',struct('format',modfor));
paty = pattern1(Nsymb,'rand',struct('format',modfor));

[sigx,normx]= digitalmod(patx,modfor,Symbolrate,'rootrc', tx);
[sigy,normy]= digitalmod(paty,modfor,Symbolrate,'rootrc', tx);
Ex = iqmodulator(Ex,sigx,struct('norm',normx));
Ey = iqmodulator(Ey,sigy,struct('norm',normy));

Esp = pbc(Ex,Ey);                           % combining the two polarizations
%% Reccuring cycle

for ns=1:Nspan
 
 E = Esp;                % Field from previous span 
 E = fiber(E,ft);          % Electric field through the  transmission fiber
 E = fiber(E,fc);          % Electric field through the compensating fiber  
 E = ampliflat(E,amp);     % amplifying the signal

 rsig = rxfrontend(E,lam,Symbolrate,rx);
 akhat = rxdsp(rsig,Symbolrate,[patx,paty],rx);
 pathat = samp2pat(akhat,modfor,rx);                   
 SNRdB(ns) = samp2snr([patx,paty],akhat,modfor);
 
 Esp=E;
disp(SNRdB);
end



%% Comparison between SNR of different powers at different spans at RDPS 0
% Graph of RDPS 0 from data generated. Single channel system
Power = [-30,-22,-18,-14,-10,-6,-2, 0, 1,2,4,6,10,15];
SNRspan1 = [2.1462,10.1853,14.1694,18.1084,21.9483,25.5618,28.6396,29.7189,30.0323,29.9125,29.3613,27.2093,20.1889,9.0902];
fig1 =plot(Power,SNRspan1, 'o-', 'linewidth',2);
hold on
SNRspan5 = [-4.9984,3.1397,7.1686,11.1771,15.1556,19.0400,22.2315,22.5850,22.1083,19.8033,18.1678,14.4201,1.2025,-4.9984];
fig2 =plot(Power,SNRspan5, 'x-', 'linewidth',2);
hold on
SNRspan10 = [-7.08218,0.2240,4.2667,8.2754,12.2619,16.1314,18.8063,18.2569,17.2529,14.1706,12.3121,8.2694,-1.3898,-7.0821];
fig3 =plot(Power,SNRspan10, '+-', 'linewidth',2);
axis([-30 15 0 32]);
hold on
grid on
legend('1st Span','5th Span', '10th Span');
xlabel('Power[dBm]');
ylabel('SNR [dB]');
