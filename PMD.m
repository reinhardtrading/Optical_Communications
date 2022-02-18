%%calculates the effect of PMD at different distances.
%vary the fiber distances.

clear
clc
rng(2020)
%% %% initializzing gstate

% Global parameters
global GSTATE
global CONSTANTS;               % CONSTANTS is a global structure variable.
%% %% Transmission Fiber Parameters
Nspan         = 1;
ft.length     = 50000;       % Length of trasmitting fiber in metres
ft.disp       = 17;           % Dispersion
ft.lambda     = 1550;
ft.ismanakov  = true;
%ft.pmdpar     = 0;
ft.alphadB    = 0.2;          % Attenuation
ft.slope      = -2*ft.disp/ft.lambda;
ft.aeff       = 80;           % Effective Area
ft.n2         = 2.5e-20;      % nonlinear index [m^2/W]
ft.dzmax      = 2E4;          % maximum SSFM step size [m]
%ft.trace      = true;         % show information on screen
%ft.nplates    = 100;           % Number of waveplates
ft.coupling   ='pol';
ft.exportpar  = true;
ft.dphimax    = 1e-4;
ft.dphi1fwm   = false;
ft.stepupd    = 'nlp';
%% %% Amplifier Parameter
amp.gain = (ft.length*ft.alphadB)*1e-3; % gain [dB]
%amp.f    = 6;      % Amplifier noise figure
%% %% Transmitter side
Nch         = 2;
%modfor      = '16qam';          % Modulation format
modfor      = 'randn';
PdBm        = 0;                % Power in dbm
tx.rolloff  = 0.01;
tx.emph     = 'asin';
lam         = 1550;             % Wavelength in nm
%deltaf      = 37.5:12.5:150;             % channel spacing [GHz]
deltaf      = 50;
Ncut = ceil(Nch/2);             % channel under test (CUT) index
%% Receiver side
rx.modformat = modfor;      % modulation format
rx.sync.type = 'da';        % time-recovery method
rx.sync.interp = 'fine';
rx.oftype = 'gauss';        % optical filter type
rx.obw = Inf;               % optical filter bandwidth normalized to symbrate
rx.eftype = 'rootrc';       % optical filter type
rx.ebw = 0.5;               % electrical filter bandwidth normalized to symbrate
rx.epar = tx.rolloff;
rx.type = 'bin';            % binary pattern
rx.adc.sps = 2;             % samples per symbol

issave = true;

% init for saving
thisfile = textread([mfilename,'.m'], '%s', 'delimiter', '\n', 'whitespace','');
try % for LinuX
    [~,thisnohup]=system(['ls -drt ',mfilename,'*nohup | tail -1']);
catch
    thisnohup = '';
end
file_name = ['Nch=',num2str(Nch),...
    '_',modfor,...
    '_PdBm=',num2str(PdBm),...
    '_Df=',num2str(deltaf),....
    '_Nspan=',num2str(Nspan),....
    ]; % save info in this file (.mat or .txt)

outdir = mfilename; % same name of this file
%outdir= './Thesis results/POGGIOLINI';
[outfile,fid] = preparefile(outdir,file_name,issave); % get output file
% fid = 1: do not write on file, otherwise write on outfile.txt


hostn = char(java.net.InetAddress.getLocalHost.getHostName); % host name



pmdvec= 0:5:50; %4:100;

for pmd= 1:length(pmdvec)
    
    ft.pmdpar=pmdvec(pmd);
    
    
    %%%% Condition on PMD  parameter value
    if ft.pmdpar == 0
        ft.nplates= 1;
    elseif ft.pmdpar ~= 0
        ft.nplates=2000;
    end
    
    lcorr         = ft.length/ft.nplates;             % length of the waveplate
    pmdpar        = ft.pmdpar;
    clight        = 299792458;
    
    
    
    ft1 = ft;
    
    
    
    
    
    for sp=1:length(deltaf)
        
        spac = deltaf(sp)*lam^2/299792458;  % channel spacing [nm]
        %% Channel
        
        Dcum= Nspan*ft.disp*(ft.length/1000);
        % properties
        Symbolrate  = 49;  % Symbol rate in Gb/s
        Nsymb_1 = (abs(Dcum)*deltaf(sp)*(lam^2)* Symbolrate* 1e-3)/clight;  % Number of symbols
        %Nsymb_min = 2.^nextpow2(Nsymb_1);
        Nsymb_1 = 2.^nextpow2(Nsymb_1);
        Nsymb_min = 2^14;
        
        
        if Nsymb_1 < 4096
            Nsymb =  Nsymb_min;
        else
            Nsymb = Nsymb_min;
        end
        
        Nt_1  = deltaf(sp)/Symbolrate;        % Number of discrete points per symbol
        Nt_1=2.^nextpow2(Nt_1);
        Nt_min = 4;
        
        if Nt_1 < Nt_min
            Nt = Nt_min;
        else
            Nt = Nt_1;
        end
        
        Nsamp       = Nsymb*Nt;
        fs          = Symbolrate*Nt;
        
        inigstate(Nsamp,fs);
        
        for m=unique([1 fid]) % if output file exists, write both on file and command window
            fprintf(m,'Start (%s):\t%s.\n\n',hostn, datestr(now));
            [~,outname] = fileparts(outfile);
            fprintf(m,'Matlab file:\t%s\nOutput dir:\t%s\n',mfilename,outdir);
            if fid ~=1, fprintf(m,'Output file:\t%s\n\n',outname); else, fprintf(m,'\n\n'); end
            fprintf(m,'-------Tx-------\n');
            fprintf(m,'Nsymb = %d. Samp x symb = %d. Sampling freq = %.1f [GHz]\n',...
                Nsymb, Nt, fs); fprintf(m,'\n\n');
            if Nch == 1, chstr='1ch'; else, chstr=sprintf('WDM: %dx%.1f=%.3f [GHz]\n',Nch,deltaf,Nch*deltaf);end
            fprintf(m,'%s @ %.1f [Gbd]\n, %s, P = %.2f [dBm]\n',modfor,Symbolrate,chstr,PdBm); fprintf(m,'\n');
        end
        
        
        
        
        
        
        % frequency in GHz
        omega        = 2*pi*GSTATE.FN(1);
        omega2       = 2*pi*GSTATE.FN(2);
        delta_omega  = omega2-omega;
        
        
        
        Plin = 10.^(PdBm/10);      % [mW]
        lamv = getlambda(lam,spac,Nch); % carrier wavelengths [nm]
        
        %%%%%Channels with different powers
        
        Plinvec=[1e-4, Plin];
        
        lamv = [lam, lam+spac];
        E = lasersource([1e-4 Plin],lamv);  % electric field. -30 -> no SPM
        [Ex,Ey] = pbs(E); % split in two orthogonal polarizations
        
        for k=1:Nch
            rng(1e5+k-ceil(Nch/2)); % ceil(..): keep seed for increasing Nch
            
            
            patx(:,k) = datapattern(Nsymb,'randn');
            paty(:,k) = datapattern(Nsymb,'randn');
            
            [elecx, normx] = digitalmod(patx(:,k),modfor,Symbolrate,'rootrc',tx);
            [elecy, normy] = digitalmod(paty(:,k),modfor,Symbolrate,'rootrc',tx);
            
            
            %% channel with different power
            Ex  = iqmodulator(Ex, elecx,struct('nch',k,'norm',normx));
            Ey  = iqmodulator(Ey, elecy,struct('nch',k,'norm',normy));
            
            
        end
        
        Esp = pbc(Ex,Ey); % combine creating a PDM signal
        
        
        % initializing the plates method
        cycles = 1;
        M = zeros(2,2,ft.nplates);
        
        M1 = zeros(2,2,ft.nplates);
        
        z= 1i*((pmdpar)*sqrt(lcorr*1e-3)*sqrt(3*pi/8)*1e-3* omega);
        z1= 1i*((pmdpar)*sqrt(lcorr*1e-3)*sqrt(3*pi/8)*1e-3*omega2);    %PMD
        
        
        E= Esp;
        
        for ns= 1:Nspan
            fprintf('Simulating Span...........%d',ns);
            Enn=E;
            for cc=1:cycles
                [E,OUT{ns}]=fiber(E,ft);
                
                E = ampliflat(E,amp);
                Ec = E;
                Et = E;
                
                
                TB=eye(2);                  % preallocating
                TB1=eye(2);
                
                for nn=ns:-1:1
                    % Reversing the polarization
                    for n =1:ft.nplates
                        
                        ft1.lin.matin(:,:,n) = OUT{nn}.lin.matin(:,:,(ft.nplates-n+1));
                        
                        ft1.lin.db0(n,:)     = -OUT{nn}.lin.db0(ft.nplates-n+1,:);
                        
                        
                        y = OUT{ns}.lin.matin(:,:,n);
                        x = OUT{ns}.lin.db0(n);
                        
                        q=sqrt(lcorr*1e-3);
                        
                        w =[exp(1i*(x*q)+z/2),exp(1i*-(x*q)-z/2)];    % creating the exponential of propagation constant plus PMD
                        
                        k = diag(w);                  % creating a diagonal matrix of W
                        
                        
                        %concatenating ft.nplates sections
                        
                        M(:,:,n)  = y*k*y';                   % M is the propagation matrix resulting from
                        
                        
                        w1=[exp(1i*(x*q)+z1/2),exp(1i*-(x*q)-z1/2 )];
                        k1 = diag(w1);
                        
                        M1(:,:,n) = y*k1*y';
                        
                        
                        
                        TB = M(:,:,n)* TB ;         % multiplying each waveplate propagation matrix
                        
                        TB1= M1(:,:,n)* TB1;          % getting M1 to be used to get the finite difference
                        
                        
                        
                    end
                    
                    ft1.disp=-ft.disp;
                    ft1.alphadB=0;
                    ft1.n2=0;
                    ft1.pmdpar=-ft.pmdpar;
                    ft1.slope  =-ft.slope;
                    
                    
                    Ec=fiber(Ec, ft1);
                end
                
                % Getting the derivative of M with respect to omega using finite
                % difference method
                
                
                % the derivative
                slope = (TB1-TB)/delta_omega;
                
                G = 1i*slope*TB';        % Where G =jdM/dw*transpose+conjugate(M)
                [V,D]= eig(G);
                
                B=abs(D(1,1)-D(2,2));     % the absolute value of the difference btw the DGD values
                e(cc) = B;
                
                
                
                %checking if the simulation is right
                averageDGD =mean(e);       % average value of experimental DGD
                
                theoriticalDGD= (pmdpar*sqrt(ft.length*1e-3));        %theoritical DGD value
                
                DGDmean=averageDGD * 1000;
                
                
                
                % Receiver
                for k=1:Nch
                    Em = sep2uni(Ec,k);
                    
                    rsig = rxfrontend(Em,lamv(k),Symbolrate,rx);
                    akhat = rxdsp(rsig,Symbolrate,[patx(:,k) paty(:,k)],rx);
                    %pathat = samp2pat(akhat,modfor,rx);
                    SNRdB(cc,k,pmd,ns) = samp2snr([patx(:,k) paty(:,k)],akhat,modfor);
                    NLIvar_ssfm(cc,k,pmd,ns) = 10*log10(Plinvec(k)) - SNRdB(cc,k,pmd,ns);
                end
                
                
                for m=unique([1 fid]) % if output file exists, write both on file and command window
                    fprintf(m,'\n...................Rx..........'); fprintf(m,'\n');
                    fprintf(m, 'PMD value being simulated......%d',pmdvec(pmd)); fprintf(m,'[ps/Sqrt[KM]]\n');
                    fprintf(m, 'Power being simulated......%d',PdBm); fprintf(m,'[dB]\n');
                    fprintf(m, 'Channel spacing being simulated......%d',deltaf(sp)); fprintf(m,'[GHZ]\n');
                    fprintf(m, 'Simulating Span......%d',ns); fprintf(m,'\n');
                    fprintf(m, 'Simulating cycle......%d',cc); fprintf(m,'\n');
                    fprintf(m,'Ch:\t\t'); for k=1:Nch, fprintf(m,'%d\t',k); end; fprintf(m,'\n');
                    fprintf(m,'SNRdB:\t\t'); for k=1:Nch, fprintf(m,'%.2f\t',SNRdB(cc,k,pmd,ns)); end; fprintf(m,' [dB] \n');
                    fprintf('\n')
                    
                    fprintf('\n')
                    fprintf('\n')
                    fprintf('\n')
                end
                
                
                
                
                %% Formula to check GNLI from Poggiolini formula (equation 40)
                alphaPG = ft.alphadB/2;                                        % Poggiolini alpha is ft.alpha/2 [dB]
                alphaPG_lin= alphaPG/4.343;
                Leff = (1-exp(-2*alphaPG_lin*ft.length))/(2*alphaPG_lin);                           % Effective length
                Leff_a= 1/(2*alphaPG_lin);                                           % Asymptotic effective length
                gamma= (ft.n2 * 2* pi)/(ft.lambda*1e-12 *ft.aeff*1e-12);           % Nonlinear coeffecient
                beta_2= abs(ft.disp)* (ft.lambda^2)/(2*pi*clight*1e-3);        % Dispersion parameter
                
                %kvalue = (Nch-1)/2;  For two channels, kvalue is 1/2 and -1/2
                k1= 1;
                k2= 0;
                
                % using channel spacing in GHz
                X1 = (pi^2)* beta_2 * Leff_a * Symbolrate *1e-6* (k1*deltaf(sp)+ Symbolrate/2);
                Y1 = (pi^2)* beta_2 * Leff_a * Symbolrate * 1e-6*(k1*deltaf(sp)- Symbolrate/2);
                X2 = (pi^2)* beta_2 * Leff_a * Symbolrate *1e-6* (k2*deltaf(sp)+ Symbolrate/2);
                Y2 = (pi^2)* beta_2 * Leff_a * Symbolrate *1e-6* (k2*deltaf(sp)- Symbolrate/2);
                
                
                XX = ((pi^2)* beta_2 * Leff_a *1e-6* Symbolrate^2)/2;
                
                %asinh_function = ((asinh(X1) - asinh(Y1)) + (asinh(X2) - asinh(Y2))) + asinh(XX);
                asinh_function = (asinh(X1) - asinh(Y1)) + asinh(XX);
                
                %(2*(1/8)* (8/9)^2)*
                NLI_Poggiolini_EQ40_lin= (((gamma^2)*(Leff^2)*((Plin/Symbolrate)^3)*(2/3)^3)/(pi*beta_2*Leff_a)) * asinh_function;  % Poggiolini NLI in linear
                NLI_Poggiolini_EQ40_dB= 10*log10(NLI_Poggiolini_EQ40_lin*Symbolrate);        % Poggiolini NLI in dB
                
                
                
                
                %%%% XPM EQ.40
                XPM_Poggiolini_EQ40_lin= gamma^2*Leff^2*Plinvec(2)^2*Plinvec(1)/Symbolrate^3*(2/3)^3/(pi*beta_2*Leff_a)*(asinh(X1) - asinh(Y1)) *Symbolrate;
                XPM_Poggiolini_EQ40_dB= 10*log10(XPM_Poggiolini_EQ40_lin);
                
                
                %%%%% SPM EQ.40
                SPM_Poggiolini_EQ40_lin= gamma^2*Leff^2*(Plinvec(2)/Symbolrate)^3*(2/3)^3/(pi*beta_2*Leff_a)* asinh(XX)*Symbolrate;
                SPM_Poggiolini_EQ40_dB= 10*log10(SPM_Poggiolini_EQ40_lin);
                
                
                
                
                %% Formula to check GNLI from Poggiolini formula (equation 39)
                
                VV1= -1i*(pi^2)*(beta_2/alphaPG_lin) * 1e-6*(k1*deltaf(sp)+ Symbolrate/2) * Symbolrate;
                VV2 = 1i*(pi^2)*(beta_2/alphaPG_lin) *1e-6* (k1*deltaf(sp)+ Symbolrate/2) * Symbolrate;
                QQ1 = -1i*(pi^2)*(beta_2/alphaPG_lin) *1e-6* (k1*deltaf(sp)- Symbolrate/2) * Symbolrate ;
                QQ2=  1i*(pi^2)*(beta_2/alphaPG_lin) *1e-6* (k1*deltaf(sp)- Symbolrate/2) * Symbolrate;
                
                
                %                 YY1= -1i*(pi^2)*(beta_2/alphaPG_lin) *1e-6* (k2*deltaf(sp)+ Symbolrate/2) * Symbolrate;
                %                 YY2 = 1i*(pi^2)*(beta_2/alphaPG_lin) * 1e-6*(k2*deltaf(sp)+ Symbolrate/2) * Symbolrate;
                %                 WW1 = -1i*(pi^2)*(beta_2/alphaPG_lin) *1e-6* (k2*deltaf(sp)- Symbolrate/2) * Symbolrate;
                %                 WW2=  1i*(pi^2)*(beta_2/alphaPG_lin) *1e-6* (k2*deltaf(sp)- Symbolrate/2) * Symbolrate;
                
                KK1= -1i*(pi^2)* beta_2 *1e-6* ((2* alphaPG_lin)^-1)*Symbolrate^2;
                KK2= 1i*(pi^2)* beta_2 * 1e-6*((2* alphaPG_lin)^-1)*Symbolrate^2;
                
                G_k1= (1i* ( dilog(1-VV1) - dilog(1-VV2)))/(2*(pi^2)*(alphaPG_lin^-1)* beta_2) - (1i* ( dilog(1-QQ1) - dilog(1-QQ2)))/(2*(pi^2)*(alphaPG_lin^-1)* beta_2);
                %G_k2= (1i* ( dilog(YY1) - dilog(YY2))/(2*(pi^2)*(alphaPG^-1)* beta_2)) - (1i* ( dilog(WW1) - dilog(WW2))/(2*(pi^2)*(alphaPG^-1)* beta_2));
                
                G_0= (1i* (dilog(1-KK1)-dilog(1-KK2)))/(2*(pi^2)*((2*alphaPG_lin)^-1)* beta_2);
                
                
                NLI_Poggiolini_EQ39_lin= (16/27)*(gamma^2)* (Leff^2)*((Plin/Symbolrate)^3)* (2*G_k1 + G_0);
                NLI_Poggiolini_EQ39_lin_dB= 10*log10(NLI_Poggiolini_EQ39_lin*Symbolrate);
                
                
                
                %%%%Nonlinear epsilon1
                lon=((pi^2)/2)*beta_2*Leff_a* 1e-6*Symbolrate^2*(Nch^2)^(Symbolrate/deltaf(sp));
                function_epsilon=1+(6/ft.length*1e-3)* Leff_a/asinh(lon);
                epsilon= (3/10)* log(function_epsilon);
                
                %%%%Nonlinear epsilon2
                xd=1-ns+ns*harmonic(ns-1);
                lon1=((pi^2)/2)*beta_2*Leff_a* 1e-6*(Symbolrate*Nch)^2;
                function_log(ns)=1+((2*Leff_a)/(ns*ft.length*1e-3)* xd/asinh(lon1));
                epsilon1(ns)= log(function_log(ns))/log(ns);
                
                %%%%Formula for Varience plot
                
                %% ALPHA2
                N = 1;             % spatial modes in SMF
                alphalin      = ft.alphadB/4.343;
                mu(pmd)         = pmdvec(pmd)*sqrt(pi/8);
%               alphadB2(pmd) = (alphalin + 2*pi*deltaf(sp)^2*mu(pmd)^2/N*1e-6);
%               alphaPG2(pmd) = alphadB2(pmd)/2;                                        % Poggiolini alpha is ft.alpha/2 [dB]
                %alphaPG_lin2(pmd)= alphaPG2(pmd); %/4.343;
                alphaPG_lin2(pmd)= (alphalin + (2*pi*deltaf(sp))^2*mu(pmd)^2/N*1e-6);
                Leff2(pmd) = (1-exp(-2*alphaPG_lin2(pmd)*ft.length))/(alphaPG_lin2(pmd));                           % Effective length
                Leff_a2(pmd)= 1/(alphaPG_lin2(pmd));                                           % Asymptotic effective length
                gamma= (ft.n2 * 2* pi)/(ft.lambda*1e-12 *ft.aeff*1e-12);           % Nonlinear coeffecient
                beta_2= abs(ft.disp)* (ft.lambda^2)/(2*pi*clight*1e-3);        % Dispersion parameter
                
                %kvalue = (Nch-1)/2;  For two channels, kvalue is 1/2 and -1/2
                k1= 1;
                k2= 0;
                
                % using channel spacing in GHz
                XX1(pmd) = (pi^2)* beta_2 * Leff_a2(pmd) * Symbolrate *1e-6* (k1*deltaf(sp)+ Symbolrate/2);
                YY1(pmd) = (pi^2)* beta_2 * Leff_a2(pmd) * Symbolrate * 1e-6*(k1*deltaf(sp)- Symbolrate/2);
                XX2(pmd) = (pi^2)* beta_2 * Leff_a2(pmd) * Symbolrate *1e-6* (k2*deltaf(sp)+ Symbolrate/2);
                YY2(pmd) = (pi^2)* beta_2 * Leff_a2(pmd) * Symbolrate *1e-6* (k2*deltaf(sp)- Symbolrate/2);
                
                
                XXX(pmd) = ((pi^2)* beta_2 * Leff_a2(pmd) *1e-6* Symbolrate^2)/2;
                
                %asinh_function = ((asinh(X1) - asinh(Y1)) + (asinh(X2) - asinh(Y2))) + asinh(XX);
                asinh_function2(pmd) = (asinh(XX1(pmd)) - asinh(YY1(pmd))) + asinh(XXX(pmd));
                
                %(2*(1/8)* (8/9)^2)*
                NLI_Poggiolini_EQ40_lin2(pmd)= (((gamma^2)*(Leff2(pmd)^2)*((Plin/Symbolrate)^3)*(2/3)^3)/(pi*beta_2*Leff_a2(pmd))) * asinh_function2(pmd);  % Poggiolini NLI in linear
                NLI_Poggiolini_EQ40_dB2(pmd)= 10*log10(NLI_Poggiolini_EQ40_lin2(pmd)*Symbolrate);        % Poggiolini NLI in dB
                
                
                
                
                %%%% XPM EQ.40
                XPM_Poggiolini_EQ40_lin2(pmd)= gamma^2*Leff2(pmd)^2*Plinvec(2)^2*Plinvec(1)/Symbolrate^3*(2/3)^3/(pi*beta_2*Leff_a2(pmd))*(asinh(XX1(pmd)) - asinh(YY1(pmd))) *Symbolrate;
                XPM_Poggiolini_EQ40_dB2(pmd)= 10*log10(XPM_Poggiolini_EQ40_lin2(pmd));
                
                
                
                
                XPM_alpha= XPM_Poggiolini_EQ40_lin/2/6;   %%divide this by 6 and 2
                XPM_alpha_dB=10*log10(XPM_alpha);
                
                
                XPM_alpha2(pmd)= XPM_Poggiolini_EQ40_lin2(pmd)/2/6;  %%divide this by 6 and 2
                
                XPM_alpha2_dB(pmd)=10*log10(XPM_alpha2(pmd));
                
                
                
                TOT_VAR(pmd) = ((2*N+1)/2*N)*((2*N+1)*XPM_alpha+ (((2*N-1)*(alphalin + 2*pi*(deltaf(sp)^2*mu(pmd)^2)/N*1e-6))/alphalin) * XPM_alpha2(pmd)); %%check the brackets
                TOTAL_VAR_dB(pmd)= 10*log10(TOT_VAR(pmd));
                
                
                %%%just check with this
                TOT_VAR_1(pmd) = ((2*N+1)/2*N)*((2*N+1)*XPM_alpha+ (2*N-1)*(alphalin + (2*pi*deltaf(sp))^2*mu(pmd)^2/N*1e-6)/alphalin * XPM_alpha2(pmd)); %%check the brackets
                TOTAL_VAR_dB_1(pmd)= 10*log10(TOT_VAR_1(pmd));
                
               %%% Chiara's
                TOT_VAR_2(pmd) = ((2*N+1)/2*N)*((2*N+1)*XPM_alpha+ (2*N-1)*(alphalin + (2*pi*deltaf(sp))^2*mu(pmd)^2/N*1e-6)/alphalin*XPM_alpha2(pmd)); %%check the brackets
                TOTAL_VAR_dB_2(pmd)= 10*log10(TOT_VAR_2(pmd));
                
                log_epsilon(ns)=(1+epsilon1(ns))*10*log10(ns);
                
                E = Enn;
                
            end
            E = Et;
            
            fprintf('\n')
        end
        
    end
    
end




save PMD_V2_50km.mat    %%% saving all variables
