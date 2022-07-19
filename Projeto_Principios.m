
% ************************************************************************%
% PROGRAMA-: Este programa faz a leitura de um sinal de m�sica stereo,    % 
%            obtem a correspondente transformada de Fourier.              %  
%            Este gera o gr�fico de um dos canais no dom�nio do tempo     %
%            e no dom�nio da frequ�ncia. E tem tamb�m a modula��o e       %
%            desmodu��o do sina, e todas as etapas tem sueus respectivos  % 
%            gr�ficos.                                                    %
% ************************************************************************%

clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%% Passo 1: Acessar o arquivo %%%%%%%%%%%%%%%%%%%%%%%

[sinal,Fs]=audioread('audio.wav'); % Recebendo o sinal 
canal = sinal(:,1);                % Canal do som
N = length(canal);                 % Tamanho do canal 

%%%%%%%%%%%%%%%%% Passo 2: Gerar gr�fico de um dos canais %%%%%%%%%%%%%%%%% 

tempo = (0:1/Fs:1/Fs*(N-1));      % Vetor com o tempo

figure(1)                         % Definindo uma tela para cada subplot
subplot(211);                     % Para deixar 2 gr�ficos na mesma figura
plot(tempo,canal);                % Plotar o Canal
title('Analise temporal do sinal m(t)'); % T�tulo do gr�fico
xlabel('Tempo (s)');                     % Texto do eixo X
ylabel('Amplitude');                     % Texto do eixo Y

%%%%%%%%%%%%%%%%%%%%%%%% Passo3: An�lise Espectral %%%%%%%%%%%%%%%%%%%%%%%%

M = fft(canal);                          % Tranformada de Fourier do sinal
w = 2*pi*Fs*(-round(N/2)+1:round(N/2))'/N; % Frequ�ncia 
y = ifft(M);                               % Inversa da tranformada

subplot(212);
plot(w,fftshift(abs(M)));                   % Plotar o espectro 
title('Espectro de Fourier do sinal m(t)'); % T�tulo do gr�fico
xlabel('Frequ�ncia (rad/s)');               % Texto do eixo X
ylabel('M�dulo |X(\omega)|');               % Texto do eixo Y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Passo 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wc = 2*pi*1187.41;                        % frequ�ncia de corte
pontoA = canal'.*2.*cos(wc .* tempo);     % Modula��o ponto A
A = fft(pontoA);                          % Tranformada Fourier do ponto A

%Plotagem An�lise temporal do sinal no ponto a
figure(2)
subplot(211);
plot(tempo,pontoA);                        
title('Analise temporal do sinal ponto A');
xlabel('Tempo (s)');
ylabel('Amplitude');

%Plotagem Espectro de Fourier do sinal no ponto a
subplot(212);
plot(w,fftshift(abs(A)));
title('Espectro de Fourier do sinal m(t) no ponto A');
xlabel('Frequ�ncia (rad/s)');
ylabel('M�dulo |X(\omega)|');

%Ponto B:
% No ponto b como a primeira parte desconsidera o ru�do
% temos que : pontoB = pontoA

%Ponto C:
pontoC = pontoA .*cos(wc * tempo);        % Modula��o ponto C
C = fft(pontoC);                          %Tranformada Fourier do ponto C

%Plotagem An�lise temporal do sinal no ponto c
figure(3)
subplot(211);
plot(tempo,pontoC);
title('Analise temporal do sinal m(t) no ponto C');
xlabel('Tempo (s)');
ylabel('Amplitude');

%Plotagem Espectro de Fourier do sinal no ponto c
subplot(212);
plot(w,fftshift(abs(C)));
title('Espectro de Fourier do sinal m(t) no ponto C');
xlabel('Frequ�ncia (rad/s)');
ylabel('M�dulo |X(\omega)|');


%%%%%%%%%%%%%%%%%%%%%%%%%%% Passo 5 e 6: Filtro %%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1.5;                         % Ganho do filtro        
n = 3;                           % Ordem do filtro
wo = 7064;                       % Frequ�ncia do filtro

aux = (wo./(i.*w+wo));           % Parte da conta do filtro
FPB = k .* pow2((aux),n);        % Fun��o tranfer�ncia do filtro

X = (fftshift(FPB))'.*C;         % O sinal demodulado passando pelo filtro
x = real(ifft(X));               % Transformada Fourier do sinal filtrado   

%Plotagem do filtro na frequ�ncia
figure(4);
plot(w,FPB);
title('Espectro de Fourier do filtro');
xlabel('Frequ�ncia (rad/s)');
ylabel('M�dulo |X(\omega)|');

%Plotagem An�lise temporal do sinal no ponto d
figure(5);
subplot(211);
plot(tempo,x);
title('An�lise temporal do sinal  m(t) ap�s o filtro');
xlabel('Tempo (s) ');
ylabel('Amplitude');

%Plotagem Espectro de Fourier do sinal no ponto d
subplot(212);
plot(w,fftshift(abs(X)));
title('Espectro de Fourier do sinal m(t) ap�s o filtro');
xlabel('Frequ�ncia (rad/s)');
ylabel('M�dulo |X(\omega)|');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Passo 7: Ru�dos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = max(canal);                    % O valor m�ximo do canal

% Gerando ru�do baixo aleat�rio
ruido_baixo = pontoA + (f/500 * randi([-100 100],1,N)); 
Y_BAIXO = fft(ruido_baixo);        % Transformada do sinal com ru�do baixo

% Gerando ru�do alto aleat�rio
ruido_alto = pontoA + (f/200 * randi([-100 100],1,N ));
Y_ALTO = fft(ruido_alto); % Transformada do sinal com ru�do alto

%Plotagem An�lise temporal do sinal com ruido baixo
figure(6);
subplot(211);
plot(tempo,ruido_baixo);
title('An�lise temporal do sinal m(t) com ru�do baixo');
xlabel('Tempo (s) ');
ylabel('Amplitude');

%Plotagem Espectro de Fourier do sinal com ruido baixo
subplot(212);
plot(w,abs(fftshift(Y_BAIXO)));
title('Espectro de Fourier do sinal m(t) com ru�do baixo');
xlabel('Frequ�ncia (rad/s)');
ylabel('M�dulo |X(\omega)|');

%Plotagem An�lise temporal do sinal com ruido alto
figure(7);
subplot(211);
plot(tempo,ruido_alto);
title('An�lise temporal do sinal m(t) com ru�do alto');
xlabel('Tempo (s) ');
ylabel('Amplitude');

%Plotagem Espectro de Fourier do sinal com ruido alto
subplot(212);
plot(w,fftshift(abs(Y_ALTO)));
title('Espectro de Fourier do sinal m(t) com ru�do alto');
xlabel('Frequ�ncia (rad/s)');
ylabel('M�dulo |X(\omega)|');

%%%%%%%%%%%%%%%% Passo 8: Sinal Demodulador com ruido alto %%%%%%%%%%%%%%%%

% Demodula��o com ru�do alto ponto B
pontoB_ruido_alto = pontoA + ruido_alto;
% Transformada do sinal com ru�do alto ponto B 
B_RUIDO_ALTO = fft(pontoB_ruido_alto);    

%Plotagem An�lise temporal do sinal com ruido alto ponto B
figure(8)
subplot(211);
plot(tempo,pontoB_ruido_alto);
title('Analise temporal do sinal m(t) no ponto B com ruido alto');
xlabel('Tempo (s)');
ylabel('Amplitude');

%Plotagem Espectro de Fourier do sinal com ruido alto ponto B
subplot(212);
plot(w,fftshift(abs(B_RUIDO_ALTO)));
title('Espectro de Fourier do sinal m(t) no ponto B com ruido alto');
xlabel('Frequ�ncia (rad/s)');
ylabel('M�dulo |X(\omega)|');

%Ponto C - RUIDO ALTO:
% Demodula��o com ru�do alto ponto C
pontoC_ruido_alto = pontoB_ruido_alto .*cos(wc * tempo); 
% Transformada do sinal com ru�do alto ponto C
C_RUIDO_ALTO = fft(pontoC_ruido_alto);

%Plotagem An�lise temporal do sinal com ruido alto ponto C
figure(9)
subplot(211);
plot(tempo,pontoC_ruido_alto);
title('Analise temporal do sinal m(t) no ponto C com ruido alto');
xlabel('Tempo (s)');
ylabel('Amplitude');

%Plotagem Espectro de Fourier do sinal com ruido alto ponto C
subplot(212);
plot(w,fftshift(abs(C_RUIDO_ALTO)));
title('Espectro de Fourier do sinal m(t) no ponto C com ruido alto');
xlabel('Frequ�ncia (rad/s)');
ylabel('M�dulo |X(\omega)|');

%%%%%%%%%%%%%%%%%%%%%%%% Filtro %%%%%%%%%%%%%%%%%%%%%%%%

% Transformada Fourier do sinal com ru�do alto filtrado 
X_RUIDO_ALTO = (fftshift(FPB))'.*C_RUIDO_ALTO; 

% O sinal original com ru�do alto passando pelo filtro
x_ruido_alto = real(ifft(X_RUIDO_ALTO));

%%%%%%%%%%%%%%%%%%%%%%%% Voltando o sinal original %%%%%%%%%%%%%%%%%%%%%%%%

%Plotagem An�lise temporal do sinal demodulado com ruido alto
figure(10);
subplot(211);
plot(tempo,x_ruido_alto);
title('An�lise temporal m(t) com ru�do alto ap�s o filtro');
xlabel('Tempo (s) ');
ylabel('Amplitude');

%Plotagem Espectro de Fourier do sinal demodulado com ruido alto
subplot(212);
plot(w,fftshift(abs(X_RUIDO_ALTO)));
title('Espectro de Fourier sinal m(t) com ruidoalto ap�s o filtro');
xlabel('Frequ�ncia (rad/s)');
ylabel('M�dulo |X(\omega)|');

%%%%%%%%%%%%%%%% Passo 8: Sinal Demodulador com ruido baixo %%%%%%%%%%%%%%%
% Demodula��o com ru�do baixo ponto B
pontoB_ruido_baixo = pontoA + ruido_baixo;
% Transformada do sinal com ru�do baixo ponto B 
B_RUIDO_BAIXO = fft(pontoB_ruido_baixo);

%Plotagem An�lise temporal do sinal com ruido baixo ponto B
figure(11)
subplot(211);
plot(tempo,pontoB_ruido_baixo);
title('Analise temporal do sinal m(t) no ponto B com ruido baixo');
xlabel('Tempo (s)');
ylabel('Amplitude');

%Plotagem Espectro de Fourier do sinal com ruido baixo ponto B
subplot(212);
plot(w,fftshift(abs(B_RUIDO_BAIXO)));
title('Espectro de Fourier do sinal m(t) no ponto B com ruido baixo');
xlabel('Frequ�ncia (rad/s)');
ylabel('M�dulo |X(\omega)|');

%Ponto C - RUIDO BAIXO:

% Demodula��o com ru�do baixo ponto C
pontoC_ruido_baixo = pontoB_ruido_baixo .*cos(wc * tempo);
% Transformada do sinal com ru�do baixo ponto C
C_RUIDO_BAIXO = fft(pontoC_ruido_baixo);

%Plotagem An�lise temporal do sinal com ruido baixo ponto C
figure(12)
subplot(211);
plot(tempo,pontoC_ruido_baixo);
title('Analise temporal do sinal m(t) no ponto C com ruido baixo');
xlabel('Tempo (s)');
ylabel('Amplitude');

%Plotagem Espectro de Fourier do sinal com ruido baixo ponto C
subplot(212);
plot(w,fftshift(abs(C_RUIDO_BAIXO)));
title('Espectro de Fourier do sinal m(t) no ponto C com ruido baixo');
xlabel('Frequ�ncia (rad/s)');
ylabel('M�dulo |X(\omega)|');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtro %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transformada Fourier do sinal com ru�do baixo filtrado 
X_RUIDO_BAIXO = (fftshift(FPB))'.*C_RUIDO_BAIXO;

% O sinal original com ru�do baixo passando pelo filtro
x_ruido_baixo = real(ifft(X_RUIDO_BAIXO));

%%%%%%%%%%%%%%%%%%%%%%%% Voltando o sinal original %%%%%%%%%%%%%%%%%%%%%%%%

%Plotagem An�lise temporal do sinal demodulado com ruido baixo
figure(13);
subplot(211);
plot(tempo,x_ruido_baixo);
title('An�lise temporal do sinal m(t) com ru�do baixo ap�s o filtro');
xlabel('Tempo (s) ');
ylabel('Amplitude');

%Plotagem Espectro de Fourier do sinal demodulado com ruido baixo
subplot(212);
plot(w,fftshift(abs(X_RUIDO_BAIXO)));
title('Espectro de Fourier do sinal m(t) com ru�do ap�s o filtro');
xlabel('Frequ�ncia (rad/s)');
ylabel('M�dulo |X(\omega)|');




%%%%%%%%%%%%%%%%%%% Reprodu��o e Salvamento dos �udios %%%%%%%%%%%%%%%%%%%%
% OBS: Para ouvir os �udios em cada momento � preciso apenas descomentar
                        % suas respectivas linhas


%sound(y,Fs);               % audio Original

%sound(x,Fs);               % audio sem ruidos

%audios com ruido 
%sound(x_ruido_alto,Fs);
%sound(x_ruido_baixo,Fs);

%Salvando os novos �udios
audiowrite('audio_sem_ruido.wav',x,Fs);
audiowrite('ruido_baixo.wav',x_ruido_baixo,Fs);
audiowrite('ruido_alto.wav',x_ruido_alto,Fs);
