
% ************************************************************************%
% PROGRAMA-: Este programa faz a leitura de um sinal de música stereo,    % 
%            obtem a correspondente transformada de Fourier.              %  
%            Este gera o gráfico de um dos canais no domínio do tempo     %
%            e no domínio da frequência. E tem também a modulação e       %
%            desmodução do sina, e todas as etapas tem sueus respectivos  % 
%            gráficos.                                                    %
% ************************************************************************%

clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%% Passo 1: Acessar o arquivo %%%%%%%%%%%%%%%%%%%%%%%

[sinal,Fs]=audioread('audio.wav'); % Recebendo o sinal 
canal = sinal(:,1);                % Canal do som
N = length(canal);                 % Tamanho do canal 

%%%%%%%%%%%%%%%%% Passo 2: Gerar gráfico de um dos canais %%%%%%%%%%%%%%%%% 

tempo = (0:1/Fs:1/Fs*(N-1));      % Vetor com o tempo

figure(1)                         % Definindo uma tela para cada subplot
subplot(211);                     % Para deixar 2 gráficos na mesma figura
plot(tempo,canal);                % Plotar o Canal
title('Analise temporal do sinal m(t)'); % Título do gráfico
xlabel('Tempo (s)');                     % Texto do eixo X
ylabel('Amplitude');                     % Texto do eixo Y

%%%%%%%%%%%%%%%%%%%%%%%% Passo3: Análise Espectral %%%%%%%%%%%%%%%%%%%%%%%%

M = fft(canal);                          % Tranformada de Fourier do sinal
w = 2*pi*Fs*(-round(N/2)+1:round(N/2))'/N; % Frequência 
y = ifft(M);                               % Inversa da tranformada

subplot(212);
plot(w,fftshift(abs(M)));                   % Plotar o espectro 
title('Espectro de Fourier do sinal m(t)'); % Título do gráfico
xlabel('Frequência (rad/s)');               % Texto do eixo X
ylabel('Módulo |X(\omega)|');               % Texto do eixo Y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Passo 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wc = 2*pi*1187.41;                        % frequência de corte
pontoA = canal'.*2.*cos(wc .* tempo);     % Modulação ponto A
A = fft(pontoA);                          % Tranformada Fourier do ponto A

%Plotagem Análise temporal do sinal no ponto a
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
xlabel('Frequência (rad/s)');
ylabel('Módulo |X(\omega)|');

%Ponto B:
% No ponto b como a primeira parte desconsidera o ruído
% temos que : pontoB = pontoA

%Ponto C:
pontoC = pontoA .*cos(wc * tempo);        % Modulação ponto C
C = fft(pontoC);                          %Tranformada Fourier do ponto C

%Plotagem Análise temporal do sinal no ponto c
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
xlabel('Frequência (rad/s)');
ylabel('Módulo |X(\omega)|');


%%%%%%%%%%%%%%%%%%%%%%%%%%% Passo 5 e 6: Filtro %%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1.5;                         % Ganho do filtro        
n = 3;                           % Ordem do filtro
wo = 7064;                       % Frequência do filtro

aux = (wo./(i.*w+wo));           % Parte da conta do filtro
FPB = k .* pow2((aux),n);        % Função tranferência do filtro

X = (fftshift(FPB))'.*C;         % O sinal demodulado passando pelo filtro
x = real(ifft(X));               % Transformada Fourier do sinal filtrado   

%Plotagem do filtro na frequência
figure(4);
plot(w,FPB);
title('Espectro de Fourier do filtro');
xlabel('Frequência (rad/s)');
ylabel('Módulo |X(\omega)|');

%Plotagem Análise temporal do sinal no ponto d
figure(5);
subplot(211);
plot(tempo,x);
title('Análise temporal do sinal  m(t) após o filtro');
xlabel('Tempo (s) ');
ylabel('Amplitude');

%Plotagem Espectro de Fourier do sinal no ponto d
subplot(212);
plot(w,fftshift(abs(X)));
title('Espectro de Fourier do sinal m(t) após o filtro');
xlabel('Frequência (rad/s)');
ylabel('Módulo |X(\omega)|');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Passo 7: Ruídos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = max(canal);                    % O valor máximo do canal

% Gerando ruído baixo aleatório
ruido_baixo = pontoA + (f/500 * randi([-100 100],1,N)); 
Y_BAIXO = fft(ruido_baixo);        % Transformada do sinal com ruído baixo

% Gerando ruído alto aleatório
ruido_alto = pontoA + (f/200 * randi([-100 100],1,N ));
Y_ALTO = fft(ruido_alto); % Transformada do sinal com ruído alto

%Plotagem Análise temporal do sinal com ruido baixo
figure(6);
subplot(211);
plot(tempo,ruido_baixo);
title('Análise temporal do sinal m(t) com ruído baixo');
xlabel('Tempo (s) ');
ylabel('Amplitude');

%Plotagem Espectro de Fourier do sinal com ruido baixo
subplot(212);
plot(w,abs(fftshift(Y_BAIXO)));
title('Espectro de Fourier do sinal m(t) com ruído baixo');
xlabel('Frequência (rad/s)');
ylabel('Módulo |X(\omega)|');

%Plotagem Análise temporal do sinal com ruido alto
figure(7);
subplot(211);
plot(tempo,ruido_alto);
title('Análise temporal do sinal m(t) com ruído alto');
xlabel('Tempo (s) ');
ylabel('Amplitude');

%Plotagem Espectro de Fourier do sinal com ruido alto
subplot(212);
plot(w,fftshift(abs(Y_ALTO)));
title('Espectro de Fourier do sinal m(t) com ruído alto');
xlabel('Frequência (rad/s)');
ylabel('Módulo |X(\omega)|');

%%%%%%%%%%%%%%%% Passo 8: Sinal Demodulador com ruido alto %%%%%%%%%%%%%%%%

% Demodulação com ruído alto ponto B
pontoB_ruido_alto = pontoA + ruido_alto;
% Transformada do sinal com ruído alto ponto B 
B_RUIDO_ALTO = fft(pontoB_ruido_alto);    

%Plotagem Análise temporal do sinal com ruido alto ponto B
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
xlabel('Frequência (rad/s)');
ylabel('Módulo |X(\omega)|');

%Ponto C - RUIDO ALTO:
% Demodulação com ruído alto ponto C
pontoC_ruido_alto = pontoB_ruido_alto .*cos(wc * tempo); 
% Transformada do sinal com ruído alto ponto C
C_RUIDO_ALTO = fft(pontoC_ruido_alto);

%Plotagem Análise temporal do sinal com ruido alto ponto C
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
xlabel('Frequência (rad/s)');
ylabel('Módulo |X(\omega)|');

%%%%%%%%%%%%%%%%%%%%%%%% Filtro %%%%%%%%%%%%%%%%%%%%%%%%

% Transformada Fourier do sinal com ruído alto filtrado 
X_RUIDO_ALTO = (fftshift(FPB))'.*C_RUIDO_ALTO; 

% O sinal original com ruído alto passando pelo filtro
x_ruido_alto = real(ifft(X_RUIDO_ALTO));

%%%%%%%%%%%%%%%%%%%%%%%% Voltando o sinal original %%%%%%%%%%%%%%%%%%%%%%%%

%Plotagem Análise temporal do sinal demodulado com ruido alto
figure(10);
subplot(211);
plot(tempo,x_ruido_alto);
title('Análise temporal m(t) com ruído alto após o filtro');
xlabel('Tempo (s) ');
ylabel('Amplitude');

%Plotagem Espectro de Fourier do sinal demodulado com ruido alto
subplot(212);
plot(w,fftshift(abs(X_RUIDO_ALTO)));
title('Espectro de Fourier sinal m(t) com ruidoalto após o filtro');
xlabel('Frequência (rad/s)');
ylabel('Módulo |X(\omega)|');

%%%%%%%%%%%%%%%% Passo 8: Sinal Demodulador com ruido baixo %%%%%%%%%%%%%%%
% Demodulação com ruído baixo ponto B
pontoB_ruido_baixo = pontoA + ruido_baixo;
% Transformada do sinal com ruído baixo ponto B 
B_RUIDO_BAIXO = fft(pontoB_ruido_baixo);

%Plotagem Análise temporal do sinal com ruido baixo ponto B
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
xlabel('Frequência (rad/s)');
ylabel('Módulo |X(\omega)|');

%Ponto C - RUIDO BAIXO:

% Demodulação com ruído baixo ponto C
pontoC_ruido_baixo = pontoB_ruido_baixo .*cos(wc * tempo);
% Transformada do sinal com ruído baixo ponto C
C_RUIDO_BAIXO = fft(pontoC_ruido_baixo);

%Plotagem Análise temporal do sinal com ruido baixo ponto C
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
xlabel('Frequência (rad/s)');
ylabel('Módulo |X(\omega)|');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtro %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transformada Fourier do sinal com ruído baixo filtrado 
X_RUIDO_BAIXO = (fftshift(FPB))'.*C_RUIDO_BAIXO;

% O sinal original com ruído baixo passando pelo filtro
x_ruido_baixo = real(ifft(X_RUIDO_BAIXO));

%%%%%%%%%%%%%%%%%%%%%%%% Voltando o sinal original %%%%%%%%%%%%%%%%%%%%%%%%

%Plotagem Análise temporal do sinal demodulado com ruido baixo
figure(13);
subplot(211);
plot(tempo,x_ruido_baixo);
title('Análise temporal do sinal m(t) com ruído baixo após o filtro');
xlabel('Tempo (s) ');
ylabel('Amplitude');

%Plotagem Espectro de Fourier do sinal demodulado com ruido baixo
subplot(212);
plot(w,fftshift(abs(X_RUIDO_BAIXO)));
title('Espectro de Fourier do sinal m(t) com ruído após o filtro');
xlabel('Frequência (rad/s)');
ylabel('Módulo |X(\omega)|');




%%%%%%%%%%%%%%%%%%% Reprodução e Salvamento dos áudios %%%%%%%%%%%%%%%%%%%%
% OBS: Para ouvir os áudios em cada momento é preciso apenas descomentar
                        % suas respectivas linhas


%sound(y,Fs);               % audio Original

%sound(x,Fs);               % audio sem ruidos

%audios com ruido 
%sound(x_ruido_alto,Fs);
%sound(x_ruido_baixo,Fs);

%Salvando os novos áudios
audiowrite('audio_sem_ruido.wav',x,Fs);
audiowrite('ruido_baixo.wav',x_ruido_baixo,Fs);
audiowrite('ruido_alto.wav',x_ruido_alto,Fs);
