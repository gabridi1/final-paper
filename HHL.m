%Solu��o de equa��es lineares baseada no algoritmo qu�ntico HHL
%Trabalho de Conclus�o de Curso - Engenharia El�trica
%Universidade de Caxias do Sul (UCS)
%Autor: Guilherme Adamatti Bridi

%---------------------------------------------------------------

format long
clc; close all; clear all

file = fopen('parametros.txt', 'w');
resultados = fopen('resultados.txt', 'w');

%--------------------Par�metros de entrada----------------------

%Sistema 2x2
%B = [4, -2;
%     -2, 4];

%p = [0.6, -0.8]';

%Sistema 4x4
B = [12, -2, 0, 0;
     -2, 12, -6, 0;
     0, -6, 12, -2;
     0, 0, -2, 12];
 
p = [0.7, -0.5, -0.5, -0.1]';

DPtol = 10; %Desvio padr�o admiss�vel na representa��o bin�ria, que pode introduzir erros te�ricos na solu��o do sistema
NQmax = 9; %N�mero de q-bits aceit�vel sendo 13 o suport�vel pela mem�ria do MATLAB

%Escrita dos par�metros de entrada
fprintf(file, '--------------------------------ARQUIVO DOS PARAMETROS---------------------------------------\n\n');
fprintf(file, '------------------------------Parametros de entrada--------------------------------------\n\n');
escreve(B, file, 'B');
escreve(p, file, 'p');
fprintf(file, '\nNQmax = %d\n', NQmax);
fprintf(file, 'DPtol = %3.3f%%\n\n', DPtol);

%--------------Prepara��o dos par�metros do circuito-----------------

%Teste de validade das dimens�es do sistema linear
if length(B) ~= length(p)
    fprintf(file, 'Dimens�es inadequadas');
    fclose(file);
    return
end

eigi = sort(eig(B)); %Autovalores da matriz B original
eigi = abs(eigi(1)); %Autovalor para complementar matrizes de ordem n�o m�ltipla de 2

%Adequa��o de sistemas de ordem n�o m�ltipla de 2
if ceil(log2(length(B))) > log2(length(B))
    for k = length(B): ceil(log2(length(B)))^2 - 1
        B(k, k) = eigi;
        p(k) = 0;
    end
end

%Normaliza��o do sistema linear
B =  B / norm(p);
p =  p / norm(p);

%Teste de hermiticidade
Bo = B';
for k = 1: length(B)
    for l = 1: length(B)
        if abs(B(k, l) - Bo(k, l)) > 0.001
            fprintf(file, 'Matriz n�o hermitiana');
            fclose(file);
            return
        end
    end
end

nb = floor(log2(length(p))); %N�mero de q-bits do registrador beta

eigSave = sort(eig(B)); %Autovalores da matriz B
[eigVectors, eigMatrix] = eig(B); %Matriz diagonalizada

%Verifica se todos os autovalores s�o positivos
for k = 1: 2^nb
    if eigSave(k) <= 0
        fprintf(file, 'Autovalor negativo ou nulo');
        fclose(file);
        return 
    end
end

s = max(eigSave) / min(eigSave); %N�mero de condi��o
nk = floor((2^(NQmax - nb - 1) - 1) / s); %Constante nk de controle que garante NQmax

%Rotina de acomoda��o dos autovalores
for k = 0: nk 
    estadosSave(1) = k + 1; %Posi��o do menor autovalor
    
    for l = 0: 2^nb - 1
        estadosSave(l + 1) = round(estadosSave(1) * eigSave(l + 1) / eigSave(1)); %Posi��o dos demais autovalores
    end

    estadosUnique = unique(estadosSave); %Elimina posi��es repetidas
    nrSave = length(estadosUnique); %N�mero de rota��es controladas
    nqSave = floor(log2(estadosUnique(nrSave) * 2)); %N�mero de q-bits no registrador alpha
    
    %Obt�m os autovalores considerados �nicos
    g = 1;
    h = 1;
    for l = 1: 2^nb
        if estadosSave(l) == estadosUnique(g)
            h = 1;
            eigUnique(g) = eigSave(l);
            g = g + 1;
        
        else
            eigUnique(g - 1) = (h * eigUnique(g - 1) + eigSave(l)) / (h + 1);
            h = h + 1;
        end
    end
    
    %C�lculo da constante t para cada autovalor
    for l = 1: 2^nb
        ti(l) = estadosSave(l) / (eigSave(l) * 2^nqSave) * 2 * pi;
    end
 
    tSave = mean(ti); %t m�dio
    
    %C�lculo do desvio padr�o de t
    dpSave = 0;
    for l = 1: 2^nb
        dpSave = dpSave + (ti(l) - tSave)^2;
    end
    dpSave = 100 * sqrt(dpSave / 2^nb) / tSave;
    
    %Verifica se dp � menor igual ao aceit�vel
    if dpSave <= DPtol 
        %Armazena os par�metros e quebra o la�o
        dp = dpSave;
        nr = nrSave;
        nq = nqSave;
        t = tSave;
        estados = estadosUnique;
        eig = eigUnique;
        break
    end
end

%Verifica se o erro est� abaixo do toler�vel
if dp > DPtol
    fprintf(file, 'Erro n�o toler�vel');
    fclose(file);
    return 
end
    
C = abs(eigSave(1)); %C definido de forma maximizar a probabilidade de leitura 1 no q-bit ancilla

bin = dec2bin(estados', nq); %Estados em bin�rio

for k = 1: nr
    R(k) = 2 * asin(C / eig(k)); %Calcula �ngulos de rota��o original
end

%Matriz auxiliar de posi��es bin�rias
eigMat = zeros(nq, nr); 
for k = 1: nq
    for l = 1: nr
        eigMat(k, l) = bin(l, k);
    end
end

%La�o para obter R, P e controle
for k = 1: nr
    comp = eigMat(:, k);
    menor = 8;
    
    for l = 1: nq
        cont = 0;
        igual = zeros(1, nr);
        
        for g = 1: nr
            if eigMat(l, g) == comp(l) && g ~= k
                    cont = cont + 1;
                    igual(g) = igual(g) - R(k);
            end
        end

        if cont < menor
             menor = cont;
             Po(k) = nq - l;
             Raux = igual;
             controle(k) = eigMat(l, k);
        end
    end
                
    R = R + Raux;
end
        
%Escrita dos par�metros do circuito
fprintf(file, '---------------------------Prepara�ao dos parametros do circuito----------------------------\n\n');
fprintf(file, 'Sistema quantico:\n\n');
escreve(B, file, 'B');
escreve(p, file, 'p');
escreve(eigSave, file, 'Autovalores');
fprintf(file, 's = %3.3f\n', s);
fprintf(file, 'nq = %d\n', nq);
fprintf(file, 'nb = %d\n\n', nb);
fprintf(file, '%d q-bits\n', nq + nb + 1);
fprintf(file, '%d bits classicos\n\n', nb + 1);
fprintf(file, 'nr = %d\n\n', nr);
escreve(eig', file, 'Autovalores unicos');
fprintf(file, 'Estados:\n');

for k = 1: nr
    fprintf(file, '%s <---> %d\n', string(bin(k, :)), estados(k));
end

fprintf(file, '\nDP%% = %3.3f%%', dp);
fprintf(file, '\nt = %4.4f', t);
fprintf(file, '\nC = %4.4f\n\n', C);

%-----------------------Algoritmo-------------------------------

fprintf(file, '--------------------------Matrizes do algoritmo-------------------------------------\n\n');

%Prepara��o do estado qu�ntico p
psy0 = [1, 0]';

for k = 1: nq
    psy0 = kron([1, 0]', psy0);
end

psy0 = kron(psy0, p);

%Transforma��o Hadamard no registrador alpha
psy1 = kron((kron(eye(2, 2), hadamard(nq))), eye(2^nb, 2^nb)) * psy0;
fprintf(file, 'Transformacao Hadamard:\n\n');
fprintf(file, 'h([0, %d])\n\n', nq - 1);

%U-controlled para a matriz e^iBt
psy2 = psy1;

fprintf(file, 'Matrizes hamiltonianas:\n');

for k = 0: nq - 1
    diag = eye(2^nb, 2^nb);
    for l = 1: 2^nb
        diag(l, l) = exp(2^k * eigMatrix(l, l) * j * t);
    end
    
    psy2 = kron(eye(2^(nq - k), 2^(nq - k)), controlled1(kron(eye(2^k, 2^k), eigVectors * diag * eigVectors'))) * psy2;
    
    fprintf(file, '\nQ-bit %d <---> %4.4f e^iBt\n\n', k, 2^(k));
    escreve(eigVectors * diag * eigVectors', file, 'U');
    
end

%Transformada de Fourier Qu�ntica inversa no registrador alpha
fprintf(file, '\nTransformada de Fourier Quantica inversa:\n\n');

for k = 0: nq/2 - 1
     fprintf(file, 'SWAP(%d, %d)\n', k, nq - k - 1);
end

for k = 0: 1: nq - 1
    for l = 0: 1: nq - 2
        if k - l > 0
            fprintf(file, 'p(-%4.4f pi)(%d, %d)\n', (2^(l - k)), k, l);
        end
    end
    
    fprintf(file, 'h(%d)\n', k);
end

psy3 = kron(kron(eye(2, 2), QFT(2^nq)'), eye(2^nb, 2^nb)) * psy2;

%Rota��o controlada no registrador ancilla
fprintf(file, '\nRotacao controlada:\n');

psy4 = psy3;

%Rota��o controlada
for k = 1: nr
    X = [0, 1; 1, 0];
    Pinv = nq - Po - 1;
    
    fprintf(file, '\nPosi�ao %d', Po(k));
    fprintf(file, '\nControle %s', string(controle(k) - 48));
    fprintf(file, '\nRotacao Ry(%5.5f)\n\n', R(k));
    escreve(Ry(R(k)), file, 'Ry');
    
    if controle(k) == '0'
        %Porta X na entrada
        psy4 = kron(kron(eye(2^(Pinv(k) + 1), 2^(Pinv(k) + 1)), X), eye(2^(nb + nq - 1 - Pinv(k)), 2^(nb + nq - 1 - Pinv(k)))) * psy4;
    end
    
    %Rota��o controlada
    psy4 = kron(controlled2(Ry(R(k)), Pinv(k) + 2), eye(2^(nb + nq - 1 - Pinv(k)), 2^(nb + nq - 1 - Pinv(k)))) * psy4;

    if controle(k) == '0'
         %Porta X na sa�da
        psy4 = kron(kron(eye(2^(Pinv(k) + 1), 2^(Pinv(k) + 1)), X), eye(2^(nb + nq - 1 - Pinv(k)), 2^(nb + nq - 1 - Pinv(k)))) * psy4;
    end
end

%Transformada de Fourier Qu�ntica no registrador alpha
fprintf(file, '\nTransformada de Fourier Quantica:\n\n');

for k = nq - 1: -1: 0
    fprintf(file, 'h(%d)\n', k);
    for l = nq - 2: -1: 0
        if k - l > 0
            fprintf(file, 'p(%4.4f pi)(%d, %d)\n', (2^(l - k)), k, l);
        end
    end
end

for k = 0: nq/2 - 1
     fprintf(file, 'SWAP(%d, %d)\n', k, nq - k - 1);
end

psy5 = kron(kron(eye(2, 2), QFT(2^nq)), eye(2^nb, 2^nb)) * psy4;
    
%U-controlled para a matriz e^-iBt
psy6 = psy5;

fprintf(file, '\nMatrizes hamiltonianas inversas:\n');

diag = eye(2^nb, 2^nb);
for k = nq - 1: -1: 0 
    for l = 1: 2^nb
        diag(l, l) = exp(-2^k * eigMatrix(l, l) * j * t);
    end
    
    psy6 = kron(eye(2^(nq - k), 2^(nq - k)), controlled1(kron(eye(2^k, 2^k), eigVectors * diag * eigVectors'))) * psy6;
    
    fprintf(file, '\nQ-bit %d <---> %4.4f e^iBt\n\n', k, 2^(k));
    escreve(eigVectors * diag * eigVectors', file, 'U');
end

%Transforma��o Hadamard no registrador alpha
psy7 = kron((kron(eye(2, 2), hadamard(nq))), eye(2^nb, 2^nb)) * psy6;
fprintf(file, '\nTransformacao Hadamard:\n\n');
fprintf(file, 'h([0, %d])', nq - 1);

%-------------------Resultados----------------------
quantica = zeros(1, 2^nb)';
histograma = zeros(1, 2^(nb + 1))';
histDec = zeros(1, 2^(nb + 1))';

probabilidade = 0;
for k = 1: 2^nb
  histograma(k, 1) = psy7(k)^2; %Localiza a solu��o na posi��o esperada
  histograma(k + 2^nb, 1) = psy7(2^(nq + nb) + k)^2; %Localiza a solu��o na posi��o esperada
  probabilidade = probabilidade + histograma(k + 2^nb, 1) * 100; %Calcula a probabilidade do estado 1, isto �, do sucesso do algoritmo, no q-bit ancilla
  quantica(k, 1) = psy7(2^(nq + nb) + k); %Localiza a solu��o na posi��o esperada
end

for k = 1: 2^(nb + 1)
  histDec(k) = k - 1; %Estados qu�nticos em decimal
end

histBin = dec2bin(histDec', nb + 1); %Estados qu�nticos em bin�rio
eX = categorical(string(histBin)'); %Eixo x do gr�fico
%grafico = bar(eX, abs(histograma)'); %Histograma
%text(1:length(abs(histograma)'), abs(histograma)',num2str(abs(histograma)), 'vert', 'bottom', 'horiz', 'center'); %Mostra valores
%title('Medida') % T�tulo do histrograma
%xlabel('Estados') % Eixo horizontal
%ylabel('Probabilidade') % Eixo vertical
%grid on

classica = linsolve(B, p); %Solu��o cl�ssica
quantica = quantica / C; %Solu��o qu�ntica
erro = 100 * mean(abs((quantica - classica) ./ classica)); %Erro te�rico
classicaN = classica / norm(classica); %Solu��o cl�ssica normalizada
quanticaN = quantica / norm(quantica); %Solu��o qu�ntica normalizada
fidelidade = abs(trace(sqrtm(sqrtm(classicaN * classicaN') * (quanticaN * quanticaN') * sqrtm(classicaN * classicaN'))))^2; %Fidelidade 

%Escrita dos resultados
fprintf(resultados, '--------------------------------ARQUIVO DOS RESULTADOS---------------------------------------\n\n');
fprintf(resultados, '----------------------Solucao---------------------------\n\n');
escreve(classica, resultados, 'Classica');
escreve(quantica, resultados, 'Quantica');
fprintf(resultados, 'Erro teorico = %3.3f%%\n', erro);
fprintf(resultados, 'Fidelidade = %5.5f\n\n', fidelidade);

fprintf(resultados, '----------------------Histograma---------------------------\n\n');
fprintf(resultados, 'p(ancilla = 1) = %3.3f%%\n\n', probabilidade);

for k = 1: 2^(nb + 1)
   fprintf(resultados, '%s <---> %4.4f\n', string(histBin(k, :)), histograma(k));
end

fprintf(resultados, '\n--------------------Vetores de estado--------------------------\n\n');
escreve(psy0, resultados, 'Estado p');
escreve(psy1, resultados, 'Hadamard');
escreve(psy2, resultados, 'U-controlled');
escreve(psy3, resultados, 'QTF-1');
escreve(psy4, resultados, 'Rotacao');
escreve(psy5, resultados, 'QTF');
escreve(psy6, resultados, 'U-controlled-1');
escreve(psy7, resultados, 'Final');

fclose(file);
fclose(resultados);

%--------------------Fun��es---------------------

%Transforma��o Hadamard
function H = hadamard(n)

H1 = [1, 1; 1, -1] / sqrt(2);

H = 1;

for k = 1 : n
    H = kron(H, H1);
end

end

%Escreve um vetor/matriz
function escreve(B, file, nome)

Dim = size(B);

fprintf(file, '%s = \n', nome);

for k = 1: Dim(1)
    fprintf(file, '          ');
    
    for l = 1: Dim(2)
        fprintf(file, ' ');
        fprintf(file, '%5.5f + %5.5f j ', real(B(k, l)), imag(B(k, l)));
    end
    
    fprintf(file, '\n');
end

end

%Operador Ry(teta)
function [R] = Ry(th)

R = [cos(th/2), - sin(th/2); sin(th/2), cos(th/2)];

end

%Operador transformada de Fourier
function m = QFT(n)

w = exp(2 * pi * j / n);
s = zeros(n);
row = 0: (n - 1);

for r = 1 : n
    m(r, :) = row *(r - 1);
end

m = w .^ m;
m = m ./ sqrt(n);

end

%Opera��o controlada tipo 1
function [C] = controlled1(U)

T = size(U);
C = eye (2 * T(1) , 2 * T(2));
k = T(1) + 1;
f = 2 * T(1);
grp = k: f;
C(grp, grp) = U;

end

%Opera��o controlada tipo 2
function [C] = controlled2(U, n)

N = 2^n;
C = eye(N, N);
sobe = 2^(n - 1);

for k = 2: 2: N/2;
    C(k, k) = U(1, 1); 
    C(k, k + sobe) = U(1, 2); 
    C(k + sobe, k) = U(2, 1);
    C(k + N/2, k + N/2) = U(2, 2);
    
end

end