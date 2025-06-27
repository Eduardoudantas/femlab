%Universidade Federal do Para - UFPa
% ITEC - Instituto de Tecnologia 
% Grupo de Dinâmica Instrumentação Processamento de sinais e imagens 
% Curso de Engenharia Civil 
% Programa didatico para Analise Linear pelo Método dos Elementos Finitos 
% Autor:  Eduardo Uchôa Dantas  
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------

% Rotina de inputs do programa para análise estática do elemento   Q4 

function Inputq4model

%       Lados da placa 
%      ------------A------------
%      |                       |
%      |                       |
%      |                       |
%      B                       B
%      |                       |
%      |                       |
%      |                       |
%      ------------A------------


% Discretização da área da placa 

% dados de entrada
numGDLno = 3;    % numero de graus de liberdade por no 
A = 2; % Dimensão A
B = 2; % Dimensão B 


nx = 2;  % Número de elementos em x 
ny = 2;  % Número de elementos em y

numNos = (nx+1)*(ny+1);      % numero de nos da estrutura
numElems = nx*ny;    % numero de elementos da estrutura

Lx = A;  
Ly = B;

dx = Lx / nx;  % Número de divisões em x
dy = Ly / ny;  % Número de divisões em y 

% Gerador das coordenadas dos nós 

% gera as coordenadas dos nos
k = 1;
for j = 0: ny
   for i = 0: nx
      coords(k,1) = dx * i;
      coords(k,2) = dy * j;
      
      k = k + 1;
   end
end
coords


% Restrições nodais 

% restricoes nodais (condicoes de apoio)
%         Dx Dy 
restrs = zeros(numNos,3); % Inicializa a matriz de restrições 

% restringe a face a esquerda da placa
for j = 1: nx + 1: (nx + 1)*ny + 1
   restrs(j,:) = [1  1  1];
end

%restringe a face a direita da placa
for j = nx+1: nx + 1: (nx + 1)*(ny + 1)
   restrs(j,:) = [0  0  1];
end

% restringe a face abaixo 
for j = 1: ny + 1
   restrs(j,:) = [1  1  1];
end

% restringe a face acima 
for j = (nx + 1)*ny + 1:(nx+1)*(ny+1)
   restrs(j,:) = [0  1  0];
end
restrs(numNos,:) = [0 1 1];
% Quaisquer outras restrições podem ser facilmente alocadas a um nó qualquer ( n ) da estrutura
% realizando a alocação de um valor em uma dada linha/coluna da matriz restrs

restrs

% Gera a conectividade dos elementos 
el = 1;
for i = 1: nx
   for j = 1: ny
      noI = i + (j - 1)*(nx+1);
      noJ = noI + 1;
      noK = noI + nx + 2;
      noL = noI + nx + 1;
      conect(el,1) = noI;
      conect(el,2) = noJ;
      conect(el,3) = noK;
      conect(el,4) = noL;
      el = el + 1;
      
   end
end
conect


% Gera as cargas nodais do elemento 

% cargas nodais
%              Fz   Mx   My
cargasNos(numNos,:) = [10 0 0];  % kN 

% Adiciona uma carga pontual em todos os nós do elemento
% Qualquer carga pontual pode facilmente ser alocada usando a matriz cargasNos


% dados dos elementos 
E = 2e6;        % Modulo de deformacao longitudinal
nu = 0.3;       % coeficiente de Poisson
t = 0.2;        % espessura do material
ordint = 2;  % ordem de integracao numerica
bz = 0;         % (kN/m2) carga distribuida por unidade de área na direcao z;

for el = 1: numElems
%                     nNos   
   dadosElem(el,:) = [E nu t ordint bz];
   tipoElem(el,:) = 'Quadrilateralq4mindlin';
end


% Realiza a chamada do arquivo principal de análise estrutural 

[matD,matPr,matS,matE] = AnalisePrincipalFEM (numGDLno,numNos,numElems,coords,restrs,cargasNos,conect,tipoElem,dadosElem );

matD   % Matriz de deslocamentos nodais 
matPr  % Matriz forças nodais, incluindo reações de apoio 
matS   % Matriz de tensões nos elementos 
matE   % Matriz de deformações nos elementos 


