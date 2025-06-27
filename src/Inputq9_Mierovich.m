%Universidade Federal do Para - UFPa
% ITEC - Instituto de Tecnologia 
% Grupo de Dinâmica Instrumentação Processamento de sinais e imagens 
% Curso de Engenharia Civil 
% Programa didatico para Analise Linear e Dinâmica  pelo Método dos Elementos Finitos 
% Autor:  Eduardo Uchôa Dantas  
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------

% Rotina de inputs do programa para análise estática do elemento   Q9 

 %Inputq9rev8x8

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







%dados de entrada 


 numGDLno = 3; % Número de graus de liberdade por nó

 A = 1.0;   % Dimensão do lado A da placa em metros 
 B = 0.05;   % Dimensão do lado B da placa em metros 
 NumElementosx =30 ;
 NumElementosy =20; 

 
  %------------- calcula o número de nós do elemento e o número de elementos 
 
 numnosx= (2*NumElementosx)+1;
 numnosy= (2*NumElementosy)+1;
 numNos =  numnosx*numnosy;   % Realizar path test desse algorítmo 
 
 
 %( Recomendação Prof Edilson ) 
 % utilizar comando unique a partir da matriz de conectividade para encontrar a quantidade de nós 
 

 numElems = NumElementosx*NumElementosy;  % numero de elementos da estrutura 


Lx=A;  % Lado x
Ly=B;  % Lado y


dx= Lx/(2*NumElementosx); % divisões em x para o elemento q9 

dy=Ly/(2*NumElementosy);  % divisões em y para o elemento q9 

% loop que gera as coordenadas dos nós 
k = 1;
for j = 0: 2*NumElementosy
   for i = 0: 2*NumElementosx
      coords(k,1) = dx * i;
      coords(k,2) = dy * j;
      
      k = k + 1;
   end
end
coordsz=zeros(numNos,1); % Vetor de coordenadas z do elemento 

%-------------------------------------RESTRIÇÕES NODAIS---------------------------------------------% 

%inicializa a matriz de restrições nodais com zeros.
% A matriz de restrições nodais indica onde estão os apoios da estrutura. 

restrs= zeros(numNos,3);  % Gdls Wz Rx Ry 

% restringe a face a esquerda da chapa 
for j = 1: (2*NumElementosx + 1): (2*NumElementosx + 1)*(2*NumElementosy + 1)     
   restrs(j,:) = [1  1  1];
end

%restringe a face a direita do da chapa
%for j = NumElementosx+1: NumElementosx + 1: (NumElementosx + 1)*(NumElementosy + 1)   % Falta Debugar 
  % restrs(j,:) = [0  0  0];
%end
% restringe a face abaixo 
%for j = 1: NumElementosy + 1  % Falta Debugar 
%   restrs(j,:) = [0  0  0];
%end

% restringe a face acima 
%for j = (NumElementosx + 1)*NumElementosy + 1:(NumElementosx+1)*(NumElementosy+1)  % Falta Debugar 
%   restrs(j,:) = [0  0  0];
%end

%Adiciona manualmente as restrições nodais na matriz restrs para posicionar as condições de contorno 

%restrs(1,:)= 1;
%restrs(6,:) = 1;
%restrs(11,:) = 1;
%restrs(16,:) = 1;
%restrs(273,:) = 1;






%------------------------------CONECTIVIDADE-----------------------------------

 % Função para calcular o índice do nó
    node_index = @(i, j) (i * (2*NumElementosx+1) + j + 1);

 %  inicializa a Matriz de conectividade
    conect = zeros(NumElementosx * NumElementosy , 9);   

 % Inicializa o loop por elemento que gera a conectividade da malha   
  element = 1;
    for i = 0:(NumElementosy-1)
        for j = 0:(NumElementosx-1)
            element_nodes = [
                node_index(2*i, 2*j),       % Nó 1
                node_index(2*i, 2*j + 2),   % Nó 2
                node_index(2*i + 2, 2*j + 2), % Nó 3
                node_index(2*i + 2, 2*j),   % Nó 4
                node_index(2*i, 2*j + 1),   % Nó 5
                node_index(2*i + 1, 2*j + 2), % Nó 6
                node_index(2*i + 2, 2*j + 1), % Nó 7
                node_index(2*i + 1, 2*j),   % Nó 8
                node_index(2*i + 1, 2*j + 1)]; % Nó 9 ( central ) 
            conect(element, :) = element_nodes;
            element = element + 1;   % Insere na matriz conectivity_matrix os nós nas colunas do elemento
        end
    end  
   
    % Essa função retorna a matriz conect com uma linha por elemento, indicando o índice de cada um dos seus nós


   

% --------------------------------------------------------------------------------------------------% 


%------------------------------CARGAS NODAIS--------------------------------------------------------%


%  Insere um carregamento constante em todos os nós 
% Essa função pode ser manipulada para que os carregamentos sejam localizados em nós específicos
% realizando uma alocação de carga pontual com ex :  cargasNos(nó , : ) = [ Fz , Mx , My ] 


% Inicializa a matriz cargas nos ( que armazenará as cargas nodais 

cargasNos = zeros(numNos,numGDLno);

% Posiciona o carregamento em todos os nós 
%                    Fz   Mx   My
%cargasNos(numNos,:) = [10 0 0];

cargasNos(:,1) = 0; % Posiciona o carregamento Fz 
cargasNos(:,2) = 0; % Posiciona o carregamento Mx 
cargasNos(:,3) = 0; % Posiciona o carregamento My 

%cargasNos(17,1) = -10;
%cargasNos(15,3) = -10
%cargasNos(289,1) = -10;
% Remove as cargas dos apoios 
%cargasNos(1,:) = 0;
%cargasNos(6,:) = 0;
%cargasNos(11,1) = 0;
%cargasNos(16,1) = 0;
%cargasNos(21,1) = 0;


%---------------------------------------------------------------------------------------------------%

%------------------------------------------------DADOS DO ELEMENTO----------------------------------% 

E = 2.1e8; % Modulo de deformação longitudinal do material 
nu=0.3;  % Coeficiente de Poisson 
t= 0.01;  % Espessura da chapa 
ordint= 4; % Ordem de integração numérica
bz=0;     % Carga por área distribuida na direção Z ( kN/m2 ) 
ro=7.850;   %(g/cm³)


% Loop que aloca os dados a todos os elementos na matriz Dados Elem e tipo Elem 

for el = 1:numElems

    dadosElem(el,:) = [E nu t ordint bz ro];    % Dados ao elemento 
    tipoElem(el,:) = 'Quadrilateralq9mindlin';  % tipo de elemento 

end    


%close all;
%figure (1);
%DesenhaMalha(numNos, numElems, coords, conect);








% Chama o arquivo que realiza a análise estrutural ( Programa geral de análise  ) 

[matD,matPr,matS,matE,matFreq,matFmodal2] = AnalisePrincipalFEM (numGDLno,numNos,numElems,coords,restrs,cargasNos,conect,tipoElem,dadosElem );


format short g

% Reorganiza as coordenadas nodais para os nós do elemento 
Coordenadasnodais = [coords,(1:numNos)'];
disp('---------Coordenadas Nodais ----------');
disp(Coordenadasnodais);




% Reorganiza a matriz de deslocamentos para apresentar os resultados 
% matD= reshape(matD,numGDLno,numNos)';
 Deslocamentos = [matD,(1:numNos)'];
 disp ('-----------DESLOCAMENTOS-------------');
 disp('        dz        rx       ry       nó   ');
disp(Deslocamentos);


% % Reorganiza a matriz de forças nodais para apresentar os resultados 
 %matPr = reshape(matPr,numGDLno,numNos)';
 Forcasnodais = [matPr,(1:numNos)'];
 disp('---------FORÇAS NODAIS----------');
 disp(Forcasnodais);

%Apresenta a matriz de frequencias naturais 

 disp('------Frequencias naturais------')
 disp(matFreq);
 
matFmodal1 = matFmodal2;


matD;   % Matriz de deslocamentos nodais 
matPr;  % Matriz forças nodais, incluindo reações de apoio 
matS;   % Matriz de tensões nos elementos 
matE;   % Matriz de deformações nos elementos 
matFreq; % Matriz de frequencias naturais 
%matFmodal; % Matriz de formas modais 
% Rotina que realiza a chamada para desenho da estrutura indeformada 




%-------------------PLOTS-----------------------------------------%

figure (1);
desenhaestinderformadaq93d(coords,coordsz,numNos,restrs,conect)

figure(3);
desenhaesdeformadaq93drev(coords,coordsz,numNos,restrs,conect,Deslocamentos)

figure(4);
desenhaformasmodais(coords,coordsz,numNos,restrs,conect,matFmodal1)

% Cria um arquivo excel com os resultados
nome_arquivo = 'Resultados_static_rev1.xlsx';


%writematrix(Deslocamentos,nome_arquivo);





