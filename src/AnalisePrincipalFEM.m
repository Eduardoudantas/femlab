%Universidade Federal do Para - UFPa
% ITEC - Instituto de Tecnologia 
% Grupo de Dinâmica Instrumentação Processamento de sinais e imagens 
% Curso de Engenharia Civil 
% Programa didatico para Analise Linear pelo Método dos Elementos Finitos 
% Autor:  Eduardo Uchôa Dantas  
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------

%Dados de Inputs 
%numGDLno = numero de graus de liberdade por nó 
%numNosEst = numero de nós da estrutura 
%numElementsEst = número de elementos da estrutura
%coords = coordenadas dos nós da estrutura 
%restrs = restrições dos nós da estrutura 
%cargasNos = cargas dos nós da estrutura 
%conect = conectividade dos elementos da estrutura 
%dadosElem = dados do elemento, nu0, E, G , etc 


%Outputs, Resultados do programa

% matD  = Matriz  de deslocamentos nodais por nó  
% matPr = Matriz de forças nodais, incluindo as reações de apoio
% matS  = Matriz de tensões nos elementos 
% matE  = Matriz de deformações nos elementos 

function [matD,matPr,matS,matE,matFreq,matFmodal1] = AnalisePrincipalFEM (numGDLno,numNos,numElems,coords,restrs,cargasNos,conect,tipoElem,dadosElem )

% Calcula o número de graus de liberdade da estrutura 
numGDLEst = numGDLno*numNos;

% Numera os graus de liberdade dos nós 
[numGDLLivre,GDLNosEst] = obtemGDLNosEst (numNos,numGDLno,restrs);

% Monta o vetor de cargas nodais 
for no = 1:numNos
    for i = 1:numGDLno
        gdl = GDLNosEst(no, i);
        P(gdl,1) = cargasNos(no, i);
    end
end

% Monta a matriz de rigidez da estrutura e o vetor de cargas nodais com zeros 

K = zeros(numGDLEst, numGDLEst);
M = zeros(numGDLEst, numGDLEst);
Peq = zeros(numGDLEst, 1);

for el = 1:numElems 

    % Determina os nos do elementos 

    numNosElem = size(conect, 2); % calcula a quantidade de nos do elemento 
    nosElem = conect(el,:); % captura o nó do elemento 

    % Captura as coordenadas dos nós do elemento 

    for i = 1:numNosElem
        coordsElem(i,:) = coords(nosElem(i),:);
    end     

   switch tipoElem(el,:)

        case 'Quadrilateralq8mindlin'
            kel = ElemQuadrilateralq8mindlin('rigidez',coordsElem,dadosElem(el,:));
            peq = ElemQuadrilateralq8mindlin('cargasequivalentes',coordsElem,dadosElem(el,:));

        case 'Quadrilateralq9mindlin'
            kel = ElemQuadrilateralq9mindlin('rigidez',coordsElem,dadosElem(el,:));
            peq = ElemQuadrilateralq9mindlin('cargasequivalentes',coordsElem,dadosElem(el,:));
            m = ElemQuadrilateralq9mindlin('massa',coordsElem,dadosElem(el,:));

        case 'Quadrilateralq4mindlin'    
            kel = ElemQuadrilateralq4mindlin('rigidez',coordsElem,dadosElem(el,:));
            peq = ElemQuadrilateralq4mindlin('cargasequivalentes',coordsElem,dadosElem(el,:));

        otherwise
            disp('funcao AnalisePrincipalFEM : Elemento inválido');
    end 

   % Determina os graus de liberdade do elemento 

   [GDLNosElem, numGDLElem ] = obtemGDLNosElem (numNosElem,nosElem,numGDLno,GDLNosEst);

   % Adiciona o vetor de cargas nodais equivalentes do elemento ao vetor de cargas da estrutura 

   for r = 1:numGDLElem
       R = GDLNosElem(r);
       Peq(R) = Peq(R) + peq(r);
   end
   P = P + Peq;


   % Realiza a montagem da matriz de rigidez e massa do elemento nas matrizes da estrutura. 

   for r = 1:numGDLElem
       R = GDLNosElem(r);
       for c = 1:numGDLElem
           C = GDLNosElem(c);
           K(R, C) = K(R, C) + kel(r, c); %Utiliza os localizadores R e C e r e c para posicionar kel em K 
           M(R, C) = M(R, C) + m(r, c);   % Utiliza os localizadores R e C e r e c para posicionar m em M 
       end
   end
end   

% Realiza o particionamento das matrizes. 

L = numGDLLivre;
T = numGDLEst;

Kll = K(1:L, 1:L);
Krl = K(L+1:T, 1:L);
Mll = M(1:L, 1:L);
Pl = P(1:L, 1);
Peqr = Peq(L+1:T, 1);

% ------------------------ ANÁLISE ESTÁTICA----------------------------------------------%

% Realiza a solução do sistema de equações. 

Dl = Kll\Pl;
Dr = zeros(T-L, 1);
D = [Dl;
     Dr];

% Calcula as reações de apoio e monta o vetor com as forças nodais 

Pr = Krl*Dl - Peqr;
P = [Pl;
     Pr];

for n = 1: numNos
    for i = 1: numGDLno
        L = GDLNosEst(n, i);
        forcasnos(n, i) = P(L);
        deslocnos(n, i) = D(L);
    end
end 
matD = deslocnos;
matPr = forcasnos;



%---------------------------ANÁLISE MODAL----------------------------------------------%
nmodos = 7; % Número de formas modais da análise 
deslocsmodais = cell(1, nmodos); % Inicializa a célula que irá armazenar os deslocamentos das formais modais 
[shape, w] = eigs(Kll,Mll,nmodos,'sm');

lambda =sqrt(w);
f = real(diag(lambda))/(2*pi);


 % Organiza e ordena as frequencias naturais e as formas modais para montar as matrizes matFreq e matFmodal 
 [freq,k] = sort(f);

 % for i = 1:nmodos 
 %     matFmodal(:,i) = shape(:,k(i));
 % end
 % 
 % matFmodal(L+1:T,:) = zeros(T-L,L);




  shape_n = real(shape(:,k));
  shape=shape_n;

  cont= 1;
  for i = 1:nmodos 
      if (freq(i)==0 || freq(i)==inf)

      else
      freqR(cont) = freq(i);
      formamodal(:, cont) = shape(:, i);
      cont = cont+1;
      matFreq = freqR;
      end
  end 

  %matFreq = freqR;

  % Inicializa a matriz de deslocamentos nodais da análise modal 
  DispNos=zeros(numNos,max(numGDLno));

  for nmodos = 1:size(freqR,2)
            Dl = formamodal(:,nmodos);  
            Dest=zeros(numGDLEst,1); % Inicializa o vetor de deslocamentos da estrutura 
            Dest(1:numGDLLivre)=Dl;  % aloca os deslocamentos livres 
            Dest((numGDLLivre+1):numGDLEst)=Dr;   % Aloca os deslocamentos restringidos 

             % %Normaliza Modos em Relação a Massa para que os deslocamentos da análise modal tenham significado físico. 
             Dest_aux=Dest;
             Mm=Dest_aux'*M*Dest_aux;
             Dest_aux2=(1/sqrt(Mm))*Dest_aux;
             Dest=Dest_aux2;    


            % % result.Numdisp(:,nmodos) = Dest;

   
           
              
            % Realiza a localização dos deslocamentos por nó e por cada grau de liberdade dos nos 
             for k = 1:numNos
                 for j = 1:max(numGDLno)
                     if GDLNosEst(k,j)>0
                        DispNos(k,j) = Dest(GDLNosEst(k,j));
                     else
                        DispNos(k,j) = 0;
                     end
         
                 end
             end
             deslocsmodais{nmodos} = DispNos;
           



           


      
        
         %Reanrranja vetor deslocamentos nodais seguindo a numeração da
         %estrutura para posterior ajuste do modelo caso necessário
         %     cont=1;
         %     for no = 1:numNos
         %         for dof=1:max(numGDLno) 
         %             if GDLNosEst(no,dof)>0
         %                 Dest(cont,nmodos)=result.shape(nmodos).disp(no,dof);
         %                cont=cont+1;
         %             end
         %         end
         %     end
         % 
         % end

 
        
  end

% Inicializa a cell que irá armazenar os deslocamentos modais
            %deslocsmodais = cell(1, nmodos);
            
            %for i = 1:nmodos
            %deslocsmodais{i} = DispNos;
            %end
  matFmodal1 = deslocsmodais{1};
  matFmodal2 = deslocsmodais{2};
  matFmodal3 = deslocsmodais{3};
  matFmodal4 = deslocsmodais{4};
  matFmodal5 = deslocsmodais{5};
  matFmodal6 = deslocsmodais{6};
  matFmodal7 = deslocsmodais{7};
  matFmodal = [matFmodal1  matFmodal2  matFmodal3  matFmodal4 ];
  



 %--------------------CALCULO DE TENSÕES E DEFORMAÇÕES---------------------------------%

 % A partir do conhecimento dos deslocamentos e das forças nodais podemos usa-los para calcular
% as tensões e as deformações nos elementos 


 for el = 1:numElems

    % Determina os nós do elemento 
    numNosElem = size(conect,2); %número de nós do elemento  ( Usar size ? ) 
    nosElem = conect(el,:);
  

    %Determina as coordenadas dos nós do elemento 
    for i = 1:numNosElem 
        coordsElem(i,:) = coords(nosElem(i),:);
    end     
    


    
    % Determina os graus de liberdade do elemento 

    %----- Estudar e implementar esta função------%

    [GDLNosElem, numGDLElem ] = obtemGDLNosElem (numNosElem,nosElem,numGDLno,GDLNosEst);



     % Obtem os deslocamentos a partir do vetor de deslocamentos da estrutura 

     for i=1:numGDLElem
         DElem(i,1)= D(GDLNosElem(i));    % estudar a função que determina os graus de liberdade 
     end     
     

    % Determina as tensões e deformações nos elementos 

    switch tipoElem(el,:)

        case 'Quadrilateralq8mindlin'
            [sigma,epslon] = ElemQuadrilateralq8mindlin('tensoesdeforms',coordsElem,dadosElem(el,:),DElem);

        case 'Quadrilateralq9mindlin'
            [sigma,epslon] = ElemQuadrilateralq9mindlin('tensoesdeforms',coordsElem,dadosElem(el,:),DElem); 

        case 'Quadrilateralq4mindlin'
            [sigma,epslon] = ElemQuadrilateralq4mindlin('tensoesdeforms',coordsElem,dadosElem(el,:),DElem);



        otherwise
            disp('AnalisePrincipalFEM: Elemento invalido');
   end

   
    
   
 end
    % Armazena e ordena as tensões e deformações 

    matE(el,:,:)=epslon;
    matS(el,:,:)=sigma;





















%------------------- FUNÇÕES DO CÓDIGO ----------------------------------------% 

% Rotina que obtem os gdl dos nós da estrutura 



function [numGDLLivre,GDLNosEst] = obtemGDLNosEst (numNos,numGDLno,restrs)

  k1=1;
  for no = 1: numNos 
      for j = 1:numGDLno
          if restrs(no,j)==0 % testa se o gl não está restringido para a linha (no) e a coluna (gdl) 
              GDLNosEst(no,j) = k1;
              k1 = k1 + 1;
          end
      end
  end

  k2=k1;

  for no = 1:numNos
      for j = 1:numGDLno
          if restrs(no,j)==1 % Testa se o gdl está restringido na linha (No) e coluna (gdl) 
              GDLNosEst(no,j) = k2;
              k2 = k2+1;
          end
      end
  end

  % Número de graus de liberdade não restringidos da estrutura 
  numGDLLivre = k1-1;
  return


  %----------------------------------------------------------------------


% Rotina que obtem os gdl dos nós do elemento 


function [GDLNosElem, numGDLElem ] = obtemGDLNosElem (numNosElem,nosElem,numGDLno,GDLNosEst)

  numGDLElem = numNosElem * numGDLno; 
  contador2 = 1;


  for i = 1:numNosElem % Loop pelos nós do elemento
      noElem = nosElem(i); % busca na matriz nosElem, que busca na matriz conect qual os nós daquele elemento 
      for j = 1:numGDLno % Loop pelos Graus de liberdade por nó do elemento 
          GDLNosElem(contador2) = GDLNosEst(noElem,j); % Checar nota de rodapé da função
          contador2 = contador2+1;
      end
  end 

return   

% Nota de rodapé : A matriz GDLNosElem, está obtendo na matriz GDLNosEst, o grau de liberdade
% correspondente ao grau de liberdade atual do nó e também correspondente a aquele nó do elemento em
%referencia ao nó da estrutura, Logo, esta função faz a interpolação entre os nós do elemento e os graus de liberdade
%locais e os nós e gdls da estrutura. 






