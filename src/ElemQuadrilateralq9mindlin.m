%Universidade Federal do Para - UFPa
% ITEC - Instituto de Tecnologia 
% Grupo de Dinâmica Instrumentação Processamento de sinais e imagens 
% Curso de Engenharia Civil 
% Programa didatico para Analise Linear pelo Método dos Elementos Finitos 
% Autor:  Eduardo Uchôa Dantas  
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
% Rotinas e funções do elemento Quadrilateral Quadrático Q9

 % [var1,var2] = ElemQuadrilateralq9mindlin(comando,coordsElem,dadosElem,DElem)
 % inputs 

 % comando = diz ao programa o que executar dentro das funções do elemento
 % coordsElem = matriz de coordenadas dos nós do elemento 
 %dadosElem = matriz que fornece os dados do elemento ( E,nu,t, ... ) 
 % Delem = Vetor de deslocamentos nodais do elemento 


 %outputs

 %  'rigidez'  = calcula a matriz de rigidez do elemento (kel)
 %  'cargasequivalentes' = calcula o vetor de cargas nodais equivalentes (peq)
 %  'tensoesdeforms' = calcula as tensões e deformações nos pontos de gauss dentro do elemento 
 %  'massa'' = calcula a matriz de massa do elemento 

 % ----------------------------------------------------------------------------------------------------
  function [var1,var2] = ElemQuadrilateralq9mindlin(comando,coordsElem,dadosElem,DElem)

  switch(comando)
      
      case 'rigidez'  %Calcula a matriz de rigidez do elemento 
          kel = MatRigidezElem(coordsElem,dadosElem);
          var1 = kel;
          var2 = 0;

      case 'cargasequivalentes' % Calcula o vetor de cargas equivalentes 
          peq = vetcargasequivalenteselem(coordsElem,dadosElem);
          var1 = peq;
          var2 = 0;

      case 'tensoesdeforms' % Calcula as tensões e deformações nos pontos de gauss no interior do elemento 
          epslon = deformselem(coordsElem,dadosElem,DElem);
          sigma = tensoeselem(dadosElem,epslon);
          var1 = sigma;
          var2 = epslon;

      case 'massa' % Calcula a matriz de massa do elemento   
          m = MatMassaElem(coordsElem, dadosElem);
          var1 = m;
         var2 = 0;

      otherwise 
          disp('ElemQuadrilateralq9mindlin: comando inválido')
  end
return






 % Funções do programa 

 
%--------------------------------MATRIZ DE RIGIDEZ DO ELEMENTO---------------------------------

% Função que realiza o cálculo da matriz de rigidez do elemento 


function kel = MatRigidezElem(coordsElem, dadosElem)

% captura as informações do elemento 

E = dadosElem(1); %Módulo de elasticidade longitudinal         
nu = dadosElem(2); %Coeficiente de Poisson 
t = dadosElem(3); %Espessura 
ordint = dadosElem(4); % Define a ordem de integração numérica por pontos de gauss
numNosElem = 9; %Número de nós do elemento

numGDLElem = 27; %Número de graus de liberdade do elemento 

% obtem a matriz constitutiva do elemento 
C = matrizconstitutivaq9mindlin(E,t,nu);

% Obtem os pontos de gauss 

[csi, eta, w] = PtsGauss2d(ordint);

%calcula o número de pontos de integração 

numpts = ordint*ordint;


% Monta a matriz de massa simples do elemento q9 
%mass=[ro*t     0         0
       %0  ro*((t^3)/12)   0 
       %0        0    ro*((t^3)/12)];

% inicializa a matriz de rigidez com zeros 

kel = zeros(numGDLElem,numGDLElem);

% Inicializa a matriz de massa do elemento com zeros 

%m = zeros(numGDLElem,numGDLElem);


for i = 1:numpts 

    %obtem a matriz das funções de forma do elemento 

    N = Calcmatrizfuncforma(csi(i),eta(i));

    % obtem as derivadas das funções de forma em relação as coordenadas naturais 

    [dNdcsi, dNdeta] = Derivfuncformnatural(csi(i), eta(i));

    % Obtem a matriz jacobiana 
    
    jac = Matrizjacobiana2d(coordsElem,dNdcsi,dNdeta);

    % Calcula o determinante da matriz jacobiana 
    
    detjac = det(jac);

    % obtem as derivadas das funções de forma em relação as coordenadas naturais 

    [dNdx, dNdy] = Derivfuncformcart(jac,dNdeta,dNdcsi);

    % monta a matriz B 

    B = CalcmatrizB(dNdy,dNdx,N);
    
    % Monta a matriz de massa do elemento 
    
    % Montar nesse formato. 
    %m = m + N'*mass*N*J*w(c);
   
   % m = m + N'*mass*N*detjac*w(i);
    
    % Calcula a matriz de rigidez através da integração numérica 

    kel = kel+ B' * C * B * detjac * w(i);       % w(i) = pesos de gauss 
    
end


return 

%----------------------------------------------------------------------------------------------%

%----------------------------- MATRIZ DE MASSA DO ELEMENTO-------------------------------------%
function m = MatMassaElem(coordsElem, dadosElem)
% captura as informações do elemento 

E = dadosElem(1); %Módulo de elasticidade longitudinal         
nu = dadosElem(2); %Coeficiente de Poisson 
t = dadosElem(3); %Espessura 
ordint = dadosElem(4); % Define a ordem de integração numérica por pontos de gauss
ro = dadosElem(6);
numNosElem = 9; %Número de nós do elemento

numGDLElem = 27; %Número de graus de liberdade do elemento 

% Obtem os pontos de gauss 

[csi, eta, w] = PtsGauss2d(ordint);

%calcula o número de pontos de integração 

numpts = ordint*ordint;


% Monta a matriz de massa simples do elemento q9 
mass=  [ro*t     0         0
        0  ro*((t^3)/12)   0 
         0        0    ro*((t^3)/12)];


% Inicializa a matriz de massa do elemento com zeros 

m = zeros(numGDLElem,numGDLElem);

for i = 1:numpts 

    %obtem a matriz das funções de forma do elemento 

    N = Calcmatrizfuncforma(csi(i),eta(i));

    % obtem as derivadas das funções de forma em relação as coordenadas naturais 

    [dNdcsi, dNdeta] = Derivfuncformnatural(csi(i), eta(i));

    % Obtem a matriz jacobiana 
    
    jac = Matrizjacobiana2d(coordsElem,dNdcsi,dNdeta);

    % Calcula o determinante da matriz jacobiana 
    
    detjac = det(jac);

    % obtem as derivadas das funções de forma em relação as coordenadas naturais 

    [dNdx, dNdy] = Derivfuncformcart(jac,dNdeta,dNdcsi);

    % monta a matriz B 

    B = CalcmatrizB(dNdy,dNdx,N);
    
    % Monta a matriz de massa do elemento
   
    m = m + N'*mass*N*detjac*w(i);
    

end
    
     


return



%----------------------------------------------------------------------------------------------%




%----------------------------VETOR DE CARGAS NODAIS EQUIVALENTES DO ELEMENTO -------------------

function  peq = vetcargasequivalenteselem(coordsElem, dadosElem)
E = dadosElem(1); % modulo de deformação longitudinal 
nu = dadosElem(2); % Coef de Poisson 
t = dadosElem(3); %Espessura da chapa 
ordint = dadosElem(4); %Ordem de integração numérica
b = zeros(3,1);
b(1,1)= dadosElem(5); % Carga distribuida em z por área kN/m2

numNosElem = 9;
numGDLElem = 27;

% Obtem a matriz constitutiva 
C = matrizconstitutivaq9mindlin(E,t,nu);

% Obtem os pontos de gauss para integração numérica 

[csi, eta, w] = PtsGauss2d(ordint);

%Calcula o numero de pontos de integração 

numpts = ordint*ordint;

%Inicializa o vetor de cargas nodais equivalentes do elemento (peq) 

peq = zeros(numGDLElem,1);


% Inicia o loop por ponto de gauss 

for i = 1:numpts

    % Obtem as funções de forma do elemento 
    N = Calcmatrizfuncforma(csi(i),eta(i));

    %Obtem as derivadas das funções de forma em coordenadas naturais 
    [dNdcsi, dNdeta] = Derivfuncformnatural(csi(i), eta(i));

    % obtem a matriz jacobiana 
    jac = Matrizjacobiana2d(coordsElem,dNdcsi,dNdeta);

    % calcula o determinante da matriz jacobiana 

    detjac = det(jac);

    %Realiza o somatório para calcular o vetor de cargas nodais equivalentes

    peq = peq + N' *b*detjac*w(i); 

end


return


%-----------------------------------------------------------------------------------------------%


%---------------------------CALCULA AS DEFORMAÇÕES NO ELEMENTO-------------------------------------% 

function epslon = deformselem(coordsElem,dadosElem,DElem)


ordInt = dadosElem(4);  % ordem de integração numérica 
numNosElem = 9;   % Numero de nós do elemento 

%Obtem os pontos de gauss para a ordem de integração definida 

[csi, eta, w] = PtsGauss2d(ordInt);

numPts = 2*ordInt;

for i = 1:numPts 

   
    % Obtem a matriz de funções de forma do elemento 

    N = Calcmatrizfuncforma(csi(i),eta(i));
    
    %obtem as derivadas das funções de forma no sistema natural. 

    [dNdcsi, dNdeta] = Derivfuncformnatural(csi(i), eta(i));

    %Obtem a matriz jacobiana de transformação 

    jac = Matrizjacobiana2d(coordsElem,dNdcsi,dNdeta);

    % Obtem as derivadas das funções de forma em relação as coordenadas cartesianas 

    [dNdx, dNdy] = Derivfuncformcart(jac,dNdeta,dNdcsi);

    % Obtem a matriz B que relaciona deslocamentos com deformações nos nós do elemento 

    B = CalcmatrizB(dNdy,dNdx,N);


    % Calcula as deformações para cada ponto de gauss e armazena no vetor epslon 

    
    epslon(:,i) = B*DElem;

end 

return





%----------------------------------------------------------------------------------------------%

%---------------------------CALCULA AS TENSÕES NO ELEMENTO-------------------------------------% 

function sigma = tensoeselem(dadosElem,epslon) 

  E = dadosElem(1);  % Modulo de deformação longitudinal
  nu = dadosElem(2); % Coef de poisson
  t = dadosElem(3);  % Espessura do elemento 
  
  % Obtem a matriz constitutiva do elemento 
  C = matrizconstitutivaq9mindlin(E,t,nu);

  % Calcula os esforços no elemento 

  sigma  = C*epslon;


return 


%----------------------------------------------------------------------------------------------%
 %---------------------------- MATRIZ DE FUNÇÕES DE FORMA---------------------------------------%   
function N = Calcmatrizfuncforma(csi,eta)
        
        N1 = csi*eta*(csi-1.0)*(eta-1)/4.0;
        N2 = csi*eta*(csi+1.0)*(eta-1)/4.0;
        N3 = csi*eta*(csi+1.0)*(eta+1)/4.0;
        N4 = csi*eta*(csi-1.0)*(eta+1)/4.0;
        N5 = eta*(1.0-csi^2)*(eta-1)/2.0;
        N6 = csi*(csi+1.0)*(1-eta^2)/2.0;
        N7 = eta*(1.0-csi^2)*(eta+1)/2.0;
        N8 = csi*(csi-1.0)*(1-eta^2)/2.0;
        N9 = (1.0-csi^2)*(1-eta^2);

  
  

 
   N = [N1 0 0 N2 0 0 N3 0 0 N4 0 0 N5 0 0 N6 0 0 N7 0 0 N8 0 0 N9 0 0;
        0 N1 0 0 N2 0 0 N3 0 0 N4 0 0 N5 0 0 N6 0 0 N7 0 0 N8 0 0 N9 0;
        0 0 N1 0 0 N2 0 0 N3 0 0 N4 0 0 N5 0 0 N6 0 0 N7 0 0 N8 0 0 N9]; 
 return

 %--------------------------------------------------------------------------------------------%

 % Define a matriz Constitutiva para o elemento pela análise de mindlin 
function C = matrizconstitutivaq9mindlin(E,t,nu)
     
     Dc=(E*t^3)/(12*(1-nu^2));
     Kappa=5/6;
     G=E/(2*(1+nu));
     C = [  Dc    Dc*nu     0        0    0; 
           Dc*nu   Dc       0        0    0;
            0     0   Dc*(1-nu)/2   0    0;
            0     0       0      Kappa*t*G  0;
            0     0       0        0   Kappa*t*G];

return

%----------------------------------------------------------------------------------------------%

 %Função que calcula as derivadas das funções de forma no sistema natural 

 % implementar esta rotina com base no codigo rev final 

 function  [dNdcsi, dNdeta] = Derivfuncformnatural(csi, eta)


  %Contém as derivadas das funções de forma em relação a csi e eta em forma
  %de duas matrizes 1x8 ( 1 linha e 9 colunas ) 

    dNdcsi9=(-2*csi)*(1-eta^2);
    dNdcsi8=-0.5*(1-eta^2)-(1/2*dNdcsi9);
	dNdcsi7=0.5*(-2*csi)*(1+eta)-(1/2*dNdcsi9);
	dNdcsi6=0.5*(1-eta^2)-(1/2*dNdcsi9);
	dNdcsi5=0.5*(-2*csi)*(1-eta)-(1/2*dNdcsi9);
	dNdcsi4=-0.25*(1+eta)-(1/2)*(dNdcsi7+dNdcsi8)-(1/4*dNdcsi9);
    dNdcsi3=0.25*(1+eta)-(1/2)*(dNdcsi6+dNdcsi7)-(1/4*dNdcsi9);
	dNdcsi2=0.25*(1-eta)-(1/2)*(dNdcsi5+dNdcsi6)-(1/4*dNdcsi9);
    dNdcsi1=-0.25*(1-eta)-(1/2)*(dNdcsi5+dNdcsi8)-(1/4*dNdcsi9);
   
   dNdeta9= (1-csi^2)*(-2*eta);
   dNdeta8= 0.5*(1-csi)*(-2*eta)-(1/2*dNdeta9);
   dNdeta7= 0.5*(1-csi^2)-(1/2*dNdeta9);
   dNdeta6= 0.5*(1+csi)*(-2*eta)-(1/2*dNdeta9);
   dNdeta5=-0.5*(1-csi^2)-(1/2*dNdeta9);
   dNdeta4=0.25*(1-csi)-(1/2)*(dNdeta7+dNdeta8)-(1/4*dNdeta9);
   dNdeta3=0.25*(1+csi)-(1/2)*(dNdeta6+dNdeta7)-(1/4*dNdeta9);
   dNdeta2=-0.25*(1+csi)-(1/2)*(dNdeta5+dNdeta6)-(1/4*dNdeta9);
   dNdeta1=-0.25*(1-csi)-(1/2)*(dNdeta5+dNdeta8)-(1/4*dNdeta9);
   
   dNdcsi=[dNdcsi1 dNdcsi2 dNdcsi3 dNdcsi4 dNdcsi5 dNdcsi6 dNdcsi7 dNdcsi8 dNdcsi9];
   dNdeta=[dNdeta1 dNdeta2 dNdeta3 dNdeta4 dNdeta5 dNdeta6 dNdeta7 dNdeta8 dNdeta9];


   % Todas as funções de forma interpolam o nó 9 pois este é o nó central do elemento Q9

  return 

  %--------------------------------------------------------------------------------------%


  % Rotina que monta a matriz jacobiana de transformação 

 
%function jac = Matrizjacobiana2d(coordsElem,dNdcsi,dNdeta)
 

 
  % monta a parte base da matriz jacobiana 

 % jacbase=[dNdcsi; % Vetor linha de derivadas das funções de forma em csi
    %       dNdeta]; % vetor linha de derivadas das funções de forma em eta 




   %monta a parte variável da matriz jacobiana 
  % jacvariavel = coordsElem;
  
  % Realiza a operação pra montagem da matriz jacobiana do elemento 
 % jac = jacbase*jacvariavel;

  %return 


  %--------------------------------------------------------------------------------------%

  % função que calcula as derivadas das funções de forma no sistema cartesiano através da matriz jacobiana 

function  [dNdx, dNdy] = Derivfuncformcart(jac,dNdeta,dNdcsi)


 % na formulação isoparamétrica, sabe-se primeiramente as derivadas em função das coordenadas naturais(dNdeta e dNdcsi) e, a partir da matriz jacobiana,
 % transforma-se essas derivadas em derivadas cartesianas (dNdx e dNdy) 

 %inicializa a matriz 2x9 que conterá as derivadas cartesianas das funções de forma 



 %monta a matriz 2x9 de derivadas das funções de forma em coordenadas naturais 

 CoordsnatN= [dNdcsi;
             dNdeta];

 %Inverte a matriz jacobiana 

 jacinv=inv(jac);

 CoordscartN = jacinv*CoordsnatN;



 dNdx = CoordscartN(1,:);
 dNdy = CoordscartN(2,:);
return 


%-----------------------------------------------------------------------------------------------------------%


% função que realiza a montagem da matriz B ( deslocamento-deformação) para o elemento de placa de mindlin 

function B = CalcmatrizB(dNdy,dNdx,N)
 


  % realiza o loop que fará a montagem da matriz B a partir das relações dos operadores diferenciais. 
  for i=1:9
      g=[0    dNdx(1,i)    0  
         0     0     dNdy(1,i)
         0  dNdy(i)  dNdx(1,i) 
         -dNdx(1,i)  N(1,i*3-2)     0  
         -dNdy(1,i)  0     N(1,i*3-2)  ];
      B(:,i*3-2:i*3)=g;
  end 

return

%----------------------------------------------------------------------------------------------------------%



