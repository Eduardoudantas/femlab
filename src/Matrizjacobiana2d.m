function jac = Matrizjacobiana2d(coordsElem,dNdcsi,dNdeta)
 

 
  % monta a parte base da matriz jacobiana 

  jacbase=[dNdcsi; % Vetor linha de derivadas das funções de forma em csi
           dNdeta]; % vetor linha de derivadas das funções de forma em eta 





   %monta a parte variável da matriz jacobiana 
   jacvariavel = coordsElem;
  
   
   
  % Realiza a operação pra montagem da matriz jacobiana do elemento 
  jac = jacbase*jacvariavel;

  return 