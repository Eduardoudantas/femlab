% Universidade Federal do Para - UFPa
% Centro Tecnologico - CT
% Departamento de Construcao Civil - DCC
% Curso de Engenharia Civil
% Disciplina: Introducao ao Metodo dos Elementos Finitos
% Programa didatico para Analise Linear por Elementos Finitos
% Autor: Remo Magalhaes de Souza       (remo@ufpa.br)
% Versão 1.0       Data: 09/2001
% arquivo para obtencao dos pontos de integracao de Gauss em uma dimensao
%---------------------------------------------------------------------------

function [csi, w] = PtsGauss1d(ordem)

  csi = zeros(ordem,1);
  w = zeros(ordem,1);

  switch ordem
    case 1
      csi(1) = 0;
      w(1) = 2;
   
    case 2
      csi(1) = -1 / sqrt(3);
      csi(2) = -csi(1);
   
      w(1) =  1;
      w(2) =  1;
   
    case 3
      csi(1) = -0.774596669241483;
      csi(2) =  0;
      csi(3) = -csi(1);
   
      w(1) = 5/9;
      w(2) = 8/9;
      w(3) = w(1);
   
    case 4
      csi(1) = -0.861136311594053;
      csi(2) = -0.339981043584856;
      csi(3) =  0.339981043584856;
      csi(4) =  0.861136311594053;
      w(1) =  0.347854845137454;
      w(2) =  0.652145154862546;
      w(3) =  0.652145154862546;
      w(4) =  0.347854845137454;
   
  end

return

