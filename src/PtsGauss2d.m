% Universidade Federal do Para - UFPa
% Centro Tecnologico - CT
% Departamento de Construcao Civil - DCC
% Curso de Engenharia Civil
% Disciplina: Introducao ao Metodo dos Elementos Finitos
% Programa didatico para Analise Linear por Elementos Finitos
% Autor: Remo Magalhaes de Souza       (remo@ufpa.br)
% Versão 1.0       Data: 09/2001
% arquivo para obtencao dos pontos de integracao de Gauss em duas dimensoes
%---------------------------------------------------------------------------

function  [csi, eta, w] = PtsGauss2d(ordem)
 
  [csi1, w1] = PtsGauss1d(ordem);
  k = 0;

  for i=1: ordem
    for j=1: ordem
      k = k + 1;
      csi(k) = csi1(j);
      eta(k) = csi1(i);
      w(k)   = w1(i) * w1(j);
    end
  end

return


