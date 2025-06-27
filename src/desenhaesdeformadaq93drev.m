% Programa para plot do elemento de placa em 3d  para o programa de análise do método dos elementos finitos 


function desenhaesdeformadaq93drev(coords,coordsz,numNos,restrs,conect,Deslocamentos)

% Unifica a matriz de coordenadas 2d e a matriz de coordenadas z para formar uma unica matriz de coordenadas 3d 
coords3d = [coords coordsz]; 

escaladesloc = 1;
escala = 0.4; % Define a escala de apresentação dos apoios no plot da malha 
% define igualdade entre os eixos x e y 
%grid on 
axis equal 
set(figure(3),'Color',[1 1 1]); 
title('Deformed Structure');

for i = 1:numNos
  coordsdef(i,1) = coords3d(i,1);
  coordsdef(i,2) = coords3d(i,2);
  coordsdef(i,3) = coords3d(i,3) + Deslocamentos(i,1)*escaladesloc;
end 



% determina as dimensões maximas e minimas da estrutura 

xmin = min(coordsdef(:,1));
xmax = max(coordsdef(:,1));

ymin = min(coordsdef(:,2));
ymax = max(coordsdef(:,2));

zmin = min(coordsdef(:,3));
zmax = max(coordsdef(:,3));



% Calcula a dimensão dx e dy 

dx = xmax - xmin; 
dy = ymax - ymin; 
dz = zmax - zmin;


lado = max([dx,dz,dy]);
delta = lado/5; 

% Ajusta os limites dos eixos do gráfico para incluir uma margem extra ( delta ) na visualização 

axis ([xmin-delta,xmax+delta,ymin-delta,ymax+delta,zmin-delta,zmax+delta]);
hold on 
set(gcf,'DefaultLineColor','red') % Define a cor default das linhas 

% Número de elementos 
numElem = size(conect(:,1));
cont=1;


% Cria a região que representa os elementos finitos quadrilaterais 
for i= 1:numElem 
    faces(cont,:) =[conect(i,1)  conect(i,2)  conect(i,3)  conect(i,4)];
    cont = cont+1;
end

patch('faces',faces,'vertices',coords3d,'FaceColor','none','Linestyle',':','EdgeColor', [0.1 0.1 0.1])
patch('faces',faces,'vertices',coordsdef,'FaceColor','none', 'EdgeColor', [1.0 0.1 0.1])


% Desenha os nós do elemento 
hold on 
for no=1: numNos
  
 %Desenha os nós na configuração indeformada
%   xNo = coords(no,1);
%   yNo = coords(no,2);
%   zNo = coords(no,3);
  %plot3(xNo,yNo,zNo, 'k.','Markersize',12);

%Desenha Apoio.

    if restrs(no,1)== 1 || restrs(no,2)== 1 || restrs(no,3)== 1
        
      
                         d = delta*escala;
                         cx = coords3d(no,1);
                         
                         cy = coords3d(no,2);
                         cz = coords3d(no,3);


                         ux =  [cx, cx ,cx-d/2,   cx-d/2, cx+d/2, cx+d/2 cx];   %Vertices do triangulo x
                         uz =  [cz, cz-d/4, cz-d/4, cz-d/2-d/4, cz-d/2-d/4, cz-d/4,cz-d/4]; %Vertices do triangulo y
                         uy =  [cy, cy, cy, cy,cy,cy,cy];

                         H=line(ux,uy,uz);
                         set(H,'color',[0,0.7,0]);
                         set(H,'LineWidth',1);


                         uy =  [cy, cy ,cy-d/2,   cy-d/2, cy+d/2, cy+d/2 cy];   %Vertices do triangulo x
                         uz =  [cz, cz-d/4, cz-d/4, cz-d/2-d/4, cz-d/2-d/4, cz-d/4,cz-d/4]; %Vertices do triangulo y
                         ux =  [cx, cx, cx, cx,cx,cx,cx];

                         H=line(ux,uy,uz);
                         set(H,'color',[0,0.7,0]);
                         set(H,'LineWidth',1);
    end

end
     
 xlabel ('X');
 ylabel ('Y');
 zlabel ('Z');
 set(gcf,'DefaultLineColor','red');
 view(3)
 return 