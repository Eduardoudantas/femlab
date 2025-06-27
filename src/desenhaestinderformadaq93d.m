% Programa para plot do elemento de placa em 3d  para o programa de análise do método dos elementos finitos 


function desenhaestinderformadaq93d(coords,coordsz,numNos,restrs,conect)

% Unifica a matriz de coordenadas 2d e a matriz de coordenadas z para formar uma unica matriz de coordenadas 3d 
coords3d = [coords coordsz]; 

escala = 0.1; % Define a escala de apresentação dos apoios no plot da malha 
% define igualdade entre os eixos x e y 
grid on 
axis equal 
set(figure(1),'Color',[1 1 1]); 
title('Underformed Structure');



figure(1);
set(gcf,'DefaultLineColor','blue') % Define a cor default das linhas 
axis equal
grid on 


% determina as dimensões maximas e minimas da estrutura 

xmin = min(coords3d(:,1));
xmax = max(coords3d(:,1));

ymin = min(coords3d(:,2));
ymax = max(coords3d(:,2));

zmin = min(coords3d(:,3));
zmax = max(coords3d(:,3));



% Calcula a dimensão dx e dy 

dx = xmax - xmin; 
dy = ymax - ymin; 
dz = zmax - zmin;


lado = max([dx,dz,dy]);
delta = lado/10; 

% Ajusta os limites dos eixos do gráfico para incluir uma margem extra ( delta ) na visualização 

axis ([xmin-delta,xmax+delta,ymin-delta,ymax+delta,zmin-delta,zmax+delta]);
%delta = lado/50;
hold on 
set(gcf,'DefaultLineColor','red') % Define a cor default das linhas 

% Número de elementos 
numElem = size(conect(:,1));


% Cria a região que representa os elementos finitos quadrilaterais 
for i= 1:numElem 
    faces(i,:) =[conect(i,1)  conect(i,5)  conect(i,2)  conect(i,6)  conect(i,3)  conect(i,7)  conect(i,4)  conect(i,8)];
end

patch('faces',faces,'vertices',coords3d,'FaceColor','c', 'EdgeColor', [0.1 0.1 0.1])
set(gcf,'DefaultLineColor','red')

% Desenha os nós do elemento 
hold on 
for no=1:numNos
    xno= coords3d(no,1);
    yno= coords3d(no,2);
    zno =coords3d(no,3);

plot3(xno,yno,zno,'r.','MarkerSize',12)
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
                         set(H,'LineWidth',[1]);


                         uy =  [cy, cy ,cy-d/2,   cy-d/2, cy+d/2, cy+d/2 cy];   %Vertices do triangulo x
                         uz =  [cz, cz-d/4, cz-d/4, cz-d/2-d/4, cz-d/2-d/4, cz-d/4,cz-d/4]; %Vertices do triangulo y
                         ux =  [cx, cx, cx, cx,cx,cx,cx];

                         H=line(ux,uy,uz);
                         set(H,'color',[0,0.7,0]);
                         set(H,'LineWidth',[1]);
    end

end     
 xlabel ('X');
 ylabel ('Y');
 zlabel ('Z');
 set(gcf,'DefaultLineColor','red');
 view(3)
 return 