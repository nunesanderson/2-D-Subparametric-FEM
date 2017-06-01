%==========================================================================
%ROTINA QUE PLOTA A MALHA E AS LINHAS EQUIPOTENCIAIS
%==========================================================================
%Monografia para conclusão do curso de pós graduação
%Aluno: Anderson Santos Nunes
%Orientador: Prof. Dr. Marcelo Grafulha Vanti
%Data: 22/08/2011
%==========================================================================

close all;

%==========================================================================
%Entrada de dados

disp('===================================================================')
disp('Arquivo pos_proc')

if problema==4
    instante = input(['Coluna da matriz de potenciais: ']);
else
    instante=1;
end

%==========================================================================
%Desenha a malha e e indica pontos com cond. de contorno

max_=max(num_Global(:,7));
min_=min(num_Global(:,7));
delta=max_-min_;

for i=1:numElem

    faces(i,1)=num_Global(i,1);
    faces(i,2)=num_Global(i,3);
    faces(i,3)=num_Global(i,5);

    tcolor(i,1)=2*(num_Global(i,7)-min_)/delta;
    tcolor(i,2)=(num_Global(i,7)-min_)/delta;
    tcolor(i,3)=.2*(num_Global(i,7)-min_)/delta;
    
    if propriedades(i,1)==1 && propriedades(i,2)==0
        tcolor(i,1)=1;
        tcolor(i,2)=1;
        tcolor(i,3)=1;
    end
end
    
figure;
patch('Faces',faces,...
      'Vertices',coordenadas,...
      'FaceVertexCData',tcolor,...
      'FaceLighting','flat',...
      'FaceColor','flat',...
      'FaceAlpha',.5,...
      'EdgeAlpha',.02);
hold on

% Insere os valores das condiçoes de coontorno em cada nó
%  for i=1:numNos
%  
%      if condContorno(i,1)==1
%         text(coordenadas(i,1),coordenadas(i,2),num2str(condContorno(i,2))); 
%      end
%  end
 
hold on
% 
% for i=1:numElem
%     x=(coordenadas(num_Global(i,1),1)+coordenadas(num_Global(i,3),1)+coordenadas(num_Global(i,5),1))/3;
%     y=(coordenadas(num_Global(i,1),2)+coordenadas(num_Global(i,3),2)+coordenadas(num_Global(i,5),2))/3;
%     text(x,y,num2str(propriedades(i,2)),'color','blue');
% end
% hold on
% for i=1:numNos
%     text(coordenadas(i,1),coordenadas(i,2),num2str((i)),'color','red')
% end

disp('Malha plotada: ok');

%==========================================================================
%Desenha linhas equipotenciais
pdeplot(coordenadas',[],num_Global(:,1:2:6)','xydata',potenciais(:,instante),'xystyle','off','contour',...
'on','levels',19,'colormap',jet,'colorbar','on')
colorbar('vert')