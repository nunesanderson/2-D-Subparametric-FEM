
hold on
for i=1:1

    %   Calcula o baricentro do elemento  
    xb=(coordenadas(num_global(i,7),1)+coordenadas(num_global(i,8),1)+coordenadas(num_global(i,9),1))/3;
    yb=(coordenadas(num_global(i,7),2)+coordenadas(num_global(i,8),2)+coordenadas(num_global(i,9),2))/3;
    text(xb,yb,num2str((i)),'color','blue');
    
    for k=1:6
        xp(k,1)=coordenadas(num_global(i,k+6),1);
        yp(k,1)=coordenadas(num_global(i,k+6),2);
        
         text(xp(k,1),yp(k,1),num2str(num_global(i,k+6)),'color','red');
    end

    for k=1:6
        if yp(k,1)==max(yp(:,1))
            ymax=yp(k,1);
            xmax=xp(k,1);
        end
    end
    ang_in=atan2(ymax-yb,xmax-xb)*180/pi;
%       Tendo o baricentro, calcula o ãngulo entre o baricentro e o 
% ponto em questão
    for k=1:6
        ang(k,1)=atan2((yp(k,1)-yb),(xp(k,1)-xb))*180/pi-ang_in;
        if ang(k,1)<0
            ang(k,1)=ang(k,1)+360;
        end
    end
  
%   Organiza a tabela de angulos em forma crescente  
      ang_sort=sort(ang,'ascend');

%Cria a nova posição da numerção, de acordo com os angulos, em forma
% crescente

    for l=1:6
        for m=1:6
            if ang(m,1)==ang_sort(l,1)
                num_Global(i,l)=num_global(i,m+6);
            end
        end
    end

    num_Global(i,7)=num_global(i,4);        %Insere o índice da physical surface

end