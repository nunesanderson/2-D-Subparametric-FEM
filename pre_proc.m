%==========================================================================
%ROTINA QUE GERA OS DADOS PARA O SOLVER
%==========================================================================
%Monografia para conclusão do curso de pós graduação
%Aluno: Anderson Santos Nunes
%Orientador: Prof. Dr. Marcelo Grafulha Vanti
%Data: 22/08/2011
%==========================================================================

Leitor_msh;

%==========================================================================

disp('===================================================================')
disp('Arquivo pre_proc')

%==========================================================================
%Obtém a quantidade de elementos e e nós e
%gera as lista de coordenadas X e Y

n=length(num_global);
numElem=n(1,1);         % Quantidade de elementos
n=length(coordenadas);
numNos=n(1,1);          %Quantidade de nós

%==========================================================================
%Organiza a numeraçao global
%         x(5)
%       /   \
%      /     \
%     x(6)    x(4)
%    /         \
%   /           \
% x(1)---x(2)----x(3)

i=0;
    
%Verifica o sentido da numeração local
i=0;
for i=1:numElem

    %   Calcula o baricentro do elemento  
    xb=(coordenadas(num_global(i,7),1)+coordenadas(num_global(i,8),1)+coordenadas(num_global(i,9),1))/3;
    yb=(coordenadas(num_global(i,7),2)+coordenadas(num_global(i,8),2)+coordenadas(num_global(i,9),2))/3;
    
    for k=1:6
        xp(k,1)=coordenadas(num_global(i,k+6),1);
        yp(k,1)=coordenadas(num_global(i,k+6),2);
    end
conta=0;
    for k=1:6
        
        if yp(k,1)==max(yp(:,1))
            conta=conta+1;
            xmax_(k,1)=xp(k,1);
        end
    end
    
    for k=1:6
        if yp(k,1)==max(yp(:,1))
            if conta==1
                ymax=yp(k,1);
                xmax=xp(k,1);
            else
                ymax=yp(k,1);
                xmax=max(xmax_(:,1));
            end
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
    
    
disp('Organização da numeração local: ok');
%==========================================================================
%Gera a tabela com com condições de contorno
%condContorno(i,1): se há cond. de contorno neste ponto
%condContorno(i,2): valor

n=length(cond_contorno);
linhasCondCont=n(1,1);%  Quantidade de linhas do arquivo com condiçoes de contorno

%Verifica quais sao as condiçoes de contorno
quaisCondCont(1,1)=cond_contorno(1,4);
cond_conta=1;
j=0;
for i=1:linhasCondCont
    nova=1;
    for j=1:cond_conta
       if cond_contorno(i,4)==quaisCondCont(1,j)
          nova=0;
        end
    end
    
    if nova==1
        cond_conta=cond_conta+1;;
        quaisCondCont(1,cond_conta)=cond_contorno(i,4);

    end
end

%Gera uma tabela como todos os nós com cond. de contorno e a sua 
% physical line correspondente
%Primeira linha
contCondCont=0;
i=1;
for j=7:9
    pontos(i,1)=cond_contorno(1,j);
    pontos(i,2)=cond_contorno(1,4);
    i=i+1;
    contCondCont=contCondCont+1;
end

%Restante
for i=2:linhasCondCont
    
    for j=7:9
        nova=1;
        
        for k=1:contCondCont
            
            if cond_contorno(i,j)==pontos(k,1)
                nova=0;
            end
        end

        if nova==1
            contCondCont=contCondCont+1;
            pontos(contCondCont,1)=cond_contorno(i,j);
            pontos(contCondCont,2)=cond_contorno(i,4);
        end
    end
end

%Entrada dos valores de contorno
valCondCont = input(['Condiçoes de contorno ' mat2str(quaisCondCont) ':']);

%Gera a matriz de condições de contorno
n=size(quaisCondCont);
condContorno=zeros(numNos,2);

n1=size(pontos);
n2=size(quaisCondCont);

for i=1:n1(1,1);
    
    for j=1:n2(1,2);
        if pontos(i,2)==quaisCondCont(1,j)
            val_=valCondCont(1,j);
        end
    end
    
    condContorno(pontos(i,1),1)=1;
    condContorno(pontos(i,1),2)=val_;
end


disp('Condições de contorno: ok');
%==========================================================================
%Propriedades 
%propriedades(i,1):índice do material (linha na tebela de materiais)
%propriedades(i,2): valor de excitação (Densidade de corrente ou carga elétrica)

%Verifica quais sao as regiões
quaisReg(1,1)=num_Global(1,7);
contReg=1;
j=0;
for i=1:numElem
    nova=1;
    for j=1:contReg
       if num_Global(i,7)==quaisReg(1,j)
          nova=0;
        end
    end
    
    if nova==1
        contReg=contReg+1;;
        quaisReg(1,contReg)=num_Global(i,7);
    end
end

%Entrada dos índices dos materiais de cada região
materiais
valMatReg = input(['Índice do material ' mat2str(quaisReg) ':']);

%Entrada dos valores de excitação
disp('Excitação: Se transiente J(t)=Jp*raiz(2)*sen(w*t)')
exctReg = input(['Excitação J[A/m2] ou Q [Q/m2] ' mat2str(quaisReg) ':']);

%Gera a matriz de condições de contorno
for i=1:numElem
    for j=1:contReg
        if quaisReg(1,j)==num_Global(i,7)
            propriedades(i,1)=valMatReg(j);
            propriedades(i,2)=exctReg(j);
            break;
        end
    end
end

disp('Propriedades: ok');

%==========================================================================
%Limpeza de variáveis não utilizadas
clear ang;
clear ang_in;
clear ang_sort;
clear conta;
clear xp;
clear yp;
clear xmax;
clear ymax;
clear xmax_;
clear xb;
clear yb;
clear cond_conta;
clear cond_contorno;
clear contCondCont;
clear contReg;
clear exctReg;
clear i;
clear j;
clear k;
clear linhasCondCont;
clear n;
clear nova;
clear num_global;
clear pontos;
clear quaisCondCont;
clear quaisReg;
% clear valCondCont;
clear valMatReg;
disp('Limpeza de variáveis: ok');