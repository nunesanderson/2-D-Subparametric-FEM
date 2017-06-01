%==========================================================================
%SOLVER TRANSIENTE DE ELEMENTOS FINTOS SUBPARAMÉTRICOS
%==========================================================================
%Monografia para conclusão do curso de pós graduação
%Aluno: Anderson Santos Nunes
%Orientador: Prof. Dr. Marcelo Grafulha Vanti
%Data: 22/08/2011
%==========================================================================

clc;
clear all;
close all;

%==========================================================================

pre_proc

%==========================================================================

disp('===================================================================')
disp('Arquivo FEMSolver_Static_2Ord')

%==========================================================================
% Carrega os arquivos de dados
% load -AscII 'coordenadas.txt'   %Coordenadas XY de cada nó (x,y)
% load -AscII 'num_Global.txt'      %Numeraçao global de cada elemento
% load -AscII 'propriedades.txt'  %Índice do material e excitaçao de cada elemento (ind. material,densi. corr.)
% load -AscII 'condContorno.txt'  %Condiçoes de contorno, em cada nós (tem? (0 ou 1),valor)

tmax = input(['Tempo de simulaçao [s]:']);
dt = input(['Intervalo de tempo [s]:']);
freq = input(['Frequencia [Hz]:']);
problema=4;
simetria=0;             %Configura a simetria planar

%==========================================================================
%Configuraçoes iniciais
%Definiçao dos pontos de integraçao numérica / triangulo de referencia
% u
% |
% |\
% | \
% |  \
% |   \
% |    \
% |(3)   \
% | |  \  \ 
% | |   \  \
% | |    \  \  
% |(1)---(2) \ 
% ------(3)-----v

%Integração numérica - 3 pontos
qtdNumInteg3=3;              %Quantidade de pontos de integração
u3=[2/3;1/6;1/6];
v3=[1/6;2/3;1/6];
w3=1/6;                      %Ponderaçao

%Integração numérica - 6 pontos
qtdNumInteg6=6;       
a=0.445948490915965;
b=0.091576213509771;
u6=[a;1-2*a;a;b;1-2*b;b];
v6=[a;a;1-2*a;b;b;1-2*b;];
w1=0.111690794839005;
w2=0.054975871827661;
w6=[w1;w1;w1;w2;w2;w2];


Integdudv=1;              %Área do elemento de referencia

gradN_prim=[-1 1 0;                 %gradN(u,v) de primeira ordem
         -1 0 1];                  



%Definiçao dos materiais
n=size(listaMat);
mu0=4*pi*1e-7;
for i=1:n(1,1)
    inv_perm(i,1)=1/(mu0*listaMat(i,3));
    cond(i,1)=listaMat(i,2);
    
end
disp('Pré loop principal: ok');

%==========================================================================
%Loop principal 

contador_tempo=1;
potenciais_=zeros(numNos,1);     %Vetor contendo os resultados
%Inicia o loop no tempo
for tempo=0:dt:tmax
    disp(['Instante: ' num2str(tempo) 'seg' ])
    
%Inicializaçao com zeros das matriz principais
    MatGlobal_esq=zeros(numNos,numNos); %Matriz global do lado esquerdo
    MatGlobal_dir=zeros(numNos,1);      %Matriz global do lado direito

%Inicia o loop para obtençao da matriz global
    for k=1:numElem
        
%Definição do valor de densidade de corrente para o elemento, em um dado
%instante de tempo
            Js=propriedades(k,2)*sin(2*pi*freq*tempo);
            
%Definiçao da propriedade do material
        Inv_perm=inv_perm(propriedades(k,1),1);
        Cond=cond(propriedades(k,1),1);
        
%Obtém os valores dos potencias nos nós deste elemento do intervalo
%de tempo anterior ao que está sendo calculado
        pot_t=zeros(1,6);
        for i=1:6
            pot_t(1,i)=potenciais_(num_Global(k,i),1);
        end


        %Definiçao da matriz Jacobiana
        coordJ=[coordenadas(num_Global(k,1),1) coordenadas(num_Global(k,1),2);
                coordenadas(num_Global(k,3),1) coordenadas(num_Global(k,3),2);
                coordenadas(num_Global(k,5),1) coordenadas(num_Global(k,5),2)];


        Jac=gradN_prim*coordJ;      %Matriz Jacobiana
        invJac=inv(Jac);            %Matriz inversa da Jacobiana
        detJac=det(Jac);

        %Integral do lado esquerdo
matLocal_esq=zeros(6,6);

%Primeiro termo - Integral (Transposta(InvJacobiana*gradN)*(InvJacobiana*gradN)*det(Jacobiana)*dudv}
%*prop.material)
        
        for pinteg=1:qtdNumInteg3
            uInteg=u3(pinteg);
            vInteg=v3(pinteg);
            tInteg=1-uInteg-vInteg;

            %gradN de secunda ordem
            gradN_sec=[(1-4*tInteg) (4*(tInteg-uInteg)) (4*uInteg-1) (4*vInteg) (0) (-4*vInteg);
                        (1-4*tInteg) (-4*uInteg) (0) (4*uInteg) (4*vInteg-1) (4*(tInteg-vInteg))]; 
             InvJacGradN=Jac\gradN_sec;
                    
            matLocal_esq=matLocal_esq+Inv_perm*(InvJacGradN)'*(InvJacGradN)*detJac*Integdudv*w3;
        end

%Segundo termo - Sigma/dt*Integral {Transposta[N]*det(Jacobiana)*dudv}  
        
        for pinteg=1:qtdNumInteg6
            uInteg=u6(pinteg);
            vInteg=v6(pinteg);
            wInteg=w6(pinteg);
            
            tInteg=1-uInteg-vInteg;
             N_sec=[(-tInteg*(1-2*tInteg)) (4*uInteg*tInteg) (-uInteg*(1-2*uInteg)) (4*uInteg*vInteg) (-vInteg*(1-2*vInteg)) (4*vInteg*tInteg)];
            matLocal_esq=(Cond/dt)*(N_sec'*N_sec)*detJac*Integdudv*wInteg+matLocal_esq;
        end
        
%Montagem da matriz global do lado esquerdo
        for im=1:6
            for jm=1:6
                MatGlobal_esq(num_Global(k,im),num_Global(k,jm))=MatGlobal_esq(num_Global(k,im),num_Global(k,jm))+matLocal_esq(im,jm);
            end
        end

%Integral do lado direito
        matLocal_dir=zeros(6,1);
% Terceiro termo -Sigma/dt*Integral(Transposta[N]*[A](t)*det[J]*dudv
        N_sec=zeros(1,6);
        for pinteg=1:qtdNumInteg6
            uInteg=u6(pinteg);
            vInteg=v6(pinteg);
            wInteg=w6(pinteg);
            
            tInteg=1-uInteg-vInteg;
             N_sec=[(-tInteg*(1-2*tInteg)) (4*uInteg*tInteg) (-uInteg*(1-2*uInteg)) (4*uInteg*vInteg) (-vInteg*(1-2*vInteg)) (4*vInteg*tInteg)];

            matLocal_dir=matLocal_dir+(Cond/dt)*N_sec'*N_sec*pot_t'*detJac*Integdudv*wInteg;
        end

%Quarto termo - Integral {Transposta[N]*Js*Det[J]*dudv}
         
         for pinteg=1:qtdNumInteg3
            uInteg=u3(pinteg);
            vInteg=v3(pinteg);
            tInteg=1-uInteg-vInteg;
            N_sec=[(-tInteg*(1-2*tInteg)) (4*uInteg*tInteg) (-uInteg*(1-2*uInteg)) (4*uInteg*vInteg) (-vInteg*(1-2*vInteg)) (4*vInteg*tInteg)];
            matLocal_dir=matLocal_dir+N_sec'*detJac*Integdudv*w3*Js;
         end
     
         
%Montagem da matriz global do lado direito
         for im=1:6
             MatGlobal_dir(num_Global(k,im),1)=MatGlobal_dir(num_Global(k,im),1)+matLocal_dir(im,1);
         end
    end

%==========================================================================
    %Aplicaçao das condiçoes de contorno
    for i=1:numNos

        if condContorno(i,1)==1
            MatGlobal_esq(i,1:numNos)=0;    %Zera a linha do nó na matriz global direita
            MatGlobal_esq(i,i)=1;           %Aplica 1 na pos. do nó na matriz global direita
            MatGlobal_dir(i,1)= condContorno(i,2);  %Aplica a cond. de contorno na pos. do nó na matriz global esquerda
        end
    end

%==========================================================================
%Resoluçao do sistema matricial através do método de Eleminaçao de Gauss
    potenciais_=zeros(numNos,1);     %Vetor contendo os resultados
    potenciais_=MatGlobal_esq\MatGlobal_dir;

%==========================================================================
%Salva os valores de A em "potenciais salvos" para os dados instantes
        potenciais(:,contador_tempo)=potenciais_;
        contador_tempo=contador_tempo+1;

end
disp('Loop no tempo: ok');
%==========================================================================
%Limpeza de variáveis não utilizadas
clear Integdudv;
clear Jac;
clear Js;
clear MatGlobal_dir;
clear MatGlobal_esq;
clear N_sec;
clear coordJ;
clear detJac;
clear gradN_prim;
clear gradN_sec;
clear i;
clear im;
clear invJac;
clear jm;
clear k;
clear matLocal_dir;
clear matLocal_esq;
clear matProp;
clear n;
clear pinteg;
%clear prop;
clear qtdNumInted;
clear r0;
clear tInteg;
clear u;
clear uInteg;
clear v;
clear vInteg;
clear w;
disp('Limpeza de variáveis: ok');

