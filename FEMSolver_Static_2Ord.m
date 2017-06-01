%==========================================================================
%SOLVER ESTÁTICO DE ELEMENTOS FINTOS SUBPARAMÉTRICOS
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
%Definições para o solver
%Entrada do tipo natureza do problema

problema = input('Eletrostática (1) Eletrocinética (2) Magnetostática (3):' );
simetria=input('Planar (0) Axissimétrico (1):');

%==========================================================================

% Carrega os arquivos de dados
% load -AscII 'coordenadas.txt'   %Coordenadas XY de cada nó (x,y)
% load -AscII 'num_Global.txt'      %Numeraçao global de cada elemento
% load -AscII 'propriedades.txt'  %Índice do material e excitaçao de cada elemento (ind. material,densi. corr.)
% load -AscII 'condContorno.txt'  %Condiçoes de contorno, em cada nós (tem? (0 ou 1),valor)
% 
% n=length(coordenadas);
% numNos=n(1,1);
% 
% n=length(num_Global);
% numElem=n(1,1);
% simetria=0;
% problema=3;
% materiais;
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
qtdNumInted=3;              %Quantidade de pontos de integração

u=[1/6;2/3;1/6];
v=[1/6;1/6;2/3];
w=1/6;                      %Ponderaçao

Integdudv=.5;              %Área do elemento de referencia

%Inicializaçao com zeros das matriz principais
MatGlobal_esq=zeros(numNos,numNos); %Matriz global do lado esquerdo
MatGlobal_dir=zeros(numNos,1);      %Matriz global do lado direito

gradN_prim=[-1 1 0;                 %gradN(u,v) de primeira ordem
         -1 0 1];
     
     r0=1; %Caso seja planar

%Definiçao dos materiais
switch problema
    case 1
        eps0=8.85*power(10,-12);
        prop=listaMat(:,1)*eps0;
    case 2
        prop=listaMat(:,2);
    case 3
        mu0=4*pi*power(10,-7);
         
 
        n=size(listaMat);
        for i=1:n(1,1)
            prop(i,1)=1/(mu0*listaMat(i,3));
        end
end
disp('Pré loop principal: ok');

%==========================================================================
%Loop principal 

for k=1:numElem
    
    %Definiçao da matriz Jacobiana
    coordJ=[coordenadas(num_Global(k,1),1) coordenadas(num_Global(k,1),2);
            coordenadas(num_Global(k,3),1) coordenadas(num_Global(k,3),2);
            coordenadas(num_Global(k,5),1) coordenadas(num_Global(k,5),2)];
        
    %Cálculo do baricentro, no caso de axissimetria
    if simetria==1
        r0=(coordJ(1,1)+coordJ(2,1)+coordJ(3,1))/3;
    end
    
    Jac=gradN_prim*coordJ;      %Matriz Jacobiana
    invJac=inv(Jac);            %Matriz inversa da Jacobiana
    detJac=det(Jac);
   
    %Definiçao da propriedade do material
    matProp=prop(propriedades(k,1),1)/r0;
    
    %Integral do lado esquerdo
    %(Transosta(InvJacobiana*gradN)*(InvJacobiana*gradN)*det(Jacobiana)*Integdudv*prop.material)
    matLocal_esq=zeros(6,6);
    
    for pinteg=1:qtdNumInted
        
        uInteg=u(pinteg,1);
        vInteg=v(pinteg,1);
        tInteg=1-uInteg-vInteg;
        
        %gradN de secunda ordem
        gradN_sec=[(1-4*tInteg) (4*(tInteg-uInteg)) (4*uInteg-1) (4*vInteg) (0) (-4*vInteg);
                    (1-4*tInteg) (-4*uInteg) (0) (4*uInteg) (4*vInteg-1) (4*(tInteg-vInteg))]; 

        matLocal_esq=matLocal_esq+matProp*(invJac*gradN_sec)'*(invJac*gradN_sec)*detJac*Integdudv*w;
    end
       
    %Montagem da matriz global do lado esquerdo
    for im=1:6
        for jm=1:6
            MatGlobal_esq(num_Global(k,im),num_Global(k,jm))=MatGlobal_esq(num_Global(k,im),num_Global(k,jm))+matLocal_esq(im,jm);
        end
    end
    %Integral do lado direito
    %Integral Trans[N]*[Js]*Det[J]*dudv
    Js=propriedades(k,2);
    matLocal_dir=zeros(6,1);
     for pinteg=1:qtdNumInted
        uInteg=u(pinteg);
        vInteg=v(pinteg);
        tInteg=1-uInteg-vInteg;
        N_sec=[-(tInteg*(1-2*tInteg)) (4*uInteg*tInteg) (-uInteg*(1-2*uInteg)) (4*uInteg*vInteg) (-vInteg*(1-2*vInteg)) (4*vInteg*tInteg)];
       
        matLocal_dir=matLocal_dir+detJac*Integdudv*w*Js*N_sec';
     end
     
     %Montagem da matriz global do lado direito
     for im=1:6
         MatGlobal_dir(num_Global(k,im),1)=MatGlobal_dir(num_Global(k,im),1)+matLocal_dir(im,1);
     end

end

disp('Loop principal: ok');
%==========================================================================
%Aplicaçao das condiçoes de contorno

for i=1:numNos
        
    if condContorno(i,1)==1
        MatGlobal_esq(i,1:numNos)=0;    %Zera a linha do nó na matriz global direita
        MatGlobal_esq(i,i)=1;           %Aplica 1 na pos. do nó na matriz global direita
        MatGlobal_dir(i,1)= condContorno(i,2);  %Aplica a cond. de contorno na pos. do nó na matriz global esquerda
    end
end

disp('Aplicação das condições de contorno: ok');

%==========================================================================
%Resoluçao do sistema matricial
potenciais=zeros(numNos,1);     %Vetor contendo os resultados
potenciais=MatGlobal_esq\MatGlobal_dir;

disp('Resoluçao do sistema matricial: ok');
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
clear prop;
clear qtdNumInted;
clear r0;
clear tInteg;
clear u;
clear uInteg;
clear v;
clear vInteg;
clear w;
disp('Limpeza de variáveis: ok');

