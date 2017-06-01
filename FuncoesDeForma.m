%==========================================================================
%VISUALIZADOR DAS FUNÇÕES DE FORMA/INTERPOLAÇÃO DE SEGUNDA ORDEM
%==========================================================================
%Monografia para conclusão do curso de pós graduação
%Aluno: Anderson Santos Nunes
%Orientador: Prof. Dr. Marcelo Grafulha Vanti
%Data: 22/08/2011
%==========================================================================

%Inicialização do Matlab
clear all;
close all;
clc;

%==========================================================================
%Quantidade de divisões nas coordenadas u ev

n=100;

%==========================================================================
%Gera as coordenadas de cada ponto da malha e o valor da função no ponto

for i=1:n+1
    for j=1:n+1
        u(i,j)=(i-1)/n;
        v(i,j)=(j-1)/n;
        t=1-u(i,j)-v(i,j);

        if u(i,j)<=1.001*(1-v(i,j))

            %Espcificação da função de forma em questão
            N(i,j)=4*t*v(i,j);

            N_=N(i,j);
            u_=u(i,j);
            v_=v(i,j);
        else
            N(i,j)=N_;
            u(i,j)=u_;
            v(i,j)=v_;
        end
    end
end

%==========================================================================
%Plota a função especificada acima

surf(u,v,N)
shading interp
colorbar('vert')
xlabel('u');
ylabel('v');
zlabel('N6');