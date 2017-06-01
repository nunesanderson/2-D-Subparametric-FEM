%==========================================================================
%VISUALIZADOR DAS FUN��ES DE FORMA/INTERPOLA��O DE SEGUNDA ORDEM
%==========================================================================
%Monografia para conclus�o do curso de p�s gradua��o
%Aluno: Anderson Santos Nunes
%Orientador: Prof. Dr. Marcelo Grafulha Vanti
%Data: 22/08/2011
%==========================================================================

%Inicializa��o do Matlab
clear all;
close all;
clc;

%==========================================================================
%Quantidade de divis�es nas coordenadas u ev

n=100;

%==========================================================================
%Gera as coordenadas de cada ponto da malha e o valor da fun��o no ponto

for i=1:n+1
    for j=1:n+1
        u(i,j)=(i-1)/n;
        v(i,j)=(j-1)/n;
        t=1-u(i,j)-v(i,j);

        if u(i,j)<=1.001*(1-v(i,j))

            %Espcifica��o da fun��o de forma em quest�o
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
%Plota a fun��o especificada acima

surf(u,v,N)
shading interp
colorbar('vert')
xlabel('u');
ylabel('v');
zlabel('N6');