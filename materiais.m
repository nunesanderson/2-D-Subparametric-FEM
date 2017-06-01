%==========================================================================
%TABELA COM OS DADOS DOS MATERIAIS
%==========================================================================
%Monografia para conclus�o do curso de p�s gradua��o
%Aluno: Anderson Santos Nunes
%Orientador: Prof. Dr. Marcelo Grafulha Vanti
%Data: 22/08/2011
%==========================================================================

%[Permissividade Condutividade Permeabilidade]
listaMat(1,:)=[1 0 1];           %ar
listaMat(2,:)=[1 58E+06 1];    %cobre
listaMat(3,:)=[1 0 5.00E+03];    %A�o Si
listaMat(4,:)=[0 5E6 300];       %A�o carbono
listaMat(5,:)=[0 1e6 50];       %A�o Condutor
listaMat(6,:)=[0 0 1000];       %A�o n�o condutor
listaMat(7,:)=[7 0 1];       %A�o n�o condutor

n=length(listaMat);
disp('�ndice  Permis. Condut. Permea.')
for i=1:n(1,1)
    disp([i listaMat(i,1) listaMat(i,2) listaMat(i,3)])
end
disp('Leitura da tabela de mateiais: ok');
