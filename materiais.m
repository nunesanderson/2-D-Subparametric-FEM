%==========================================================================
%TABELA COM OS DADOS DOS MATERIAIS
%==========================================================================
%Monografia para conclusão do curso de pós graduação
%Aluno: Anderson Santos Nunes
%Orientador: Prof. Dr. Marcelo Grafulha Vanti
%Data: 22/08/2011
%==========================================================================

%[Permissividade Condutividade Permeabilidade]
listaMat(1,:)=[1 0 1];           %ar
listaMat(2,:)=[1 58E+06 1];    %cobre
listaMat(3,:)=[1 0 5.00E+03];    %Aço Si
listaMat(4,:)=[0 5E6 300];       %Aço carbono
listaMat(5,:)=[0 1e6 50];       %Aço Condutor
listaMat(6,:)=[0 0 1000];       %Aço não condutor
listaMat(7,:)=[7 0 1];       %Aço não condutor

n=length(listaMat);
disp('Índice  Permis. Condut. Permea.')
for i=1:n(1,1)
    disp([i listaMat(i,1) listaMat(i,2) listaMat(i,3)])
end
disp('Leitura da tabela de mateiais: ok');
