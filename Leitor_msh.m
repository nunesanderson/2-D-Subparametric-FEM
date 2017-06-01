%==========================================================================
%ROTINA PARA LEITURA DOS DADOS DAS MALHAS GERADAS PELO PROGRAMA GMSH
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

disp('===================================================================')
disp('Arquivo Leitor_msh')


%==========================================================================
%Leitura do arquivo
%Entrada do nome do arquivo
arquivo = input('Nome do aquivo (com extensão):','s')


% Abre o arquivo
fid = fopen(arquivo, 'r');
qtd_=0; %Quantidade de $

linhas_coordenadas_nos=1;       %Controle de linha da leitura das coordenadas
linha_cond_contorno=1;          %Controle de linha da leitura das condiçoes de contorno
linha_num_global=1;             %Controle de linha da leitura da numeração global

%Loop que le cada linha do arquivo
while 1
tline = fgetl(fid);
linha=cellstr(char(tline));
linha=char(linha);

A = sscanf(linha,'%f');

% Controla a linha com $
    if linha(1,1)=='$'
        
        qtd_=qtd_+1;
        
    end

%Cria a tabela com as coordenadas
    if qtd_==3 & length(A)>1
        
        coordenadas(linhas_coordenadas_nos,1)=A(2);
        coordenadas(linhas_coordenadas_nos,2)=A(3);
        linhas_coordenadas_nos=linhas_coordenadas_nos+1;
        
    end

% Cria a tabela com os dados das condiçoes de contorno
    if qtd_==5 & length(A)>1 & A(2)==8
        
        col_cond_contorno=1;
        
        for j=1:9
            
            cond_contorno(linha_cond_contorno,col_cond_contorno)=A(j);
            col_cond_contorno=col_cond_contorno+1;
            
        end

        linha_cond_contorno=linha_cond_contorno+1;
    end
    
% Cria a tabela com a numeraçao global
    if qtd_==5 & length(A)>1 & A(2)==9
        
        col_num_global=1;
        
        for j=1:12
            
            num_global(linha_num_global,col_num_global)=A(j);
            col_num_global=col_num_global+1;
            
        end

        linha_num_global=linha_num_global+1;
    end

%Controla o fim da leitura do arquivo
    
    if tline==-1
        break;
    end
    
end
fclose(fid);
disp('Leitura do aquivo .msh: ok');

%==========================================================================
%Limpeza de variáveis não utilizadas
clear ans;
clear A;
clear arquivo;
clear col_cond_contorno;
clear col_num_global;
clear fid;
clear j;
clear linha;
clear linha_cond_contorno;
clear linha_num_global;
clear linha_num_global;
clear qtd_;
clear tline;
clear linhas_coordenadas_nos;
disp('Limpeza de variáveis: ok');

