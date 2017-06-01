
gradN_prim=[-1 1 0;                 %gradN(u,v) de primeira ordem
         -1 0 1];          
for k=1:numElem
    
%Coordenadas para a matriz Jacobiana
    coordJ=[coordenadas(num_Global(k,1),1) coordenadas(num_Global(k,1),2);
            coordenadas(num_Global(k,3),1) coordenadas(num_Global(k,3),2);
            coordenadas(num_Global(k,5),1) coordenadas(num_Global(k,5),2)];
%Matriz Jacobiana
    Jac=gradN_prim*coordJ;

% Converte o ponto do elemento real para elemento de referência
    ptUV=inv(Jac')*(ptXY-[coordJ(1,1);coordJ(1,2)]);
    u=ptUV(1);
    v=ptUV(2);
%N de primeira ordem
    N_prim=[1-u-v u v ];
    t=1-u-v;
 
%Identifica o elemento que contém o ponto(x,y) dado e calcula o
%gradPotencial neste ponto
    if N_prim(1)>=0 && N_prim(1)<=1 && N_prim(2)>=0 && N_prim(2)<=1 && N_prim(3)>=0 && N_prim(3)<=1
          N_sec=[-(t*(1-2*t)) (4*u*t) (-u*(1-2*u)) (4*u*v) (-v*(1-2*v)) (4*v*t)];
          gradN_sec=[(1-4*t) (4*(t-u)) (4*u-1) (4*v) (0) (-4*v);
                    (1-4*t) (-4*u) (0) (4*u) (4*v-1) (4*(t-v))]; 
                
%Cálculo do baricentro, no caso de axissimetria
        r0=1;
        if simetria==1
            r0=(coordJ(1,1)+coordJ(2,1)+coordJ(3,1))/3;
        end 
        
        for elem=1:6
             pot(elem,1)=potenciais(num_Global(k,elem),instante);
        end
        potXY=N_sec*pot;
        gradPotXY=inv(Jac)*(gradN_sec)*pot*(1/r0);
    end
end