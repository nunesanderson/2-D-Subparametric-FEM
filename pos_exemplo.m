
if problema==4
    instante = input(['Coluna da matriz de potenciais: ']);
else
    instante=1;
end

hold on

linha=1;
mm=power(10,-3);
xlinha=[0.09;0.11];
ylinha=[0.05;0.05];
npontos=150;

dx=(xlinha(2)-xlinha(1))/npontos;
dy=(ylinha(2)-ylinha(1))/npontos;

pt=0;
for pt=1:npontos+1;

    %Gera os pontos da reta
    xmat(pt)=(pt-1)*dx+xlinha(1);
    ymat(pt)=(pt-1)*dy+ylinha(1);
    xpt(pt)=(pt-1)*dx;
    ypt(pt)=(pt-1)*dy;
    comp=sqrt(xpt(pt)^2+ypt(pt)^2);
    ptXY=[xmat(pt);ymat(pt)];
    
    %Chama o arquivo que calcula o pontencial e os gradientes no ponto
    dadoXY
    
    if problema==3|problema==4
        bx(pt,1)=comp;
        bx(pt,2)=gradPotXY(2);
        by(pt,1)=comp;
        by(pt,2)=-gradPotXY(1);
        b(pt,1)=comp;
        b(pt,2)=sqrt(bx(pt,2)^2+by(pt,2)^2);
        a(pt,1)=comp;
        a(pt,2)=potXY;
    end

    text(xmat(pt),ymat(pt),'-' ); 
end
disp('Curva: ok')


plot(a(:,1),a(:,2))