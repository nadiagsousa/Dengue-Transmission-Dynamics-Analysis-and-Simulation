//Comparing vector–host and SIR models for dengue transmission
//Modelo Vector–host model
//Nádia Guimarães Sousa
//20/05/2021

clc
clear
mode(-1)

//Leitura de dados
sheets = readxls('Dados.xls')
dados=sheets(6)
[mm,nn]=size(dados)
dados_2020=dados(3:mm,7)

function dxdt=SIR(t,x)
    HS=x(1)
    HI=x(2)
    HR=x(3)
    vS=x(4)//VS/V----> retirando a dependencia de V
    vI=x(5)//VI/V
    H=HS+HI+HR
   // V=VS+VI
    Bh=miH*H
    //Bv=miv*V
    dxdt(1)=Bh-bethaaH*vI*HS-miH*HS //bethaaH=m*c*bethah
    dxdt(2)=bethaaH*vI*HS-gamaH*HI-miH*HI
    dxdt(3)=gamaH*HI-miH*HR
    dxdt(4)=miv-bethaav*HI/H*vS-miv*vS//bethaav=c*bethav
    dxdt(5)=bethaav*HI/H*vS-miv*vI
endfunction

//Programa Principal
//Entrada de dados

//bethaaH=0.072*365/12//0.0686*365/12//meses^-1 taxa de transmissão mosquito-para-humano
//bethaav=0.46*365/12//0.4307*365/12//meses^-1 taxa de transmissão humano-para-mosquito
gamaH=1/2.8*365/12//1/2.28*365/12//meses^-1 taxa de recuperação (2 a 7 dias)
miH=1/(76*12)//1/(76*365)//meses^-1 taxa de mortalidade
miv=0.058*365/12//0.0605*365/12//meses^-1


H0=21119536//46808000
HI0=3990
vI0=0.001//0.0009
vs0=1-vI0
xinicial=[H0-HI0;HI0;0;vs0;vI0]
t=0:0.1:11
//x=ode(xinicial,t(1),t,SIR)
x=xinicial
t0=3.1
betha0H=0.05*365/12//mês^-1 mosquito para humano
betha_minH=0.038*365/12//mês^-1 mosquito para humano
betha0V=0.38*365/12//mês^-1  humano para mosquito
betha_minV=0.2*365/12//mês^-1 humano para mosquito
r=0.8//mês^-1

for i=2:length(t)
    if t(i)<=t0 then
        bethaaH=betha0H
        bethaav=betha0V
        R0=bethaaH*bethaav/(miv*(miH+gamaH))
        disp('R0')
        disp(R0)
        disp('--------')
    else
        BetaH=betha_minH+(betha0H-betha_minH)*exp(-r*(t-t0))
        bethaaH=BetaH(i)
        BetaV=betha_minV+(betha0V-betha_minV)*exp(-r*(t-t0))
        bethaaV=BetaV(i)
        R0=bethaaH*bethaav/(miv*(miH+gamaH))
        disp('R0')
        disp(R0)
        disp('--------')
    end
    x(1:5,i)=ode(xinicial,t(i-1),t(i),SIR)
    xinicial=x(1:5,i)

end

//saída de dados
scf(0)
clf()
subplot(3,1,1),plot(t,x(1,:),'b')
xtitle('','$Meses$','$Indivíduos \ suscetíveis \ (H_S)$')
subplot(3,1,2),plot(t,x(2,:),'b')
xtitle('','$Meses$','$Indivíduos \ Infectados \ (H_I)$')
subplot(3,1,3),plot(t,x(3,:),'b')
xtitle('','$Meses$','$Indivíduos \ Recuperados\ (H_R)$')

//Infectados acumulado
tt=0:1:11
Pos=[]
for i=1:length(t)
    for j=1:length(tt)
        if t(i)==tt(j) then
            pos=i
            Pos=[Pos;pos]
        end
   end
end
calc=[42929,x(2,Pos)]
for i=1:length(tt)
   Calc(i)=sum(calc(1:i))
end

Dados_2020=[42929;dados_2020]
dados_2020=[HI0;dados_2020]

scf(1)
clf
subplot(2,1,1),plot(t,x(2,:),'r')
subplot(2,1,1),plot(tt,dados_2020','ko')
xtitle('Mensal','t[meses]','HI')
subplot(2,1,2),plot(tt,Calc','r')
xtitle('Acumulado','t[meses]','HI')
for i=1:length(Dados_2020)
   subplot(2,1,2),plot(tt(i),sum(Dados_2020(1:i)),'ko') 
   D_2020(i)=sum(Dados_2020(1:i))
end

scf(2)
clf()
subplot(2,1,1),plot(t,x(4,:),'r')
xtitle('Mosquitos','t[meses]','Vs/V[]')
subplot(2,1,2),plot(t,x(5,:),'r')
xtitle('','t[meses]','VI/V[]')

com_mosquitos=[t',x']
fprintfMat('com_mosquitos.xls',com_mosquitos)
acumulado_com_mosquitos=[tt',Calc,D_2020]
fprintfMat('acumulado_com_mosquitos.xls',acumulado_com_mosquitos)

//scf(3)
//clf()
//plot([0:1:11],dados_2011','r')
//plot([0:1:11],dados_2012','b')
//plot([0:1:11],dados_2013','m')
//plot([0:1:11],dados_2014','k')
//plot([0:1:11],dados_2015','g')
//plot([0:1:11],dados_2016','y-.')
//plot([0:1:11],dados_2017(2:$)','cyan')
//plot([0:1:11],dados_2018','y')
//plot([0:1:11],dados_2019','b-.')
//legend(['2011','2012','2013','2014','2015','2016','2017','2018','2019'])

scf(4)
clf()
plot(tt,Calc','b')
xtitle('','$Meses$','$Total \ de \ casos \ (H_I)$')
for i=1:length(Dados_2020)
   plot(tt(i),sum(Dados_2020(1:i)),'ko') 
   D_2020(i)=sum(Dados_2020(1:i))
end

ybar=mean(dados_2020)
SQtot=sum((dados_2020-ybar)^2)
SQexp=sum((calc-ybar)^2)
SQres=sum((dados_2020'-calc(2:$))^2)
R2C=1-SQres/SQtot
//R2=SQexp/SQtot
disp('R2 - número de casos')
disp(R2C)

ybar=mean(D_2020)
SQtot=sum((D_2020-ybar)^2)
SQexp=sum((Calc-ybar)^2)
SQres=sum((D_2020-Calc)^2)
R2CA=1-SQres/SQtot
//R2=SQexp/SQtot
disp('R2 - número de casos acumulados')
disp(R2CA)


disp('Fim')




