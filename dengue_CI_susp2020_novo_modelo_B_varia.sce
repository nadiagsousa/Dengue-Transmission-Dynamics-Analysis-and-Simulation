//Comparing vector–host and SIR models for dengue transmission
//Modelo SIR
//Nádia Guimarães Sousa
//Minas Gerais 2020
//20/05/2021

clc
clear
mode(-1)

//Leitura de dados
sheets = readxls('Dados.xls')
dados=sheets(6)
[mm,nn]=size(dados)
dados_2020=dados(3:mm,7)

function f=ee(x)
    HS=x(1)
    HI=x(2)
    HR=x(3)
    H=HS+HI+HR
    Bh=miH*H
    f(1)=Bh-betha*HI/H*HS-miH*HS
    f(2)=betha*HI/H*HS-gamaH*HI-miH*HI
    f(3)=gamaH*HI-miH*HR
endfunction

function dxdt=SIR(t,x)
    HS=x(1)
    HI=x(2)
    HR=x(3)
    H=HS+HI+HR
    Bh=miH*H
    dxdt(1)=Bh-betha*HI/H*HS-miH*HS
    dxdt(2)=betha*HI/H*HS-gamaH*HI-miH*HI
    dxdt(3)=gamaH*HI-miH*HR
endfunction

//Programa Principal
//Entrada de dados
betha_aux=0.513*365/12//0.5//mês^-1 taxa de transmissão humano-para-humano//0.511*365/12
gamaH=1/2.55*365/12//1/2.49*365/12//mês^-1 taxa de recuperação (2 a 7 dias)//1/2.26*365/12
miH=1/(76*12)//mês^-1 taxa de mortalidade
//R0=betha/(miH+gamaH)
//disp(R0)

//Processamento
//xchute=[15674842;0.4321982;5444693.5]
//[sol,v,info]=fsolve(xchute,ee)
//if info==1 then
//    disp(sol)
//else
//    disp('tente novamente')
//end


H0=21119536
HI0=3990//dados_2020(1)//42929//1000

xinicial=[H0-HI0;HI0;0]
t=0:0.1:11

x=xinicial
t0=0.695
//t1=6
betha0=0.399*365/12//0.471*365/12//mês^-1
betha_min=0.197*365/12//0.3*365/12//mês^-1
r=0.095//0.0099//mês^-1
for i=2:length(t)
    if t(i)<=t0 then
        betha=betha0
        R0(i)=betha/(miH+gamaH)
        disp('R0')
        disp(R0(i))
        disp('--------')
    else
        Beta=betha_min+(betha0-betha_min)*exp(-r*(t-t0))
        betha=Beta(i)
        R0(i)=betha/(miH+gamaH)
        disp('R0')
        disp(R0(i))
        disp('--------')
    end
    x(1:3,i)=ode(xinicial,t(i-1),t(i),SIR)
    xinicial=x(1:3,i)

end

//saída de dados
scf(0)
clf()
subplot(3,1,1),plot(t,x(1,:),'r')
//xtitle('','t[meses]','HS')
subplot(3,1,2),plot(t,x(2,:),'r')
//xtitle('','t[meses]','HI')
subplot(3,1,3),plot(t,x(3,:),'r')
//xtitle('','t[meses]','HR')

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
subplot(2,1,1),plot(t,x(2,:),'r-.')
subplot(2,1,1),plot(tt,dados_2020','ko')
xtitle('Mensal','t[meses]','HI')
subplot(2,1,2),plot(tt,Calc','r-.')
xtitle('Acumulado','t[meses]','HI')
for i=1:length(Dados_2020)
   subplot(2,1,2),plot(tt(i),sum(Dados_2020(1:i)),'ko') 
   D_2020(i)=sum(Dados_2020(1:i))
end

sem_mosquitos=[t',x']
fprintfMat('sem_mosquitos.xls',sem_mosquitos)



acumulado_sem_mosquitos=[tt',Calc,D_2020]
fprintfMat('acumulado_sem_mosquitos.xls',acumulado_sem_mosquitos)

scf(4)
clf()
plot(tt,Calc','r')
//xtitle('Acumulado','t[meses]','HI')
//for i=1:length(dados_2017)
//   plot(tt(i),sum(dados_2017(1:i)),'ko') 
//   D_2017(i)=sum(dados_2017(1:i))
//end

//scf(2)
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

scf(5)
clf()
R0(1)=R0(2)
plot(t,R0,'r')
xtitle('','$Meses$','$Reprodutividade\ basal - R_0$')

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




