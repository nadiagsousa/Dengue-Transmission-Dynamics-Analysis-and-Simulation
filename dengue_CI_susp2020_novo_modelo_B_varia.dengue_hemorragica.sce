//Comparing vector–host and SIR models for dengue transmission
//Modelo SIR
//Nádia Guimarães Sousa
//08/09/2021

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
    HC=x(4)
    H=HS+HI+HR+HC
    Bh=miH*H
    dxdt(1)=Bh-betha*HI/H*HS-miH*HS
    dxdt(2)=betha*HI/H*HS-gamaH*HI-miH*HI
    dxdt(3)=gamaH*HI-miH*HR
    dxdt(4)=p*betha*(HI/H)*HS
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
p=0.0022

H0=21119536
HI0=3990//dados_2020(1)//42929//1000

xinicial=[H0-HI0;HI0;0;0]
t=0:0.1:11

x=xinicial
t0=0.695//1
//t1=6
betha0=0.399*365/12//0.42*365/12//0.471*365/12//mês^-1
betha_min=0.197*365/12//0.299*365/12//0.3*365/12//mês^-1
r=0.095//0.1//0.0099//mês^-1
for i=2:length(t)
    if t(i)<=t0 then
        betha=betha0
        R0=betha/(miH+gamaH)
        disp('R0')
        disp(R0)
        disp('--------')
    else
        Beta=betha_min+(betha0-betha_min)*exp(-r*(t-t0))
        betha=Beta(i)
        R0=betha/(miH+gamaH)
        disp('R0')
        disp(R0)
        disp('--------')
    end
    x(1:4,i)=ode(xinicial,t(i-1),t(i),SIR)
    xinicial=x(1:4,i)

end

//saída de dados
scf(0)
clf()
subplot(2,2,1),plot(t,x(1,:),'r-.')
xtitle('','t[meses]','HS')
subplot(2,2,2),plot(t,x(2,:),'b-.')
xtitle('','t[meses]','HI')
subplot(2,2,3),plot(t,x(3,:),'b-.')
xtitle('','t[meses]','HR')
subplot(2,2,4),plot(t,x(4,:),'r')
xtitle('','t[meses]','HC')

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
for i=1:length(dados_2020)
   subplot(2,1,2),plot(tt(i),sum(dados_2020(1:i)),'ko') 
   D_2020(i)=sum(dados_2020(1:i))
end


scf(2)
clf
plot(t,x(4,:),'r')
xtitle('','$Meses$','$Total \ de \ casos \ (H_C)$')

H_sem_mosquitos=[x(4,:)']
fprintfMat('H_sem_mosquitos.xls',H_sem_mosquitos)

aaa
ybar=mean(dados_2017)
SQtot=sum((dados_2017-ybar)^2)
SQexp=sum((calc-ybar)^2)
SQres=sum((dados_2017'-calc)^2)
R2C=1-SQres/SQtot
//R2=SQexp/SQtot
disp('R2 - número de casos')
disp(R2C)

ybar=mean(D_2017)
SQtot=sum((D_2017-ybar)^2)
SQexp=sum((Calc-ybar)^2)
SQres=sum((D_2017-Calc)^2)
R2CA=1-SQres/SQtot
//R2=SQexp/SQtot
disp('R2 - número de casos acumulados')
disp(R2CA)


disp('Fim')




