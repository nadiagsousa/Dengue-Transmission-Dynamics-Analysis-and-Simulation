//Comparing vector–host and SIR models for dengue transmission
//Modelo Vector–host model
//Modificação - modelo completo para fases de mosquitos
//Controle com Larvicida
//Nádia Guimarães Sousa
//05/01/2021

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
    e=x(4)//E/V
    l=x(5)//L/V
    p=x(6)//P/V
    vS=x(7)//VS/V----> retirando a dependencia de V
    vI=x(8)//VI/V
    H=HS+HI+HR
    //V=VS+VI
    Bh=miH*H
    //Bv=miv*V
    
    dxdt(1)=Bh-bethaaH*vI*HS-miH*HS //bethaaH=m*c*bethah
    dxdt(2)=bethaaH*vI*HS-gamaH*HI-miH*HI
    dxdt(3)=gamaH*HI-miH*HR
    dxdt(4)=phi*(1-e*(vS+vI)/CL)*f-(sigmae-mie)*e
    dxdt(5)=sigmae*e-(sigmaL+miL+miLL)*l
    dxdt(6)=sigmaL*l-(sigmap+mip+mipL)*p
    dxdt(7)=sigmap*p-bethaav*HI/H*vS-(miv+mivL)*vS//bethaav=c*bethav
    dxdt(8)=bethaav*HI/H*vS-(miv+mivL)*vI

endfunction

//Programa Principal
//Entrada de dados
phi=0.74*365/12//taxa percapta de oviposição
CL=2000 //quatidade de criadouros
f=0.4 //taxa da população de mosquitos que são fêmeas
//Período desfavorável
//sigmae=1/5*365/12
//sigmaL=1/15*365/12
//sigmap=1/11*365/12
//Período favorável
sigmae=1/3.3*365/12 //taxa de transição da fase ovo para larva (sigmae^-1 é o período médio de eclosão)
sigmaL=1/8*365/12//taxa de transição da fase larva para pupa 
sigmap=1/3.1*365/12//taxa de transição da fase pupa para mosquito adulto
//-------------------------------------
mie=1/100*365/12//taxa de mortalidade natural
miL=1/2*365/12
mip=1/60*365/12
//Larvicida
//miLL=0
//mipL=0
// Inseticida
mivL=0
//bethaaH=0.072*365/12//0.0686*365/12//meses^-1 taxa de transmissão mosquito-para-humano
//bethaav=0.46*365/12//0.4307*365/12//meses^-1 taxa de transmissão humano-para-mosquito
gamaH=1/3*365/12//1/2.28*365/12//meses^-1 taxa de recuperação (2 a 7 dias)
miH=1/(76*12)//1/(76*365)//meses^-1 taxa de mortalidade
miv=0.058*365/12//0.0605*365/12//meses^-1


H0=21119536//46808000
HI0=3990
v0=50
vI0=0.001//0.0009
vs0=1-vI0//1-vI0
e0=0
l0=0
p0=0
xinicial=[H0-HI0;HI0;0;e0;l0;p0;vs0;vI0]
t=0:0.1:11
//x=ode(xinicial,t(1),t,SIR)
x=xinicial
t0=3.1//3.1
betha0H=0.05*365/12//0.062*365/12//mês^-1 mosquito para humano
betha_minH=0.038*365/12//0.032*365/12//mês^-1 mosquito para humano
betha0V=0.38*365/12//0.38*365/12//mês^-1  humano para mosquito
betha_minV=0.2*365/12//0.2*365/12//mês^-1 humano para mosquito
r=0.8//0.8//mês^-1

for i=2:length(t)
    if t(i)<=t0 then
        bethaaH=betha0H
        bethaav=betha0V
       // R0(i)=bethaaH*bethaav/(miv*(miH+gamaH))
      //  disp('R0')
      //  disp(R0(i))
     //   disp('--------')
    else
        BetaH=betha_minH+(betha0H-betha_minH)*exp(-r*(t-t0))
        bethaaH=BetaH(i)
        BetaV=betha_minV+(betha0V-betha_minV)*exp(-r*(t-t0))
        bethaaV=BetaV(i)
     //   R0(i)=bethaaH*bethaav/(miv*(miH+gamaH))
    //    disp('R0')
    //    disp(R0(i))
    //    disp('--------')
    end
  // Larvicida 
    if t(i)>=60/30 & t(i)<=110/30 then
      miLL=1/1.8*365/12//1/2*365/12 //milL=2 
      mipL=1/55*365/12//1/60*365/12

  else
      miLL=0 
      mipL=0
    end  
    x(1:8,i)=ode(xinicial,t(i-1),t(i),SIR)
    xinicial=x(1:8,i)

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
   D_2020(i)=sum(dados_2020(1:i))
end

scf(2)
clf()
subplot(2,1,1),plot(t,x(7,:),'r')
xtitle('Mosquitos','t[meses]','Vs/V[]')
subplot(2,1,2),plot(t,x(8,:),'r')
xtitle('','t[meses]','VI/V[]')

scf(3)
clf()
subplot(3,1,1),plot(t,x(7,:)+x(8,:),'k')
xtitle('','$Meses$','$v[-]$')
subplot(3,1,2),plot(t,x(7,:)+x(8,:),'k')
xtitle('','$Meses$','$v[-]$')
subplot(3,1,3),plot(t,x(7,:)+x(8,:),'k')
xtitle('','$Meses$','$v[-]$')

scf(4)
clf()
plot(tt,Calc','k')
xtitle('','$Meses$','$Total \ de \ casos \ (H_I)$')
for i=1:length(Dados_2020)
   plot(tt(i),sum(Dados_2020(1:i)),'ko') 
   D_2020(i)=sum(Dados_2020(1:i))
end

controle_larvicida=[t',x']
fprintfMat('controle_larvicida.xls',controle_larvicida)
acumulado_larvicida=[tt',Calc,D_2020]
fprintfMat('acumulado_larvicida.xls',acumulado_larvicida)

//scf(5)
//clf()
//R0(1)=R0(2)
//plot(t,R0,'b-.')
//xtitle('','$Meses$','$Reprodutividade\ basal - R_0$')


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




