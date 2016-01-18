clear all
SDfile=10;
numfiles=1;
maxN=4.60387725965*10^21;
b=[-0.05:-0.02: -0.12];
lb=length(b)
%ls=lb;
sigmf=[500:100:1200];
ls=lb;
ft=14;
ctN=19;
myEffStress=cell(SDfile,numfiles);
myMeanStresses=cell(SDfile,numfiles);
myStresses=cell(SDfile,numfiles);
myPrinc=cell(SDfile,numfiles);
myStiff=cell(1,SDfile);
myDefl=cell(SDfile,numfiles);
Nfun=cell(1,ls);
for k=1:SDfile
    myStiffName=sprintf('Concrete_E_%d.txt',k-1);
    myStiff{k}=importdata(myStiffName);
    StiffM(:,k)= myStiff{k};
    StiffMS=sort(StiffM,1);
  l=1
        myStressName=sprintf('RailMaxStress_%d_%d.txt',k-1,l-1);
        myPrincName=sprintf('RailMaxPrinc_%d%d.txt',k-1,l-1);
        myDeflname=sprintf('RailMaxDisp_%d%d.txt',k-1,l-1);
        myEffname=sprintf('Eff_Stress%d.txt',k-1);
        myMeanname=sprintf('Mean_Stress%d.txt',k-1);
        sigmae{k}=importdata(myEffname);
        sigmam{k}=importdata(myMeanname);
        myStresses{k,l}=importdata(myStressName);
        myDefl{k,l}=importdata(myDeflname);
        myPrinc{k,l}=importdata(myPrincName);
       DeflD(:,k)=myDefl{k,l}
        L=size(myStresses{k,l});
        stressM((l-1)*L(1)+1:l*L(1),k)=myStresses{k,l};
        princM((l-1)*L(1)+1:l*L(1),k)=myPrinc{k,l}(:,1);
      for j=1:ls
          if ls==length(sigmf)
    Nfun{j}(:,k)=(1.0/2.0*(sigmae{k}/sigmf(j)).^(1.0/b(6))).*((1.0-sigmam{k}/sigmf(j)).^(-1.0/b(6)));
          else
    Nfun{j}(:,k)=(1.0/2.0*(sigmae{k}/sigmf(end-2)).^(1.0/b(j))).*((1.0-sigmam{k}/sigmf(end-2)).^(-1.0/b(j))); 
          end
      end
    stressAv(k,1)=mean(stressM(:,k));
    princAv(k,1)=mean(princM(:,k));
    StressMax(k)=max(stressM(:,k));
    princMax(k,1)=max(princM(:,k));
    StiffSD(k)=std(myStiff{k},1);
    Deflmax(k,1)=max(DeflD(:,k));
    DeflAv(k,1)=mean(DeflD(:,k));
    StiffAv(k,1)=mean(myStiff{k});
    Stiff_ND(:,k)=1.0/(StiffSD(k)*sqrt(2*pi))*exp(-((StiffMS(:,k)-StiffAv(k)).^2)./(2*(StiffSD(k))^2));
    for j=1:ls
    Ncymax(k,j)=max(Nfun{j}(:,k));
    Ncymin(k,j)=min(Nfun{j}(:,k));    
    end
end 
for j=1:ls
Nmx(j)=max(max(Nfun{j}));
Nfunmax(:,j)=Ncymax(:,j)/Nmx(j);
Nfunmin(:,j)=Ncymin(:,j)/Nmx(j);
NfunminN(:,j)=Ncymin(:,j)./Ncymax(:,j);
end
save cycles
save StiffAv_150.txt -ascii -double StiffAv
save StressMax_150.txt -ascii -double StressMax
save princMax_150.txt -ascii -double princMax
save StiffSD_150.txt -ascii -double StiffSD
save stressAv_150.txt -ascii -double stressAv
save princAv_150.txt -ascii -double princAv
save StiffMS_150.txt -ascii -double StiffMS
save Stiff_ND_150.txt -ascii -double Stiff_ND
save Defldata.txt -ascii -double DeflD

figure(1)
axes1=axes('Parent',figure(1),'FontSize',16,'FontName','Times New Roman');
% xlim(axes1,[0 Rl]);
hold(axes1,'on');
box(axes1,'on');
scatter(StiffAv,Deflmax,'d','filled')
xlabel('Average of track moduli (MPa)','FontSize',16)
ylabel('Maximum deflections (mm)','FontSize',16)

figure(2)
axes2=axes('Parent',figure(2),'FontSize',16,'FontName','Times New Roman');
% xlim(axes1,[0 Rl]);
hold(axes2,'on');
box(axes2,'on');
scatter(StiffAv,princMax,'d','filled')
xlabel('Average of track  moduli (MPa)','FontSize',16)
ylabel('Maximum principal stress (MPa)','FontSize',16)

figure(3)
axes3=axes('Parent',figure(3),'FontSize',16,'FontName','Times New Roman');
% xlim(axes1,[0 Rl]);
hold(axes3,'on');
box(axes3,'on');
scatter(StiffAv,princAv,'s','filled')
xlabel('Average of track moduli (MPa)','FontSize',16)
ylabel('Average of maximum principal stress (MPa)','FontSize',16)

figure(4)
axes4=axes('Parent',figure(4),'FontSize',16,'FontName','Times New Roman');
% xlim(axes1,[0 Rl]);
hold(axes4,'on');
box(axes4,'on');
scatter(StiffAv,DeflAv,'s','filled')
xlabel('Average of track  moduli (MPa)','FontSize',16)
ylabel('Average of maximum deflections (mm)','FontSize',16)

for k=1:ls
   figure(k+5)
axe{k+5}=axes('Parent',figure(k+5),'FontSize',ft,'FontName','Times New Roman');
xlim(axes1,[0 600]);
hold(axe{k+5},'on');
box(axe{k+5},'on');
%scatter(StiffAv,Nfunmax,'s','filled')
scatter(StiffAv,Nfunmin(:,k),'o','filled','MarkerEdgeColor','k')
scatter(StiffAv,NfunminN(:,k),'s','filled','MarkerEdgeColor','b')
set(gca,'yscale','log');
xlabel('Average of track moduli (MPa)','FontSize',ft)
ylabel('Normalized load cycles to  fatigue life','FontSize',ft)
legend('minimum fatigue normalized to highests maximum fatigue up to 600Mpa modulus','minimum fatigue normalized to highest maximum fatigue for each average of moduli')
end

for k=1:1
   figure(k+20)
axe{k+20}=axes('Parent',figure(k+20),'FontSize',ft,'FontName','Times New Roman');
xlim(axes1,[0 600]);
hold(axe{k+20},'on');
box(axe{k+20},'on');
%scatter(StiffAv,Nfunmax,'s','filled')
linestyle={'o','s','d','*','+','<','>','v'};
linecolor={'k','b','r','g','k','b','r','g'};

for j=1:ls
    if ls==length(sigmf)
        leg{j}=['{\it\sigma_f} = ',num2str(sigmf(j)),' MPa'];
    else
      leg{j}=['{\itb} = ',num2str(b(j))];
    end
scatter(StiffAv,Nfunmin(:,j),60,linestyle{j},'filled','DisplayName',leg{j},'MarkerEdgeColor',linecolor{j});
% scatter(StiffAv,NfunminN(:,2),'s','filled','MarkerEdgeColor','b')
% scatter(StiffAv,NfunminN(:,3),'*','filled','MarkerEdgeColor','r')
% scatter(StiffAv,NfunminN(:,4),'+','filled','MarkerEdgeColor','g')
% scatter(StiffAv,NfunminN(:,5),'<','filled','MarkerEdgeColor','k')
% scatter(StiffAv,NfunminN(:,6),'>','filled','MarkerEdgeColor','b')
% scatter(StiffAv,NfunminN(:,7),'v','filled','MarkerEdgeColor','r')
set(gca,'yscale','log');
end
xlabel('Average of track moduli (MPa)','interpreter','Latex','FontSize',ft);
ylabel('Normalized  load cycles  to  fatigue life','FontSize',ft);
% if ls==8
% legend('b=-0.05','b=-0.06','b=-0.07','b=-0.08','b=-0.09','b=-0.1','b=-0.11','b=-0.12');
% else
legend('show');
% end
end
for k=1:1
    m=30
figure(m)
axe{m}=axes('Parent',figure(m),'FontSize',ft,'FontName','Times New Roman');
xlim(axes1,[0 600]);
hold(axe{m},'on');
box(axe{m},'on');
%scatter(StiffAv,Nfunmax,'s','filled')
semilogy(StiffAv,(NfunminN(:,1)),'MarkerEdgeColor','k')
semilogy(StiffAv,(NfunminN(:,2)),'MarkerEdgeColor','b')
semilogy(StiffAv,(NfunminN(:,3)),'MarkerEdgeColor','r')
semilogy(StiffAv,(NfunminN(:,4)),'g')
semilogy(StiffAv,(NfunminN(:,5)),'k')
semilogy(StiffAv,(NfunminN(:,6)),'b')
semilogy(StiffAv,(NfunminN(:,7)),'r')
set(gca,'yscale','log');
xlabel(' Average  of track  moduli  (MPa) ' ,'FontSize',ft);
ylabel('Normalized  load  cycles  to  fatigue  life ' ,'FontSize',ft);
legend('b=-0.05','b=-0.06','b=-0.07','b=-0.08','b=-0.09','b=-0.1','b=-0.11','b=-0.12');
end
