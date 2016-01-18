ft=16
load 'C:\Users\celestink\Documents\TestingStation\ChngStiffreducedAv\folders200\datas\StressMax_200.txt'
load 'C:\Users\celestink\Documents\TestingStation\ChngStiffreducedAv\folders200\datas\StiffSD_200.txt'
load 'C:\Users\celestink\Documents\TestingStation\ChngStiffreducedAv\folders200\datas\StiffAv_200.txt'
load 'C:\Users\celestink\Documents\TestingStation\ChngStiffreducedAv\folders200\datas\stressAv_200.txt'
load 'C:\Users\celestink\Documents\TestingStation\ChngStiffreducedAv\folders200\datas\princMax_200.txt'
load 'C:\Users\celestink\Documents\TestingStation\ChngStiffreducedAv\folders200\datas\princAv_200.txt'

% load 'C:\Users\celestink\Documents\TestingStation\InclModels\ChangingSDx10_3\datas\StressMax_3.txt'
% load 'C:\Users\celestink\Documents\TestingStation\InclModels\ChangingSDx10_3\datas\StiffSD_3.txt'
% load 'C:\Users\celestink\Documents\TestingStation\InclModels\ChangingSDx10_3\datas\StiffAv_3.txt'
% load 'C:\Users\celestink\Documents\TestingStation\InclModels\ChangingSDx10_3\datas\stressAv_3.txt'
% 
% load 'C:\Users\celestink\Documents\TestingStation\InclModels\ChangingSDx10_4\datas\StressMax_4.txt'
% load 'C:\Users\celestink\Documents\TestingStation\InclModels\ChangingSDx10_4\datas\StiffSD_4.txt'
% load 'C:\Users\celestink\Documents\TestingStation\InclModels\ChangingSDx10_4\datas\StiffAv_4.txt'
% load 'C:\Users\celestink\Documents\TestingStation\InclModels\ChangingSDx10_4\datas\stressAv_4.txt'

figure(1)

axes1=axes('Parent',figure(1),'FontSize',ft,'FontName','Times New Roman'); 
hold(axes1,'on');
  scatter(StiffSD_200,StressMax_200,40,'d','MarkerEdgeColor','k'),hold on;
%   scatter(StiffSD_3,StressMax_3,40,'o','MarkerEdgeColor','b'),hold on; 
%   scatter(StiffSD_4,StressMax_4,40,'+','MarkerEdgeColor','r'),hold on;
%   plot(StiffSD_3,StressMax_3,'b'),hold on; 
  plot(StiffSD_200,StressMax_200,'k'),hold on;
%   plot(StiffSD_4,StressMax_4,'r')
  xlabel('SD of track modulus (MPa)','FontName','Times New Roman', 'FontSize',ft);
  ylabel('Maximum longitudinal stress (MPa)','FontName','Times New Roman',  'FontSize',ft);
legend('mean stiffness=200 MPa')
 figure(2)
 axes2=axes('Parent',figure(2),'FontSize',ft,'FontName','Times New Roman'); 
 hold(axes2,'on');
  scatter(StiffSD_200,stressAv_200,40,'d','MarkerEdgeColor','k'),hold on;
%   scatter(StiffSD_3,stressAv_3,40,'o','MarkerEdgeColor','b'),hold on; 
%   scatter(StiffSD_4,stressAv_4,40,'+','MarkerEdgeColor','r'),hold on;
%  plot(StiffSD_3,stressAv_3,'b'),hold on; 
  plot(StiffSD_200,stressAv_200,'k'),hold on;
 % plot(StiffSD_4,stressAv_4,'r');
  xlabel('SD of track modulus (MPa)','FontName','Times New Roman', 'FontSize',ft);
  ylabel('Average maximum longitudinal stress (MPa)','FontName','Times New Roman',  'FontSize',ft);
legend('mean stiffness=200 MPa')

 figure(3)
 axes3=axes('Parent',figure(3),'FontSize',ft,'FontName','Times New Roman');
 hold(axes3,'on');
  scatter(StiffSD_200,princAv_200,40,'d','MarkerEdgeColor','k'),hold on;
%   scatter(StiffSD_3,stressAv_3,40,'o','MarkerEdgeColor','b'),hold on; 
%   scatter(StiffSD_4,stressAv_4,40,'+','MarkerEdgeColor','r'),hold on;
%  plot(StiffSD_3,stressAv_3,'b'),hold on; 
  plot(StiffSD_200,princAv_200,'k'),hold on;
 % plot(StiffSD_4,stressAv_4,'r');
  xlabel('SD of track modulus (MPa)','FontName','Times New Roman', 'FontSize',ft);
  ylabel('Average maximum principal stress (MPa)','FontName','Times New Roman',  'FontSize',ft);
legend('mean stiffness=20 0MPa')
figure(4)

axes4=axes('Parent',figure(4),'FontSize',ft,'FontName','Times New Roman');  
hold(axes4,'on');
  scatter(StiffSD_200,princMax_200,40,'d','MarkerEdgeColor','k'),hold on;
%   scatter(StiffSD_3,StressMax_3,40,'o','MarkerEdgeColor','b'),hold on; 
%   scatter(StiffSD_4,StressMax_4,40,'+','MarkerEdgeColor','r'),hold on;
%   plot(StiffSD_3,StressMax_3,'b'),hold on; 
  plot(StiffSD_200,princMax_200,'k'),hold on;
%   plot(StiffSD_4,StressMax_4,'r')
  xlabel('SD of track modulus (MPa)','FontName','Times New Roman', 'FontSize',ft);
  ylabel('Maximum principal stress (MPa)','FontName','Times New Roman',  'FontSize',ft);
legend('mean stiffness=200 MPa')

figure(5)

axes5=axes('Parent',figure(5),'FontSize',ft,'FontName','Times New Roman');  
hold(axes5,'on');
  scatter(StiffSD_200,StiffAv_200,40,'d','MarkerEdgeColor','k'),hold on;
%   scatter(StiffSD_3,StressMax_3,40,'o','MarkerEdgeColor','b'),hold on; 
%   scatter(StiffSD_4,StressMax_4,40,'+','MarkerEdgeColor','r'),hold on;
%   plot(StiffSD_3,StressMax_3,'b'),hold on; 
  plot(StiffSD_200,StiffAv_200,'k'),hold on;
%   plot(StiffSD_4,StressMax_4,'r')
  xlabel('SD of track modulus (MPa)','FontName','Times New Roman', 'FontSize',ft);
  ylabel('Average of track modulus (MPa)','FontName','Times New Roman',  'FontSize',ft);
legend('mean stiffness=200 MPa')
