clear all;close all; clc;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This routine plot the validation of Aures' tonality code for pure tone signals (85 dBSPL @ 1 kHz)
%  and varying emergence level [0 10 20 30 40 50 60 70 80] dBSPL over a
%  bandpass (white noise) in the critical band centered around the tone 
%   
%   The reference data is is taken:
%   Aaron Hastings,a) Kyoung Hoon Lee,a) Patricia Davies,a) and Aimée M. Surprenantb)
%  "Measurement of the attributes of complex tonal components commonly found in product sound" 
%  Noise Control Eng. J. 51 (4), 2003 Jul–Aug
%- fig 1
%
%   Gil Felix Greco, Braunschweig, 22.03.2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% save settings

save_figs=0;

%% reference data 
%Source: Aaron Hastings,a) Kyoung Hoon Lee,a) Patricia Davies,a) and Aimée M. Surprenantb)
%  "Measurement of the attributes of complex tonal components commonly found in product sound" 
%  Noise Control Eng. J. 51 (4), 2003 Jul–Aug
%- fig 1

ref=[
0  0
10  0.07533719859035015
20  0.2106359366048477
30  0.39928132412994843
40  0.5923422380983235
50  0.7627954783966964
60  0.884892892884656
70  0.9592512633054631
80  0.9871305523642211
     ];

%% load signals

tag_str={'0','10','20','30','40','50','60','70','80'};

for i=1:9
    
[x(:,i),fs]=audioread(['MakeToneEmergence\1Bark_tone_prominence_' sprintf('%s',char(tag_str(i))) 'dB_fc_1khz_44khz_64bit.wav']);

end

%% compute tonality

for i=1:9
    
T{i} = Tonality_Aures1985(x(:,i),fs,0,0,true);

end

%%

for i=1:9
    
results(i)= [T{i}.Kmean]; % create vector with time-averaged tonality values

end

%% plot

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
 
yyaxis left

% reference results
a=plot( ref(:,1),ref(:,2),'k*-','MarkerSize',8);hold on;

% SQAT
b=plot(0:10:80,results,'ko:','MarkerSize',8);

% dummy plot to make x-axis larger
plot(1,0.1);hold on;
plot(5,0.1);hold on;

ylim([0 1.1]);
xlim([0 80]);

ylabel('Aures tonality, $K$ (t.u.)','Interpreter','Latex');
xlabel(sprintf('Signal to noise ratio in the critical band\n centered about the tone (dBSPL)'),'Interpreter','Latex'); 
grid off

legend([a,b],{ 'Literature',...
        'SQAT'},...
        'Location','SE','Interpreter','Latex');

legend boxoff

yyaxis right

plot(0:10:80,results'-ref(:,2),'o-','MarkerSize',4,'MarkerFaceColor',[0.85,0.33,0.10],'HandleVisibility','off');hold on;

ylim([-0.2 0.2]);

ylabel('SQAT minus Literature (t.u.)','Interpreter','Latex');

ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues =  0:10:80;

set(ax,'XTick',[0 20 40 60 80]);     
set(ax,'YTick',[0 .2 .4 .6 .8 1]);

yticks([-0.2 -0.1 0 0.1 0.2])

ax.YAxis(1).Color = 'k';

set(gcf,'color','w');

%%

if save_figs==1
figuresdir = 'figs\'; 
saveas(gcf,strcat(figuresdir,'tonality_validation_SNR_tone_85dBSPL_1khz'), 'fig');
saveas(gcf,strcat(figuresdir,'tonality_validation_SNR_tone_85dBSPL_1khz'), 'pdf');
saveas(gcf,strcat(figuresdir,'tonality_validation_SNR_tone_85dBSPL_1khz'), 'png');
else
end
