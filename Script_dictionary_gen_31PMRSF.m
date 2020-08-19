%% Input experimental setup
clear all
close all
clc
b1Scale = 1; %b1 scale, actual flip angle = designed flip angle * b1Scale
gamma = 17.235; %gyromagnetic, unit: MHz/T
b0 = 9.4; %main field strength, unit: T

load('MRSF_Seq_#0.mat');
addpath('Dictionary function_folders _pcodes'); %-> pcode folder for the functions 
%% Generate model variables for dictionary simulation
clear Variable
%input the range and step size of variables here
%------------- Set T1, T2, M0, off-res for each compartment ---------------
Variable.pool = struct('T1',[],'T2',[],'M0',[],'OffRes',[],...
    'rOffResPool',[],'rOffRes',[],'OutFluxRate',[]);
nPool = 2;
%for example, pool 1 is assigned for PCr
Variable.pool(1).T1 = 3.3; %T1,unit: s
Variable.pool(1).T2 = 0.12; %T2,unit: s
Variable.pool(1).M0 = 3; %equilibrium magnetization, a.u.
Variable.pool(1).OffRes = 0; %apparent off resonance in one voxel, relative to rotating frame, unit: Hz
Variable.pool(1).OutFluxRate = cell(1,nPool);
Variable.pool(1).OutFluxRate{2} = [0:0.05:0.5]; %rate of MT from PCr to ATP
%for example, pool 2 is assigned for ATP 
Variable.pool(2).T1 = 0.8;
Variable.pool(2).T2 = 0.016; 
Variable.pool(2).M0 = 1;
Variable.pool(2).OffRes = []; %the off res will be a fixed value relative another pool
Variable.pool(2).rOffResPool = 1; %offset will be relative to pool1
Variable.pool(2).rOffRes = gamma*b0*2.4; %offset value (Hz)
Variable.pool(2).OutFluxRate = cell(1,nPool);
Variable.pool(2).OutFluxRate{1} = 'Equil'; %flux of MT from ATP to PCr will be in equilibrium with PCr to ATP 

%-------- Set off-resonance within voxel shared by all compartments -------
Variable.nSpins = 51; %number of spins in one voxel 
Variable.vOffResMax = [-75; +75];%maximal off res relative to the apparent off res in one voxel, unit: Hz
lwfwhm = 15; %prepare parameter for Lorentzian linewidth, 1 value at a time, unit: Hz
Variable.wVOffRes = func_lineshape('Lorentzian',Variable.nSpins,lwfwhm); %the probability of each off res value 
Variable.wVOffRes = Variable.wVOffRes./sum(Variable.wVOffRes);

%% Perform Bloch simulation
nSubDic = 1;
nDic = 1;
dummy = 0;

%seqDMpar: Number of dictionary entry, Number of Variables
%seqDMdicFID: 3*nPool(Mx1,My1,Mz1,...), adcNpt, NR, Number of Variables  
tic
DictionaryFID = func_dictionary_generation(Sequence, Variable, b1Scale, gamma, b0, dummy, nSubDic, nDic);
toc

%FID should be sum of Mxy in PCr and ATP 
fidAll = cat(4,DictionaryFID.Entry); %extract all dictionary entries
entryFIDsum = squeeze(sum(fidAll(1:3:end,:,:,:),1)+1i*sum(fidAll(2:3:end,:,:,:),1)); %Mx+1i*My
NR = size(entryFIDsum,2);
% do FFT along time domain to get spectrum
entrySpec = squeeze(fftshift(fft(ifftshift(entryFIDsum,1),[],1),1));
pcrSignal = reshape(multiaa(squeeze(entrySpec(36,:,:))),NR,[]);
atpSignal = reshape(multiaa(squeeze(entrySpec(33,:,:))),NR,[]);
% organize the data according to PCr/ATP acquisition frame
fingerprint = zeros(NR,length(DictionaryFID));
for ii = 1:16
    fingerprint((1:10)+ii*20-20,:) = pcrSignal((1:10)+ii*20-20,:);
    fingerprint((1:10)+ii*20-10,:) = atpSignal((1:10)+ii*20-10,:); 
end
% assign fingerprint to dictionary structure
Dictionary = struct('T1',{DictionaryFID.T1},'T2',{DictionaryFID.T2},...
    'M0',{DictionaryFID.M0},'OffRes',{DictionaryFID.OffRes},...
    'MtRate',{DictionaryFID.MtRate},...
    'Entry',mat2cell(fingerprint,NR,ones(1,length(DictionaryFID))));
%if you need all T1 values in a matrix, use T1 = cat(1,Dictionary.T1), T1 has size: #Entry x #Compartment
%if you need all fingerprints in a matrix, use dic = cat(2,Dictionary.Entry).', dic has size: #Entry x #TR
%%
kf = Variable.pool(1).OutFluxRate{1,2};
c = 1; 
cc = 10; 
figure('name','Simulated fingerprints')
hold on
plot(real(pcrSignal(:,c)),'k','LineWidth',1.5)
plot(real(pcrSignal(:,cc)),'r','LineWidth',1.5)
plot(real(atpSignal(:,c)-60),'k','LineWidth',1.5)
plot(real(atpSignal(:,cc)-60),'r','LineWidth',1.5)
xlim([0 320])
ylim([-80 120])
set(gca,'XTick',0:80:320)
ylabel('Signal (a.u.)')
legend(sprintf('Kf = %g',kf(c)),sprintf('Kf = %g',kf(cc)))
xlabel('#TR')
set(gca,'FontSize',14,'FontWeight','bold','LineWidth',1.5,'box','on')