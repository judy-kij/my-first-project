function run_e205b_paired_policy_validation()
% E20.5B | INDEPENDENT PAIRED-NOISE POLICY VALIDATION
% ----------------------------------------------------------------------
% Purpose: compare FROZEN E20.5A online feedback against frozen
% E20.4/v4 fixed-low, fixed-high, and E20.3B latched feedback, with
% IDENTICAL stochastic noise realizations within each paired scenario.
%
% The code DOES NOT optimize or retune any policy with validation data.
% This new experiment covers a predeclared, REDUCED boundary-focused
% grid; it does not repeat all 30 previous fault conditions.
%
% Important scientific restrictions:
%   - True dR is used only by PLANT measurement and physical update.
%   - Common random numbers are an evaluation technique, not information
%     provided to online controllers.
%   - Resource caps are checked PATHWISE, not using means alone.
%   - Reused sequential threshold is fixed, not recalibrated here.
%   - Bonferroni CP bounds certify individual policy PD/FAR; paired
%     exact McNemar tests compare detections on the same noise paths.
%   - Paired tests do not establish global optimality or strict physical
%     coupling beyond the archived model.
%   - UNAVAILABLE means unresolved at this session's stop/deadline;
%     it does NOT certify a universally impossible diagnosis.
%   - No Simulink, solver, hardware, or thermal residual calibration.
%
% Save next to run_e205a_online_three_state.m in project ROOT.
% Run: run_e205b_paired_policy_validation
% Output: evidence_e205b_paired_policy_validation[_run_TIMESTAMP]
% ----------------------------------------------------------------------
clc;close all; timer=tic;
root=fileparts(mfilename('fullpath'));
assert(isfolder(root),'Place this .m file in the E11 project root.');
srcA=fullfile(root,'evidence_e205a_online_three_state','E205A_results.mat');
srcV4=fullfile(root,'evidence_e204_finite_time_boundary_v4','E204_v4_results.mat');
assert(isfile(srcA),'Missing E20.5A MAT: %s',srcA);
assert(isfile(srcV4),'Missing E20.4 v4 MAT: %s',srcV4);
x=load(srcA,'Result');y=load(srcV4,'S');
assert(isfield(x,'Result') && isfield(y,'S'),'Expected Result and S structs.');
A=x.Result;V=y.S; P=A.P; D=A.D;
assert(abs(P.R0-5e-4)<1e-12 && abs(P.Rupper-8e-4)<1e-12,...
    'Archived resistance parameters unexpected.');
assert(abs(double(A.threshold)-double(V.threshold))<1e-12,...
    'The archived sequential thresholds differ.');
assert(P.Ihigh<=P.Icap+1e-12 && P.Ibase>0 && P.Ihigh>P.Ibase);
assert(P.Nswitch==10 && P.sigmaRaw>0 && P.sigmaShared>=0);
assert(P.TempCap>0 && P.Rth>0 && P.Cth>0 && P.alpha>=0);
assert(D.forecastSteps>=1 && D.minProjectedGainZ>0);
assert(D.designDeltaR_mOhm>0 && D.triggerLowerFraction>=0);
assert(numel(V.C.currentPolicies)==5 && numel(D.DeltaR_mOhm)==30);
assert(isfield(A,'noTrueDeltaRInController') && A.noTrueDeltaRInController);

% Predeclare a small but scientifically relevant, never retuned grid.
C=struct();
C.EnergyBudget_J=[1.0;2.0;3.0];
C.Time_s=[1.0 1.5 2.0 3.0 4.0 5.0];
C.DeltaR_mOhm=[0.05;0.10;0.12;0.14;0.15;0.16;0.18;0.20;0.22;0.25];
C.policyNames={'Fixed low','Fixed high','Archived feedback','E20.5A online'};
C.nHealthy=16000;
C.nFault=8000;
C.seed=2050201;             % Fresh data; NOT E20.5A or E20.4 seeds.
C.alphaFamily=0.05; % Global familywise level across all reported inferential claims.
C.alphaCPGroup=0.025;
C.alphaMcNemarGroup=0.0125;
C.alphaDeltaGroup=0.0125;
C.Nsteps=round(C.Time_s/P.Ts);
assert(all(abs(C.Nsteps*P.Ts-C.Time_s)<1e-10));
assert(C.Nsteps(end)<=P.Nmax);
assert(max(C.DeltaR_mOhm)*1e-3+P.R0<=P.Rupper+1e-12);
for jj=1:numel(C.DeltaR_mOhm)
    assert(any(abs(double(D.DeltaR_mOhm(:))-C.DeltaR_mOhm(jj))<1e-9),...
        'A selected dR value is outside the archived E20.5A grid.');
end
assert(~ismember(C.seed,[D.seed,V.C.seed]));
Pthr=double(A.threshold);
nE=numel(C.EnergyBudget_J);nT=numel(C.Time_s);
nR=numel(C.DeltaR_mOhm);nP=numel(C.policyNames);
C.nCPClaims=nP*nE*nT*(1+2*nR); % healthy FAR upper; fault PD lower/upper
C.alphaCP=C.alphaCPGroup/C.nCPClaims;
C.nPairedClaims=(nP-1)*nE*nT*nR;
C.alphaPair=C.alphaMcNemarGroup/C.nPairedClaims;
% Global familywise 0.05: disjoint allocations 0.025 CP + 0.0125
% McNemar + 0.0125 paired Hoeffding bounds (union bound).
C.alphaScope=['Global FWER <= 0.05 via three allocated families: ',...
    'CP 0.025; paired exact McNemar 0.0125; paired effect CI 0.0125.'];
C.nHealthyPerCondition=C.nHealthy;
C.nFaultPerCondition=C.nFault;

outDir=fullfile(root,'evidence_e205b_paired_policy_validation');
if isfolder(outDir)
    listing=dir(outDir);
    listing=listing(~ismember({listing.name},{'.','..'}));
    if ~isempty(listing)
        outDir=fullfile(root,['evidence_e205b_paired_policy_validation_run_',...
            datestr(now,'yyyymmdd_HHMMSSFFF')]);
    end
end
if ~isfolder(outDir)
    mkdir(outDir);
else
    remaining=dir(outDir);
    remaining=remaining(~ismember({remaining.name},{'.','..'}));
    assert(isempty(remaining),'Output folder collision.');
end
fprintf('\nE20.5B INDEPENDENT PAIRED VALIDATION\n');
fprintf('Inputs: %s\n        %s\n',srcA,srcV4);
fprintf('New seed: %d | healthy: %d | per fault: %d\n',...
    C.seed,C.nHealthy,C.nFault);
fprintf('Budgets %d, horizons %d, selected dR %d, policies %d\n',...
    nE,nT,nR,nP);
fprintf('CP claims %d; paired claims %d\n',...
    C.nCPClaims,C.nPairedClaims);
fprintf('E20.5A controller and E20.3B feedback are frozen; no tuning.\n');

rng(C.seed,'twister');initialRNG=rng;
% Four policies observe the SAME Gaussian raw noise and SAME shared bias
% in each healthy or fault condition. The different current selections
% produce different measurements/trajectories (as they should).
farCount=zeros(nE,nT,nP);
meanHenergy=zeros(nE,nT,nP);meanHtemp=zeros(nE,nT,nP);
healthyViolation=zeros(nE,nP);healthyPairDiscord=zeros(nE,nT,nP-1,2);
[B,N]=localNoise(C.nHealthy,C.Nsteps(end),P);
for ie=1:nE
    Q=cell(1,nP);
    for ip=1:nP
        Q{ip}=localBatch(0,ip,B,N,C.EnergyBudget_J(ie),P,D,Pthr,C.Nsteps);
        healthyViolation(ie,ip)=nnz(Q{ip}.violation);
    end
    for it=1:nT
        h=C.Nsteps(it);
        isNew=Q{4}.alarmK<=h;
        for ip=1:nP
            alarms=Q{ip}.alarmK<=h;
            farCount(ie,it,ip)=nnz(alarms);
            meanHenergy(ie,it,ip)=mean(Q{ip}.energy(:,it));
            meanHtemp(ie,it,ip)=mean(Q{ip}.temp(:,it));
            if ip<4
                healthyPairDiscord(ie,it,ip,1)=nnz(isNew & ~alarms);
                healthyPairDiscord(ie,it,ip,2)=nnz(~isNew & alarms);
            end
        end
    end
end
clear B N Q;
assert(~any(healthyViolation(:)),'Healthy-path physical constraint violation.');

nCount=zeros(nE,nT,nR,nP);
meanEnergy=zeros(nE,nT,nR,nP);meanTemp=meanEnergy;
meanN=meanEnergy;boostRate=meanEnergy;meanBoost=meanEnergy;
resourceRate=meanEnergy;deadlineRate=meanEnergy;
physicalViolations=zeros(nE,nR,nP);
newOnly=zeros(nE,nT,nR,nP-1);
baseOnly=zeros(nE,nT,nR,nP-1);
pairedDeltaEnergy=zeros(nE,nT,nR,nP-1);
pairedEnergySE=zeros(nE,nT,nR,nP-1);
pairedDeltaObs=zeros(nE,nT,nR,nP-1);
for ir=1:nR
    dR=C.DeltaR_mOhm(ir);
    fprintf('Paired fault: dR %.3f mOhm (%d/%d)\n',dR,ir,nR);
    [B,N]=localNoise(C.nFault,C.Nsteps(end),P);
    for ie=1:nE
        Q=cell(1,nP);
        for ip=1:nP
            Q{ip}=localBatch(dR,ip,B,N,C.EnergyBudget_J(ie),P,D,Pthr,C.Nsteps);
            physicalViolations(ie,ir,ip)=nnz(Q{ip}.violation);
        end
        for it=1:nT
            h=C.Nsteps(it);
            newDetected=Q{4}.alarmK<=h;
            for ip=1:nP
                baseDetected=Q{ip}.alarmK<=h;
                nCount(ie,it,ir,ip)=nnz(baseDetected);
                meanEnergy(ie,it,ir,ip)=mean(Q{ip}.energy(:,it));
                meanTemp(ie,it,ir,ip)=mean(Q{ip}.temp(:,it));
                meanN(ie,it,ir,ip)=mean(Q{ip}.observations(:,it));
                boostRate(ie,it,ir,ip)=mean(Q{ip}.boostSteps(:,it)>0);
                meanBoost(ie,it,ir,ip)=mean(Q{ip}.boostSteps(:,it));
                stopped=Q{ip}.resourceStopK<=h & ~baseDetected;
                resourceRate(ie,it,ir,ip)=mean(stopped);
                deadlineRate(ie,it,ir,ip)=mean(~baseDetected & ~stopped);
                assert(abs(nCount(ie,it,ir,ip)/C.nFault+...
                    resourceRate(ie,it,ir,ip)+...
                    deadlineRate(ie,it,ir,ip)-1)<2e-12,...
                    'Terminal state partition is inconsistent.');
                if ip<4
                    newOnly(ie,it,ir,ip)=nnz(newDetected & ~baseDetected);
                    baseOnly(ie,it,ir,ip)=nnz(~newDetected & baseDetected);
                    eDiff=Q{4}.energy(:,it)-Q{ip}.energy(:,it);
                    pairedDeltaEnergy(ie,it,ir,ip)=mean(eDiff);
                    pairedEnergySE(ie,it,ir,ip)=std(eDiff,0)/sqrt(C.nFault);
                    pairedDeltaObs(ie,it,ir,ip)=mean( ...
                        Q{4}.observations(:,it)-Q{ip}.observations(:,it));
                end
            end
        end
    end
end
assert(~any(physicalViolations(:)),'Fault-path physical constraint violation.');
assert(all(diff(nCount,1,2)>=0,'all'),...
    'Alarms must be non-decreasing with deadline.');
for ie=1:nE
    e=reshape(meanEnergy(ie,:,:,:),[],1);
    assert(all(e<=C.EnergyBudget_J(ie)+1e-8));
end
assert(all(meanTemp(:)<=P.TempCap+1e-8));

% CP confidence bounds and new joint certificates (common declared family).
FAR=farCount/C.nHealthy;
FARup=zeros(size(FAR));
for j=1:numel(FARup)
    FARup(j)=localCPUpper(farCount(j),C.nHealthy,C.alphaCP);
end
PDEst=nCount/C.nFault;
PDlo=zeros(size(PDEst));PDhi=zeros(size(PDEst));
for j=1:numel(PDEst)
    PDlo(j)=localCPLower(nCount(j),C.nFault,C.alphaCP);
    PDhi(j)=localCPUpper(nCount(j),C.nFault,C.alphaCP);
end
Cert=false(size(PDEst));
for ip=1:nP
    for ir=1:nR
        Cert(:,:,ir,ip)=PDlo(:,:,ir,ip)>=P.PDtarget &...
            FARup(:,:,ip)<=P.FARtarget;
    end
end
Bound=nan(nE,nT,nP);
for ie=1:nE
    for it=1:nT
        for ip=1:nP
            ix=find(reshape(Cert(ie,it,:,ip),[],1),1,'first');
            if ~isempty(ix),Bound(ie,it,ip)=C.DeltaR_mOhm(ix);end
        end
    end
end

% Matched-pair estimands: DeltaPD = PD(new online) - PD(baseline).
% Exact McNemar: conditional, two-sided under equal marginal PD.
% A simultaneous magnitude bound also follows from Hoeffding since each
% paired difference is in [-1,1] (often conservative).
delta=zeros(nE,nT,nR,nP-1);pM=ones(size(delta));
loDelta=zeros(size(delta));hiDelta=zeros(size(delta));
halfwidth=sqrt(2*log(2*C.nPairedClaims/C.alphaDeltaGroup)/C.nFault);
for ip=1:nP-1
    more=newOnly(:,:,:,ip);less=baseOnly(:,:,:,ip);
    delta(:,:,:,ip)=(more-less)/C.nFault;
    loDelta(:,:,:,ip)=max(-1,delta(:,:,:,ip)-halfwidth);
    hiDelta(:,:,:,ip)=min(1,delta(:,:,:,ip)+halfwidth);
    for j=1:numel(more)
        pM(j+(ip-1)*numel(more))=localMcNemar(more(j),less(j));
    end
end
sigBetter=pM<=C.alphaPair & delta>0;
sigWorse=pM<=C.alphaPair & delta<0;

% Tables with all important denominators and column definitions.
HealthyTable=table();FaultTable=table();PairedTable=table();BoundaryTable=table();
for ie=1:nE
    for it=1:nT
        for ip=1:nP
            ht=table(C.EnergyBudget_J(ie),C.Time_s(it),string(C.policyNames{ip}),...
                C.nHealthy,farCount(ie,it,ip),FAR(ie,it,ip),...
                FARup(ie,it,ip),FARup(ie,it,ip)<=P.FARtarget,...
                meanHenergy(ie,it,ip),meanHtemp(ie,it,ip),...
                'VariableNames',{'EnergyBudget_J','Horizon_s','Policy',...
                'Trials','FalseAlarms','EmpiricalFAR','SimultaneousFARUpper',...
                'FARCertified','MeanHealthyEnergy_J','MeanHealthyRise_K'});
            HealthyTable=[HealthyTable;ht]; %#ok<AGROW>
            Bt=table(C.EnergyBudget_J(ie),C.Time_s(it),...
                string(C.policyNames{ip}),Bound(ie,it,ip),...
                string(localStatus(Bound(ie,it,ip))),...
                'VariableNames',{'EnergyBudget_J','Horizon_s','Policy',...
                'FirstCertifiedDeltaR_mOhm','Status'});
            BoundaryTable=[BoundaryTable;Bt]; %#ok<AGROW>
            for ir=1:nR
                ft=table(C.EnergyBudget_J(ie),C.Time_s(it),...
                    C.DeltaR_mOhm(ir),string(C.policyNames{ip}),...
                    C.nFault,nCount(ie,it,ir,ip),PDEst(ie,it,ir,ip),...
                    PDlo(ie,it,ir,ip),PDhi(ie,it,ir,ip),...
                    FARup(ie,it,ip),Cert(ie,it,ir,ip),...
                    meanEnergy(ie,it,ir,ip),meanTemp(ie,it,ir,ip),...
                    meanN(ie,it,ir,ip),boostRate(ie,it,ir,ip),...
                    meanBoost(ie,it,ir,ip),resourceRate(ie,it,ir,ip),...
                    deadlineRate(ie,it,ir,ip),...
                    'VariableNames',{'EnergyBudget_J','Horizon_s',...
                    'TrueDeltaR_AnalysisOnly_mOhm','Policy','Trials',...
                    'Detections','EmpiricalPD','SimultaneousPDLower',...
                    'SimultaneousPDUpper','SimultaneousFARUpper',...
                    'JointCertified','MeanContactEnergy_J','MeanTempRise_K',...
                    'MeanObservedSamples','HighCurrentTrialRate','MeanHighCurrentSteps',...
                    'ResourceUnavailableRate','DeadlineUnavailableRate'});
                FaultTable=[FaultTable;ft]; %#ok<AGROW>
            end
        end
        for ip=1:nP-1
            for ir=1:nR
                pt=table(C.EnergyBudget_J(ie),C.Time_s(it),...
                    C.DeltaR_mOhm(ir),string(C.policyNames{ip}),...
                    string(C.policyNames{4}),C.nFault,...
                    newOnly(ie,it,ir,ip),baseOnly(ie,it,ir,ip),...
                    delta(ie,it,ir,ip),loDelta(ie,it,ir,ip),...
                    hiDelta(ie,it,ir,ip),pM(ie,it,ir,ip),...
                    sigBetter(ie,it,ir,ip),sigWorse(ie,it,ir,ip),...
                    pairedDeltaEnergy(ie,it,ir,ip),...
                    pairedEnergySE(ie,it,ir,ip),...
                    pairedDeltaObs(ie,it,ir,ip),...
                    'VariableNames',{'EnergyBudget_J','Horizon_s',...
                    'TrueDeltaR_AnalysisOnly_mOhm','BaselinePolicy','NewPolicy',...
                    'PairedTrials','NewOnlyDetections','BaselineOnlyDetections',...
                    'DeltaPD_NewMinusBaseline','SimultaneousDeltaLower',...
                    'SimultaneousDeltaUpper','McNemarExactTwoSidedP',...
                    'NewBetterBonferroni','NewWorseBonferroni',...
                    'MeanDeltaEnergy_J_NewMinusBaseline',...
                    'DescriptiveSE_DeltaEnergy_J',...
                    'MeanDeltaObservedSamples_NewMinusBaseline'});
                PairedTable=[PairedTable;pt]; %#ok<AGROW>
            end
        end
    end
end

writetable(HealthyTable,fullfile(outDir,'E205B_healthy.csv'));
writetable(FaultTable,fullfile(outDir,'E205B_fault.csv'));
writetable(PairedTable,fullfile(outDir,'E205B_paired.csv'));
writetable(BoundaryTable,fullfile(outDir,'E205B_policy_boundary.csv'));

Results=struct();
Results.C=C;Results.P=P;Results.D=D;Results.threshold=Pthr;
Results.initialRNG=initialRNG;
Results.healthyTable=HealthyTable;Results.faultTable=FaultTable;
Results.pairedTable=PairedTable;Results.boundaryTable=BoundaryTable;
Results.PD=PDEst;Results.PDLower=PDlo;Results.PDUpper=PDhi;
Results.FAR=FAR;Results.FARUpper=FARup;Results.JointCertified=Cert;
Results.Boundary=Bound;Results.MeanEnergy_J=meanEnergy;
Results.MeanTempRise_K=meanTemp;Results.MeanObservedSamples=meanN;
Results.BoostRate=boostRate;Results.MeanBoostSteps=meanBoost;
Results.UnavailableResourceRate=resourceRate;
Results.UnavailableDeadlineRate=deadlineRate;
Results.NewOnlyDetections=newOnly;Results.BaselineOnlyDetections=baseOnly;
Results.PairedDeltaPD=delta;Results.PairedMcNemarP=pM;
Results.PairedMeanDeltaEnergy_J=pairedDeltaEnergy;
Results.PairedDescriptiveSEEnergy_J=pairedEnergySE;
Results.PairedMeanDeltaObservedSamples=pairedDeltaObs;
Results.PairedDeltaSimultaneousLower=loDelta;
Results.PairedDeltaSimultaneousUpper=hiDelta;
Results.PathViolationFault=physicalViolations;
Results.PathViolationHealthy=healthyViolation;
Results.scientificScope='Four frozen policies on a reduced, matched-noise grid, within original residual/thermal model';
Results.notSimulinkValidated=true;Results.notGloballyOptimal=true;
Results.sourceE205A=srcA;Results.sourceE204V4=srcV4;
Results.runtimeSeconds=toc(timer);
try,Results.matlabVersion=version;catch,Results.matlabVersion='unknown';end
save(fullfile(outDir,'E205B_results.mat'),'Results','-v7.3');
copyfile([mfilename('fullpath') '.m'],...
    fullfile(outDir,'run_e205b_paired_policy_validation_source_snapshot.m'));

localFigures(outDir,C,P,PDEst,PDlo,delta,loDelta,hiDelta,...
    meanEnergy,pairedDeltaEnergy,Bound);
reportPath=fullfile(outDir,'E205B_report.txt');
fid=fopen(reportPath,'w');assert(fid>0,'Report cannot be created.');
cleanupObj=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'E20.5B FROZEN POLICY, PAIRED-NOISE VALIDATION\n');
fprintf(fid,'==========================================\n');
fprintf(fid,'Run time %.2f s | random seed %d\n',toc(timer),C.seed);
fprintf(fid,'New simulation: healthy %d, faults per dR %d, dR points %d\n',...
    C.nHealthy,C.nFault,nR);
fprintf(fid,'CP claims %d: per-claim alpha %.9g\n',C.nCPClaims,C.alphaCP);
fprintf(fid,'Paired hypotheses %d: Bonferroni alpha %.9g\n',...
    C.nPairedClaims,C.alphaPair);
fprintf(fid,'Global FWER <=0.05: CP allocation %.4f, McNemar allocation %.4f, paired CI allocation %.4f.\n',...
    C.alphaCPGroup,C.alphaMcNemarGroup,C.alphaDeltaGroup);
fprintf(fid,'All simulated paths passed physical caps: %d healthy, %d fault violations.\n',...
    nnz(healthyViolation),nnz(physicalViolations));
fprintf(fid,'FAR certified %d/%d policy/budget/time conditions.\n',...
    nnz(FARup<=P.FARtarget),numel(FARup));
fprintf(fid,'Joint certified fault-policy conditions %d/%d.\n',...
    nnz(Cert),numel(Cert));
fprintf(fid,'New better (paired exact, Bonferroni) %d/%d; worse %d/%d.\n',...
    nnz(sigBetter),numel(sigBetter),nnz(sigWorse),numel(sigWorse));
fprintf(fid,'Paired Hoeffding half-width %.6g (conservative)\n',halfwidth);
fprintf(fid,'\nPOLICY BOUNDARIES AT E=3.0 J (subset grid only):\n');
ie=find(abs(C.EnergyBudget_J-3)<1e-12,1);
for it=1:nT
    fprintf(fid,'T=%.1f s | ',C.Time_s(it));
    for ip=1:nP
        fprintf(fid,'%s %.3f | ',C.policyNames{ip},Bound(ie,it,ip));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\nPRIMARY PREDECLARED CONDITION E=3J T=3s dR=0.15mOhm\n');
it=find(abs(C.Time_s-3)<1e-12,1);
ir=find(abs(C.DeltaR_mOhm-0.15)<1e-12,1);
assert(~isempty(ir));
for ip=1:nP
    fprintf(fid,'%s: PD %.6f, lower %.6f, energy %.5f J, rise %.6f K, obs %.3f samples\n',...
        C.policyNames{ip},PDEst(ie,it,ir,ip),PDlo(ie,it,ir,ip),...
        meanEnergy(ie,it,ir,ip),meanTemp(ie,it,ir,ip),...
        meanN(ie,it,ir,ip));
end
for ip=1:nP-1
    fprintf(fid,'New minus %s: DeltaPD %.6f, NewOnly %d, BaseOnly %d, exact P %.9g, better=%d, worse=%d; paired DeltaEnergy %.6f J (descriptive SE %.6f)\n',...
        C.policyNames{ip},delta(ie,it,ir,ip),...
        newOnly(ie,it,ir,ip),baseOnly(ie,it,ir,ip),...
        pM(ie,it,ir,ip),sigBetter(ie,it,ir,ip),sigWorse(ie,it,ir,ip),...
        pairedDeltaEnergy(ie,it,ir,ip),pairedEnergySE(ie,it,ir,ip));
end
fprintf(fid,'\nRESTRICTIONS\n');
fprintf(fid,'* Matching random noise across policies allows paired tests, not shared controller information.\n');
fprintf(fid,'* No tuning based on the present test noise or observed results.\n');
fprintf(fid,'* First certified dR uses a sparse, predeclared grid; do not compare numerically as a continuous optimum.\n');
fprintf(fid,'* Mean resource results do not prove pathwise safety on their own; simulation asserts path checks.\n');
fprintf(fid,'* High-current usage means commanded Ihigh; fixed-high is always high, not an online boost event.\n');
fprintf(fid,'* Paired McNemar tests compare detection events by each fixed deadline. Energy differences/SE are DESCRIPTIVE, not simultaneous certificates.\n');
fprintf(fid,'* Both Gaussian residual mean and thermal guard retain historical model assumptions.\n');
fprintf(fid,'* Unavailable means unresolved, NOT no fault or proven impossibility.\n');
fprintf(fid,'* No Simulink, solver convergence, hardware, or temperature-dependent measurement validation.\n');

fprintf('\n==================================================\n');
fprintf('E20.5B PAIRED POLICY VALIDATION COMPLETE (MODEL ONLY)\n');
fprintf('Healthy FAR certified: %d/%d\n',nnz(FARup<=P.FARtarget),numel(FARup));
fprintf('Joint fault certified: %d/%d\n',nnz(Cert),numel(Cert));
fprintf('Online better / worse paired tests: %d / %d (of %d)\n',...
    nnz(sigBetter),nnz(sigWorse),numel(sigBetter));
fprintf('Pathwise constraint violations: %d\n',...
    nnz(physicalViolations)+nnz(healthyViolation));
fprintf('Primary E=3 J, T=3 s, dR=.15 mOhm: new PD %.4f vs fixed-high %.4f\n',...
    PDEst(ie,it,ir,4),PDEst(ie,it,ir,2));
fprintf('Elapsed %.1f s\nSaved: %s\n',toc(timer),outDir);
fprintf('==================================================\n');
end

function [B,N]=localNoise(nTrials,nSteps,P)
B=P.sigmaShared*randn(nTrials,1);
N=P.sigmaRaw*randn(nTrials,nSteps);
end

function Q=localBatch(dR,mode,B,Noise,Ecap,P,D,thr,checkSteps)
% Freeze mode 3 to the original E20.3B/v4 latched feedback;
% freeze mode 4 to the E20.5A score/resource online controller.
% Only the PLANT uses dR (measurement mean and physical resistance).
m=numel(B);nMax=checkSteps(end);nT=numel(checkSteps);
assert(iscolumn(B) && isequal(size(Noise),[m,nMax]));
Rplant=P.R0+dR*1e-3;gap=max(dR-P.robustOffset,0);
assert(Rplant>0 && Rplant<=P.Rupper+1e-12);
Ttrue=zeros(m,1);Etrue=zeros(m,1);
Tguard=zeros(m,1);Eguard=zeros(m,1);
sumY=zeros(m,1);score=zeros(m,1);
active=true(m,1);latched=false(m,1);
alarmK=inf(m,1);resourceStopK=inf(m,1);
obs=zeros(m,1);nBoost=zeros(m,1);violation=false(m,1);
Ecp=zeros(m,nT);Tcp=zeros(m,nT);Ncp=zeros(m,nT);Bcp=zeros(m,nT);
for k=1:nMax
    if any(active)
        current=P.Ibase*ones(m,1);
        [tLow,eLow]=localStep(Tguard,Eguard,P.Ibase,P.Rupper,P);
        [tHigh,eHigh]=localStep(Tguard,Eguard,P.Ihigh,P.Rupper,P);
        canLow=tLow<=P.TempCap+1e-12 & eLow<=Ecap+1e-12;
        canHigh=tHigh<=P.TempCap+1e-12 & eHigh<=Ecap+1e-12;
        if mode==2
            current(active)=P.Ihigh;
        elseif mode==3
            if k==P.Nswitch+1
                latched(active)=score(active)<P.switchFraction*thr;
            end
            % Fall back from high to low PERMANENTLY if guard blocks high.
            blocked=active & latched & ~canHigh;
            latched(blocked)=false;
            current(active & latched)=P.Ihigh;
        elseif mode==4 && k>P.Nswitch
            look=min(D.forecastSteps,nMax-k+1);
            nextN=k-1+look;
            denom=sqrt(nextN*P.sigmaRaw^2+nextN^2*P.sigmaShared^2);
            designGap=max(D.designDeltaR_mOhm-P.robustOffset,0);
            predLow=(sumY+look*designGap*P.Ibase)/denom;
            predHigh=(sumY+look*designGap*P.Ihigh)/denom;
            boost=active & canHigh & ...
                score>=D.triggerLowerFraction*thr & score<thr & ...
                predLow<thr & (predHigh-predLow)>=D.minProjectedGainZ;
            current(boost)=P.Ihigh;
        end
        useHigh=current>P.Ibase+1e-10;
        admitted=active & canLow;
        if mode==2
            admitted=active & canHigh;
        end
        resourceNow=active & ~admitted;
        resourceStopK(resourceNow)=k;
        active(resourceNow)=false;
        if any(admitted)
            Tnext=tLow;Enext=eLow;
            chosenHigh=admitted & useHigh;
            Tnext(chosenHigh)=tHigh(chosenHigh);
            Enext(chosenHigh)=eHigh(chosenHigh);
            assert(all(Tnext(admitted)<=P.TempCap+1e-10 &...
                Enext(admitted)<=Ecap+1e-10));
            Tguard(admitted)=Tnext(admitted);
            Eguard(admitted)=Enext(admitted);
            [tp,ep]=localStep(Ttrue(admitted),Etrue(admitted),...
                current(admitted),Rplant,P);
            Ttrue(admitted)=tp;Etrue(admitted)=ep;
            violation(admitted)=violation(admitted) | ...
                tp>P.TempCap+1e-9 | ep>Ecap+1e-9 | ...
                tp<-1e-9 | ep<-1e-9;
            obs(admitted)=obs(admitted)+1;
            nBoost(chosenHigh)=nBoost(chosenHigh)+1;
            measurement=gap*current(admitted)+B(admitted)+Noise(admitted,k);
            sumY(admitted)=sumY(admitted)+measurement;
            den=sqrt(k*P.sigmaRaw^2+k^2*P.sigmaShared^2);
            score(admitted)=sumY(admitted)/den;
            detected=admitted & score>=thr;
            alarmK(detected)=k;active(detected)=false;
        end
    end
    ix=find(checkSteps==k,1);
    if ~isempty(ix)
        Ecp(:,ix)=Etrue;Tcp(:,ix)=Ttrue;
        Ncp(:,ix)=obs;Bcp(:,ix)=nBoost;
    end
end
Q.alarmK=alarmK;Q.resourceStopK=resourceStopK;
Q.energy=Ecp;Q.temp=Tcp;Q.observations=Ncp;
Q.boostSteps=Bcp;Q.violation=violation;
end

function [Tnew,Enew]=localStep(T,E,I,R,P)
tau=P.Rth*P.Cth;
gain=I.^2.*R.*P.alpha.*P.Rth;
assert(all(gain(:)<1),'Electrothermal instability; parameters are invalid.');
rate=(1-gain)./tau;
Tinf=P.Rth.*I.^2.*R./(1-gain);
Tnew=Tinf+(T-Tinf).*exp(-rate*P.Ts);
Tint=Tinf*P.Ts+(T-Tinf).*(-expm1(-rate*P.Ts))./rate;
Enew=E+I.^2.*R.*(P.Ts+P.alpha.*Tint);
end

function lo=localCPLower(x,n,a)
if x==0,lo=0;else,lo=betaincinv(a,x,n-x+1);end
end
function hi=localCPUpper(x,n,a)
if x==n,hi=1;else,hi=betaincinv(1-a,x+1,n-x);end
end
function p=localMcNemar(n10,n01)
% 2-sided exact conditional binomial test; k=min discordant counts.
d=n10+n01;
if d==0,p=1;return;end
k=min(n10,n01);
if k==d,p=1;return;end
p=min(1,2*betainc(0.5,d-k,k+1));
end
function s=localStatus(v)
if isfinite(v),s='CERTIFIED_IN_SELECTED_GRID';
else,s='NOT_CERTIFIED_NOT_PROVEN_INFEASIBLE';end
end

function localFigures(outDir,C,P,PD,PL,delta,lo,hi,energy,dE,boundary)
[~,ie]=min(abs(C.EnergyBudget_J-3));
[~,it]=min(abs(C.Time_s-3));
R=C.DeltaR_mOhm;T=C.Time_s;nP=numel(C.policyNames);
colors=[.23 .46 .69;.83 .39 .22;.55 .40 .65;.19 .52 .38];
f=figure('Color','w','Position',[90 90 850 485]);ax=axes('Parent',f);hold(ax,'on');
for p=1:nP
    plot(ax,R,reshape(PD(ie,it,:,p),[],1),'-o','Color',colors(p,:),...
        'LineWidth',1.6,'MarkerSize',4,'DisplayName',C.policyNames{p});
end
plot(ax,R,reshape(PL(ie,it,:,4),[],1),'--',...
    'Color',colors(4,:),'LineWidth',1.35,...
    'DisplayName','Online: simultaneous PD lower');
plot(ax,[R(1) R(end)],[P.PDtarget P.PDtarget],':',...
    'Color',[.25 .25 .25],'LineWidth',1.3,'DisplayName','PD target');
xlabel(ax,'True degradation (m\Omega) [analysis only]');
ylabel(ax,'Detection probability / lower confidence bound');
title(ax,'Paired-noise benchmark, Emax=3 J, T=3 s');
legend(ax,'Location','southeast','Box','off');ylim(ax,[0 1.025]);
localStyle(ax);localExport(f,outDir,'E205B_paired_PD');close(f);

f=figure('Color','w','Position',[90 90 850 485]);ax=axes('Parent',f);hold(ax,'on');
for p=1:3
    plot(ax,R,reshape(delta(ie,it,:,p),[],1),'-o',...
        'Color',colors(p,:),'LineWidth',1.7,'MarkerSize',5,...
        'DisplayName',['Online minus ' C.policyNames{p}]);
end
plot(ax,[R(1) R(end)],[0 0],':','Color',[.2 .2 .2],...
    'LineWidth',1.15,'DisplayName','Zero difference');
xlabel(ax,'True degradation (m\Omega) [analysis only]');
ylabel(ax,'Paired detection-rate difference');
title(ax,'Paired effect at Emax=3 J, T=3 s');
legend(ax,'Location','best','Box','off');localStyle(ax);
localExport(f,outDir,'E205B_paired_deltaPD');close(f);

f=figure('Color','w','Position',[90 90 810 485]);ax=axes('Parent',f);hold(ax,'on');
for p=1:nP
    plot(ax,R,reshape(energy(ie,it,:,p),[],1),'-o',...
        'Color',colors(p,:),'LineWidth',1.55,'MarkerSize',4,...
        'DisplayName',C.policyNames{p});
end
xlabel(ax,'True degradation (m\Omega) [analysis only]');
ylabel(ax,'Mean integrated contact energy (J)');
title(ax,'Energy usage including early detection, Emax=3 J, T=3 s');
legend(ax,'Location','best','Box','off');localStyle(ax);
localExport(f,outDir,'E205B_energy_comparison');close(f);

% Descriptive paired energy differences: positive = new uses more energy.
f=figure('Color','w','Position',[90 90 820 470]);ax=axes('Parent',f);hold(ax,'on');
for p=1:3
    plot(ax,R,reshape(dE(ie,it,:,p),[],1),'-o',...
        'Color',colors(p,:),'LineWidth',1.6,'MarkerSize',4,...
        'DisplayName',['Online minus ' C.policyNames{p}]);
end
plot(ax,[R(1) R(end)],[0 0],':','Color',[.25 .25 .25],...
    'DisplayName','Zero difference');
xlabel(ax,'True degradation (m\Omega) [analysis only]');
ylabel(ax,'Paired mean contact-energy difference (J)');
title(ax,'Descriptive paired energy contrast, Emax=3 J, T=3 s');
legend(ax,'Location','best','Box','off');localStyle(ax);
localExport(f,outDir,'E205B_paired_energy_difference');close(f);

f=figure('Color','w','Position',[90 90 810 485]);ax=axes('Parent',f);hold(ax,'on');
styles={'-o','--s','-.^','-v'};
for p=1:nP
    plot(ax,T,reshape(boundary(ie,:,p),[],1),styles{p},...
        'Color',colors(p,:),'LineWidth',1.7,'MarkerSize',6,...
        'DisplayName',C.policyNames{p});
end
xlabel(ax,'Requested observation horizon (s)');
ylabel(ax,'First certified degradation (m\Omega)');
title(ax,'Selected-grid policy boundaries, Emax=3 J');
legend(ax,'Location','bestoutside','Box','off');localStyle(ax);
localExport(f,outDir,'E205B_policy_boundaries');close(f);
end
function localStyle(ax)
set(ax,'FontName','Times New Roman','FontSize',11,...
    'LineWidth',0.9,'TickDir','out','Box','off');
grid(ax,'on');
try,ax.GridAlpha=0.13;catch,end
try,ax.Toolbar.Visible='off';catch,end
end
function localExport(f,outDir,name)
png=fullfile(outDir,[name '.png']);
pdf=fullfile(outDir,[name '.pdf']);
try
    exportgraphics(f,png,'Resolution',300);
    exportgraphics(f,pdf,'ContentType','vector');
catch ME
    warning('Figure export fallback: %s',ME.message);
    print(f,png,'-dpng','-r300');print(f,pdf,'-dpdf','-painters');
end
end
