%% Poldrack

%%% 1st Level Categorical Design for all sessions/subjects
clear;
spm('defaults','FMRI')

global defaults

mainpath = dir('Data path');  % Working directory
                   
TR = 1.98;

for Sub = 1:20    % number of subjects
    behavdir = dir(strcat('Onset files','\*.mat'));

    for RunNum = 1:5          % Session ID: 1->3   
        
            
      for trial = 1:120 % 60
        SPM.nscan          = 206;
        SPM.xBF.name = 'hrf';       % name of basis function
        SPM.xBF.length = 32.0513;   % length in seconds of basis
        SPM.xBF.order = 1;          % order of basis set
        SPM.xBF.T = 32;             % number of subdivisions of TR
        SPM.xBF.T0 = 16;            % first time bin (see slice timing)
        SPM.xBF.UNITS = 'scans';    % options: 'scans'|'secs' for onsets
        SPM.xBF.Volterra = 1;       % order of convolution
        currDir = ['Data path\',mainpath(Sub+2).name];
        cd(currDir);
        if RunNum==1
        outputdir = (strcat('Output path',num2str(Sub)));
        mkdir(strcat('Output path',num2str(Sub)));
        end
        
        onsetDir = ('Onset files');
        cd(onsetDir);
        load(behavdir((Sub-1)*5+RunNum).name);    % Load the .mat file for stimulus onset timings
        dir{RunNum} = strcat(currDir,'\run',num2str(RunNum));
        tmp{RunNum} = spm_select('fplist',dir{RunNum},'^swars.*\.img');
        trialdir = strcat('Output path',num2str(Sub),'\Trial',num2str((RunNum-1)*120+trial));
        mkdir(strcat('Output path',num2str(Sub),'\Trial',num2str((RunNum-1)*120+trial)));
        

        Beta(1:20)=Onset(1:20,2);
        Beta(21:40)=Onset(1:20,3);
        Beta(41:60)=Onset(1:20,4);
        
%         Beta2 = [Beta; Beta + 0.7576];  % split each trial to 2
%         Beta2 = reshape(Beta2, 1, []);
        
        DurTime(1:20)=1.5152;
        DurTime(21:40)=1.5152;
        DurTime(41:60)=1.5152;

%         DurTime(1:120)=0.7576; % 0.75TR 

        
        BetaEx=Beta2;
        BetaEx(trial)=[];
        DurEx=DurTime;
        DurEx(trial)=[];
         
        for c = 1:1
            SPM.Sess(1).U(c).name = {strcat('trial',num2str(trial))};   
            SPM.Sess(1).U(c).ons = Beta2(trial)';
            SPM.Sess(1).U(c).dur = DurTime(trial);
            SPM.Sess(1).U(c).P(1).name = 'none';       % Parametric Modulation; 'none' for now
        end
        
         for c = 2:2         
            SPM.Sess(1).U(c).name = {'Other'};    
            SPM.Sess(1).U(c).ons = BetaEx;
            SPM.Sess(1).U(c).dur =DurEx;
            SPM.Sess(1).U(c).P(1).name = 'none';       % Parametric Modulation; 'none' for now
         end
        
        
        cd(strcat(currDir,'\run',num2str(RunNum)));
        HeadMotionDir=strcat(currDir,'\run',num2str(RunNum));
        rnam = {'X','Y','Z','x','y','z'};
        fn = spm_select('list',HeadMotionDir,'^rp.*\.txt');
        [r1,r2,r3,r4,r5,r6] = textread(fn,'%f%f%f%f%f%f');

        SPM.Sess(1).C.C = [r1 r2 r3 r4 r5 r6];
        SPM.Sess(1).C.name = rnam;


     
    %-------------------------------------------------------------
%         SPM.xGX.iGXcalc = 'none';
        SPM.xGX.iGXcalc = 'Scaling';

    %-------------------------------------------------------------
        SPM.xX.K(trial).HParam = 128;
        
    %-------------------------------------------------------------
        SPM.xVi.form       = 'AR(1)';
        
        SPM.xY.P = cat(1,tmp{:});
        SPM.xY.RT = TR;
        SPM = spm_fmri_spm_ui(SPM);
        cd(trialdir);
        SPM = spm_spm(SPM);
  
        clear SPM ;
        clear dir  
    
        
     end
     clear tmp
    end
    clear behavdir
end





