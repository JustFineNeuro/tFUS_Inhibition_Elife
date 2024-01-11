for Group=1:2
    %% Rject whack electrodes
    
    if Group==1
        
        subj=4:24
        rejchan{1}={'TP10','O2'}; %starts at subject 4
        rejchan{2}={''};
        rejchan{3}={''}; %Go back and remove (T7, TP7,TP8,T8), may require more epoch clearing
        rejchan{4}={'T8'};
        rejchan{5}={'CP2'};
        rejchan{6}={'T7','FT7'};
        rejchan{7}={'TP10'};
        rejchan{8}={''};
        rejchan{9}={'AF8'};
        rejchan{10}={''};
        rejchan{11}={'C5','FC5'};
        rejchan{12}={'T7'};
        rejchan{13}={'AF8','FP2'};
        rejchan{14}={''};
        rejchan{15}={''};
        rejchan{16}={''};
        rejchan{17}={'AF8'};
        rejchan{18}={'T7','P8'};
        rejchan{19}={''};
        rejchan{20}={''};
        rejchan{21}={'AF7','T8'};
    elseif Group==2
        subj=2:20
        %Rject whack electrodes
        rejchan{1}={''}; %starts at subject 2
        rejchan{2}={'TP7','TP10'};%DONE
        rejchan{3}={'F3','F4','AF7','AF3'}; %DONE
        rejchan{4}={'AF8','TP8','FP2'};%DONE
        rejchan{5}={''};%DONE
        rejchan{6}={'PO3'};%DONE
        rejchan{7}={''};%DONE
        rejchan{8}={'TP10'};%DONE
        rejchan{9}={'TP8'};%DONE
        rejchan{10}={'FP1','FP2','AF7','AF8'};%DONE
        rejchan{11}={'T7'}; %DONE
        rejchan{12}={'T7','T8'};%DONE
        rejchan{13}={'F5','AF4'};%DONE
        rejchan{14}={'AF7','TP10'};%DONE
        rejchan{15}={'AF8','PO9'};%DONE
        rejchan{16}={'P8'};%DONE
        rejchan{17}={'T8','P7','F6'};%DONE
        rejchan{18}={'PO10','TP7'};%DONE
        rejchan{19}={'F8','P8','TP10'};%DONE
    end
    
    
    %%Set directories
    eeglab_path='C:\Users\Justin\Documents\MATLAB\TOOLBOXES\NEURAL\eeglab_current\eeglab2019_1';
    
    dd=strcat('D:\TFUS_MODEL_STUDY\EEG\Group',num2str(Group),'\Long_EEG\')
    load(strcat('ChanKeep_Group',num2str(Group),'.mat')) %Electrodes to keep before ICA
    load('D:\TFUS_MODEL_STUDY\EXTERNAL_FILES\chans.mat')
    save_path=strcat('D:\TFUS_MODEL_STUDY\EEG\Group',num2str(Group),'\');
    chan_path='D:\\TFUS_MODEL_STUDY\EXTERNAL_FILES\';
    addpath('D:\TFUS_MODEL_STUDY\EXTERNAL_CODE\structfind')
    
    addpath(chan_path)
    addpath(eeglab_path)
    eeglab
    
    %% STEP 1.  clean before first ICA--- (DONE: 8/10/20- )
    for s=1:length(subj)
        if Group ==1
            
            %load dataset used for AMICA
            EEG = pop_loadset('filename',strcat('\S',num2str(subj(s)),'_step_2.set'),'filepath','D:\TFUS_MODEL_STUDY\EEG\Group1');
            
            oglocs=EEG.chanlocs;
            %Collect datafile and originally rejected epochs
            exdat=EEG.externaldat;
            rjman=EEG.rejected_manual;
            
            %Find trials still existing after 2nd rejection on epoched data
            tnum=structfind(EEG.event,'type','TRIAL_NUMBER');
            for i=1:length(tnum)
                tnumm{i}=EEG.event(tnum(i)).code;
            end
            
            for i=1:length(tnumm)
                temp=strsplit(tnumm{i},'_');
                tnum(i)=str2num(temp{2});
            end
            
            
            dd='D:\TFUS_MODEL_STUDY\EEG\Group1\Long_EEG';
            
            %Load continuous dataset
            EEG = pop_loadset('filename',strcat('S',num2str(subj(s)),'_long_EEG.set'),'filepath',dd);
            EEG = pop_interp(EEG,chanlocs)
            
            
            
            
            
            
            % Remove data-windows containing rejected epochs  using rjman
            GO_EV=structfind(EEG.event,'type','S4'); %Location of GO event
            templat=[EEG.event(GO_EV(find(rjman))).latency];
            GO_EV=GO_EV(find(rjman));
            for i=1:length(templat)
                [e,idx]=min(abs(EEG.times-templat(i)));
                rg(i,:)=[(idx-(750/(1000/EEG.srate))),(idx+(750/(1000/EEG.srate)))];
            end
            
            for i=1:length(GO_EV)
                tp=EEG.event(GO_EV(i));
                EEG.event(end+1).type=strcat('REJECT_EVENT_',num2str(i));
                EEG.event(end).latency= tp.latency+2;
                EEG.event(end).code='REJECT_EVENT_A';
                EEG.event(end).duration=1;
                
            end
            
            [~,index] = sortrows([EEG.event.latency].'); EEG.event = EEG.event(index(1:end)); clear index
            
            for i=1:length(GO_EV)
                EEG = pop_rmdat( EEG, {strcat('REJECT_EVENT_',num2str(i))},[-1.25 1.5] ,1);
            end
            
            
            %Then add pseudo-trial #s for 2nd manual rejection
            GO_EV=structfind(EEG.event,'type','S4'); %Location of GO event
            
            for i=1:length(GO_EV)
                tp=EEG.event(GO_EV(i));
                EEG.event(end+1).type=strcat('TRIAL_',num2str(i));
                EEG.event(end).latency= tp.latency+2;
                EEG.event(end).code='REJECT_EVENT_A';
                EEG.event(end).duration=1;
                
            end
            [~,index] = sortrows([EEG.event.latency].'); EEG.event = EEG.event(index(1:end)); clear index
            
            
            %Compare GO_EV and tnum for existing trials, then remove
            
            rm=setdiff(1:length(GO_EV),tnum);
            for i=1:length(rm)
                EEG = pop_rmdat( EEG, {strcat('TRIAL_',num2str(rm(i)))},[-1.25 1.5] ,1);
            end
            
            %remove inserted boundary events
            removeIdxArray = [];
            for b = 1:length(EEG.event)
                
                if strcmp(EEG.event(b).type, 'boundary')
                    removeIdxArray = [ removeIdxArray b ];
                end
                
            end
            EEG.event( removeIdxArray ) = [];
            
            EEG = fullRankAveRef(EEG);    %Reref
            
            clocs=EEG.chanlocs;
            EEG = clean_drifts(EEG,[0.25 0.75]);
            EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',35);
            
            EEG = clean_artifacts(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off','BurstCriterion',10,'WindowCriterion',0.25,'BurstRejection','on');
            EEG = fullRankAveRef(EEG);
            
            %%  Find electrodes to reject
            keepchan=[];
            for i=1:length(locsKeep{s})
                keepchan{i}=locsKeep{s}(i).labels;
            end
            
            %% remove electreodes from long_eeg before ICA
            EEG = pop_select(EEG, 'channel', keepchan )
            
            
            %% Run crap
            EEG = pop_continuousartdet( EEG , 'ampth',  100, 'chanArray',  1:EEG.nbchan, 'colorseg', [ 1 0.9765 0.5294], 'forder',  100, 'numChanThreshold',  1, 'stepms',  250, 'threshType', 'peak-to-peak', 'winms',  500 ); % GUI: 13-Jul-2020 20:51:35
            
            pop_saveset(EEG,'filepath',strcat('D:\TFUS_MODEL_STUDY\EEG\Group1\Long_EEG\S',num2str(subj(s)),'_long_step_2.set'));%save
            
        elseif Group==2 %This may not be right and I may need to grab the proper electrodes to keep.... but who cvares
            clearvars -except subj s rejchan dd chan_path eeglab_path
            
            dd='D:\TFUS_MODEL_STUDY\EEG\Group2\Long_EEG';
            
            %Load continuous dataset
            EEG = pop_loadset('filename',strcat('S',num2str(subj(s)),'_long_EEG.set'),'filepath',dd);
            
            
            %% interpolate the elecrodes pre-ASR
            
            
            EEG=pop_select(EEG,'nochannel',rejchan{s});
            
            EEG.preinterplocs=EEG.chanlocs;
            
            %% load channel file
            load(strcat(chan_path,'chans.mat'));
            EEG=interpol_chan(EEG,chanlocs);
            
            %% average refeence while correting for rank
            EEG = fullRankAveRef(EEG);    %Reref
            
            clocs=EEG.chanlocs;
            EEG = clean_drifts(EEG,[0.25 0.75]);
            EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',40);
            
            EEG = clean_artifacts(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off','BurstCriterion',10,'WindowCriterion',0.25,'BurstRejection','on');
            EEG = fullRankAveRef(EEG);
            
            %%  remove electrodes again as we will run ICA only good channels only
            for i=1:length(EEG.preinterplocs)
                keep{i}=EEG.preinterplocs(i).labels;
            end
            EEG=pop_select(EEG,'channel',keep)
            
            
            %% Run crap
            EEG = pop_continuousartdet( EEG , 'ampth',  200, 'chanArray',  1:EEG.nbchan, 'colorseg', [ 1 0.9765 0.5294], 'forder',  100, 'numChanThreshold',  1, 'stepms',  250, 'threshType', 'peak-to-peak', 'winms',  500 ); % GUI: 13-Jul-2020 20:51:35
            
            pop_saveset(EEG,'filepath',strcat(dd,'\S',num2str(subj(s)),'_long_step_2.set'));%save
        end
    end
    
end

%% REMOVE BOUNDARY EVENTS (DONE- 8/10/20 - JMF)
for Group=1:2
    if Group==1
        subj=4:24
    elseif Group==2
        subj=2:20
    end
    
    for s=1:length(subj)
        EEG=pop_loadset('filepath',strcat('D:\TFUS_MODEL_STUDY\EEG\Group',num2str(Group),'\Long_EEG\S',num2str(subj(s)),'_long_step_2.set'));%save
        
        EEG.event(structfind(EEG.event,'type','boundary'))=[];
        %% REMOVE ALL DATA BEFORE FIRST EVET
        
        [~,idx]=min(abs(EEG.times-EEG.event(1).latency));
        EEG = eeg_eegrej( EEG, [1 idx-25] );
        EEG.event(structfind(EEG.event,'type','boundary'))=[];
        
        pop_saveset(EEG,'filepath',strcat('D:\TFUS_MODEL_STUDY\EEG\Group',num2str(Group),'\Long_EEG\'),'filename',strcat('S',num2str(subj(s)),'_long_step_2.set'));%save
    end
end



%% ICA for cleaning only
for Group=1:2
    if Group==1
        subj=4:24
    elseif Group==2
        
        subj=2:20
    end
    for s=1:length(subj)
        EEG=pop_loadset('filepath',strcat('D:\TFUS_MODEL_STUDY\EEG\Group',num2str(Group),'\Long_EEG\'),'filename',strcat('S',num2str(subj(s)),'_long_step_2.set'));%save
        
        %% Run ICA Extended infomax
        varargin{1}='lrate';
        varargin{2}=[1.0000e-03];
        varargin{3}='extended';
        varargin{4}=1;
        data=EEG.data;
        data = reshape(data, size(data,1), size(data,2)*size(data,3));
        [nchans,nframes] = size(data);
        [weights_epoch,sphere_epoch]=mycudaica(data,['S' num2str(subj(s))],varargin{1},varargin{2},varargin{3},varargin{4});
        %Reorganze by variance
        data = sphere_epoch * data;
        winv = inv(weights_epoch*sphere_epoch);
        meanvar = sum(winv.^2).*sum((data').^2)/((nchans*nframes)-1);
        [~, windex] = sort(meanvar);
        ncomps=nchans;
        windex = windex(ncomps:-1:1); % order from large to small
        weights_epoch = weights_epoch(windex,:);% reorder the weight matrix
        EEG.icaact=[];
        EEG.icawinv=[];
        EEG.icasphere=sphere_epoch;
        EEG.icaweights=weights_epoch;
        EEG.icachansind=1:EEG.nbchan;
        save(strcat('D:\TFUS_MODEL_STUDY\EEG\Group',num2str(Group),'\Long_EEG\ICA\EINFOMAX\S',num2str(subj(s)),'_weights_long.mat'),'weights_epoch')
        save(strcat('D:\TFUS_MODEL_STUDY\EEG\Group',num2str(Group),'\Long_EEG\ICA\EINFOMAX\S',num2str(subj(s)),'_sphere_long.mat'),'sphere_epoch')
        
        EEG=eeg_checkset(EEG,'ica')
        pop_saveset(EEG,'filepath',strcat('D:\TFUS_MODEL_STUDY\EEG\Group',num2str(Group),'\Long_EEG\ICA\EINFOMAX\S',num2str(subj(s)),'_long_step_2_ICA_1.set'));%save
    end
end



%% RUN AMICA THRU SERVER
for Group=1:2
    for models=1:3
        if Group==1
            subj=4:24;
        elseif Group==2
            subj=2:20
        end
        mkdir(strcat('D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA'));
        %Make directories
        %% COPY DATAFILES
        
        for s=subj
            mkdir(strcat('D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA\S',num2str(s),'G',num2str(Group),'M',num2str(models)));
        end
        
        for s=1:length(subj)
            fid = fopen('D:\TFUS_MODEL_STUDY\EEG\Group1\Long_EEG\ICA\AMICA\NSG_AMICA.m','r');
            i = 1;
            tline = fgetl(fid);
            A{i} = tline;
            while ischar(tline)
                i = i+1;
                tline = fgetl(fid);
                A{i} = tline;
            end
            fclose(fid);
            A{14}=strcat('num_models=',num2str(models));
            A{6}=strcat('filename = ''S',num2str(subj(s)),'_long_step_2.set''')
            A{end-1}='';
            %     A{end-1}=strcat('pop_saveset(EEG,''filename'',''S',num2str(subj(s)),'_long_step_2_AMICA.set'')');
            % Change cell A
            % Write cell A into txt
            fid = fopen(strcat('D:\TFUS_MODEL_STUDY\EEG\\Long_EEG\ICA\AMICA\S',num2str(subj(s)),'G',num2str(Group),'M',num2str(models),'\NSG_AMICA.m'),'w+');
            for i = 1:numel(A)
                if A{i+1} == -1
                    fprintf(fid,'%s', A{i});
                    break
                else
                    fprintf(fid,'%s\n', A{i});
                end
            end
            fclose(fid)
            copyfile(strcat('D:\TFUS_MODEL_STUDY\EEG\Group',num2str(Group),'\Long_EEG\S',num2str(subj(s)),'_long_step_2.set'),strcat('D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA\S',num2str(subj(s)),'G',num2str(Group),'M',num2str(models)))
            copyfile(strcat('D:\TFUS_MODEL_STUDY\EEG\Group',num2str(Group),'\Long_EEG\S',num2str(subj(s)),'_long_step_2.fdt'),strcat('D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA\S',num2str(subj(s)),'G',num2str(Group),'M',num2str(models)))
        end
        
        
       
    end
end


%% Submit AMICA jobs to server
for Group=1:2
    for models=1:3
        if Group==1
            subj=4:24;
        elseif Group==2
            subj=2:20;
        end
        
        for i=1:length(subj)
            job_url{i}=nsg_run(strcat('D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA\S',num2str(subj(i)),'G',num2str(Group),'M',num2str(models)),'jobid',strcat('AMICA_LONG_S',num2str(subj(i)),'G',num2str(Group),'M',num2str(models)),'outfile',strcat('NSG_AMICA_','G',num2str(Group),'S',num2str(subj(i)),'M',num2str(models)),'runtime',48,'filename','NSG_AMICA.m','subdirname','','statusemail','true','emailaddress','')
        end
    end
end

%% Download AMICA JOBS

tmpalljobs    = nsg_jobs;
for i=1:length((tmpalljobs.joblist.jobs.jobstatus))
    nsg_download(tmpalljobs.joblist.jobs.jobstatus{i}.selfUri.url)
end




%% Load AMICA and sphere into long_EEG.set 
%% THIS ROUTINE:

%1. builds design matrix, cleans data, removes weird or bad events and
%savee a file that is just GO trials for the GLM
%2. DO like original, except load IC and dipole informaton from :        pop_loadset(EEG,'filepath',strcat('D:\TFUS_MODEL_STUDY\EEG\Group',num2str(Group),'\Long_EEG\ICA\AMICA\S',num2str(subj(s)),'_long_step_3_POST_AMICA_ALL.set'));%save
%3. Calculate smoothed model probabilities over time. 
%4. Calculate ICAACCT with and without model probability weighting. 
%5. DO GLM: ERP, BETA - mean and burst

%% if model probabilities pre-weighting doesn't work well, try calculating model probabilities (evoked) ahead of time and weighting...
%% % If that fails, just use 1 model IC and call it a day.


for Group=1:2
    if Group==1
        subj=4:24;
    elseif Group==2
        subj=2:20;
    end
    
    parfor s=1:length(subj)
        make_go_step5(s,subj,Group)
    end
end



%% Do stop trials


for Group=1:2
    if Group==1
        subj=4:24;
    elseif Group==2
        subj=2:20;
    end
    
    parfor s=1:length(subj)
        make_stop_step5(s,subj,Group)
    end
end
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','E_mail','justfineneuro@gmail.com');
setpref('Internet','SMTP_Username','justfineneuro');
setpref('Internet','SMTP_Password','12BUSter34');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
sendmail('justfineneuro@gmail.com','Data made') ;


%% Do dipfitting, eloreta, and classificaiton per ica in each ICA model set
clear all
addpath('C:\Users\Justin\Documents\MATLAB\TOOLBOXES\NEURAL\eeglab_current\eeglab2019_1')
prefix='D:\\TFUS_MODEL_STUDY\EEG\Long_EEG';
addpath('C:\Users\Justin\Documents\MATLAB\TOOLBOXES\NEURAL\eeglab\plugins\Fieldtrip-lite20210128')
% addpath('C:\Users\Justin\Documents\MATLAB\TOOLBOXES\NEURAL\eeglab')
addpath(genpath('C:\Users\Justin\Documents\MATLAB\TOOLBOXES\NEURAL\eeglab\plugins\dipfit3.7'))
addpath 'D:\TFUS_MODEL_STUDY\EXTERNAL_CODE'
eeglab_path='C:\Users\Justin\Documents\MATLAB\TOOLBOXES\NEURAL\eeglab_current\eeglab2019_1';
eeglab
load('D:\TFUS_MODEL_STUDY\BEHAVIORAL_OUTPUT\BEHAVE.mat')


for Group=1:2
    if Group==1
        subj=4:24;
    elseif Group==2
        subj=2:20;
    end
    parfor s=1:length(subj)
        AMICA_DIP_CLASSIFY(subj,Group,s)
    end
end




for Group=1:2
    if Group==1
        subj=4:24;
    elseif Group==2
        subj=2:20;
    end
    for s=1:length(subj)
        for cond={'GO','STOP'}
            try
                EEG = pop_loadset('filename',strcat('S',num2str(subj(s)),'G',num2str(Group),'_long_step_6_AMICA_',cond{1},'.set'),'filepath','D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA');

                for models=1:2
                    EEG = pop_loadmodout(EEG,strcat('D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA\NSG_AMICA_G',num2str(Group),'S',num2str(subj(s)),'M',num2str(models),'\','S',num2str(subj(s)),'G',num2str(Group),'M',num2str(models),'\amicaouttmp'))
                    EEG.amica(models).amica=EEG.etc.amica;
                    EEG.icasphere=[];
                    EEG.icawinv=[];
                    EEG.icaweights=[];
                    EEG.icachansind=[];
                    %Make an EEG.icaact for model set 1 (go from EEG.icaact
                    %to EEG.amica(models).icaact);
                    
                    
                end
                
                %% reject the undesired ICs manually
                %% compile one largeass EEG.icawinv,EEG.icaapshere (Columns are electrodes, so leave all of those)
                %% Make a largess EEG.icaact == EEG.icawinv * EEG.data
                %% make a set probability weighted using the computed v_smooth in EEG.amica().v_smooth
                %% Compute TF power, alpha, beta (reg and burst: count and amplitude<-- will require some resoting of beta_burst calculation output to put back into EEG.icaact)
                pop_saveset(EEG,'filepath',strcat('D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA\S',num2str(subj(s)),'G',num2str(Group),'_long_step_6_AMICA_',cond{1},'.set'));%save
            end
        end
    end
end

setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','E_mail','justfineneuro@gmail.com');
setpref('Internet','SMTP_Username','justfineneuro');
setpref('Internet','SMTP_Password','12BUSter34');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
sendmail('justfineneuro@gmail.com','Dipfit Done') ;

% Make all meta data
brainprob=0.50;
rvthresh=0.18;

for Group=1:2
    if Group==1
        subj=4:24;
    elseif Group==2
        subj=2:20;
    end
    parfor s=1:length(subj)
        MakeClustMetaData(subj,s,Group,rvthresh,brainprob)
    end
end




%% Gather ICA activatons, compute time frequency mean power and beta burst_ also gets information for clustering
brainprob=0.50;
rvthresh=0.18;

backproj=0;
for Group=1:2
    if Group==1
        subj=4:24;
    elseif Group==2
        subj=2:20;
    end
    parfor s=1:length(subj)
        Step_7_Evoked_TimeFrequency_GlM_Prep(subj,s,Group,rvthresh,brainprob,backproj)
    end
end
% 
% setpref('Internet','SMTP_Server','smtp.gmail.com');
% setpref('Internet','E_mail','justfineneuro@gmail.com');
% setpref('Internet','SMTP_Username','justfineneuro');
% setpref('Internet','SMTP_Password','');
% props = java.lang.System.getProperties;
% props.setProperty('mail.smtp.auth','true');
% props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
% props.setProperty('mail.smtp.socketFactory.port','465');
% sendmail('justfineneuro@gmail.com','Time frequency dnone') ;
% 

%% Do beta burst separately...must use continuous filter because we will still need to reject these bitches later

for Group=1:2
    if Group==1
        subj=4:24;
    elseif Group==2
        subj=2:20;
    end
    parfor s=1:length(subj)
       Parallel_BetaBurstIC(subj,s,Group)
    end
end


%% Calculate burst rate
addpath('C:\Users\Justin\Documents\MATLAB\TOOLBOXES\NEURAL\eeglab_current\eeglab2019_1')
eeglab
for Group=1:2
    if Group==1
        subj=4:24;
    elseif Group==2
        subj=2:20;
    end
    for s=1:length(subj)
        for cond={'GO','STOP'}
            try
            EEG = pop_loadset('filename',strcat('S',num2str(subj(s)),'G',num2str(Group),'_long_step_8_AMICA_',cond{1},'.set'),'filepath','D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA');
            binsize=0.01; %10 ms bin size
            kernel.taxis = -3:binsize:3;  % make a time axis of 1000 ms  %%TODO change to sectons and multiple by sampling rate (smaples/second)
            gauss_SD=100/1000/binsize;
            kernel.kernel = normpdf(kernel.taxis, 0, gauss_SD);
            kernel.kernel = kernel.kernel ./ sum(kernel.kernel);
            
            
            btburst=EEG.ICBetaBursts.binary';
            [r,c]=size(btburst);
            btburst=btburst(:);
            BstRate=reshape(conv(btburst,kernel.kernel/binsize,'same'),r,c)';
            EEG.ICBetaBursts.rate=BstRate;
            pop_saveset(EEG,'filepath',strcat('D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA\S',num2str(subj(s)),'G',num2str(Group),'_long_step_8_AMICA_',cond{1},'.set'));%save
            movefile(strcat('D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA\S',num2str(subj(s)),'G',num2str(Group),'_long_step_8_AMICA_',cond{1},'.set'),strcat('D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA\S',num2str(subj(s)),'G',num2str(Group),'_GLM_READY_',cond{1},'.set'))
            movefile(strcat('D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA\S',num2str(subj(s)),'G',num2str(Group),'_long_step_8_AMICA_',cond{1},'.fdt'),strcat('D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA\S',num2str(subj(s)),'G',num2str(Group),'_GLM_READY_',cond{1},'.fdt'))
            end
        end
    end
end

% 
% setpref('Internet','SMTP_Server','smtp.gmail.com');
% setpref('Internet','E_mail','justfineneuro@gmail.com');
% setpref('Internet','SMTP_Username','justfineneuro');
% setpref('Internet','SMTP_Password','');
% props = java.lang.System.getProperties;
% props.setProperty('mail.smtp.auth','true');
% props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
% props.setProperty('mail.smtp.socketFactory.port','465');
% sendmail('justfineneuro@gmail.com','Finished beta burst') ;


%% Make a separate iclisting file columns: (1) amica #models  (2)amica model #, (3) IC in original list, (4) computation: ERP, time frequency etc... 
% cfg.group(1).subj=[4:24];
% cfg.group(2).subj=[2:20];
% cfg.Groups=1:2;
% 
% for g=cfg.Groups
%     for s=1:length(cfg.group(g).subj)
%         EEG = pop_loadset('filename',strcat('S',num2str(cfg.group(g).subj(s)),'G',num2str(g),'_GLM_READY_STOP.set'),'filepath','D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA');
%         iclist=EEG.iclisting(find(EEG.iclisting(:,end)==1),:)
%         save(strcat('iclist_','S',num2str(cfg.group(g).subj(s)),'G',num2str(g),'.mat'),'iclist')
%     end
% end

%% DO GLM: ERP, ALPHA, BETA

clear all
addpath('C:\Users\Justin\Documents\MATLAB\TOOLBOXES\NEURAL\eeglab_current\eeglab2019_1')
cfg.par_fold='D:\TFUS_MODEL_STUDY\';
cfg.output_path=[cfg.par_fold,'EEG_OUTPUT\AllGroups\']
cfg.tf_splines=40;
cfg.erp_basis=20;
cfg.do_burst=0;
cfg.erp_times=[-0.4 1.25];
cfg.tf_times=[-0.5 1.5];
cfg.do_bootstrap=[];
cfg.backproj=0;
cfg.do_poisson=0;
cfg.do_logistic=0;
cfg.use_gpu=0;

for amica_num=1
    cfg.amica_num=amica_num;
    for Group=1:2
        if Group==1
            subj=4:24;
            
        elseif Group==2
            subj=2:20;
            
        end
        

        for cond={'STOP'}
            
%             gpuDevice(1)
%             delete(gcp('nocreate'))
            %         parpool(2)
%             sjsj=waitbar(0,'subjects')
            parfor s=1:length(subj)
                GLM_MAIN(subj,s,cond,Group,cfg)
%                 waitbar(s/length(subj),sjsj)
            end
        end
    end
end







setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','E_mail','justfineneuro@gmail.com');
setpref('Internet','SMTP_Username','justfineneuro');
setpref('Internet','SMTP_Password','12BUSter34');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
sendmail('justfineneuro@gmail.com','Finished GLM') ;



%% START CLUSTERING, GROUPIG, ANALYSIS ICs 

% %% CFG file settings

% 
% 
% 
% %%% These settings matter for cfg.cluster_type='template_matching' (not
% %%% done)
% % %
% % % cfg.matching_StaticVars(1).props.type='dippos';
% % % cfg.matching_StaticVars(1).props.template=[-40 -24 50];
% % % cfg.matching_StaticVars(1).props.name='lm1';
% % % cfg.matching_StaticVars(2).props.type='dippos'
% % % cfg.matching_StaticVars(2).props.template=[44 46 12];
% % % cfg.matching_StaticVars(2).props.name='rifg';
% 
% % % cfg.matching_StaticVars(3).name='ScalpWinv';
% 
% 
% 
% %%% These settings matter for cfg.cluster_type='cluster'
% 
% cfg.cluster_StaticVars(1).name='dippos';
% cfg.cluster_StaticVars(2).name='dipmom';
% cfg.cluster_StaticVars(3).name='ScalpWinv';
% cfg.cluster_DynamicVars=[];
% %
% % cfg.cluster_StaticVars(1).name='dippos';
% % cfg.cluster_StaticVars(2).name='ScalpWinv';
% %
% %
% % %
% % cfg.cluster_DynamicVars(1).type='ERP';
% % cfg.cluster_DynamicVars(1).subname=[];
% % cfg.cluster_DynamicVars(1).cond='STOP';
% % cfg.cluster_DynamicVars(1).GLM_MODEL='A';
% % cfg.cluster_DynamicVars(1).beta_cols=[1]; % SS stopping only
% % cfg.cluster_DynamicVars(1).twin=[0 0.5];
% % %
% %
% % cfg.cluster_DynamicVars(2).type='BetaBurst';
% % cfg.cluster_DynamicVars(2).subname='logexp'; %Other is .logexp
% % cfg.cluster_DynamicVars(2).cond='STOP'
% % cfg.cluster_DynamicVars(2).GLM_MODEL='ONLYSTOP';
% % cfg.cluster_DynamicVars(2).beta_cols=[1];  % SS stopping only
% % cfg.cluster_DynamicVars(2).twin=[0 0.5];
% %
% 
% %
% % cfg.cluster_DynamicVars(3).type='BetaBurst';
% % cfg.cluster_DynamicVars(3).subname='logexp'; %Other is .logexp
% % cfg.cluster_DynamicVars(3).cond='SSRT'
% % cfg.cluster_DynamicVars(3).GLM_MODEL='ONLYSSRT';
% % cfg.cluster_DynamicVars(3).beta_cols=[1];  % SS stopping only
% % cfg.cluster_DynamicVars(3).twin=[-0.10 -0.01];
% %
% % cfg.bootstrap_DynamicVars(1).type='BetaBurst';
% % cfg.bootstrap_DynamicVars(1).subname='logexp'; %Other is .logexp
% % cfg.bootstrap_DynamicVars(1).cond='SSRT'
% % cfg.bootstrap_DynamicVars(1).GLM_MODEL='ONLYSSRT';
% % cfg.bootstrap_DynamicVars(1).beta_cols=[1]
% % cfg.bootstrap_DynamicVars(1).twin=[-0.10 -0.01];
% 


%%% Setup cluster settings as cfg.
clear all
cfg=[]
cfg.group(1).subj=[4:15,17:18,20:24];
cfg.group(2).subj=[2:15,17 19:20];
cfg.par_fold='D:\TFUS_MODEL_STUDY\';
cfg.plot=1;


% %% CFG file settings
% 
cfg.Groups=1:2;
cfg.cluster_StaticVars(1).name='dippos';
cfg.cluster_StaticVars(2).name='dipmom';
cfg.cluster_StaticVars(3).name='ScalpWinv';
cfg.cluster_DynamicVars=[];
%%%% Best settings so far %%%%%%
cfg.amica_num=2;

if cfg.amica_num==1
    cfg.clust_file_prefix='ClustData_'
         cfg.algo.options.clusters=11;

elseif cfg.amica_num==2
    cfg.clust_file_prefix='ClustData_2MOD_'
     cfg.algo.options.clusters=13;

end
cfg.clust_dat_folder='D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA';
cfg.glm_dat_folder='D:\TFUS_MODEL_STUDY\EEG_OUTPUT\AllGroups\';
cfg.cluster_type='cluster'
%Algorithm options
cfg.algo.name='kmedoids';
cfg.algo.input_type='dissimilarity';

%This will change pending the algorithm used

cfg.algo.weights=[10 4 4]
 cfg.algo.weights=[15 6 2]./sum([15 6 2]);

% % reducing data down more
cfg.reduce_ICS_brain=.75;
cfg.cut_ICS_eye=.7; %Not implamented
cfg.save_scalp=1;
[CLUST_OUT]=ClusteringICs(cfg)
save(['CLUSTER_OUTPUT_',date,'.mat'],'CLUST_OUT')

cfg_cluster=cfg;
%% Gather variables to plot using cluster assignment + if amica_num==2, and a subject has multiple ICs within a cluster, summing necessary to peform averge.


%%%% MOST IMPORTANT CONSIDERATION IS HOW IC's are weighted and selected...
%%%% Easiest is just take all of them and average --- mod_prob==zeros
cfg=[];


cfg.par_fold='D:\TFUS_MODEL_STUDY\';
cfg.save_fold='D:\TFUS_MODEL_STUDY\EEG_OUTPUT\';
cfg.save=1;
cfg.edit_master_cfg=1;


cfg.group(1).subj=[4:15,17:18,20:24];
cfg.group(2).subj=[2:15,17 19:20];
load('D:\TFUS_MODEL_STUDY\EXTERNAL_FILES\chans.mat')
cfg.chanlocs=chanlocs;
cfg.amica_num=2;
cfg.reduce_ICS_brain=.8;

if cfg.amica_num==1
    cfg.clust_file_prefix='ClustData_'
elseif cfg.amica_num==2
    cfg.clust_file_prefix='ClustData_2MOD_'
end
cfg.clust_dat_folder='D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA';
cfg.glm_dat_folder='D:\TFUS_MODEL_STUDY\EEG_OUTPUT\AllGroups\';
cfg.measures={'ERP','ERP','ERP','TF','TF','TF','TF','TF','TF','TF','TF','BetaBurst','BetaBurst','BetaBurst','BetaBurst'};
cfg.conditions={'STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP'}
cfg.IC_BACKPROJ={'IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC'};
cfg.BACKPROJ_elec={'','','','','','','','','','','','','','',''};
cfg.BACKPROJ_operation={'','','','','','','','','','','','','','',''};
cfg.IC_SELECTOR.type={'ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL'}; %Can be 1 closest to dipole center, can be highest power at other electrode, etc...
cfg.IC_SELECTOR.subinfo={'ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL'};
cfg.stats_model={'A','D','E','A','B','C','D','A','B','C','D','A','B','A','B','A','D'};
cfg.subname={'','','','LowBeta','LowBeta','LowBeta','LowBeta','alpha','alpha','alpha','alpha','normal','normal','logit_raw','logit_raw'};
cfg.mod_prob_weight=ones(1,length(cfg.stats_model)); %% todo:  Add options to select IC with largest value at a speciric electrode...
cfg.gather_type={'cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster'}
cfg.CLUST_OUT=CLUST_OUT;
cfg.prefix='weighted_amica_weight_ALL_2MOD'

[STUDY_GLM,cfg]=gather_data(cfg)




cfg.par_fold='D:\TFUS_MODEL_STUDY\';
cfg.save_fold='D:\TFUS_MODEL_STUDY\EEG_OUTPUT\';
cfg.save=1;
cfg.edit_master_cfg=1;


cfg.group(1).subj=[4:15,17:18,20:24];
cfg.group(2).subj=[2:15,17 19:20];
load('D:\TFUS_MODEL_STUDY\EXTERNAL_FILES\chans.mat')
cfg.chanlocs=chanlocs;
cfg.amica_num=2;
cfg.reduce_ICS_brain=.8;

if cfg.amica_num==1
    cfg.clust_file_prefix='ClustData_'
elseif cfg.amica_num==2
    cfg.clust_file_prefix='ClustData_2MOD_'
end
cfg.clust_dat_folder='D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA';
cfg.glm_dat_folder='D:\TFUS_MODEL_STUDY\EEG_OUTPUT\AllGroups\';
cfg.measures={'ERP','ERP','ERP','TF','TF','TF','TF','TF','TF','TF','TF','BetaBurst','BetaBurst','BetaBurst','BetaBurst'};
cfg.conditions={'STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP'}
cfg.IC_BACKPROJ={'IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC'};
cfg.BACKPROJ_elec={'','','','','','','','','','','','','','',''};
cfg.BACKPROJ_operation={'','','','','','','','','','','','','','',''};
cfg.IC_SELECTOR.type={'dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist'}; %Can be 1 closest to dipole center, can be highest power at other electrode, etc...
cfg.IC_SELECTOR.subinfo={'dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist','dipdist'};
cfg.stats_model={'A','D','E','A','B','C','D','A','B','C','D','A','B','A','B','A','D'};
cfg.subname={'','','','LowBeta','LowBeta','LowBeta','LowBeta','alpha','alpha','alpha','alpha','normal','normal','logit_raw','logit_raw'};
cfg.mod_prob_weight=ones(1,length(cfg.stats_model)); %% todo:  Add options to select IC with largest value at a speciric electrode...
cfg.gather_type={'cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster'}
cfg.CLUST_OUT=CLUST_OUT;
cfg.prefix='weighted_amica_weight_dipdist_2MOD'

[STUDY_GLM,cfg]=gather_data(cfg)

% %IF cfg.gather_type{1}==template_matching
% load('D:\TFUS_MODEL_STUDY\EXTERNAL_FILES\SCALP_TEMPLATES.mat')
% cfg.IC_SELECTOR.type={'dipdist','maxelec','scalpwinv'};
% cfg.IC_SELECTOR.inputs.dipdist={{-40 -24 50},{0 10 28},{44 46 12}};
% cfg.IC_SELECTOR.inputs.maxelec={{'C3','CP3'},{'FCz'},{'F4','F6','F8'}};
% cfg.IC_SELECTOR.inputs.scalpwinv={{SCALP_TEMPLATES.SCALPWINV.LM1},{SCALP_TEMPLATES.SCALPWINV.PSMA},{SCALP_TEMPLATES.SCALPWINV.RIFG}};
% 

%%%%%%% JUST ERPs
cfg=[];


cfg.par_fold='D:\TFUS_MODEL_STUDY\';
cfg.save_fold='D:\TFUS_MODEL_STUDY\EEG_OUTPUT\STUDY_GLMS\';
cfg.save=1;
cfg.edit_master_cfg=1;


cfg.group(1).subj=[4:15,17:18,20:24];
cfg.group(2).subj=[2:15,17 19:20];
load('D:\TFUS_MODEL_STUDY\EXTERNAL_FILES\chans.mat')
load('D:\TFUS_MODEL_STUDY\EXTERNAL_FILES\CLUSTER_OUTPUT_2MOD_BEST_11-Mar-2021.mat')
cfg.chanlocs=chanlocs;
cfg.amica_num=2;
cfg.reduce_ICS_brain=.65;

if cfg.amica_num==1
    cfg.clust_file_prefix='ClustData_'
elseif cfg.amica_num==2
    cfg.clust_file_prefix='ClustData_2MOD_'
end
cfg.clust_dat_folder='D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA';
cfg.glm_dat_folder='D:\TFUS_MODEL_STUDY\EEG_OUTPUT\AllGroups\';
% cfg.measures={'ERP','ERP','ERP','ERP'};
% cfg.conditions={'STOP','STOP','STOP','STOP'}
% cfg.IC_BACKPROJ={'IC','IC','IC','IC'};
% % cfg.BACKPROJ_elec={'','','',''};
% % cfg.BACKPROJ_operation={'','','mean','mean'};
% cfg.IC_SELECTOR.type={'ALL','ALL','ALL','ALL'}; %Can be 1 closest to dipole center, can be highest power at other electrode, etc...
% cfg.IC_SELECTOR.subinfo={'ALL','ALL','ALL','ALL'};
% cfg.stats_model={'A','B','C','D'};
% cfg.subname={'','','',''};
% cfg.mod_prob_weight=ones(1,length(cfg.stats_model)); %% todo:  Add options to select IC with largest value at a speciric electrode...
% cfg.gather_type={'cluster','cluster','cluster','cluster'}
cfg.CLUST_OUT=CLUST_OUT;
cfg.prefix='weighted_amica_weight_dipdist_2MOD'

[STUDY_GLM,cfg]=gather_data(cfg)
cfg.mod_prob_weight=zeros(1,length(cfg.stats_model)); %% todo:  Add options to select IC with largest value at a speciric electrode...

cfg.prefix='equal_amica_weight_2MOD'

[STUDY_GLM,cfg]=gather_data(cfg)


load('D:\TFUS_MODEL_STUDY\EXTERNAL_FILES\CLUSTER_OUTPUT_1MOD_BEST_13-Mar-2021.mat')
cfg.CLUST_OUT=CLUST_OUT;

cfg.amica_num=1;
cfg.reduce_ICS_brain=.65;

if cfg.amica_num==1
    cfg.clust_file_prefix='ClustData_'
elseif cfg.amica_num==2
    cfg.clust_file_prefix='ClustData_2MOD_'
end
cfg.prefix='equal_amica_weight_1MOD'

[STUDY_GLM,cfg]=gather_data(cfg)


%% Gather data Model 1

load('CLUSTER_OUTPUT_1MOD_03-May-2021.mat')
cfg=[];
cfg.par_fold='D:\TFUS_MODEL_STUDY\';
cfg.save_fold='D:\TFUS_MODEL_STUDY\EEG_OUTPUT\';
cfg.save=1;
cfg.edit_master_cfg=1;


cfg.group(1).subj=[4:15,17:18,20:24];
cfg.group(2).subj=[2:15,17 19:20];
load('D:\TFUS_MODEL_STUDY\EXTERNAL_FILES\chans.mat')
cfg.chanlocs=chanlocs;
cfg.amica_num=1;
cfg.reduce_ICS_brain=.8;

if cfg.amica_num==1
    cfg.clust_file_prefix='ClustData_'
elseif cfg.amica_num==2
    cfg.clust_file_prefix='ClustData_2MOD_'
end
cfg.clust_dat_folder='D:\TFUS_MODEL_STUDY\EEG\Long_EEG\ICA\AMICA';
cfg.glm_dat_folder='D:\TFUS_MODEL_STUDY\EEG_OUTPUT\AllGroups\';
cfg.measures={'ERP','ERP','ERP','TF','TF','TF','TF','TF','TF','TF','TF','BetaBurst','BetaBurst','BetaBurst','BetaBurst','TF','TF','TF','TF','TF'};
cfg.conditions={'STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP','STOP'}
cfg.IC_BACKPROJ={'IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC','IC'};
cfg.BACKPROJ_elec={'','','','','','','','','','','','','','','','','','','',''};
cfg.BACKPROJ_operation={'','','','','','','','','','','','','','','''','','','',''};
cfg.IC_SELECTOR.type={'ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL'}; %Can be 1 closest to dipole center, can be highest power at other electrode, etc...
cfg.IC_SELECTOR.subinfo={'ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL','ALL'};
cfg.stats_model={'A','D','E','A','B','C','D','A','B','C','D','A','B','A','B','A','D','A','B','C','D','E'};
cfg.subname={'','','','LowBeta','LowBeta','LowBeta','LowBeta','alpha','alpha','alpha','alpha','normal','normal','logit_raw','logit_raw','theta','theta','theta','theta','theta'};
cfg.mod_prob_weight=ones(1,length(cfg.stats_model)); %% todo:  Add options to select IC with largest value at a speciric electrode...
cfg.gather_type={'cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster','cluster'}
cfg.CLUST_OUT=CLUST_OUT;
cfg.prefix='weighted_ALL_1MOD_V2'

[STUDY_GLM,cfg]=gather_data(cfg)



%% plot all ICs
figure
for i=1:15
subplot(5,5,i); topoplot(CLUST_OUT.scalpwinv.average{i,1},chanlocs)
end
%% STats and Plotting from STUDY_GLM


cfg.operations.plot={'plot_Group_mean','plot_Grand_mean'}

%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  STATS BELOW - manual right now                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.glmfilename='STUDY_GLM_ICS_weighted_amica_weight_dipdist_2MOD_15-Mar-2021.mat';
cfg.data_file=['D:\TFUS_MODEL_STUDY\EEG_OUTPUT\STUDY_GLMS\',cfg.glmfilename]
cfg.glmvars_file='D:\TFUS_MODEL_STUDY\EEG_OUTPUT\STUDY_GLMS\MetaData\GLM_VARS.mat';
cfg.glm_model_name='A'
cfg.glm_model_DV={'ERP','TF'}

cfg.ref_cols{1,1}={'intercept','stop'} %This reference is also successful stopping
cfg.ref_cols{2,1}={'ssd_upper','stop'} %This reference is also successful stopping
%%% Get ref columns
ref_cols=[];
for rfidx=1:size(cfg.ref_cols,1)
ref_cols=[ref_cols;intersect(structfind(eval(['GLM_VARS.MODELS.',cfg.event{evidx},'.',cfg.glm_model_DV{dvidx},'.',cfg.glm_model_name]),'name',cfg.ref_cols{rfidx}{1,1}),structfind(eval(['GLM_VARS.MODELS.',cfg.event{evidx},'.',cfg.glm_model_DV{dvidx},'.',cfg.glm_model_name]),'event',cfg.ref_cols{rfidx}{1,2}))];
end

cfg.stats_models='SSvUS + Stim + Group + SSvUS*Stim + SSvUS*Stim*Group'
cfg.reg_model_type='bayes','ols','glmnet','lmm'
cfg.clust_analyze=
cfg.design_col.cluster=1;
cfg.design_col.group=2;
cfg.design_col.subjects=3;
cfg.smooth_time=gausswin(7);
cfg.timepoints=[-0.2 0.5] <-- example
cfg.times=GLM_VARS.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  ENCODIG/DECODIG - BACKPROJECT ONLY USED ICs        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



