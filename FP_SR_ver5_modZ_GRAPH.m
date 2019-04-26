%% ModZ + subtract

%Channel (background transparent)

close all
%% All Trial

for currentTrial= 1:SessionData.nTrials
    TimeStamp_AllTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
    TimStamp_AllTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
    TimeStamp_AllTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
    TimeStamp_AllTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
    
end

TimeStamp_AllTrial.CenterLitOn(isnan(TimeStamp_AllTrial.CenterLitOn))=[];
TimeStamp_AllTrial.Initiation(isnan(TimeStamp_AllTrial.Initiation))=[];
TimeStamp_AllTrial.Choice(isnan(TimeStamp_AllTrial.Choice))=[];
TimeStamp_AllTrial.RewardRetrieval(isnan(TimeStamp_AllTrial.RewardRetrieval))=[];


% Rewarded Trial

for currentTrial= 1:SessionData.nTrials
    if SessionData.Reward(currentTrial) == 1
        Reward_Trial(currentTrial) = currentTrial;
    else
        Reward_Trial(currentTrial) = nan;
    end
    
end
Reward_Trial(isnan(Reward_Trial))=[];




for currentTrial= 1:SessionData.nTrials
    if SessionData.Reward(currentTrial) == 1
        TimeStamp_RewardedTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
        TimeStamp_RewardedTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
        TimeStamp_RewardedTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
        TimeStamp_RewardedTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
    else
        TimeStamp_RewardedTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_RewardedTrial.Initiation(currentTrial) = nan;
        TimeStamp_RewardedTrial.Choice(currentTrial) = nan;
        TimeStamp_RewardedTrial.RewardRetrieval(currentTrial) = nan;
    end
end

TimeStamp_RewardedTrial.CenterLitOn(isnan(TimeStamp_RewardedTrial.CenterLitOn))=[];
TimeStamp_RewardedTrial.Initiation(isnan(TimeStamp_RewardedTrial.Initiation))=[];
TimeStamp_RewardedTrial.Choice(isnan(TimeStamp_RewardedTrial.Choice))=[];
TimeStamp_RewardedTrial.RewardRetrieval(isnan(TimeStamp_RewardedTrial.RewardRetrieval))=[];

% Unrewarded Trial
for currentTrial= 1:SessionData.nTrials
    if SessionData.Reward(currentTrial) == -1
        Unreward_Trial(currentTrial) = currentTrial;
    else
        Unreward_Trial(currentTrial) = nan;
    end
    
end
Unreward_Trial(isnan(Unreward_Trial))=[];


for currentTrial= 1:SessionData.nTrials
    if SessionData.Reward(currentTrial) == -1
        TimeStamp_UnrewardedTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
        TimeStamp_UnrewardedTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
        TimeStamp_UnrewardedTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
        TimeStamp_UnrewardedTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
    else
        TimeStamp_UnrewardedTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_UnrewardedTrial.Initiation(currentTrial) = nan;
        TimeStamp_UnrewardedTrial.Choice(currentTrial) = nan;
        TimeStamp_UnrewardedTrial.RewardRetrieval(currentTrial) = nan;
    end
end

TimeStamp_UnrewardedTrial.CenterLitOn(isnan(TimeStamp_UnrewardedTrial.CenterLitOn))=[];
TimeStamp_UnrewardedTrial.Initiation(isnan(TimeStamp_UnrewardedTrial.Initiation))=[];
TimeStamp_UnrewardedTrial.Choice(isnan(TimeStamp_UnrewardedTrial.Choice))=[];
TimeStamp_UnrewardedTrial.RewardRetrieval(isnan(TimeStamp_UnrewardedTrial.RewardRetrieval))=[];

% Previous Reward history

for currentTrial= 2:SessionData.nTrials
    if SessionData.Reward(currentTrial-1) == 1
        TimeStamp_RewardAfterTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
        TimeStamp_RewardAfterTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
        TimeStamp_RewardAfterTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
        TimeStamp_RewardAfterTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
    else
        TimeStamp_RewardAfterTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_RewardAfterTrial.Initiation(currentTrial) = nan;
        TimeStamp_RewardAfterTrial.Choice(currentTrial) = nan;
        TimeStamp_RewardAfterTrial.RewardRetrieval(currentTrial) = nan;
    end
end

TimeStamp_RewardAfterTrial.CenterLitOn(isnan(TimeStamp_RewardAfterTrial.CenterLitOn))=[];
TimeStamp_RewardAfterTrial.Initiation(isnan(TimeStamp_RewardAfterTrial.Initiation))=[];
TimeStamp_RewardAfterTrial.Choice(isnan(TimeStamp_RewardAfterTrial.Choice))=[];
TimeStamp_RewardAfterTrial.RewardRetrieval(isnan(TimeStamp_RewardAfterTrial.RewardRetrieval))=[];

% Previous Unreward history

for currentTrial= 2:SessionData.nTrials
    if ~(SessionData.Reward(currentTrial-1) == 1)
        
        TimeStamp_UnrewardAfterTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
        TimeStamp_UnrewardAfterTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
        TimeStamp_UnrewardAfterTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
        TimeStamp_UnrewardAfterTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
        
    else
        
        TimeStamp_UnrewardAfterTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_UnrewardAfterTrial.Initiation(currentTrial) = nan;
        TimeStamp_UnrewardAfterTrial.Choice(currentTrial) = nan;
        TimeStamp_UnrewardAfterTrial.RewardRetrieval(currentTrial) = nan;
        
    end
end

TimeStamp_UnrewardAfterTrial.CenterLitOn(isnan(TimeStamp_UnrewardAfterTrial.CenterLitOn))=[];
TimeStamp_UnrewardAfterTrial.Initiation(isnan(TimeStamp_UnrewardAfterTrial.Initiation))=[];
TimeStamp_UnrewardAfterTrial.Choice(isnan(TimeStamp_UnrewardAfterTrial.Choice))=[];
TimeStamp_UnrewardAfterTrial.RewardRetrieval(isnan(TimeStamp_UnrewardAfterTrial.RewardRetrieval))=[];

% Unrewarded(preRewarded)

for currentTrial= 1:SessionData.nTrials
    if currentTrial == 1
        TrialID_RewardAfterUnrewardTrial(currentTrial) = nan;
    else
        if SessionData.Reward(currentTrial-1) == 1 && ~(SessionData.Reward(currentTrial) == 1)
            TrialID_RewardAfterUnrewardTrial(currentTrial) = currentTrial;
        else
            TrialID_RewardAfterUnrewardTrial(currentTrial) = nan;
        end
    end
    
end
TrialID_RewardAfterUnrewardTrial(isnan(TrialID_RewardAfterUnrewardTrial))=[];

% Trial 1

for currentTrial= 1:SessionData.nTrials
    if currentTrial == 1
        TimeStamp_RewardAfterUnrewardTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_RewardAfterUnrewardTrial.Initiation(currentTrial) = nan;
        TimeStamp_RewardAfterUnrewardTrial.Choice(currentTrial) = nan;
        TimeStamp_RewardAfterUnrewardTrial.RewardRetrieval(currentTrial) = nan;
    else
        if SessionData.Reward(currentTrial-1) == 1 && ~(SessionData.Reward(currentTrial) == 1)
            TimeStamp_RewardAfterUnrewardTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
            TimeStamp_RewardAfterUnrewardTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
            TimeStamp_RewardAfterUnrewardTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
            TimeStamp_RewardAfterUnrewardTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
        else
            TimeStamp_RewardAfterUnrewardTrial.CenterLitOn(currentTrial) = nan;
            TimeStamp_RewardAfterUnrewardTrial.Initiation(currentTrial) = nan;
            TimeStamp_RewardAfterUnrewardTrial.Choice(currentTrial) = nan;
            TimeStamp_RewardAfterUnrewardTrial.RewardRetrieval(currentTrial) = nan;
        end
    end
end

TimeStamp_RewardAfterUnrewardTrial.CenterLitOn(isnan(TimeStamp_RewardAfterUnrewardTrial.CenterLitOn))=[];
TimeStamp_RewardAfterUnrewardTrial.Initiation(isnan(TimeStamp_RewardAfterUnrewardTrial.Initiation))=[];
TimeStamp_RewardAfterUnrewardTrial.Choice(isnan(TimeStamp_RewardAfterUnrewardTrial.Choice))=[];
TimeStamp_RewardAfterUnrewardTrial.RewardRetrieval(isnan(TimeStamp_RewardAfterUnrewardTrial.RewardRetrieval))=[];

% Rewarded (previously unrewarded)
for currentTrial= 1:SessionData.nTrials
    if currentTrial == 1
        TrialID_RewTrialPreUnrewarded(currentTrial) = nan;
    else
        if ~(SessionData.Reward(currentTrial-1) == 1) && SessionData.Reward(currentTrial) == 1
            TrialID_RewTrialPreUnrewarded(currentTrial) = currentTrial;
        else
            TrialID_RewTrialPreUnrewarded(currentTrial) = nan;
        end
    end
    
end
TrialID_RewTrialPreUnrewarded(isnan(TrialID_RewTrialPreUnrewarded))=[];


for currentTrial= 1:SessionData.nTrials
    if currentTrial == 1
        TimeStamp_UnrewardAfterRewardTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_UnrewardAfterRewardTrial.Initiation(currentTrial) = nan;
        TimeStamp_UnrewardAfterRewardTrial.Choice(currentTrial) = nan;
        TimeStamp_UnrewardAfterRewardTrial.RewardRetrieval(currentTrial) = nan;
    else
        if ~(SessionData.Reward(currentTrial-1) == 1) && SessionData.Reward(currentTrial) == 1
            TimeStamp_UnrewardAfterRewardTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
            TimeStamp_UnrewardAfterRewardTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
            TimeStamp_UnrewardAfterRewardTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
            TimeStamp_UnrewardAfterRewardTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
        else
            TimeStamp_UnrewardAfterRewardTrial.CenterLitOn(currentTrial) = nan;
            TimeStamp_UnrewardAfterRewardTrial.Initiation(currentTrial) = nan;
            TimeStamp_UnrewardAfterRewardTrial.Choice(currentTrial) = nan;
            TimeStamp_UnrewardAfterRewardTrial.RewardRetrieval(currentTrial) = nan;
        end
    end
end

TimeStamp_UnrewardAfterRewardTrial.CenterLitOn(isnan(TimeStamp_UnrewardAfterRewardTrial.CenterLitOn))=[];
TimeStamp_UnrewardAfterRewardTrial.Initiation(isnan(TimeStamp_UnrewardAfterRewardTrial.Initiation))=[];
TimeStamp_UnrewardAfterRewardTrial.Choice(isnan(TimeStamp_UnrewardAfterRewardTrial.Choice))=[];
TimeStamp_UnrewardAfterRewardTrial.RewardRetrieval(isnan(TimeStamp_UnrewardAfterRewardTrial.RewardRetrieval))=[];

% Rewarded Trial (Previously Rewarded)(TrialLookup = 8)

for currentTrial= 1:SessionData.nTrials
    if currentTrial == 1
        TimeStamp_RewardAfterRewardTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_RewardAfterRewardTrial.Initiation(currentTrial) = nan;
        TimeStamp_RewardAfterRewardTrial.Choice(currentTrial) = nan;
        TimeStamp_RewardAfterRewardTrial.RewardRetrieval(currentTrial) = nan;
    else
        if SessionData.Reward(currentTrial-1) == 1 && SessionData.Reward(currentTrial) == 1
            TimeStamp_RewardAfterRewardTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
            TimeStamp_RewardAfterRewardTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
            TimeStamp_RewardAfterRewardTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
            TimeStamp_RewardAfterRewardTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
        else
            TimeStamp_RewardAfterRewardTrial.CenterLitOn(currentTrial) = nan;
            TimeStamp_RewardAfterRewardTrial.Initiation(currentTrial) = nan;
            TimeStamp_RewardAfterRewardTrial.Choice(currentTrial) = nan;
            TimeStamp_RewardAfterRewardTrial.RewardRetrieval(currentTrial) = nan;
        end
    end
end

TimeStamp_RewardAfterRewardTrial.CenterLitOn(isnan(TimeStamp_RewardAfterRewardTrial.CenterLitOn))=[];
TimeStamp_RewardAfterRewardTrial.Initiation(isnan(TimeStamp_RewardAfterRewardTrial.Initiation))=[];
TimeStamp_RewardAfterRewardTrial.Choice(isnan(TimeStamp_RewardAfterRewardTrial.Choice))=[];
TimeStamp_RewardAfterRewardTrial.RewardRetrieval(isnan(TimeStamp_RewardAfterRewardTrial.RewardRetrieval))=[];

% Unrewarded Trial (Previously Unrewarded)(TrialLookup = 9)

for currentTrial= 1:SessionData.nTrials
    if currentTrial == 1
        TimeStamp_UnrewardAfterUnrewardTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_UnrewardAfterUnrewardTrial.Initiation(currentTrial) = nan;
        TimeStamp_UnrewardAfterUnrewardTrial.Choice(currentTrial) = nan;
        TimeStamp_UnrewardAfterUnrewardTrial.RewardRetrieval(currentTrial) = nan;
    else
        if SessionData.Reward(currentTrial-1) == -1 && SessionData.Reward(currentTrial) == -1
            TimeStamp_UnrewardAfterUnrewardTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
            TimeStamp_UnrewardAfterUnrewardTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
            TimeStamp_UnrewardAfterUnrewardTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
            TimeStamp_UnrewardAfterUnrewardTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
        else
            TimeStamp_UnrewardAfterUnrewardTrial.CenterLitOn(currentTrial) = nan;
            TimeStamp_UnrewardAfterUnrewardTrial.Initiation(currentTrial) = nan;
            TimeStamp_UnrewardAfterUnrewardTrial.Choice(currentTrial) = nan;
            TimeStamp_UnrewardAfterUnrewardTrial.RewardRetrieval(currentTrial) = nan;
        end
    end
end

TimeStamp_UnrewardAfterUnrewardTrial.CenterLitOn(isnan(TimeStamp_UnrewardAfterUnrewardTrial.CenterLitOn))=[];
TimeStamp_UnrewardAfterUnrewardTrial.Initiation(isnan(TimeStamp_UnrewardAfterUnrewardTrial.Initiation))=[];
TimeStamp_UnrewardAfterUnrewardTrial.Choice(isnan(TimeStamp_UnrewardAfterUnrewardTrial.Choice))=[];
TimeStamp_UnrewardAfterUnrewardTrial.RewardRetrieval(isnan(TimeStamp_UnrewardAfterUnrewardTrial.RewardRetrieval))=[];

% BlockEnd Trial

for currentTrial= 1:SessionData.nTrials
    if SessionData.Tag.Trial_BlockEnd(currentTrial) == 1
        TimeStamp_BlockEndTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
        TimeStamp_BlockEndTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
        TimeStamp_BlockEndTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
        TimeStamp_BlockEndTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
    else
        TimeStamp_BlockEndTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_BlockEndTrial.Initiation(currentTrial) = nan;
        TimeStamp_BlockEndTrial.Choice(currentTrial) = nan;
        TimeStamp_BlockEndTrial.RewardRetrieval(currentTrial) = nan;
    end
end

TimeStamp_BlockEndTrial.CenterLitOn(isnan(TimeStamp_BlockEndTrial.CenterLitOn))=[];
TimeStamp_BlockEndTrial.Initiation(isnan(TimeStamp_BlockEndTrial.Initiation))=[];
TimeStamp_BlockEndTrial.Choice(isnan(TimeStamp_BlockEndTrial.Choice))=[];
TimeStamp_BlockEndTrial.RewardRetrieval(isnan(TimeStamp_BlockEndTrial.RewardRetrieval))=[];

% BlockStart Trial
for currentTrial= 1:SessionData.nTrials
    if SessionData.Tag.Trial_BlockStart_1(currentTrial) == 1
        TimeStamp_BlockStartTrial_1.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
        TimeStamp_BlockStartTrial_1.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
        TimeStamp_BlockStartTrial_1.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
        TimeStamp_BlockStartTrial_1.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
    else
        TimeStamp_BlockStartTrial_1.CenterLitOn(currentTrial) = nan;
        TimeStamp_BlockStartTrial_1.Initiation(currentTrial) = nan;
        TimeStamp_BlockStartTrial_1.Choice(currentTrial) = nan;
        TimeStamp_BlockStartTrial_1.RewardRetrieval(currentTrial) = nan;
    end
end

TimeStamp_BlockStartTrial_1.CenterLitOn(isnan(TimeStamp_BlockStartTrial_1.CenterLitOn))=[];
TimeStamp_BlockStartTrial_1.Initiation(isnan(TimeStamp_BlockStartTrial_1.Initiation))=[];
TimeStamp_BlockStartTrial_1.Choice(isnan(TimeStamp_BlockStartTrial_1.Choice))=[];
TimeStamp_BlockStartTrial_1.RewardRetrieval(isnan(TimeStamp_BlockStartTrial_1.RewardRetrieval))=[];
% BlockStart Trial_2
for currentTrial= 1:SessionData.nTrials
    if SessionData.Tag.Trial_BlockStart_2(currentTrial) == 1
        TimeStamp_BlockStartTrial_2.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
        TimeStamp_BlockStartTrial_2.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
        TimeStamp_BlockStartTrial_2.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
        TimeStamp_BlockStartTrial_2.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
    else
        TimeStamp_BlockStartTrial_2.CenterLitOn(currentTrial) = nan;
        TimeStamp_BlockStartTrial_2.Initiation(currentTrial) = nan;
        TimeStamp_BlockStartTrial_2.Choice(currentTrial) = nan;
        TimeStamp_BlockStartTrial_2.RewardRetrieval(currentTrial) = nan;
    end
end

TimeStamp_BlockStartTrial_2.CenterLitOn(isnan(TimeStamp_BlockStartTrial_2.CenterLitOn))=[];
TimeStamp_BlockStartTrial_2.Initiation(isnan(TimeStamp_BlockStartTrial_2.Initiation))=[];
TimeStamp_BlockStartTrial_2.Choice(isnan(TimeStamp_BlockStartTrial_2.Choice))=[];
TimeStamp_BlockStartTrial_2.RewardRetrieval(isnan(TimeStamp_BlockStartTrial_2.RewardRetrieval))=[];

% Manual input

for currentTrial= 1:SessionData.nTrials
    if sum(currentTrial == Manual_Trial) == 1
        TimeStamp_TEST.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
        TimeStamp_TEST.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
        TimeStamp_TEST.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
        TimeStamp_TEST.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
    else
        TimeStamp_TEST.CenterLitOn(currentTrial) = nan;
        TimeStamp_TEST.Initiation(currentTrial) = nan;
        TimeStamp_TEST.Choice(currentTrial) = nan;
        TimeStamp_TEST.RewardRetrieval(currentTrial) = nan;
    end
end

TimeStamp_TEST.CenterLitOn(isnan(TimeStamp_TEST.CenterLitOn))=[];
TimeStamp_TEST.Initiation(isnan(TimeStamp_TEST.Initiation))=[];
TimeStamp_TEST.Choice(isnan(TimeStamp_TEST.Choice))=[];
TimeStamp_TEST.RewardRetrieval(isnan(TimeStamp_TEST.RewardRetrieval))=[];


%% Trial Selection for PSTH
Trial_List = {'All' 'Rewarded','Unrewarded', ...
    'Trial_preRew', 'Trial_preUnrew', ...
    'UnrewTrial_preRew', 'RewaTrial_preUnrewd',...
    'RewTrial_preRew','UnrewaTrial_preUnrew',...
    'BlockEndTrial', 'BlockStartTrial_1', 'BlockStartTrial_2', 'TEST'
    };

% make  arrays based on behavior timestamps, then average rows into a mean vector for plotting

[TrialLookup,tf] = listdlg('ListString',Trial_List);

if TrialLookup == 1
    TimeStamp_Choice_Analysis = TimeStamp_AllTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_AllTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_AllTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_AllTrial.RewardRetrieval;
elseif TrialLookup == 2
    TimeStamp_Choice_Analysis = TimeStamp_RewardedTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_RewardedTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_RewardedTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_RewardedTrial.RewardRetrieval;
elseif TrialLookup == 3
    TimeStamp_Choice_Analysis = TimeStamp_UnrewardedTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_UnrewardedTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_UnrewardedTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_UnrewardedTrial.RewardRetrieval;
elseif TrialLookup == 4
    TimeStamp_Choice_Analysis = TimeStamp_RewardAfterTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_RewardAfterTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_RewardAfterTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_RewardAfterTrial.RewardRetrieval;
elseif TrialLookup == 5
    TimeStamp_Choice_Analysis = TimeStamp_UnrewardAfterTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_UnrewardAfterTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_UnrewardAfterTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_UnrewardAfterTrial.RewardRetrieval;
elseif TrialLookup == 6
    TimeStamp_Choice_Analysis = TimeStamp_RewardAfterUnrewardTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_RewardAfterUnrewardTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_RewardAfterUnrewardTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_RewardAfterUnrewardTrial.RewardRetrieval;
elseif TrialLookup == 7
    TimeStamp_Choice_Analysis = TimeStamp_UnrewardAfterRewardTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_UnrewardAfterRewardTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_UnrewardAfterRewardTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_UnrewardAfterRewardTrial.RewardRetrieval;
elseif TrialLookup == 8
    TimeStamp_Choice_Analysis = TimeStamp_RewardAfterRewardTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_RewardAfterRewardTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_RewardAfterRewardTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_RewardAfterRewardTrial.RewardRetrieval;
elseif TrialLookup ==9
    TimeStamp_Choice_Analysis = TimeStamp_UnrewardAfterUnrewardTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_UnrewardAfterUnrewardTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_UnrewardAfterUnrewardTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_UnrewardAfterUnrewardTrial.RewardRetrieval;
elseif TrialLookup == 10
    TimeStamp_Choice_Analysis = TimeStamp_BlockEndTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_BlockEndTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_BlockEndTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_BlockEndTrial.RewardRetrieval;
elseif TrialLookup == 11
    TimeStamp_Choice_Analysis = TimeStamp_BlockStartTrial_1.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_BlockStartTrial_1.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_BlockStartTrial_1.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_BlockStartTrial_1.RewardRetrieval;
elseif TrialLookup == 12
    TimeStamp_Choice_Analysis = TimeStamp_BlockStartTrial_2.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_BlockStartTrial_2.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_BlockStartTrial_2.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_BlockStartTrial_2.RewardRetrieval;
elseif TrialLookup == 13
    TimeStamp_Choice_Analysis = TimeStamp_TEST.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_TEST.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_TEST.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_TEST.RewardRetrieval;
end

%% Graph setting - Together

% x-Axis Analysis window
nSecPrev = 10; %change to make different window
nSecPost = 10;

% x-Axis View window
nSecPrevWindow = 1; %change to make different window
nSecPostWindow = 1;

% y-axis
ymin = -1.5;
ymax = 1.5;

% convert seconds to TDT timestamps
nTsPrev = round (nSecPrev * samplingRate);
nTsPost = round (nSecPost * samplingRate);

totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;

% make raw_PSTH for Correct Press
nCHOICE = length(TimeStamp_Choice_Analysis);
PsthArray_CHOICE_ch1 = NaN(nCHOICE,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_CHOICE_ch2 = NaN(nCHOICE,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nCHOICE
    thisTime = TimeStamp_Choice_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_CHOICE_ch1(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_CHOICE_ch2(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch2, thisIndex, nTsPrev, nTsPost);
end

FIGURE_CHOICE.Err_CHOICE_ch1 = (nanstd(PsthArray_CHOICE_ch1))/sqrt(size(PsthArray_CHOICE_ch1,1));
FIGURE_CHOICE.Psth_CHOICE_ch1 = nanmean(PsthArray_CHOICE_ch1);
FIGURE_CHOICE.Err_Positive_ch1 = FIGURE_CHOICE.Psth_CHOICE_ch1 + FIGURE_CHOICE.Err_CHOICE_ch1;
FIGURE_CHOICE.Err_Negative_ch1 = FIGURE_CHOICE.Psth_CHOICE_ch1 - FIGURE_CHOICE.Err_CHOICE_ch1;

FIGURE_CHOICE.Err_CHOICE_ch2 = (nanstd(PsthArray_CHOICE_ch2))/sqrt(size(PsthArray_CHOICE_ch2,1));
FIGURE_CHOICE.Psth_CHOICE_ch2 = nanmean(PsthArray_CHOICE_ch2);
FIGURE_CHOICE.Err_Positive_ch2 = FIGURE_CHOICE.Psth_CHOICE_ch2 + FIGURE_CHOICE.Err_CHOICE_ch2;
FIGURE_CHOICE.Err_Negative_ch2 = FIGURE_CHOICE.Psth_CHOICE_ch2 - FIGURE_CHOICE.Err_CHOICE_ch2;


F(11) = figure('color', [1 1 1]); subplot(1,5,3)
hold on

o = 0.5;
r1 = 255;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);
fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch1, fliplr(FIGURE_CHOICE.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch2, fliplr(FIGURE_CHOICE.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('Choice', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);


legend BOXOFF;
xlim ([(-1 * nSecPrevWindow),nSecPostWindow]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
line('XData', [(-1 * nSecPrevWindow),nSecPostWindow], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
legend ([h1 h2], 'pDMS', 'aDMS');
legend('off')

set(gca, 'layer', 'top');
set(gcf, 'Position', [0 600 1600 400])


% make raw_PSTH for Initiation
nINITIATION = length(TimeStamp_Initiation_Analysis);
PsthArray_INITIATION_ch1 = NaN(nINITIATION,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_INITIATION_ch2 = NaN(nINITIATION,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nINITIATION
    thisTime = TimeStamp_Initiation_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_INITIATION_ch1(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_INITIATION_ch2(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch2, thisIndex, nTsPrev, nTsPost);
end

FIGURE_INITIATION.Err_INITIATION_ch1 = (nanstd(PsthArray_INITIATION_ch1))/sqrt(size(PsthArray_INITIATION_ch1,1));
FIGURE_INITIATION.Psth_INITIATION_ch1 = nanmean(PsthArray_INITIATION_ch1);
FIGURE_INITIATION.Err_Positive_ch1 = FIGURE_INITIATION.Psth_INITIATION_ch1 + FIGURE_INITIATION.Err_INITIATION_ch1;
FIGURE_INITIATION.Err_Negative_ch1 = FIGURE_INITIATION.Psth_INITIATION_ch1 - FIGURE_INITIATION.Err_INITIATION_ch1;

FIGURE_INITIATION.Err_INITIATION_ch2 = (nanstd(PsthArray_INITIATION_ch2))/sqrt(size(PsthArray_INITIATION_ch2,1));
FIGURE_INITIATION.Psth_INITIATION_ch2 = nanmean(PsthArray_INITIATION_ch2);
FIGURE_INITIATION.Err_Positive_ch2 = FIGURE_INITIATION.Psth_INITIATION_ch2 + FIGURE_INITIATION.Err_INITIATION_ch2;
FIGURE_INITIATION.Err_Negative_ch2 = FIGURE_INITIATION.Psth_INITIATION_ch2 - FIGURE_INITIATION.Err_INITIATION_ch2;


subplot(1,5,2);
hold on

o = 0.5;
r1 = 255;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);
fill([timeAxis, fliplr(timeAxis)],[FIGURE_INITIATION.Err_Positive_ch1, fliplr(FIGURE_INITIATION.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,FIGURE_INITIATION.Psth_INITIATION_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);


r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[FIGURE_INITIATION.Err_Positive_ch2, fliplr(FIGURE_INITIATION.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,FIGURE_INITIATION.Psth_INITIATION_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('Initiation', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([(-1 * nSecPrevWindow),nSecPostWindow]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
line('XData', [(-1 * nSecPrevWindow),nSecPostWindow], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
legend ([h1 h2], 'pDMS', 'aDMS');
legend('off')

set(gca, 'layer', 'top');
%set(gcf, 'Position', [500 600 2000 400])


% make raw_PSTH for Liton
nCENTERLITON = length(TimeStamp_CenterLitOn_Analysis);
PsthArray_CENTERLITON_ch1 = NaN(nCENTERLITON,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_CENTERLITON_ch2 = NaN(nCENTERLITON,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nCENTERLITON
    thisTime = TimeStamp_CenterLitOn_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_CENTERLITON_ch1(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_CENTERLITON_ch2(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch2, thisIndex, nTsPrev, nTsPost);
end


FIGURE_CENTERLITON.Err_CENTERLITON_ch1 = (nanstd(PsthArray_CENTERLITON_ch1))/sqrt(size(PsthArray_CENTERLITON_ch1,1));
FIGURE_CENTERLITON.Psth_CENTERLITON_ch1 = nanmean(PsthArray_CENTERLITON_ch1);
FIGURE_CENTERLITON.Err_Positive_ch1 = FIGURE_CENTERLITON.Psth_CENTERLITON_ch1 + FIGURE_CENTERLITON.Err_CENTERLITON_ch1;
FIGURE_CENTERLITON.Err_Negative_ch1 = FIGURE_CENTERLITON.Psth_CENTERLITON_ch1 - FIGURE_CENTERLITON.Err_CENTERLITON_ch1;

FIGURE_CENTERLITON.Err_CENTERLITON_ch2 = (nanstd(PsthArray_CENTERLITON_ch2))/sqrt(size(PsthArray_CENTERLITON_ch2,1));
FIGURE_CENTERLITON.Psth_CENTERLITON_ch2 = nanmean(PsthArray_CENTERLITON_ch2);
FIGURE_CENTERLITON.Err_Positive_ch2 = FIGURE_CENTERLITON.Psth_CENTERLITON_ch2 + FIGURE_CENTERLITON.Err_CENTERLITON_ch2;
FIGURE_CENTERLITON.Err_Negative_ch2 = FIGURE_CENTERLITON.Psth_CENTERLITON_ch2 - FIGURE_CENTERLITON.Err_CENTERLITON_ch2;



%F(13) = figure('color', [1 1 1]);
subplot(1,5,1)
hold on

o = 0.5;
r1 = 255;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);
fill([timeAxis, fliplr(timeAxis)],[FIGURE_CENTERLITON.Err_Positive_ch1, fliplr(FIGURE_CENTERLITON.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,FIGURE_CENTERLITON.Psth_CENTERLITON_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[FIGURE_CENTERLITON.Err_Positive_ch2, fliplr(FIGURE_CENTERLITON.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,FIGURE_CENTERLITON.Psth_CENTERLITON_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('CenterLightON', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([(-1 * nSecPrevWindow),nSecPostWindow]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
line('XData', [(-1 * nSecPrevWindow),nSecPostWindow], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
legend ([h1 h2], 'pDMS', 'aDMS');
legend('off')

set(gca, 'layer', 'top');
%set(gcf, 'Position', [0 600 500 400])

nREWARD = length(TimeStamp_RewardRetrieval_Analysis);
PsthArray_REWARD_ch1 = NaN(nREWARD,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_REWARD_ch2 = NaN(nREWARD,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nREWARD
    thisTime = TimeStamp_RewardRetrieval_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_REWARD_ch1(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_REWARD_ch2(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch2, thisIndex, nTsPrev, nTsPost);
end


FIGURE_REWARD.Err_REWARD_ch1 = (nanstd(PsthArray_REWARD_ch1))/sqrt(size(PsthArray_REWARD_ch1,1));
FIGURE_REWARD.Psth_REWARD_ch1 = nanmean(PsthArray_REWARD_ch1);
FIGURE_REWARD.Err_Positive_ch1 = FIGURE_REWARD.Psth_REWARD_ch1 + FIGURE_REWARD.Err_REWARD_ch1;
FIGURE_REWARD.Err_Negative_ch1 = FIGURE_REWARD.Psth_REWARD_ch1 - FIGURE_REWARD.Err_REWARD_ch1;

FIGURE_REWARD.Err_REWARD_ch2 = (nanstd(PsthArray_REWARD_ch2))/sqrt(size(PsthArray_REWARD_ch2,1));
FIGURE_REWARD.Psth_REWARD_ch2 = nanmean(PsthArray_REWARD_ch2);
FIGURE_REWARD.Err_Positive_ch2 = FIGURE_REWARD.Psth_REWARD_ch2 + FIGURE_REWARD.Err_REWARD_ch2;
FIGURE_REWARD.Err_Negative_ch2 = FIGURE_REWARD.Psth_REWARD_ch2 - FIGURE_REWARD.Err_REWARD_ch2;

%F(14) = figure('color', [1 1 1]);
subplot(1,5,4)
hold on

o = 0.5;
r1 = 255;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);
fill([timeAxis, fliplr(timeAxis)],[FIGURE_REWARD.Err_Positive_ch1, fliplr(FIGURE_REWARD.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,FIGURE_REWARD.Psth_REWARD_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[FIGURE_REWARD.Err_Positive_ch2, fliplr(FIGURE_REWARD.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,FIGURE_REWARD.Psth_REWARD_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('RewardRetrieval', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([(-1 * nSecPrevWindow),nSecPostWindow]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
line('XData', [(-1 * nSecPrevWindow),nSecPostWindow], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
legend ([h1 h2], 'pDMS', 'aDMS');
legend('off')

set(gca, 'layer', 'top');
%set(gcf, 'Position', [500 0 500 400])

% make raw_PSTH for next ISI
nCHOICE = length(TimeStamp_Choice_Analysis);
PsthArray_CHOICE_ch1 = NaN(nCHOICE,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_CHOICE_ch2 = NaN(nCHOICE,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nCHOICE
    thisTime = TimeStamp_Choice_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_CHOICE_ch1(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_CHOICE_ch2(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch2, thisIndex, nTsPrev, nTsPost);
end

FIGURE_CHOICE.Err_CHOICE_ch1 = (nanstd(PsthArray_CHOICE_ch1))/sqrt(size(PsthArray_CHOICE_ch1,1));
FIGURE_CHOICE.Psth_CHOICE_ch1 = nanmean(PsthArray_CHOICE_ch1);
FIGURE_CHOICE.Err_Positive_ch1 = FIGURE_CHOICE.Psth_CHOICE_ch1 + FIGURE_CHOICE.Err_CHOICE_ch1;
FIGURE_CHOICE.Err_Negative_ch1 = FIGURE_CHOICE.Psth_CHOICE_ch1 - FIGURE_CHOICE.Err_CHOICE_ch1;

FIGURE_CHOICE.Err_CHOICE_ch2 = (nanstd(PsthArray_CHOICE_ch2))/sqrt(size(PsthArray_CHOICE_ch2,1));
FIGURE_CHOICE.Psth_CHOICE_ch2 = nanmean(PsthArray_CHOICE_ch2);
FIGURE_CHOICE.Err_Positive_ch2 = FIGURE_CHOICE.Psth_CHOICE_ch2 + FIGURE_CHOICE.Err_CHOICE_ch2;
FIGURE_CHOICE.Err_Negative_ch2 = FIGURE_CHOICE.Psth_CHOICE_ch2 - FIGURE_CHOICE.Err_CHOICE_ch2;

subplot(1,5,5)
hold on

o = 0.5;
r1 = 255;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);
fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch1, fliplr(FIGURE_CHOICE.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch2, fliplr(FIGURE_CHOICE.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('ISI', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([0,8]);
xticks([0 3 8]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
line('XData', [0,8], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [3 3], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')

legend ([h1 h2], 'pDMS', 'aDMS');

set(gca, 'layer', 'top');
%set(gcf, 'Position', [0 0 2000 400])

%% Graph setting - Ch1

% make raw_PSTH for Correct Press
nCHOICE = length(TimeStamp_Choice_Analysis);
PsthArray_CHOICE_ch1 = NaN(nCHOICE,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_CHOICE_ch2 = NaN(nCHOICE,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nCHOICE
    thisTime = TimeStamp_Choice_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_CHOICE_ch1(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_CHOICE_ch2(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch2, thisIndex, nTsPrev, nTsPost);
end

FIGURE_CHOICE.Err_CHOICE_ch1 = (nanstd(PsthArray_CHOICE_ch1))/sqrt(size(PsthArray_CHOICE_ch1,1));
FIGURE_CHOICE.Psth_CHOICE_ch1 = nanmean(PsthArray_CHOICE_ch1);
FIGURE_CHOICE.Err_Positive_ch1 = FIGURE_CHOICE.Psth_CHOICE_ch1 + FIGURE_CHOICE.Err_CHOICE_ch1;
FIGURE_CHOICE.Err_Negative_ch1 = FIGURE_CHOICE.Psth_CHOICE_ch1 - FIGURE_CHOICE.Err_CHOICE_ch1;

FIGURE_CHOICE.Err_CHOICE_ch2 = (nanstd(PsthArray_CHOICE_ch2))/sqrt(size(PsthArray_CHOICE_ch2,1));
FIGURE_CHOICE.Psth_CHOICE_ch2 = nanmean(PsthArray_CHOICE_ch2);
FIGURE_CHOICE.Err_Positive_ch2 = FIGURE_CHOICE.Psth_CHOICE_ch2 + FIGURE_CHOICE.Err_CHOICE_ch2;
FIGURE_CHOICE.Err_Negative_ch2 = FIGURE_CHOICE.Psth_CHOICE_ch2 - FIGURE_CHOICE.Err_CHOICE_ch2;


F(11) = figure('color', [1 1 1]); subplot(1,5,3)
hold on

o = 0.5;
r1 = 255;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);
fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch1, fliplr(FIGURE_CHOICE.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

% r2 = 0;
% g2 = 0;
% b2 = 255;
% rgb2_o = opacity (o, r2, g2, b2);
%
% fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch2, fliplr(FIGURE_CHOICE.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
%h2= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('Choice', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([(-1 * nSecPrevWindow),nSecPostWindow]);
ylim ([ymin, ymax]);

set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
set(gca, 'Color', 'none')
line('XData', [(-1 * nSecPrevWindow),nSecPostWindow], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')

%legend ([h1 h2], 'pDMS', 'aDMS');
legend('off')

set(gca, 'layer', 'top');
set(gcf, 'Position', [100 400 1600 400])

%
% h1f=fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch1, fliplr(FIGURE_CHOICE.Err_Negative_ch1)], 'g', 'EdgeColor', 'none');
% set(h1f,'facealpha',.5);
% h1= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch1,'Color','r', 'Linewidth', 3);

% PsthArray_CORRECTPRESS_Block.Block1 = PsthArray_CORRECTPRESS(1:50, :);
% PsthArray_CORRECTPRESS_Block.Block2 = PsthArray_CORRECTPRESS(51:100, :);
% PsthArray_CORRECTPRESS_Block.Block3 = PsthArray_CORRECTPRESS(101:150, :);
% %PsthArray_CORRECTPRESS_Block.Block4 = PsthArray_CORRECTPRESS(151:200, :);



% make raw_PSTH for Initiation
nINITIATION = length(TimeStamp_Initiation_Analysis);
PsthArray_INITIATION_ch1 = NaN(nINITIATION,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_INITIATION_ch2 = NaN(nINITIATION,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nINITIATION
    thisTime = TimeStamp_Initiation_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_INITIATION_ch1(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_INITIATION_ch2(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch2, thisIndex, nTsPrev, nTsPost);
end

FIGURE_INITIATION.Err_INITIATION_ch1 = (nanstd(PsthArray_INITIATION_ch1))/sqrt(size(PsthArray_INITIATION_ch1,1));
FIGURE_INITIATION.Psth_INITIATION_ch1 = nanmean(PsthArray_INITIATION_ch1);
FIGURE_INITIATION.Err_Positive_ch1 = FIGURE_INITIATION.Psth_INITIATION_ch1 + FIGURE_INITIATION.Err_INITIATION_ch1;
FIGURE_INITIATION.Err_Negative_ch1 = FIGURE_INITIATION.Psth_INITIATION_ch1 - FIGURE_INITIATION.Err_INITIATION_ch1;

FIGURE_INITIATION.Err_INITIATION_ch2 = (nanstd(PsthArray_INITIATION_ch2))/sqrt(size(PsthArray_INITIATION_ch2,1));
FIGURE_INITIATION.Psth_INITIATION_ch2 = nanmean(PsthArray_INITIATION_ch2);
FIGURE_INITIATION.Err_Positive_ch2 = FIGURE_INITIATION.Psth_INITIATION_ch2 + FIGURE_INITIATION.Err_INITIATION_ch2;
FIGURE_INITIATION.Err_Negative_ch2 = FIGURE_INITIATION.Psth_INITIATION_ch2 - FIGURE_INITIATION.Err_INITIATION_ch2;




%F(12) = figure('color', [1 1 1]);
subplot(1,5,2);
hold on

o = 0.5;
r1 = 255;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);
fill([timeAxis, fliplr(timeAxis)],[FIGURE_INITIATION.Err_Positive_ch1, fliplr(FIGURE_INITIATION.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,FIGURE_INITIATION.Psth_INITIATION_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);


% r2 = 0;
% g2 = 0;
% b2 = 255;
% rgb2_o = opacity (o, r2, g2, b2);
%
% fill([timeAxis, fliplr(timeAxis)],[FIGURE_INITIATION.Err_Positive_ch2, fliplr(FIGURE_INITIATION.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
% h2= plot (timeAxis,FIGURE_INITIATION.Psth_INITIATION_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('Initiation', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([(-1 * nSecPrevWindow),nSecPostWindow]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
set(gca, 'Color', 'none')
line('XData', [(-1 * nSecPrevWindow),nSecPostWindow], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
%legend ([h1 h2], 'pDMS', 'aDMS');
legend('off')

set(gca, 'layer', 'top');
%set(gcf, 'Position', [500 600 2000 400])


% PsthArray_INITIATION_Block.Block1 = PsthArray_INITIATION(1:50, :);
% PsthArray_INITIATION_Block.Block2 = PsthArray_INITIATION(51:100, :);
% PsthArray_INITIATION_Block.Block3 = PsthArray_INITIATION(101:150, :);
% PsthArray_INITIATION_Block.Block4 = PsthArray_INITIATION(151:200, :);

% make raw_PSTH for Liton
nCENTERLITON = length(TimeStamp_CenterLitOn_Analysis);
PsthArray_CENTERLITON_ch1 = NaN(nCENTERLITON,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_CENTERLITON_ch2 = NaN(nCENTERLITON,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nCENTERLITON
    thisTime = TimeStamp_CenterLitOn_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_CENTERLITON_ch1(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_CENTERLITON_ch2(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch2, thisIndex, nTsPrev, nTsPost);
end


FIGURE_CENTERLITON.Err_CENTERLITON_ch1 = (nanstd(PsthArray_CENTERLITON_ch1))/sqrt(size(PsthArray_CENTERLITON_ch1,1));
FIGURE_CENTERLITON.Psth_CENTERLITON_ch1 = nanmean(PsthArray_CENTERLITON_ch1);
FIGURE_CENTERLITON.Err_Positive_ch1 = FIGURE_CENTERLITON.Psth_CENTERLITON_ch1 + FIGURE_CENTERLITON.Err_CENTERLITON_ch1;
FIGURE_CENTERLITON.Err_Negative_ch1 = FIGURE_CENTERLITON.Psth_CENTERLITON_ch1 - FIGURE_CENTERLITON.Err_CENTERLITON_ch1;

FIGURE_CENTERLITON.Err_CENTERLITON_ch2 = (nanstd(PsthArray_CENTERLITON_ch2))/sqrt(size(PsthArray_CENTERLITON_ch2,1));
FIGURE_CENTERLITON.Psth_CENTERLITON_ch2 = nanmean(PsthArray_CENTERLITON_ch2);
FIGURE_CENTERLITON.Err_Positive_ch2 = FIGURE_CENTERLITON.Psth_CENTERLITON_ch2 + FIGURE_CENTERLITON.Err_CENTERLITON_ch2;
FIGURE_CENTERLITON.Err_Negative_ch2 = FIGURE_CENTERLITON.Psth_CENTERLITON_ch2 - FIGURE_CENTERLITON.Err_CENTERLITON_ch2;



%F(13) = figure('color', [1 1 1]);
subplot(1,5,1)
hold on

o = 0.5;
r1 = 255;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);
fill([timeAxis, fliplr(timeAxis)],[FIGURE_CENTERLITON.Err_Positive_ch1, fliplr(FIGURE_CENTERLITON.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,FIGURE_CENTERLITON.Psth_CENTERLITON_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

% r2 = 0;
% g2 = 0;
% b2 = 255;
% rgb2_o = opacity (o, r2, g2, b2);
%
% fill([timeAxis, fliplr(timeAxis)],[FIGURE_CENTERLITON.Err_Positive_ch2, fliplr(FIGURE_CENTERLITON.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
% h2= plot (timeAxis,FIGURE_CENTERLITON.Psth_CENTERLITON_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('CenterLightON', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([(-1 * nSecPrevWindow),nSecPostWindow]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
set(gca, 'Color', 'none')
line('XData', [(-1 * nSecPrevWindow),nSecPostWindow], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
%legend ([h1 h2], 'pDMS', 'aDMS');
legend('off')

set(gca, 'layer', 'top');
%set(gcf, 'Position', [0 600 500 400])


% PsthArray_LITON_Block.Block1 = PsthArray_CENTERLITON(1:50, :);
% PsthArray_LITON_Block.Block2 = PsthArray_CENTERLITON(51:100, :);
% PsthArray_LITON_Block.Block3 = PsthArray_CENTERLITON(101:150, :);
% PsthArray_LITON_Block.Block4 = PsthArray_CENTERLITON(151:200, :);

nREWARD = length(TimeStamp_RewardRetrieval_Analysis);
PsthArray_REWARD_ch1 = NaN(nREWARD,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_REWARD_ch2 = NaN(nREWARD,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nREWARD
    thisTime = TimeStamp_RewardRetrieval_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_REWARD_ch1(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_REWARD_ch2(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch2, thisIndex, nTsPrev, nTsPost);
end


FIGURE_REWARD.Err_REWARD_ch1 = (nanstd(PsthArray_REWARD_ch1))/sqrt(size(PsthArray_REWARD_ch1,1));
FIGURE_REWARD.Psth_REWARD_ch1 = nanmean(PsthArray_REWARD_ch1);
FIGURE_REWARD.Err_Positive_ch1 = FIGURE_REWARD.Psth_REWARD_ch1 + FIGURE_REWARD.Err_REWARD_ch1;
FIGURE_REWARD.Err_Negative_ch1 = FIGURE_REWARD.Psth_REWARD_ch1 - FIGURE_REWARD.Err_REWARD_ch1;

FIGURE_REWARD.Err_REWARD_ch2 = (nanstd(PsthArray_REWARD_ch2))/sqrt(size(PsthArray_REWARD_ch2,1));
FIGURE_REWARD.Psth_REWARD_ch2 = nanmean(PsthArray_REWARD_ch2);
FIGURE_REWARD.Err_Positive_ch2 = FIGURE_REWARD.Psth_REWARD_ch2 + FIGURE_REWARD.Err_REWARD_ch2;
FIGURE_REWARD.Err_Negative_ch2 = FIGURE_REWARD.Psth_REWARD_ch2 - FIGURE_REWARD.Err_REWARD_ch2;

%F(14) = figure('color', [1 1 1]);
subplot(1,5,4)
hold on

o = 0.5;
r1 = 255;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);
fill([timeAxis, fliplr(timeAxis)],[FIGURE_REWARD.Err_Positive_ch1, fliplr(FIGURE_REWARD.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,FIGURE_REWARD.Psth_REWARD_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

% r2 = 0;
% g2 = 0;
% b2 = 255;
% rgb2_o = opacity (o, r2, g2, b2);
%
% fill([timeAxis, fliplr(timeAxis)],[FIGURE_REWARD.Err_Positive_ch2, fliplr(FIGURE_REWARD.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
% h2= plot (timeAxis,FIGURE_REWARD.Psth_REWARD_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('RewardRetrieval', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([(-1 * nSecPrevWindow),nSecPostWindow]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
set(gca, 'Color', 'none')
line('XData', [(-1 * nSecPrevWindow),nSecPostWindow], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
%legend ([h1 h2], 'pDMS', 'aDMS');
legend('off')

set(gca, 'layer', 'top');
%set(gcf, 'Position', [500 0 500 400])
% PsthArray_REWARD_Block.Block1 = PsthArray_REWARD(1:50, :);
% PsthArray_REWARD_Block.Block2 = PsthArray_REWARD(51:100, :);
% PsthArray_REWARD_Block.Block3 = PsthArray_REWARD(101:150, :);
% PsthArray_REWARD_Block.Block4 = PsthArray_REWARD(151:200, :);

% make raw_PSTH for next ISI
nCHOICE = length(TimeStamp_Choice_Analysis);
PsthArray_CHOICE_ch1 = NaN(nCHOICE,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_CHOICE_ch2 = NaN(nCHOICE,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nCHOICE
    thisTime = TimeStamp_Choice_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_CHOICE_ch1(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_CHOICE_ch2(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch2, thisIndex, nTsPrev, nTsPost);
end

FIGURE_CHOICE.Err_CHOICE_ch1 = (nanstd(PsthArray_CHOICE_ch1))/sqrt(size(PsthArray_CHOICE_ch1,1));
FIGURE_CHOICE.Psth_CHOICE_ch1 = nanmean(PsthArray_CHOICE_ch1);
FIGURE_CHOICE.Err_Positive_ch1 = FIGURE_CHOICE.Psth_CHOICE_ch1 + FIGURE_CHOICE.Err_CHOICE_ch1;
FIGURE_CHOICE.Err_Negative_ch1 = FIGURE_CHOICE.Psth_CHOICE_ch1 - FIGURE_CHOICE.Err_CHOICE_ch1;

FIGURE_CHOICE.Err_CHOICE_ch2 = (nanstd(PsthArray_CHOICE_ch2))/sqrt(size(PsthArray_CHOICE_ch2,1));
FIGURE_CHOICE.Psth_CHOICE_ch2 = nanmean(PsthArray_CHOICE_ch2);
FIGURE_CHOICE.Err_Positive_ch2 = FIGURE_CHOICE.Psth_CHOICE_ch2 + FIGURE_CHOICE.Err_CHOICE_ch2;
FIGURE_CHOICE.Err_Negative_ch2 = FIGURE_CHOICE.Psth_CHOICE_ch2 - FIGURE_CHOICE.Err_CHOICE_ch2;


%F(11) = figure('color', [1 1 1]);
subplot(1,5,5)
hold on

o = 0.5;
r1 = 255;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);
fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch1, fliplr(FIGURE_CHOICE.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

% r2 = 0;
% g2 = 0;
% b2 = 255;
% rgb2_o = opacity (o, r2, g2, b2);
%
% fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch2, fliplr(FIGURE_CHOICE.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
% h2= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('ISI', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([0,8]);
xticks([0 3 8]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
set(gca, 'Color', 'none')
line('XData', [0,8], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [3 3], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')

legend ([h1], 'pDMS');

set(gca, 'layer', 'top');
%set(gcf, 'Position', [0 0 2000 400])


%% Graph setting - ch2

% x-Axis Analysis window
nSecPrev = 10; %change to make different window
nSecPost = 10;

% x-Axis View window
nSecPrevWindow = 1; %change to make different window
nSecPostWindow = 1;

% y-axis
ymin = -1.5;
ymax = 1.5;

% convert seconds to TDT timestamps
nTsPrev = round (nSecPrev * samplingRate);
nTsPost = round (nSecPost * samplingRate);

totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;

% make raw_PSTH for Correct Press
nCHOICE = length(TimeStamp_Choice_Analysis);
PsthArray_CHOICE_ch1 = NaN(nCHOICE,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_CHOICE_ch2 = NaN(nCHOICE,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nCHOICE
    thisTime = TimeStamp_Choice_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_CHOICE_ch1(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_CHOICE_ch2(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch2, thisIndex, nTsPrev, nTsPost);
end

FIGURE_CHOICE.Err_CHOICE_ch1 = (nanstd(PsthArray_CHOICE_ch1))/sqrt(size(PsthArray_CHOICE_ch1,1));
FIGURE_CHOICE.Psth_CHOICE_ch1 = nanmean(PsthArray_CHOICE_ch1);
FIGURE_CHOICE.Err_Positive_ch1 = FIGURE_CHOICE.Psth_CHOICE_ch1 + FIGURE_CHOICE.Err_CHOICE_ch1;
FIGURE_CHOICE.Err_Negative_ch1 = FIGURE_CHOICE.Psth_CHOICE_ch1 - FIGURE_CHOICE.Err_CHOICE_ch1;

FIGURE_CHOICE.Err_CHOICE_ch2 = (nanstd(PsthArray_CHOICE_ch2))/sqrt(size(PsthArray_CHOICE_ch2,1));
FIGURE_CHOICE.Psth_CHOICE_ch2 = nanmean(PsthArray_CHOICE_ch2);
FIGURE_CHOICE.Err_Positive_ch2 = FIGURE_CHOICE.Psth_CHOICE_ch2 + FIGURE_CHOICE.Err_CHOICE_ch2;
FIGURE_CHOICE.Err_Negative_ch2 = FIGURE_CHOICE.Psth_CHOICE_ch2 - FIGURE_CHOICE.Err_CHOICE_ch2;


F(11) = figure('color', [1 1 1]); subplot(1,5,3)
hold on

% o = 0.5;
% r1 = 255;
% g1 = 0;
% b1 = 0;
% rgb1_o = opacity (o, r1, g1, b1);
% fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch1, fliplr(FIGURE_CHOICE.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
% h1= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch2, fliplr(FIGURE_CHOICE.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('Choice', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([(-1 * nSecPrevWindow),nSecPostWindow]);
ylim ([ymin, ymax]);

set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
set(gca,'Color','none')
line('XData', [(-1 * nSecPrevWindow),nSecPostWindow], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')

%legend ([h1 h2], 'pDMS', 'aDMS');
legend('off')

set(gca, 'layer', 'top');
set(gcf, 'Position', [400 400 1600 400])

%
% h1f=fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch1, fliplr(FIGURE_CHOICE.Err_Negative_ch1)], 'g', 'EdgeColor', 'none');
% set(h1f,'facealpha',.5);
% h1= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch1,'Color','r', 'Linewidth', 3);

% PsthArray_CORRECTPRESS_Block.Block1 = PsthArray_CORRECTPRESS(1:50, :);
% PsthArray_CORRECTPRESS_Block.Block2 = PsthArray_CORRECTPRESS(51:100, :);
% PsthArray_CORRECTPRESS_Block.Block3 = PsthArray_CORRECTPRESS(101:150, :);
% %PsthArray_CORRECTPRESS_Block.Block4 = PsthArray_CORRECTPRESS(151:200, :);



% make raw_PSTH for Initiation
nINITIATION = length(TimeStamp_Initiation_Analysis);
PsthArray_INITIATION_ch1 = NaN(nINITIATION,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_INITIATION_ch2 = NaN(nINITIATION,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nINITIATION
    thisTime = TimeStamp_Initiation_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_INITIATION_ch1(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_INITIATION_ch2(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch2, thisIndex, nTsPrev, nTsPost);
end

FIGURE_INITIATION.Err_INITIATION_ch1 = (nanstd(PsthArray_INITIATION_ch1))/sqrt(size(PsthArray_INITIATION_ch1,1));
FIGURE_INITIATION.Psth_INITIATION_ch1 = nanmean(PsthArray_INITIATION_ch1);
FIGURE_INITIATION.Err_Positive_ch1 = FIGURE_INITIATION.Psth_INITIATION_ch1 + FIGURE_INITIATION.Err_INITIATION_ch1;
FIGURE_INITIATION.Err_Negative_ch1 = FIGURE_INITIATION.Psth_INITIATION_ch1 - FIGURE_INITIATION.Err_INITIATION_ch1;

FIGURE_INITIATION.Err_INITIATION_ch2 = (nanstd(PsthArray_INITIATION_ch2))/sqrt(size(PsthArray_INITIATION_ch2,1));
FIGURE_INITIATION.Psth_INITIATION_ch2 = nanmean(PsthArray_INITIATION_ch2);
FIGURE_INITIATION.Err_Positive_ch2 = FIGURE_INITIATION.Psth_INITIATION_ch2 + FIGURE_INITIATION.Err_INITIATION_ch2;
FIGURE_INITIATION.Err_Negative_ch2 = FIGURE_INITIATION.Psth_INITIATION_ch2 - FIGURE_INITIATION.Err_INITIATION_ch2;




%F(12) = figure('color', [1 1 1]);
subplot(1,5,2);
hold on

% o = 0.5;
% r1 = 255;
% g1 = 0;
% b1 = 0;
% rgb1_o = opacity (o, r1, g1, b1);
% fill([timeAxis, fliplr(timeAxis)],[FIGURE_INITIATION.Err_Positive_ch1, fliplr(FIGURE_INITIATION.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
% h1= plot (timeAxis,FIGURE_INITIATION.Psth_INITIATION_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);


r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[FIGURE_INITIATION.Err_Positive_ch2, fliplr(FIGURE_INITIATION.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,FIGURE_INITIATION.Psth_INITIATION_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('Initiation', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([(-1 * nSecPrevWindow),nSecPostWindow]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
set(gca,'Color','none')
line('XData', [(-1 * nSecPrevWindow),nSecPostWindow], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
%legend ([h1 h2], 'pDMS', 'aDMS');
legend('off')

set(gca, 'layer', 'top');
%set(gcf, 'Position', [500 600 2000 400])


% PsthArray_INITIATION_Block.Block1 = PsthArray_INITIATION(1:50, :);
% PsthArray_INITIATION_Block.Block2 = PsthArray_INITIATION(51:100, :);
% PsthArray_INITIATION_Block.Block3 = PsthArray_INITIATION(101:150, :);
% PsthArray_INITIATION_Block.Block4 = PsthArray_INITIATION(151:200, :);

% make raw_PSTH for Liton
nCENTERLITON = length(TimeStamp_CenterLitOn_Analysis);
PsthArray_CENTERLITON_ch1 = NaN(nCENTERLITON,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_CENTERLITON_ch2 = NaN(nCENTERLITON,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nCENTERLITON
    thisTime = TimeStamp_CenterLitOn_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_CENTERLITON_ch1(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_CENTERLITON_ch2(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch2, thisIndex, nTsPrev, nTsPost);
end


FIGURE_CENTERLITON.Err_CENTERLITON_ch1 = (nanstd(PsthArray_CENTERLITON_ch1))/sqrt(size(PsthArray_CENTERLITON_ch1,1));
FIGURE_CENTERLITON.Psth_CENTERLITON_ch1 = nanmean(PsthArray_CENTERLITON_ch1);
FIGURE_CENTERLITON.Err_Positive_ch1 = FIGURE_CENTERLITON.Psth_CENTERLITON_ch1 + FIGURE_CENTERLITON.Err_CENTERLITON_ch1;
FIGURE_CENTERLITON.Err_Negative_ch1 = FIGURE_CENTERLITON.Psth_CENTERLITON_ch1 - FIGURE_CENTERLITON.Err_CENTERLITON_ch1;

FIGURE_CENTERLITON.Err_CENTERLITON_ch2 = (nanstd(PsthArray_CENTERLITON_ch2))/sqrt(size(PsthArray_CENTERLITON_ch2,1));
FIGURE_CENTERLITON.Psth_CENTERLITON_ch2 = nanmean(PsthArray_CENTERLITON_ch2);
FIGURE_CENTERLITON.Err_Positive_ch2 = FIGURE_CENTERLITON.Psth_CENTERLITON_ch2 + FIGURE_CENTERLITON.Err_CENTERLITON_ch2;
FIGURE_CENTERLITON.Err_Negative_ch2 = FIGURE_CENTERLITON.Psth_CENTERLITON_ch2 - FIGURE_CENTERLITON.Err_CENTERLITON_ch2;



%F(13) = figure('color', [1 1 1]);
subplot(1,5,1)
hold on

% o = 0.5;
% r1 = 255;
% g1 = 0;
% b1 = 0;
% rgb1_o = opacity (o, r1, g1, b1);
% fill([timeAxis, fliplr(timeAxis)],[FIGURE_CENTERLITON.Err_Positive_ch1, fliplr(FIGURE_CENTERLITON.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
% h1= plot (timeAxis,FIGURE_CENTERLITON.Psth_CENTERLITON_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[FIGURE_CENTERLITON.Err_Positive_ch2, fliplr(FIGURE_CENTERLITON.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,FIGURE_CENTERLITON.Psth_CENTERLITON_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('CenterLightON', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([(-1 * nSecPrevWindow),nSecPostWindow]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
set(gca,'Color','none')
line('XData', [(-1 * nSecPrevWindow),nSecPostWindow], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
%legend ([h1 h2], 'pDMS', 'aDMS');
legend('off')

set(gca, 'layer', 'top');
%set(gcf, 'Position', [0 600 500 400])


% PsthArray_LITON_Block.Block1 = PsthArray_CENTERLITON(1:50, :);
% PsthArray_LITON_Block.Block2 = PsthArray_CENTERLITON(51:100, :);
% PsthArray_LITON_Block.Block3 = PsthArray_CENTERLITON(101:150, :);
% PsthArray_LITON_Block.Block4 = PsthArray_CENTERLITON(151:200, :);

nREWARD = length(TimeStamp_RewardRetrieval_Analysis);
PsthArray_REWARD_ch1 = NaN(nREWARD,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_REWARD_ch2 = NaN(nREWARD,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nREWARD
    thisTime = TimeStamp_RewardRetrieval_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_REWARD_ch1(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_REWARD_ch2(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch2, thisIndex, nTsPrev, nTsPost);
end


FIGURE_REWARD.Err_REWARD_ch1 = (nanstd(PsthArray_REWARD_ch1))/sqrt(size(PsthArray_REWARD_ch1,1));
FIGURE_REWARD.Psth_REWARD_ch1 = nanmean(PsthArray_REWARD_ch1);
FIGURE_REWARD.Err_Positive_ch1 = FIGURE_REWARD.Psth_REWARD_ch1 + FIGURE_REWARD.Err_REWARD_ch1;
FIGURE_REWARD.Err_Negative_ch1 = FIGURE_REWARD.Psth_REWARD_ch1 - FIGURE_REWARD.Err_REWARD_ch1;

FIGURE_REWARD.Err_REWARD_ch2 = (nanstd(PsthArray_REWARD_ch2))/sqrt(size(PsthArray_REWARD_ch2,1));
FIGURE_REWARD.Psth_REWARD_ch2 = nanmean(PsthArray_REWARD_ch2);
FIGURE_REWARD.Err_Positive_ch2 = FIGURE_REWARD.Psth_REWARD_ch2 + FIGURE_REWARD.Err_REWARD_ch2;
FIGURE_REWARD.Err_Negative_ch2 = FIGURE_REWARD.Psth_REWARD_ch2 - FIGURE_REWARD.Err_REWARD_ch2;

%F(14) = figure('color', [1 1 1]);
subplot(1,5,4)
hold on

% o = 0.5;
% r1 = 255;
% g1 = 0;
% b1 = 0;
% rgb1_o = opacity (o, r1, g1, b1);
% fill([timeAxis, fliplr(timeAxis)],[FIGURE_REWARD.Err_Positive_ch1, fliplr(FIGURE_REWARD.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
% h1= plot (timeAxis,FIGURE_REWARD.Psth_REWARD_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[FIGURE_REWARD.Err_Positive_ch2, fliplr(FIGURE_REWARD.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,FIGURE_REWARD.Psth_REWARD_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('RewardRetrieval', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([(-1 * nSecPrevWindow),nSecPostWindow]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
set(gca,'Color','none')
line('XData', [(-1 * nSecPrevWindow),nSecPostWindow], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
%legend ([h1 h2], 'pDMS', 'aDMS');
legend('off')

set(gca, 'layer', 'top');
%set(gcf, 'Position', [500 0 500 400])
% PsthArray_REWARD_Block.Block1 = PsthArray_REWARD(1:50, :);
% PsthArray_REWARD_Block.Block2 = PsthArray_REWARD(51:100, :);
% PsthArray_REWARD_Block.Block3 = PsthArray_REWARD(101:150, :);
% PsthArray_REWARD_Block.Block4 = PsthArray_REWARD(151:200, :);

% make raw_PSTH for next ISI
nCHOICE = length(TimeStamp_Choice_Analysis);
PsthArray_CHOICE_ch1 = NaN(nCHOICE,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_CHOICE_ch2 = NaN(nCHOICE,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nCHOICE
    thisTime = TimeStamp_Choice_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_CHOICE_ch1(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_CHOICE_ch2(i,:) = processPhotDataRow_normDat(normDat_DB_modZ_ch2, thisIndex, nTsPrev, nTsPost);
end

FIGURE_CHOICE.Err_CHOICE_ch1 = (nanstd(PsthArray_CHOICE_ch1))/sqrt(size(PsthArray_CHOICE_ch1,1));
FIGURE_CHOICE.Psth_CHOICE_ch1 = nanmean(PsthArray_CHOICE_ch1);
FIGURE_CHOICE.Err_Positive_ch1 = FIGURE_CHOICE.Psth_CHOICE_ch1 + FIGURE_CHOICE.Err_CHOICE_ch1;
FIGURE_CHOICE.Err_Negative_ch1 = FIGURE_CHOICE.Psth_CHOICE_ch1 - FIGURE_CHOICE.Err_CHOICE_ch1;

FIGURE_CHOICE.Err_CHOICE_ch2 = (nanstd(PsthArray_CHOICE_ch2))/sqrt(size(PsthArray_CHOICE_ch2,1));
FIGURE_CHOICE.Psth_CHOICE_ch2 = nanmean(PsthArray_CHOICE_ch2);
FIGURE_CHOICE.Err_Positive_ch2 = FIGURE_CHOICE.Psth_CHOICE_ch2 + FIGURE_CHOICE.Err_CHOICE_ch2;
FIGURE_CHOICE.Err_Negative_ch2 = FIGURE_CHOICE.Psth_CHOICE_ch2 - FIGURE_CHOICE.Err_CHOICE_ch2;


%F(11) = figure('color', [1 1 1]);
subplot(1,5,5)
hold on

% o = 0.5;
% r1 = 255;
% g1 = 0;
% b1 = 0;
% rgb1_o = opacity (o, r1, g1, b1);
% fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch1, fliplr(FIGURE_CHOICE.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
% h1= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch2, fliplr(FIGURE_CHOICE.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('ISI', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([0,8]);
xticks([0 3 8]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
set(gca,'Color','none')
line('XData', [0,8], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [3 3], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')

legend ([h2], 'aDMS');

set(gca, 'layer', 'top');
%set(gcf, 'Position', [0 0 2000 400])


