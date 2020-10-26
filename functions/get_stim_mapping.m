function [] = get_stim_mapping(subjs)
%saves the mapping of stimulus to node for each subject

save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
nNode = 10;

for s = 1:numel(subjs)
    subj = subjs{s};
    HUP_ID = load([save_dir, subj, '/header_clean.mat'], 'HUP_ID');
    if isstruct(HUP_ID)
       HUP_ID = HUP_ID.HUP_ID; 
    end
    try
        stim_map = readtable(['/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavioral_data_raw/',HUP_ID, '_typingTask/subj', subj, '_log_prac.csv']);
    catch
        stim_map = readtable(['/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavioral_data_raw/',HUP_ID, '_typingTask/subj', subj, '_log_prac1.csv']);
    end
    stim_map = stim_map.resp;
    
    order = zeros(1,nNode);
    for n = 1:nNode
        if stim_map{n} == 'q'
            order(n) = 1;
        elseif stim_map{n} == 'w'
            order(n) = 2;
        elseif stim_map{n} == 'e'
            order(n) = 3;
        elseif stim_map{n} == 'r'
            order(n) = 4;
        elseif stim_map{n} == 'v'
            order(n) = 5;
        elseif stim_map{n} == 'b'
            order(n) = 6;
        elseif stim_map{n} == 'u'
            order(n) = 7;
        elseif stim_map{n} == 'i'
            order(n) = 8;
        elseif stim_map{n} == 'o'
            order(n) = 9;
        elseif stim_map{n} == 'p'
            order(n) = 10;
        end
    end
    save([save_dir, subj, '/order.mat'],'order')
end
end

