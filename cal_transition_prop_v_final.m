
    %% To calculate the transition matrix and dwelling time based on cluster indices 'k' for each frame of the Relative Phase Movie 
    %
    % The shape of 'IDX' : time x 1
    % 'unit' need to calibrate the output's unit (unit should be in second).
    %
    % In the transition matrix, (i,j) means the transition from mode j to mode i.
    %
    % made by Youngjai Park, November 21th, 2023.

    function properties = cal_transition_prop_v_final(IDX, K, unit)

    
    properties = struct;
    properties.switches = zeros([1,K]);
    properties.dwell_time = zeros([1,K]);
    properties.dwell_dist = cell([1,K]);
    properties.occurrence = zeros([1,K]);
    properties.trans_freq_mat = zeros([K,K]);
    properties.trans_prob_mat = zeros([K,K]);
    properties.K = K;
    properties.unit = unit;
    
    for k = 1:K
        IDX_1 = 1.*(IDX==k); % binarize the given time series
        % So, IDX_1 array has occuring transition event when an element of IDX_1 is 1 or -1.
        
        % calculate the number of transition from 'mode k' per a second. (Hertz)
        properties.switches(k) = sum(diff(IDX_1)==-1)./(unit*(length(IDX)-1));
        
        % calculate the proportion of 'mode k'. (ratio)
        properties.occurrence(k) = sum(IDX_1)./length(IDX);
    
        % calculate the dwelling time. (Second)
        IDX_2 = diff([0; IDX_1; 0]); % add '0' at the first and the last
        trans_id = find(IDX_2~=0); % find indexes when transition occurs.
        dwell_tmp = diff(trans_id);
        properties.dwell_dist{k} = dwell_tmp(1:2:end).*unit;
        properties.dwell_time(k) = mean(properties.dwell_dist{k});
    end
    
    % calculate a transition matrix per a second. (Hertz)
    IDX_trans = IDX(1:end-1) - IDX(2:end).^K;
    trans_mat = zeros(K);
    for i = 1:K
        for j = 1:K
            % to make all cases as different values, we use IDX_trans
            % (i,j) means the transition from mode j to mode i 
            trans_mat(i,j) = sum(IDX_trans==(j-i^K));
        end
    end
    properties.trans_freq_mat = trans_mat./(unit*(length(IDX_trans)));
    properties.trans_prob_mat = trans_mat./sum(trans_mat,1);

end
