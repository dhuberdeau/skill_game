function [ind, dist] = knnsearch_unique(X, Y)

% function [ind, dist] = knnsearch_unique(X, Y)
%
% Input:
%   X = the index list
%   Y = the search or master list
%
% Finds the index of values in X that are best-matched to values in Y,
% forbidding unique pairings between X and Y. I.e. only unique indicies of
% X values are paired to Y values. The algorithm uses knnsearch with a
% greedy iterative algorithm on top that searches for pairs from the most
% minimal and in ascending order when pair conflicts arrise.
%
% Output:
%   ind = the indicies of values from X closest-paired with Y values
%   dist = the distance between values of Y and the corresponding closest
%   paired value from X, given the uniqueness constraint.
%
% David Huberdeau
% Copyright 2018

if ~(isempty(X) || isempty(Y))
    [indx_, distx_] = knnsearch(X, Y, 'K', length(X));

    % For each row (i.e. each sample in Y), choose the index with
    % minimum value that hasn't been previously chosen, then move
    % to next row and continue. If a value has been chosen, go down
    % the columns until you find a unique value.

    unique_distX = nan(length(Y), 1);
    unique_indX = nan(length(Y), 1);

    rows_left = 1:length(Y);
    row_counter = 0;
    while ~isempty(rows_left)
        temp_distx = distx_(rows_left, :);
        temp_indx = indx_(rows_left, :);
        i_col = 1;
        unique_found = 0;
        while ~unique_found && i_col <= size(temp_distx,2)
            [dist_min, index_min] = min(temp_distx(:, i_col));
            this_min_indx = temp_indx(index_min, i_col);
            if ~ismember(this_min_indx, unique_indX)
%                 unique_indX(sum(~isnan(unique_indX)) + 1) = this_min_indx;
%                 unique_distX(sum(~isnan(unique_distX)) + 1) = dist_min;
                unique_indX(rows_left(index_min)) = this_min_indx;
                unique_distX(rows_left(index_min)) = dist_min;
                unique_found = 1;
                rows_left = setdiff(rows_left, rows_left(index_min));
            else
                i_col = i_col + 1;
            end
        end
        row_counter = row_counter + 1;
        if row_counter > length(Y)
            warning('No matched pair found for at least one value')
            break
        end
    end

    ind = unique_indX;
    dist = unique_distX;
    
else
    warning('No values to match');
    ind = nan;
    dist = nan;
end





%% unused earlier drafts:

            % PROBLEM: knnsearch could end up having large number of
            % repeated pairings between X & Y, thus highly biasing the
            % measure of corresponding difference in success rates. Thus,
            % changing code to make pairings unique only. 
%             [indx, distx] = knnsearch(X, Y);
%             unique_indX = nan(size(Y));
%             unique_distX = nan(size(Y));
            % get indicies and corresponding distances of matches:
%             [indx_, distx_] = knnsearch(X(:), Y(:), 'K', length(X));
            % conservative algorithm:
            % Sort the maximum distances of each row, exclude any above the 
            % threshold, then pick the indicies in order of highest to 
            % lowest distances.
%             distx_cutoff = distx_;
%             distx_cutoff(distx_cutoff > TH_DIST) = nan;
%             [val_max, ind_max] = sort(max(distx_cutoff,[],2));
%             distx_maxsort = distx_(ind_max, end:-1:1);
%             indx_maxsort = indx_(ind_max, end:-1:1);
%             unique_inds_left = 1;
%             kth_col = 1; kth_row = 1;
%             used_inds = [];
%             remaining_inds = 1:length(X);
%             kth_ind = 1;
%             while unique_inds_left
% %                 [this_dist, k_col] = min(distx_maxsort(:, kth_col));
%                 this_dist = distx_maxsort(kth_row, kth_col);
%                 this_ind = indx_maxsort(kth_row, kth_col);
%                 if ~ismember(this_ind, used_inds)
%                     used_inds(length(used_inds) + 1) = this_ind;
%                     unique_distX(kth_ind) = this_dist;
%                     unique_indX(kth_ind) = this_ind;
%                     kth_ind = kth_ind + 1;
%                     remaining_inds = setdiff(remaining_inds, this_ind);
%                 end
%                 if length(setdiff(indx_maxsort(:), remaining_inds)) > 0
%                     unique_inds_left = 1;
%                 else
%                     unique_inds_left = 0;
%                 end
%                 kth_col = kth_col + 1;
%                 if kth_col > size(distx_maxsort,2)
%                     kth_col = 1;
%                     kth_row = kth_row + 1;
%                 end
%                 if kth_row > size(distx_maxsort,1)
%                     unique_inds_left = 0;
%                 end
%             end
            
%             % Greedy algorithm:
            % For each subsequent column of distances, pick the minimum and
            % exclude that index from further selection. Continue until you
            % reach the end or a threshold distance:
%             below_TH = min(distx_(:,1)) < TH_DIST;
%             kth_col = 1;
%             while below_TH
%                 [dist_min, index_min] = min(distx_(:, kth_col));
%                 unique_distX(kth_col) = dist_min;
%                 unique_indX(kth_col) = indx_(index_min);
%             end
            
            % Only algorithm that makes sense:
            % For each row (i.e. each sample in Y), choose the index with
            % minimum value that hasn't been previously chosen, then move
            % to next row and continue. If a value has been chosen, go down
            % the columns until you find a unique value.
            
%             unique_distX = nan(length(Y), 1);
%             unique_indX = nan(length(Y), 1);
%             
%             rows_left = 1:length(Y);
%             row_counter = 0;
%             while ~isempty(rows_left)
%                 temp_distx = distx_(rows_left, :);
%                 temp_indx = indx_(rows_left, :);
%                 i_col = 1;
%                 unique_found = 0;
%                 while ~unique_found && i_col <= size(temp_distx,2)
%                     [dist_min, index_min] = min(temp_distx(:, i_col));
%                     this_min_indx = temp_indx(index_min, i_col);
%                     if ~ismember(this_min_indx, unique_indX)
%                         unique_indX(sum(~isnan(unique_indX)) + 1) = this_min_indx;
%                         unique_distX(sum(~isnan(unique_distX)) + 1) = dist_min;
%                         unique_found = 1;
%                         rows_left = setdiff(rows_left, rows_left(index_min));
%                     else
%                         i_col = i_col + 1;
%                     end
%                 end
%                 row_counter = row_counter + 1;
%                 if row_counter > length(Y)
%                     warning('No matched pair found for at least one value')
%                     break
%                 end
%             end
% 
%             indx = unique_indX;
%             distx = unique_distX;