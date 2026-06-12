%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to apply percolation analysis on a weighted Graph as described in Bardella et al. 2016 and Bardella et al. 2020
% 
% [input] 
%  Matrix: symmetric matrix 
%  
%
% [output] 
%  perc_threshold: first threshold at which #of connected components (#CC) > 1
%  perc_threshold_step: step at which perc_threshold is found 
%  Perc_Matrix: cleaned matrix with links < perc_threshold set to 0
%  Sort_value: link values in ascending order
%  n_com_size: #CC for each step
%
%
%  If used for your research project please cite: 
%  Bardella et al. 2016 (https://doi.org/10.1038/srep32060)
%  Bardella et al. 2020 (https://doi.org/10.1016/j.neuroimage.2019.116354)
%  Bardella et al. 2024 (https://doi.org/10.1162/netn_a_00365)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [perc_threshold,Perc_Matrix, Sort_value,n_com_size,perc_threshold_step] = Percolation (Matrix)

Sort_value = unique(sort(Matrix));
dont_enter=0;

        matrix_corr = Matrix;
        matrix = Matrix;
        
        for t = 1:length(Sort_value)
            matrix(matrix == Sort_value(t)) = 0;
            [n_com, bin_size] = conncomp(graph(matrix,'upper'));
            n_com_size(t,1)=length(bin_size);
           
            if max(n_com)>1 && dont_enter==0
                dont_enter=1;
                matrix2 = matrix_corr;
                corr_sort = Sort_value(t-1);
                x1 = matrix2 <= corr_sort;
                matrix2(x1) = 0;
                perc_threshold = corr_sort;
                Perc_Matrix = matrix2;
                perc_threshold_step=t-1;

            end
        end
end