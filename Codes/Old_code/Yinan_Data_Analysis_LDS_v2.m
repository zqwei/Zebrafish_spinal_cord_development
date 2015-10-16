function   Yinan_Data_Analysis_LDS_v2 (starting_point, xDim, dff, idx, nfit)

    idxLR        = (idx == 1) | (idx ==2);
    Y            = dff(idxLR, starting_point:starting_point + 5000); % Suggested by Yinan
    nLength      = 2000;
    nTrial       = 1000;
    rand_Y       = randomCropSingleTrial(Y, nLength, nTrial);
    mean_type    = 'Constant_mean';
    tol          = 1e-5;
    cyc          = 10000;
%     [yDim, T, K] = size(rand_Y);
    is_suc = false;

    while ~is_suc
        try
            Ph      = lds_uni_latent(rand_Y, xDim, 'mean_type',mean_type,'tol',tol,'cyc',cyc);
            is_suc  = true;
        catch
            is_suc  = false;
        end
    end
    
    Ph              = rmfield(Ph, {'Y', 'Xk_t'});

    LT                   = size(Y,2);
    [err, y_est, ~] = loo (Y, Ph, [0, LT]); %#ok<ASGLU>

    correct = 1- err; %#ok<NASGU>
    
    save(['Fitted_result/Latent_system_',num2str(starting_point,'%05d'),...
         '_',num2str(xDim,'%02d'),'_',num2str(nfit,'%02d'),'.mat'],'correct',...
         'y_est','Ph');