function [] = o_deconvolution_Batch()

global SETTINGS

SETTINGS.PLOT_GRAPHS = false;

arr_ex_num = {'1','2','3'};
arr_noise = {1e-6,1e-8,1e-10,1e-12};
arr_bool_preproc = {true, false};


for i1 = 1:1:length(arr_ex_num)
    ex_num = arr_ex_num{i1};
    for i2 = 1:1:length(arr_noise)
        emin = arr_noise{i2};
        for i3 = 1:1:length(arr_bool_preproc)
            bool_preproc = arr_bool_preproc{i3};
            
            Test_Deconvolution(ex_num,emin,bool_preproc)
            
        end
    end
end


end