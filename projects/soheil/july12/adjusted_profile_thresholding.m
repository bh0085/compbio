function [pf_th,th_used]=adjusted_profile_thresholding(pf,num_max)

% it threshold this profile to have at most num_max rows

row_size=size(pf,1);
if row_size<=num_max
    pf_th=pf;
    th_used=0;
else
    th=0.3;
    stop=0;
    pf_th=profile_thresholding(pf,th);
    if size(pf_th,1)<=num_max
        stop=1;
        th_used=th;
    else
       while stop==0 
        th=th+0.01;
        pf_th=profile_thresholding(pf,th);
        if size(pf_th,1)<=num_max
        stop=1;
        th_used=th;
        end
       end
        
    end
    
    
end

