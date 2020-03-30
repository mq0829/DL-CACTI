function  [im_denoised]     =  FFD_Net_Denoiserfast(input, imageNoiseSigma,net)



inputNoiseSigma   =   imageNoiseSigma;

format compact;
global sigmas;



    
sigmas = inputNoiseSigma; 
    
 res    = vl_simplenn(net, single(input),[],[],'conserveMemory',true,'mode','test'); % matconvnet default % use this if you did  not install matconvnet; very slow
    
    im_denoised = res(end).x;
    



end

