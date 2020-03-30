function  [im_denoised]     =  FFD_Net_DenoiserGPU(input, imageNoiseSigma,net)



inputNoiseSigma   =   imageNoiseSigma;

net = vl_simplenn_move(net, 'gpu') ;

format compact;
global sigmas;

input = gpuArray(single(input));

    
sigmas = inputNoiseSigma; 
    
 res    = vl_simplenn(net, input,[],[],'conserveMemory',true,'mode','test'); % matconvnet default % use this if you did  not install matconvnet; very slow
 
 
    output = res(end).x;
   im_denoised = gather(output);
     %   input  = gather(input);  
%if useGPU
        
   % end


end

