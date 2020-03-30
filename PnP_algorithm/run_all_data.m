filenames = {'duomino_cr_10','duomino_cr_20','duomino_cr_30','duomino_cr_40','duomino_cr_50',...
    'hand_cr_10','hand_cr_20','hand_cr_30','hand_cr_40','hand_cr_50',...
    'pendulumBall_cr_10','pendulumBall_cr_20','pendulumBall_cr_30','pendulumBall_cr_40','pendulumBall_cr_50',...
    'pingpang_cr_10','pingpang_cr_20','pingpang_cr_30','pingpang_cr_40','pingpang_cr_50',...
    'waterBalloon_cr_10','waterBalloon_cr_20','waterBalloon_cr_30','waterBalloon_cr_40','waterBalloon_cr_50',...
    'big_long_cr_10','big_long_cr_20','big_long_cr_30','big_long_cr_40','big_long_cr_50'};
   
crs = [10 20 30 40 50 10 20 30 40 50 10 20 30 40 50 10 20 30 40 50 10 20 30 40 50 10 20 30 40 50];

sigmas = [20 25 45 60 60 20 25 45 60 60 20 25 45 60 60 20 25 45 60 60 20 25 45 60 60 20 25 45 60 60];

maxiters = [30 50 150 150 150 30 50 150 150 150 30 50 150 150 150 30 50 150 150 150 30 50 150 150 150 30 50 150 150 150];

for i_file=26:30
    cr = crs(i_file);    % overall cr
    
    filename = filenames{i_file};
        
    sigma = sigmas(i_file);
    
    maxiter = maxiters(i_file);
    
    main_real_cacti_20191010_GPU_single_recon;
end