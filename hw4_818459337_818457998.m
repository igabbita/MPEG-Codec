close all;
clear all;
clc;

Video_Frame = aviread('C:\Users\Iswarya Gabbita\Documents\MATLAB\2nd Sem\walk_qcif.avi');
MBSize= 8; % Macroblock size
No_Frames= 5; % Number of frames to be worked on in the video

Width = size(Video_Frame(1).cdata, 1); % Width of frame
Height = size(Video_Frame(1).cdata, 2); %Height of frame
MBCountX = Width/MBSize; % Number of macroblocks in a row
MBCountY = Height/MBSize; % Number of macroblocks in a column

E_Y = zeros(Width, Height, No_Frames); % Array of encoded Y frames
E_Cb = zeros(Width/2, Height/2, No_Frames); % Array of encoded Y frames
E_Cr = zeros(Width/2, Height/2, No_Frames); % Array of encoded Y frames

MV_Luma = zeros(2, Width*Height/MBSize.^2, No_Frames); % Stores Motion Vectors of luma component of all frames 
MV_Chroma = zeros(2, Width*Height/MBSize.^2, No_Frames); % Stores Motion Vectors of chroma component of all frames

%Encode the I frame
[E_Y(:,:,1) E_Cb(:,:,1) E_Cr(:,:,1)] = encodeImage(rgb2ycbcr(Video_Frame(1).cdata));

initialStepSize = 8; % Initial step size for the logarithmic search

dec_I = (predictRefImage(E_Y(:,:,1), E_Cb(:,:,1), E_Cr(:,:,1)));
figure; imshow(ycbcr2rgb(dec_I));

Ref_Y = double(dec_I(:,:,1)); % Get Y component for Reference frame
Ref_Cb = double(dec_I(:,:,2));
Ref_Cr = double(dec_I(:,:,3));

%%Loop over frames
for n = 1:4
    TargetFrame = double(rgb2ycbcr(video_frames(n+1).cdata));
    
    %Use 2D Logarithmic Search to find motion vectors
    Target_Y = double(TargetFrame(:,:,1)); %Get Y component for Target frame

    
    %Estimate motion vectors for the Y-component frame using 2D logarithmic search
    [MotionVect, Computations] = MotionEstLog(Target_Y, Ref_Y, MBSize, initialStepSize);
    MV_Luma(:,:,n+1) = MotionVect;
    
    %Display motion vectors using quiver
    [x, y] = meshgrid(1:1:MBCountY, 1:1:MBCountX);
    MBCount = 1;
    for i = 1:MBCountX
        for j = 1:MBCountY
            mvx(i, j) = MotionVect(1, MBCount);
            mvy(i, j) = MotionVect(2, MBCount);
            MBCount = MBCount + 1;
        end
    end
    
    figure
    quiver(x, y, mvx, mvy, 1.5)
    title(sprintf('2D Logarithmic Motion vectors for frame %i',n+1))
    
 %Scale down the motion vectors of Y component to get motion vectors for
    %the chroma components of the video frame
    MotionVector_Chroma = MotionVect/2;
    MV_Chroma(:,:,n+1) = MotionVector_Chroma;
   
    
    % Compute motion compensated Y-component using motion vectors
    ImgComp_Y = MotionComp(ref_Y, motionVect, mbSize);
    
    % Subsample Cb and Cr bands using 4:2:0 
    ref_Cb_SS = ref_Cb(1:2:size(ref_Cb, 1), 1:2:size(ref_Cb, 2));
    ref_Cr_SS = ref_Cr(1:2:size(ref_Cr, 1), 1:2:size(ref_Cr, 2));
    
    % Compute motion compensated Cb,Cr-components using scaled motion vectors
    ImgComp_Cb = MotionComp(ref_Cb_SS, MotionVect_chroma, mbSize/2);
    ImgComp_Cr = MotionComp(ref_Cr_SS, MotionVect_chroma, mbSize/2);
    
    %Upsample using Linear Interpolation and display Cb and Cr bands
    ImgComp_Cb_US = upsample(ImgComp_Cb);
    ImgComp_Cr_US = upsample(ImgComp_Cr);
    
    %Convert YCbCr upsampled image to RGB and display
    Comp_YCbCr_US = (zeros(144, 176, 3));
    Comp_YCbCr_US(:,:,1) =  (ImgComp_Y);
    Comp_YCbCr_US(:,:,2) =  (ImgComp_Cb_US);
    Comp_YCbCr_US(:,:,3) =  (ImgComp_Cr_US);

    %Get the difference of the target frame and its estimate
    diff = Target_Frame - Comp_YCbCr_US;
   %Encode the difference
    [Enc_Y(:,:,n+1) Enc_Cb(:,:,n+1) Enc_Cr(:,:,n+1)] = encodeImage(diff);
    
    %Predict/Decode the error frame and add it to the motion compensated image
    %so it can be used as reference in the next iteration
    Decode_diff = predictRefImage(Enc_Y(:,:,n+1), Enc_Cb(:,:,n+1), Enc_Cr(:,:,n+1));

    ref_frame = double(Decode_diff) + Comp_YCbCr_US;
    %figure; imshow(ycbcr2rgb(uint8(ref_frame)))
    ref_Y = double(ref_frame(:,:,1));
    ref_Cb = double(ref_frame(:,:,2));
    ref_Cr = double(ref_frame(:,:,3));
end
