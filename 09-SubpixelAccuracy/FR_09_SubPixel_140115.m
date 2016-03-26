clc;
clear all;
close all;
%%
%=========User-Defined Varibles============
%C:\\Users\\Chayut\\Dropbox\\Academic File 3.2\\Intern (img)
FNAME = '..\\..\\..\\Resource\\Foreman Sequences\\CIF BMP\\foreman_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Waterfall Sequence\\waterfall_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Tempete Sequence\\tempete_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Stefan Sequence\\stefan_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Foreman Sequences\\CIF BMP\\foreman_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Akiyo Sequence\\akiyo_qcif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Coastguard Sequence\\cif\\coastguard_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Flower Sequence\\flower_cif_%d.bmp';
RNAME= 'Foreman_';

StartFrame = 1; %Frame No. of the sequence to start
StopFrame = 40;

iFrameHeight =288; 
iFrameWidth =  352;
Nshift = 11;
DFDFactor = 1; 
iBlkSize = 8;

iZeroFillTreashold = 4;

%%
%============================================
TIME = datestr(now,'yymmddHHMM');
NFrames =StopFrame - StartFrame; %No of Frames to be loaded
iVerticalrun =iFrameHeight/iBlkSize;
iHorizontalrun =iFrameWidth/iBlkSize;
NRows =2*Nshift+1;
NPixels = iFrameHeight*iFrameWidth;
iTotalBlk =iHorizontalrun*iVerticalrun;
MADResult = zeros(NRows,NRows,16);%results in each run
ResultVectorArrayPreviousStackint = int8(zeros(iVerticalrun,iHorizontalrun,3));%array of vector for Sub-pixel
ResultVectorArrayPreviousdouble = zeros(iVerticalrun,iHorizontalrun,2);%array of vector for Sub-pixel
ResultVectorArrayPreviousHalfStackint = int8(zeros(iVerticalrun,iHorizontalrun,3));%array of vector for Sub-pixel
ResultVectorArrayPreviousHalfdouble = zeros(iVerticalrun,iHorizontalrun,2);%array of vector for Sub-pixel
ResultVectorArrayPreviousFullPxint = zeros(iVerticalrun,iHorizontalrun,2); %array of vector for Full-pixel
OrigLumui8 = cell(1,NFrames);
DFStack = uint8( zeros(iFrameHeight,iFrameWidth,3));
iPadLength = 1; 
%======================================================

%%
hWB = waitbar(0,'Loading Source...');
%Load Source
for i=StartFrame:StopFrame
    OrigLumui8{i} = rgb2gray(imread(sprintf(FNAME,i))); % read file with conversion to RGB
end


%%
%============================================
for m=StartFrame+1:StopFrame
    OrigFrameCurrentui8 = OrigLumui8{m};
    OrigFramePrevisouui8 = OrigLumui8{m-1};

    FramePreviousSbpxPadui8 = padarray(OrigFramePrevisouui8,[iPadLength iPadLength],'replicate','post'); %pad for sub-px
    
    % iFrame2Height = 2* iFrameHeight;
    % iFrame2Width = 2*iFrameWidth;
    % QuadFrameui8 = imresize(OrigFrameui8,[iFrame2Height iFrame2Width],'Bilinear');

    FramePreviousStack = uint8(  zeros(iFrameHeight,iFrameWidth,16));
    %%
    waitbar(0,hWB,'Matching...'); 
    %===========Generate 16 Frames=================================
    for i = 1:iFrameHeight
        for j = 1:iFrameWidth

            
            %obtain 4 neghbouring pixels
            MatB  =double(FramePreviousSbpxPadui8(i:i+1,j:j+1));
            for p =1:4
                for q = 1:4
                    Stackindex = (p-1)*4+q;
                    X = (p-1)/4;
                    Y = (q-1)/4;
                    MatA = [1-X,X];
                    MatC = [1-Y; Y];
                    FramePreviousStack(i,j,Stackindex) = uint8(MatA*MatB*MatC);
                    
                end
            end
        end
        
        waitbar(i/iFrameHeight,hWB,sprintf('Preparing Frame... Shift:%d, Frame:%d/%d, Blk:%d/%d ',Nshift,m,StopFrame,i,iFrameHeight)); 
    end
    
    FramePreviousStackPad = padarray(FramePreviousStack,[Nshift Nshift],'replicate','both');

    
    
    %%
    %==================Vector-map ====================
    for SourceBlkNoV = 1:1:iVerticalrun
        for SourceBlkNoH  = 1:1:iHorizontalrun
            Vsourceblkoffset = (SourceBlkNoV-1) * iBlkSize;
            Hsourceblkoffset = (SourceBlkNoH-1) * iBlkSize;
            %Obtain FrameN
            BlkOrigCurrent=   OrigFrameCurrentui8(...
                Vsourceblkoffset+1:...
                Vsourceblkoffset+iBlkSize,...
                Hsourceblkoffset+1:...
                Hsourceblkoffset+iBlkSize);

            BlkCurrentStack = repmat(BlkOrigCurrent,[1 1 16]);
            
            %---------------N-1---------------------
            %Obtain FrameN-1 and difference (with shift offset)
            for i =-Nshift:1:Nshift 
                for j =-Nshift:1:Nshift

                    BlkPreviousStack =  FramePreviousStackPad(...
                         Nshift+i +Vsourceblkoffset  +1:...
                         Nshift+i +Vsourceblkoffset  +iBlkSize,...
                         Nshift+j +Hsourceblkoffset  +1:...
                         Nshift+j +Hsourceblkoffset + +iBlkSize,:); 

                    %Pixel-wise matching criterion operation 

                    iMADStacki16 = abs(int16(BlkCurrentStack)...
                        -int16(BlkPreviousStack)) ;

                    MADResult(i+Nshift+1,j+Nshift+1,:) =...
                        sum(sum(iMADStacki16));%sumation of all pixel
                end
            end

            %========collect result of quater-px
            [R1,I1] = min(MADResult);
            [R2,I2] = min(R1);
            [R3,I3] = min(R2);

            VerticalInd=  I1(:,I2(I3),I3);
            HorizontalInd=I2(I3);
            FrameInd = I3;
            
            VSbpxOffset = double(floor((FrameInd-1)/4))*0.25;
            HSbpxOffset = rem(FrameInd-1,4)*0.25;
            
            ResultVectorArrayPreviousStackint(SourceBlkNoV,SourceBlkNoH,:)...
                = [VerticalInd-Nshift-1,HorizontalInd-Nshift-1,FrameInd];
            
            ResultVectorArrayPreviousdouble(SourceBlkNoV,SourceBlkNoH,:)= [VerticalInd-Nshift-1+VSbpxOffset,HorizontalInd-Nshift-1+HSbpxOffset];
            %========collect result of half-px
            MADresultHalf = cat(3,MADResult(:,:,1),MADResult(:,:,3),MADResult(:,:,9),MADResult(:,:,11));
            
            [R1,I1] = min(MADresultHalf);
            [R2,I2] = min(R1);
            [R3,I3] = min(R2);

            VerticalInd=  I1(:,I2(I3),I3);
            HorizontalInd=I2(I3);
            FrameInd = I3;
            
            VSbpxOffset = double(floor((FrameInd-1)/2))*0.5;
            HSbpxOffset = rem(FrameInd-1,2)*0.5;
            
            ResultVectorArrayPreviousHalfStackint(SourceBlkNoV,SourceBlkNoH,:)...
                = [VerticalInd-Nshift-1,HorizontalInd-Nshift-1,FrameInd];
            
            ResultVectorArrayPreviousHalfdouble(SourceBlkNoV,SourceBlkNoH,:)= [VerticalInd-Nshift-1+VSbpxOffset,HorizontalInd-Nshift-1+HSbpxOffset];
            
            
            %========collect result of full-px
            [R1, I1] = min(MADResult(:,:,1));
            [R2, I2] = min(R1);
            VerticalLoc=I1(I2);
            HorizontalLoc=I2;
            ResultVectorArrayPreviousFullPxint(SourceBlkNoV,SourceBlkNoH,:) = [VerticalLoc-Nshift-1,HorizontalLoc-Nshift-1];
            
        end
        iBlkNo = ((SourceBlkNoV-1)*iHorizontalrun  + SourceBlkNoH) ;
        waitbar(iBlkNo/iTotalBlk,hWB,sprintf('Matching... Shift:%d, Frame:%d/%d, Blk:%d/%d ',Nshift,m,StopFrame,iBlkNo,iTotalBlk)); 
       end

    
   %%
   %-------------------Reconstruction----------------------
   
   waitbar(1,hWB,'Reconstructing...');
   for SourceBlkNoV = 1:1:iVerticalrun
        for SourceBlkNoH  = 1:1:iHorizontalrun
           
            Vsourceblkoffset = (SourceBlkNoV-1) * iBlkSize;
            Hsourceblkoffset = (SourceBlkNoH-1) * iBlkSize;
            %---------------------
           FramePreviousReconPad = padarray(OrigFramePrevisouui8,[Nshift+iPadLength Nshift+iPadLength],'replicate','both');
            
           
           
           
           
           %----------------------
           VerticalShift = ResultVectorArrayPreviousFullPxint(SourceBlkNoV,SourceBlkNoH,1) ;
           HorizontalShift = ResultVectorArrayPreviousFullPxint(SourceBlkNoV,SourceBlkNoH,2) ;
            
           BlkPreviousFull = FramePreviousReconPad(...
                 Nshift+iPadLength+VerticalShift +Vsourceblkoffset  +1:...
                 Nshift+iPadLength+VerticalShift +Vsourceblkoffset  +iBlkSize,...
                 Nshift+iPadLength+HorizontalShift +Hsourceblkoffset +1:...
                 Nshift+iPadLength+HorizontalShift +Hsourceblkoffset +iBlkSize);
   
           %----------------------
           VerticalShift=  ResultVectorArrayPreviousHalfdouble(SourceBlkNoV,SourceBlkNoH,1);
           HorizontalShift =  ResultVectorArrayPreviousHalfdouble(SourceBlkNoV,SourceBlkNoH,2);   
           
           VPxShiftd = floor(VerticalShift);
           HPxShiftd = floor(HorizontalShift);
           VPxShift = int8(VPxShiftd );
           HPxShift = int8(HPxShiftd );
           VSbpxShift = VerticalShift - VPxShiftd;
           HSbpxShift = HorizontalShift - HPxShiftd;
           
          BlkHalfSbpx = uint8(zeros(iBlkSize,iBlkSize));
          for i = 1:iBlkSize
            for j = 1:iBlkSize
                 OrigPxMat = double(FramePreviousReconPad(...
                 Nshift+iPadLength+VPxShiftd +Vsourceblkoffset +i  :...
                 Nshift+iPadLength+VPxShiftd +Vsourceblkoffset +i +1,...
                 Nshift+iPadLength+HPxShiftd +Hsourceblkoffset +j :...
                 Nshift+iPadLength+HPxShiftd +Hsourceblkoffset +j +1));
                
                 MatA = [1-VSbpxShift,VSbpxShift];
                 MatC = [1-HSbpxShift; HSbpxShift];
                 BlkHalfSbpx(i,j) = uint8(MatA*OrigPxMat*MatC);  
            end
          end
          
          %------------------------------------------
           VerticalShift=  ResultVectorArrayPreviousdouble(SourceBlkNoV,SourceBlkNoH,1);
           HorizontalShift =  ResultVectorArrayPreviousdouble(SourceBlkNoV,SourceBlkNoH,2);   
           
           VPxShiftd = floor(VerticalShift);
           HPxShiftd = floor(HorizontalShift);
           VPxShift = int8(VPxShiftd );
           HPxShift = int8(HPxShiftd );
           VSbpxShift = VerticalShift - VPxShiftd;
           HSbpxShift = HorizontalShift - HPxShiftd;
           
          BlkSbpx = uint8(zeros(iBlkSize,iBlkSize));
          for i = 1:iBlkSize
            for j = 1:iBlkSize
                 OrigPxMat = double(FramePreviousReconPad(...
                 Nshift+iPadLength+VPxShiftd +Vsourceblkoffset +i  :...
                 Nshift+iPadLength+VPxShiftd +Vsourceblkoffset +i +1,...
                 Nshift+iPadLength+HPxShiftd +Hsourceblkoffset +j :...
                 Nshift+iPadLength+HPxShiftd +Hsourceblkoffset +j +1));
                
                 MatA = [1-VSbpxShift,VSbpxShift];
                 MatC = [1-HSbpxShift; HSbpxShift];
                 BlkSbpx(i,j) = uint8(MatA*OrigPxMat*MatC);  
            end
          end
          %-------------------------------------------
          
          
          
          
          DFStack(...
                 Vsourceblkoffset  +1:...
                 Vsourceblkoffset  +iBlkSize,...
                 Hsourceblkoffset  +1:...
                 Hsourceblkoffset + iBlkSize,:)...
            = cat(3,BlkPreviousFull,BlkHalfSbpx,BlkSbpx);
          
          
                    
          FrameStack = cat(3, OrigFramePrevisouui8,DFStack);
        end
   end
    
    
    %%
    %-------Performance Evaluation------------------
    
    DFDStack =  int16(FrameStack)- repmat(int16(OrigFrameCurrentui8),[1 1 4]);
    MSEStack = permute(sum(sum(double(DFDStack).^2))./NPixels,[2 3 1]) ;
    MSEprofile(m-1,:)= MSEStack;
   
    
    %===Display===
%     AvgMSEprofile = mean(MSEprofile,1);
%     figure(1);
%     plot(1:m-1,MSEprofile(:,1),'-o',1:m-1,MSEprofile(:,2),'-x',1:m-1,MSEprofile(:,3),'-*',1:m-1,MSEprofile(:,4),'-+');
%     legend('FD','Full-Pixel','Half-Pixel','Quater-Pixel','Location','NorthWest');
%     title('MSE');
%     hold on 
%     plot(1:m-1,ones(1,m-1)*AvgMSEprofile(1),'--',...
%         1:m-1,ones(1,m-1)*AvgMSEprofile(2),'--',...
%         1:m-1,ones(1,m-1)*AvgMSEprofile(3),'--',...
%         1:m-1,ones(1,m-1)*AvgMSEprofile(4),'--');
%     hold off

    % with frame diff
    waitbar(1,hWB,'Plotting...');
   
    AvgMSEprofile = mean(MSEprofile,1);
    figure(1);
    plot(StartFrame+1:m,MSEprofile(:,2),'-x',...
        StartFrame+1:m,MSEprofile(:,3),'-*',...
        StartFrame+1:m,MSEprofile(:,4),'-+');
    legend('Full-Pixel','Half-Pixel','Quater-Pixel','Location','NorthWest');
    title('MSE');
    hold on 
    plot(StartFrame+1:m,ones(1,m-1)*AvgMSEprofile(2),'--',...
        StartFrame+1:m,ones(1,m-1)*AvgMSEprofile(3),'--',...
        StartFrame+1:m,ones(1,m-1)*AvgMSEprofile(4),'--');
    hold off
    
    %----------------------------
    DisDFDStack = uint8(DFDStack+127);
    GridOut1 = [FrameStack(:,:,2),FrameStack(:,:,3),FrameStack(:,:,4)];
    GridOut2 = [DisDFDStack(:,:,2),DisDFDStack(:,:,3),DisDFDStack(:,:,4)];
    GridOut = cat(1,GridOut1,GridOut2);
    strText = 'Full-Pixel, Half-Pixel, Quater-Pixel';
    RGB = insertText( GridOut, [10,10],strText, 'FontSize', 18, 'FontSize', 18);
    figure(2);
    imshow(RGB);

end

%%


toc



