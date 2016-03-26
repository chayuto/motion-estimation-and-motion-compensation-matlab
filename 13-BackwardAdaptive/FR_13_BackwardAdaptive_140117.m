
%%
clear all;
clc;
close all;
%%
FNAME = '..\\..\\..\\Resource\\Foreman Sequences\\CIF BMP\\foreman_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Waterfall Sequence\\waterfall_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Tempete Sequence\\tempete_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Stefan Sequence\\stefan_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Foreman Sequences\\CIF BMP\\foreman_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Akiyo Sequence\\akiyo_qcif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Coastguard Sequence\\cif\\coastguard_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Flower Sequence\\flower_cif_%d.bmp';
RNAME= 'Foreman_';
TIME = datestr(now,'yymmddHHMM');

iFrameHeight =288; 
iFrameWidth =  352;

StartFrame = 1; %Frame No. of the sequence to start
StopFrame = 40;


DFDFactor = 3;
Nshift = 11;
iBlkSize =8;
iEdgeSize = 2;
iPadLength =1;
%%
%============================================
NFrames =StopFrame - StartFrame; %No of Frames to be loaded
iVertivalrun =iFrameHeight/iBlkSize;
iHorizontalrun =iFrameWidth/iBlkSize;
iTotalBlk =iHorizontalrun*iVertivalrun;
NRows =2*Nshift+1;
NPixels = iFrameHeight*iFrameWidth;
MADResultPrevious = zeros(NRows,NRows);%results in each run

OrigLumui8 = cell(1,NFrames);

%======================================================

%%
%Load Source
for i=StartFrame:StopFrame
    OrigLumui8{i} = rgb2gray(imread(sprintf(FNAME,i))); % read file with conversion to RGB
end

hWB = waitbar(0,'Matching...');
%%
%============================================
for  m=StartFrame+1:StopFrame-1 % loop for all frame-pair
    
    OrigFrameCurrentui8 = OrigLumui8{m};
    OrigFramePrevisouui8 = OrigLumui8{m-1};
    
    FramePreviousSbpxPadui8 = padarray(OrigFramePrevisouui8,[iPadLength iPadLength],'replicate','post'); %pad for sub-px
    %%
    FramePreviousStack = uint8(  zeros(iFrameHeight,iFrameWidth,16));
    %%
    %===========Generate 16 Frames=================================
    for i = 1:iFrameHeight
        for j = 1:iFrameWidth

            MatB  =double( FramePreviousSbpxPadui8(i:i+1,j:j+1));
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
        waitbar(i/iFrameHeight,hWB,sprintf('Preparing Frame Stack.... Frame:%d/%d',m,StopFrame));
            
    end
    
    FramePreviousStackPad = padarray(FramePreviousStack,[Nshift Nshift],'replicate','both');
     %--------------------------------------------
   
   for SourceBlkNoV = 2:1:iVertivalrun
        for SourceBlkNoH  = 2:1:iHorizontalrun
            Vsourceblkoffset = (SourceBlkNoV-1) * iBlkSize;
            Hsourceblkoffset = (SourceBlkNoH-1) * iBlkSize;
            
            BlkA = OrigFrameCurrentui8(...
                Vsourceblkoffset+1-iEdgeSize:Vsourceblkoffset,...
                Hsourceblkoffset+1-iEdgeSize:Hsourceblkoffset+iBlkSize);
            BlkB = OrigFrameCurrentui8(...
                Vsourceblkoffset+1:Vsourceblkoffset+iBlkSize,...
                Hsourceblkoffset+1-iEdgeSize:Hsourceblkoffset);
            BlkOrig = OrigFrameCurrentui8(...
                Vsourceblkoffset+1:Vsourceblkoffset+iBlkSize,...
                 Hsourceblkoffset+1:Hsourceblkoffset+iBlkSize);

             BlkAStack = repmat(BlkA,[1 1 16]);
             BlkBStack = repmat(BlkB,[1 1 16]);
             BlkOrigStack = repmat(BlkOrig,[1 1 16]);
            %%
            %Previous Prediction
            
             for i =-Nshift:1:Nshift 
                for j =-Nshift:1:Nshift

                    BlkASearch =   FramePreviousStackPad(...
                        Nshift+i+Vsourceblkoffset+1-iEdgeSize:...
                        Nshift+i+Vsourceblkoffset,...
                        Nshift+j+Hsourceblkoffset+1-iEdgeSize:...
                        Nshift+j+Hsourceblkoffset+iBlkSize,:); 

                    BlkBSearch =  FramePreviousStackPad(...
                        Nshift+i+Vsourceblkoffset+1:...
                        Nshift+i+Vsourceblkoffset+iBlkSize,...
                        Nshift+j+Hsourceblkoffset+1-iEdgeSize:...
                        Nshift+j+Hsourceblkoffset,:);

                    BlkOrigSearch =   FramePreviousStackPad(...
                                 Nshift+i +Vsourceblkoffset  +1:...
                                 Nshift+i +Vsourceblkoffset  +iBlkSize,...
                                 Nshift+j +Hsourceblkoffset  +1:...
                                 Nshift+j +Hsourceblkoffset + +iBlkSize,:); 

                    iMADBlkAi16 =  sum(sum(abs(int16(BlkASearch)-int16(BlkAStack))),2);
                    iMADBlkBi16 =  sum(sum(abs(int16(BlkBSearch)-int16(BlkBStack))),2);
                    iMADBlkOrigi16 =  sum(sum(abs(int16(BlkOrigSearch)-int16(BlkOrigStack))),2);

                    
                    aMADBlkOrig(i+Nshift+1,j+Nshift+1,:) = iMADBlkOrigi16;
                    aMADEdge(i+Nshift+1,j+Nshift+1,:) = iMADBlkAi16 + iMADBlkBi16;

                end
             end
                
             %%
             %========collect result of sub-px edge
            [R1,I1] = min(aMADEdge);
            [R2,I2] = min(R1);
            [R3,I3] = min(R2);

            VerticalInd=  I1(:,I2(I3),I3);
            HorizontalInd=I2(I3);
            FrameInd = I3;
            
            VSbpxOffset = double(floor((FrameInd-1)/4))*0.25;
            HSbpxOffset = rem(FrameInd-1,4)*0.25;
            
         
            VectorArrayEdgeSbpxdouble(SourceBlkNoV,SourceBlkNoH,:)= [VerticalInd-Nshift-1+VSbpxOffset,HorizontalInd-Nshift-1+HSbpxOffset];
             
            %========collect result of full-px edge
            [R1, I1] = min(aMADEdge(:,:,1));
            [R2, I2] = min(R1);
            VerticalLoc=I1(I2);
            HorizontalLoc=I2;
            VectorArrayEdgedouble(SourceBlkNoV,SourceBlkNoH,:) = [VerticalLoc-Nshift-1,HorizontalLoc-Nshift-1];
           
           
            %========collect result of full-px blk
            [R1, I1] = min(aMADBlkOrig(:,:,1));
            [R2, I2] = min(R1);
            VerticalLoc=I1(I2);
            HorizontalLoc=I2;
            VectorArrayBlkdouble(SourceBlkNoV,SourceBlkNoH,:) = [VerticalLoc-Nshift-1,HorizontalLoc-Nshift-1];
           
          
        end
         iBlkNo = ((SourceBlkNoV-1)*iHorizontalrun  + SourceBlkNoH) ;
         waitbar(iBlkNo/iTotalBlk,hWB,sprintf('Matching.... Frame:%d/%d, Blk:%d/%d ',m,StopFrame,iBlkNo,iTotalBlk));
   end
   
   
   
     
     %==================DF  ============================

    for SourceBlkNoV = 1:1:iVertivalrun
        for SourceBlkNoH  = 1:1:iHorizontalrun

           Vsourceblkoffset = (SourceBlkNoV-1) * iBlkSize;
           Hsourceblkoffset = (SourceBlkNoH-1) * iBlkSize;
            
           %---------------------
           FramePreviousReconPad = padarray(OrigFramePrevisouui8,[Nshift+iPadLength Nshift+iPadLength],'replicate','both');
           
            
            %%
           %======== obtain edge sbpx blk  ============
           VerticalShift=  VectorArrayEdgeSbpxdouble(SourceBlkNoV,SourceBlkNoH,1);
           HorizontalShift =  VectorArrayEdgeSbpxdouble(SourceBlkNoV,SourceBlkNoH,2);   
           
           VPxShiftd = floor(VerticalShift);
           HPxShiftd = floor(HorizontalShift);
           VPxShift = int8(VPxShiftd );
           HPxShift = int8(HPxShiftd );
           VSbpxShift = VerticalShift - VPxShiftd;
           HSbpxShift = HorizontalShift - HPxShiftd;
           
          BlkEdgeSbpx = uint8(zeros(iBlkSize,iBlkSize));
          for i = 1:iBlkSize
            for j = 1:iBlkSize
                 OrigPxMat = double(FramePreviousReconPad(...
                 Nshift+iPadLength+VPxShiftd +Vsourceblkoffset +i  :...
                 Nshift+iPadLength+VPxShiftd +Vsourceblkoffset +i +1,...
                 Nshift+iPadLength+HPxShiftd +Hsourceblkoffset +j :...
                 Nshift+iPadLength+HPxShiftd +Hsourceblkoffset +j +1));
                
                 MatA = [1-VSbpxShift,VSbpxShift];
                 MatC = [1-HSbpxShift; HSbpxShift];
                 BlkEdgeSbpx(i,j) = uint8(MatA*OrigPxMat*MatC);  
            end
          end
          %================ Orig Edge Blk ========================
            
           VerticalShift = VectorArrayEdgedouble(SourceBlkNoV,SourceBlkNoH,1) ;
           HorizontalShift = VectorArrayEdgedouble(SourceBlkNoV,SourceBlkNoH,2) ;
            
           BlkEdge = FramePreviousReconPad(...
                 Nshift+iPadLength+VerticalShift +Vsourceblkoffset  +1:...
                 Nshift+iPadLength+VerticalShift +Vsourceblkoffset  +iBlkSize,...
                 Nshift+iPadLength+HorizontalShift +Hsourceblkoffset +1:...
                 Nshift+iPadLength+HorizontalShift +Hsourceblkoffset +iBlkSize);
          
          %================ Orig Blk  ========================
            
           VerticalShift = VectorArrayBlkdouble(SourceBlkNoV,SourceBlkNoH,1) ;
           HorizontalShift = VectorArrayBlkdouble(SourceBlkNoV,SourceBlkNoH,2) ;
            
           BlkOrigBlk = FramePreviousReconPad(...
                 Nshift+iPadLength+VerticalShift +Vsourceblkoffset  +1:...
                 Nshift+iPadLength+VerticalShift +Vsourceblkoffset  +iBlkSize,...
                 Nshift+iPadLength+HorizontalShift +Hsourceblkoffset +1:...
                 Nshift+iPadLength+HorizontalShift +Hsourceblkoffset +iBlkSize);
      
             
             %%
            %Place DF
            DisplacedFrameEdgeSbpxui8(...
                         Vsourceblkoffset  +1:...
                         Vsourceblkoffset  +iBlkSize,...
                         Hsourceblkoffset  +1:...
                         Hsourceblkoffset  +iBlkSize)...
                =  BlkEdgeSbpx;
            
            DisplacedFrameEdgeui8(...
                         Vsourceblkoffset  +1:...
                         Vsourceblkoffset  +iBlkSize,...
                         Hsourceblkoffset  +1:...
                         Hsourceblkoffset  +iBlkSize)...
                =  BlkEdge;

            DisplacedFrameBlkui8(...
                         Vsourceblkoffset  +1:...
                         Vsourceblkoffset  +iBlkSize,...
                         Hsourceblkoffset  +1:...
                         Hsourceblkoffset  +iBlkSize)...
                =  BlkOrigBlk;  
            
        end
     iBlkNo = ((SourceBlkNoV-1)*iHorizontalrun  + SourceBlkNoH) ;
     waitbar(iBlkNo/iTotalBlk,hWB,sprintf('Placing Block.... Frame:%d/%d, Blk:%d/%d ',m,StopFrame,iBlkNo,iTotalBlk));
    end
     
%      figure(1);
%      subplot(1,2,1);
%      imshow(DisplacedFrameEdgeui8);
%      
%      subplot(1,2,2);
%      imshow(DisplacedFrameBlkui8);

     
     FrameDiffEdgeSbpxi16 = int16(OrigFrameCurrentui8)-int16(DisplacedFrameEdgeSbpxui8); % get array of difference
     FrameDiffEdgei16 = int16(OrigFrameCurrentui8)-int16(DisplacedFrameEdgeui8); % get array of difference
     FrameDiffBlki16 = int16(OrigFrameCurrentui8)-int16(DisplacedFrameBlkui8); % get array of difference
     
    
     FrameDiffEdgeSbpxd = double(FrameDiffEdgeSbpxi16);
     FrameDiffEdged = double(FrameDiffEdgei16);
     FrameDiffBlkd = double(FrameDiffBlki16);
     
    
     FDEdgeSbpxMSE = sum(sum(FrameDiffEdgeSbpxi16.^2),2) /NPixels
     MSEprofile(m-1,2) =  FDEdgeSbpxMSE;
     FDEdgeMSE = sum(sum(FrameDiffEdgei16.^2),2) /NPixels
     MSEprofile(m-1,1) =  FDEdgeMSE;
     FDBlkMSE = sum(sum(FrameDiffBlki16.^2),2) /NPixels
     MSEprofile(m-1,3) =  FDBlkMSE;
     
    AvgMSEprofile = mean(MSEprofile,1);
    figure(1);
    plot(StartFrame+1:m,MSEprofile(:,1),'-o',...
        StartFrame+1:m,MSEprofile(:,2),'-x',...
        StartFrame+1:m,MSEprofile(:,3),'-+');
    legend('BAP','BAP with Sub-px','FAP','Location','SouthEast');
     title('MSE');
    hold on 
    plot(StartFrame+1:m,ones(1,m-1)*AvgMSEprofile(1),'--',...
        StartFrame+1:m,ones(1,m-1)*AvgMSEprofile(2),'--',...
        StartFrame+1:m,ones(1,m-1)*AvgMSEprofile(3),'--');
    hold off
    
     
   
    %save(strcat(RNAME,'InterSbpx_MSE_',num2str(iBlkSize),'_',num2str(iEdgeSize),'_', TIME,'_',num2str(StartFrame),'_',num2str(StopFrame)),'MSEprofile');
end
close(hWB);







