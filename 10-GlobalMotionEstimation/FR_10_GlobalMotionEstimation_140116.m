
clear all;
close all;
clc;
tic;

%%
%FNAME = '..\\..\\..\\Resource\\Foreman Sequences\\CIF BMP\\foreman_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Waterfall Sequence\\waterfall_cif_%d.bmp';
FNAME = '..\\..\\..\\Resource\\Tempete Sequence\\tempete_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Stefan Sequence\\stefan_cif_%d.bmp';
RNAME= 'Waterfall_';
iFrameHeight =288; 
iFrameWidth =  352;

StartFrame = 1; %Frame No. of the sequence to start
StopFrame = 89;

Nshift = 11;
iBlkSize =8;

iPadLength2 = 5;
iPadLength = 1; 
%%
%============================================
NFrames =StopFrame - StartFrame; %No of Frames to be loaded
iVerticalrun =iFrameHeight/iBlkSize;
iHorizontalrun =iFrameWidth/iBlkSize;
iTotalBlk =iHorizontalrun*iVerticalrun;
NRows =2*Nshift+1;
NPixels = iFrameHeight*iFrameWidth;
MADResult = zeros(NRows,NRows,16);%results in each run
ResultVectorArrayPreviousStackint = int8(zeros(iVerticalrun,iHorizontalrun,3));%array of vector for Sub-pixel

ResultVectorArrayPreviousdouble = zeros(iVerticalrun,iHorizontalrun,2);%array of vector for Sub-pixel
ResultVectorArrayPreviousFullPxint = zeros(iVerticalrun,iHorizontalrun,2); %array of vector for Full-pixel
OrigLumui8 = cell(1,NFrames);
MSEprofile = zeros(NFrames,2);
VectorVprofile = zeros(NFrames,iVerticalrun);
VectorHprofile = zeros(NFrames,iHorizontalrun);

DFStack = uint8( zeros(iFrameHeight,iFrameWidth,2));

NewAccumVGlobalMV =  zeros(iFrameHeight,iFrameWidth);
NewAccumHGlobalMV =  zeros(iFrameHeight,iFrameWidth);
AccumVGlobalMV =  zeros(iFrameHeight,iFrameWidth);
AccumHGlobalMV =  zeros(iFrameHeight,iFrameWidth);

%======================================================

%%
%Load Source
hWB = waitbar(0,'Loading Source...');
for i=StartFrame:StopFrame
    OrigLumui8{i} = rgb2gray(imread(sprintf(FNAME,i))); % read file with conversion to RGB
end




%%
%============================================
for m=StartFrame+1:StopFrame
    OrigFrameCurrentui8 = OrigLumui8{m};
    OrigFramePreviousui8 = OrigLumui8{m-1};
  
    FramePreviousPadui8 = padarray(OrigFramePreviousui8,[Nshift Nshift],'replicate','both'); % Padding of Frame N-1
   %--------------------------------------------
   waitbar(0,hWB,'Matching...'); 
   for SourceBlkNoV = 1:1:iVerticalrun
        for SourceBlkNoH  = 1:1:iHorizontalrun
            Vsourceblkoffset = (SourceBlkNoV-1) * iBlkSize;
            Hsourceblkoffset = (SourceBlkNoH-1) * iBlkSize;
           
            BlkOrig = OrigFrameCurrentui8(...
                Vsourceblkoffset+1:Vsourceblkoffset+iBlkSize,...
                 Hsourceblkoffset+1:Hsourceblkoffset+iBlkSize);

             for i =-Nshift:1:Nshift 
                for j =-Nshift:1:Nshift

                    

                    BlkOrigSearch =  FramePreviousPadui8(...
                                 Nshift+i +Vsourceblkoffset  +1:...
                                 Nshift+i +Vsourceblkoffset  +iBlkSize,...
                                 Nshift+j +Hsourceblkoffset  +1:...
                                 Nshift+j +Hsourceblkoffset + +iBlkSize); 

                    
                    iMADBlkOrigi16 =  sum(sum(abs(int16(BlkOrigSearch)-int16(BlkOrig))),2);

                    aMADBlkOrig(i+Nshift+1,j+Nshift+1) = iMADBlkOrigi16;
             

                end
             end
               % find minimum
             [R1,I1] = min(aMADBlkOrig);
             [~,I2] = min(R1);
             
             if SourceBlkNoH ~= 1 
                 % if I1(I2),I2
                 iMADLowest = aMADBlkOrig(I1(I2),I2);
                 
                 iShiftPreviousV = VectorArray(SourceBlkNoV,SourceBlkNoH-1,1);
                 iShiftPreviousH = VectorArray(SourceBlkNoV,SourceBlkNoH-1,2);
                 
                 iMADPrevious =  aMADBlkOrig(iShiftPreviousV+Nshift+1,iShiftPreviousH+Nshift+1);
                 iMADCenter = aMADBlkOrig(Nshift+1,Nshift+1);
                 iMADTresh = iBlkSize * 4;
                 iMADDiff = abs(iMADLowest - iMADPrevious);

                 if  iMADDiff <= iMADTresh

                 VectorArray(SourceBlkNoV,SourceBlkNoH,:) = [iShiftPreviousV,iShiftPreviousH]; %reset to zero shift

                 else

                 VectorArray(SourceBlkNoV,SourceBlkNoH,:) = [I1(I2)-Nshift-1,I2-Nshift-1];

                 end
             else% first block of the row
                 
                 VectorArray(SourceBlkNoV,SourceBlkNoH,:) = [I1(I2)-Nshift-1,I2-Nshift-1]; 
                 
             end
             
             VectorArray2(SourceBlkNoV,SourceBlkNoH,:)= [I1(I2)-Nshift-1,I2-Nshift-1]; %without matching error elimination
             aQuiverX (SourceBlkNoV,SourceBlkNoH) = Vsourceblkoffset +1;
             aQuiverY (SourceBlkNoV,SourceBlkNoH) = Hsourceblkoffset +1;
        end
        
        iBlkNo = ((SourceBlkNoV-1)*iHorizontalrun  + SourceBlkNoH) ;
        waitbar(iBlkNo/iTotalBlk,hWB,sprintf('Matching... Shift:%d, Frame:%d/%d, Blk:%d/%d ',Nshift,m,StopFrame,iBlkNo,iTotalBlk)); 
       end
   
    
    
    
    %%=====with========

    VectorVprofile(m,:) = rot90(sum(VectorArray(:,:,1)/iHorizontalrun,2));
    VectorHprofile(m,:) = sum(VectorArray(:,:,2)/iVerticalrun,1);

    %%
    %pixel wise matrix
    [c1,c0 ] = polyfit(1:iBlkSize:iFrameHeight,VectorVprofile(m,:),1);
    X = 1:iFrameHeight;
    Y = X*(c1(1))+c1(2);
    regressVMat(m,:) = Y;
    figure(1)
    plot(1:iBlkSize:iFrameHeight,VectorVprofile(m,:));
    hold on
    plot(1:iFrameHeight,regressVMat(m,:),'--');
    title('Vertical Motion Vectors');
    hold off
    
    [c1,c0 ] = polyfit(1:iBlkSize:iFrameWidth,VectorHprofile(m,:),1); 
    X = 1:iFrameWidth;
    Y = X*(c1(1))+c1(2);
    regressHMat(m,:) = Y;
    figure(2)
    plot(1:iBlkSize:iFrameWidth,VectorHprofile(m,:));
    hold on
    plot(1:iFrameWidth,regressHMat(m,:),'--');
    title('Horizontal Motion Vectors');
    hold off
    

    %%
    
    Temp = padarray(OrigFramePreviousui8,[iPadLength iPadLength],'replicate','post'); %pad for sub-px
    FramePreviousSbpxPadui8 = padarray(Temp,[Nshift+iPadLength2 Nshift+iPadLength2],'replicate','both'); %pad for sub-px
    
    GlobalVMotion =repmat(rot90(regressVMat(m,:),3),1,iFrameWidth);
    GlobalHMotion = repmat(regressHMat(m,:),iFrameHeight,1);
    
    for i = 1:iFrameHeight 
        for j = 1:iFrameWidth

           VerticalShift=  GlobalVMotion(i,j);
           HorizontalShift =  GlobalHMotion(i,j);   
           
           %calculate pixel shift and sub-pixel shift
           VPxShiftd = floor(VerticalShift);
           HPxShiftd = floor(HorizontalShift);
           
           VPxShift = int16(VPxShiftd );
           HPxShift = int16(HPxShiftd );
           VSbpxShift = VerticalShift - VPxShiftd;
           HSbpxShift = HorizontalShift - HPxShiftd;
           
           %calculate pixel index
           pixelVIndex = VPxShift+Nshift+iPadLength2 +i;
           pixelHIndex = HPxShift+Nshift+iPadLength2 +j;
           
            %arrang 4 orignal neghbouring pixels in matrix
            MatB = double(FramePreviousSbpxPadui8(...
                 pixelVIndex  :...
                 pixelVIndex +1,...
                 pixelHIndex :...
                 pixelHIndex +1));
            
           
            
                    X = VSbpxShift;
                    Y = HSbpxShift;
                    MatA = [1-X,X];
                    MatC = [1-Y; Y];
           %calculate pixel value with bilinear interpolation         
           OutFrameGMZero(i,j) = uint8(MatA*MatB*MatC); 
      
        end
    end
    
    %===Displace Frame Local Vector===
    
    FramePreviousPadui8 = padarray(OrigFramePreviousui8,[Nshift Nshift],'replicate','both'); % Padding of Frame N-1
    for SourceBlkNoV = 1:1:iVerticalrun
        for SourceBlkNoH  = 1:1:iHorizontalrun

            Vsourceblkoffset = (SourceBlkNoV-1) * iBlkSize;
            Hsourceblkoffset = (SourceBlkNoH-1) * iBlkSize;
            
            %Retrieve Vector map of Previous Frame
            iVerticalVectorPrevious = VectorArray(SourceBlkNoV,SourceBlkNoH,1);
            iHorizontalVectorPrevous = VectorArray(SourceBlkNoV,SourceBlkNoH,2);
            
            %===Obtaining Blocks===
            %obtain Previous blocks from vector map
            BlkFramePrevious =  FramePreviousPadui8(...
                         Nshift+iVerticalVectorPrevious +Vsourceblkoffset +1:...
                         Nshift+iVerticalVectorPrevious +Vsourceblkoffset +iBlkSize,...
                         Nshift+iHorizontalVectorPrevous+Hsourceblkoffset +1:...
                         Nshift+iHorizontalVectorPrevous+Hsourceblkoffset +iBlkSize);
            %-place DF N-1
            FrameDisplacedFramePreviousui8(...
                         Vsourceblkoffset  +1:...
                         Vsourceblkoffset  +iBlkSize,...
                         Hsourceblkoffset  +1:...
                         Hsourceblkoffset  +iBlkSize)...
                =  BlkFramePrevious ;
        end
    end
    
    
    
    %%
    
    
    DFDGM = int16(OutFrameGMZero)-int16(OrigFrameCurrentui8);
    DFDLoc = int16(FrameDisplacedFramePreviousui8)-int16(OrigFrameCurrentui8);
    FD = int16(OrigFramePreviousui8)-int16(OrigFrameCurrentui8);
    
    Out1 =[OrigFrameCurrentui8,OutFrameGMZero,FrameDisplacedFramePreviousui8];
    Out2 =[uint8(FD+127),uint8(DFDGM+127),uint8(DFDLoc+127)];
    
    
    OutCompare = cat(1,Out1,Out2);
    
    figure(3)
    imshow(OutCompare);
%     %imwrite(OutCompare,strcat(RNAME,num2str(m),'_.bmp'),'bmp');
%     
%     
%     
%     
%     if m == StartFrame+1;
% 		imwrite(OutCompare,strcat('Z_',RNAME,'.gif'),'gif','LoopCount',Inf,'DelayTime',0.07);
% 	else
% 		imwrite(OutCompare,strcat('Z_',RNAME,'.gif'),'gif','WriteMode','append','DelayTime',0.07);
%     end
    
end




%%

toc