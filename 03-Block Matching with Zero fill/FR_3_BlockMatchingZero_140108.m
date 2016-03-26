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

iFrameHeight =288; 
iFrameWidth =  352;
StartFrame = 1; %Frame No. of the sequence to start
StopFrame = 40;
DFDFactor = 3; 
Nshift = 11;
iBlkSize =8;
iZeroFillTreashold = 4;

%%
%============================================
TIME = datestr(now,'yymmddHHMM');
iPixPerBlk = iBlkSize^2;
NFrames =StopFrame - StartFrame; %No of Frames to be loaded
iVerticalrun =iFrameHeight/iBlkSize;
iHorizontalrun =iFrameWidth/iBlkSize;
iTotalBlk =iHorizontalrun*iVerticalrun;
NRows =2*Nshift+1;
NPixels = iFrameHeight*iFrameWidth;
OrigLumui8 = cell(1,NFrames);

%%


%======================================================

%%
hWB = waitbar(0,'Matching...');
%Load Source
for i=StartFrame:StopFrame
    OrigLumui8{i} = rgb2gray(imread(sprintf(FNAME,i))); % read file with conversion to RGB
end





%============================================
for  m=StartFrame+1:StopFrame-1 % loop for all frame-pair
    
    OrigFrameCurrentui8 = OrigLumui8{m};
    OrigFramePrevisouui8 = OrigLumui8{m-1};
  
    FramePreviousPadui8 = padarray(OrigFramePrevisouui8,[Nshift Nshift],'replicate','both'); % Padding of Frame N-1
   %--------------------------------------------
   
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
                 iMADTresh = iBlkSize * iZeroFillTreashold;
                 iMADDiff = abs(iMADLowest - iMADPrevious);

                 if  iMADDiff <= iMADTresh

                 VectorArray(SourceBlkNoV,SourceBlkNoH,:) = [iShiftPreviousV,iShiftPreviousH]; %set vector to value of previous block

                 else

                 VectorArray(SourceBlkNoV,SourceBlkNoH,:) = [I1(I2)-Nshift-1,I2-Nshift-1];

                 end
             else% first block of the row
                 
                 VectorArray(SourceBlkNoV,SourceBlkNoH,:) = [I1(I2)-Nshift-1,I2-Nshift-1]; 
                 
             end
             VectorArray2(SourceBlkNoV,SourceBlkNoH,:)= [I1(I2)-Nshift-1,I2-Nshift-1]; %without matching error elimination
        
             aQuiverX (SourceBlkNoV,SourceBlkNoH) = Vsourceblkoffset +1;%for quiver plotting
             aQuiverY (SourceBlkNoV,SourceBlkNoH) = Hsourceblkoffset +1;%for quiver plotting
        end
        
        iBlkNo = ((SourceBlkNoV-1)*iHorizontalrun  + SourceBlkNoH) ;
        waitbar(iBlkNo/iTotalBlk,hWB,sprintf('Matching... Shift:%d, Frame:%d/%d, Blk:%d/%d ',Nshift,m,StopFrame,iBlkNo,iTotalBlk)); 
   end
   %%
   FigHandle = figure(1);
   set(FigHandle, 'Position', [100, 100, 1049, 895]);
   imshow(OrigFrameCurrentui8);
   hold on
   quiver(aQuiverY,aQuiverX,-VectorArray(:,:,2),-VectorArray(:,:,1));
   hold off

   FrameCaptureStruct = getframe;  
   CapturedFrame = imresize(FrameCaptureStruct.cdata,[iFrameHeight*2,iFrameWidth*2]);
   
   
   imwrite(CapturedFrame,strcat('BM_',RNAME,num2str(iBlkSize),'_',num2str(m),'.bmp'),'bmp');
   %%
   [A,map] = rgb2ind(CapturedFrame,256); 
   if m == StartFrame+1;
		imwrite(A,map,strcat('BM_',RNAME,num2str(iBlkSize),'.gif'),'gif','LoopCount',Inf,'DelayTime',0.12);
   else
		imwrite(A,map,strcat('BM_',RNAME,num2str(iBlkSize),'.gif'),'gif','WriteMode','append','DelayTime',0.12);
   end
   
end
close(hWB);
