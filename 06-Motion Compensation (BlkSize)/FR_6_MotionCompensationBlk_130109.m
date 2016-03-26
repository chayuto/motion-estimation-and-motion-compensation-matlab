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


iZeroFillTreashold = 4;

%%
%============================================
TIME = datestr(now,'yymmddHHMM');
NFrames =StopFrame - StartFrame; %No of Frames to be loaded
OrigLumui8 = cell(1,NFrames);

BlkSizeStack = [16,8,4];



%%


%======================================================

%%
hWB = waitbar(0,'Matching...');
%Load Source
for i=StartFrame:StopFrame
    OrigLumui8{i} = rgb2gray(imread(sprintf(FNAME,i))); % read file with conversion to RGB
end





%============================================

for  m=StartFrame+1:StopFrame % loop for all frame-pair
    for r = 1:3
    
        iBlkSize = BlkSizeStack(r);
        NRows =2*Nshift+1;
        iPixPerBlk = iBlkSize^2;
        iVerticalrun =iFrameHeight/iBlkSize;
        iHorizontalrun =iFrameWidth/iBlkSize;
        iTotalBlk =iHorizontalrun*iVerticalrun;
        NPixelsFrame = iFrameHeight*iFrameWidth;
        
        %%
        %============================================
        VectorArray= zeros(iVerticalrun,iHorizontalrun,2);
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

                 aMADBlkOrig = zeros(NRows,NRows);
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

                 VectorArray(SourceBlkNoV,SourceBlkNoH,:)= [I1(I2)-Nshift-1,I2-Nshift-1]; %without matching error elimination


            end

            iBlkNo = ((SourceBlkNoV-1)*iHorizontalrun  + SourceBlkNoH) ;
            waitbar(iBlkNo/iTotalBlk,hWB,sprintf('Matching... Shift:%d, Frame:%d/%d, Blk:%d/%d ',Nshift,m,StopFrame,iBlkNo,iTotalBlk)); 
       end
       %%

       for SourceBlkNoV = 1:1:iVerticalrun
            for SourceBlkNoH  = 1:1:iHorizontalrun
                FramePreviousPadui8 = padarray(OrigFramePrevisouui8,[Nshift Nshift],'replicate','both'); % Padding of Frame N-1
      
                Vsourceblkoffset = (SourceBlkNoV-1) * iBlkSize;
                Hsourceblkoffset = (SourceBlkNoH-1) * iBlkSize;


                iVerticalVectorPrevious = VectorArray(SourceBlkNoV,SourceBlkNoH,1);
                iHorizontalVectorPrevous = VectorArray(SourceBlkNoV,SourceBlkNoH,2);

               
                
                %obtain blocks from vector map
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
       
        
       FrameDiffPreviousi16 = int16(OrigFrameCurrentui8)-int16(OrigFramePrevisouui8); % get array of difference
       FrameDiffPreviousd = double(FrameDiffPreviousi16);
       FDPreviousMSE = sum(sum(FrameDiffPreviousd.^2),2) /NPixelsFrame;
       MSEprofile(m-1,1) =  FDPreviousMSE;

       FrameDFDPreviousi16 = int16(OrigFrameCurrentui8)-int16( FrameDisplacedFramePreviousui8); % get array of difference
       FrameDFDPreviousd= double(FrameDFDPreviousi16 );
       DFDMSEPrevious = sum(sum(FrameDFDPreviousd.^2),2) /NPixelsFrame;
       MSEprofile(m-1,r+1) =  DFDMSEPrevious;
       
       FrameStack(:,:,r) = FrameDisplacedFramePreviousui8;
       DFDFrameStack(:,:,r) = uint8(FrameDFDPreviousi16*DFDFactor+127);
    end

   
    AvgMSEprofile = mean(MSEprofile,1);
    figure(1);
    plot(StartFrame+1:m,MSEprofile(:,1),'-o',...
        StartFrame+1:m,MSEprofile(:,2),'-x',...
        StartFrame+1:m,MSEprofile(:,3),'-*',...
        StartFrame+1:m,MSEprofile(:,4),'-+');
    legend('FD','DFD(16)','DFD(8)','DFD(4)','Location','NorthWest');
    title('MSE');
    hold on 
    plot(StartFrame+1:m,ones(1,m-1)*AvgMSEprofile(1),'--',...
        StartFrame+1:m,ones(1,m-1)*AvgMSEprofile(2),'--',...
        StartFrame+1:m,ones(1,m-1)*AvgMSEprofile(3),'--',...
        StartFrame+1:m,ones(1,m-1)*AvgMSEprofile(4),'--');
    hold off
    FrameCaptureStruct = getframe;  
    CapturedFrame = imresize(FrameCaptureStruct.cdata,[iFrameHeight,iFrameWidth]);
   
   
   GridOut1 = [FrameStack(:,:,1),FrameStack(:,:,2),FrameStack(:,:,3)];
   GridOut2 = [DFDFrameStack(:,:,1),DFDFrameStack(:,:,2),DFDFrameStack(:,:,3)];
   GridOut = cat(1,GridOut1,GridOut2);
   strText = 'Block Size: 4,8 and 16';
   RGB = insertText( GridOut, [10,10],strText, 'FontSize', 18, 'FontSize', 18);
   figure(2);
   imshow(RGB);
   
   imwrite(RGB,strcat('MC2_',RNAME,num2str(Nshift),'_',num2str(m),'.bmp'),'bmp');
  
   [img,idx] = rgb2ind(RGB,255);
   if m == StartFrame+1;
		imwrite(img,idx,strcat('MC2_',RNAME,num2str(Nshift),'.gif'),'gif','LoopCount',Inf,'DelayTime',0.13);
   else
		imwrite(img,idx,strcat('MC2_',RNAME,num2str(Nshift),'.gif'),'gif','WriteMode','append','DelayTime',0.13);
   end
   
 end

close(hWB);
