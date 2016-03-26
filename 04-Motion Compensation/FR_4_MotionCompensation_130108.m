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
DFDFactor = 1; 
Nshift = 7;
iBlkSize =16;
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
NPixelsFrame = iFrameHeight*iFrameWidth;
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
for  m=StartFrame+1:StopFrame % loop for all frame-pair
    
    OrigFrameCurrentui8 = OrigLumui8{m};
    OrigFramePreviousui8 = OrigLumui8{m-1};
  
    FramePreviousPadui8 = padarray(OrigFramePreviousui8,[Nshift Nshift],'replicate','both'); % Padding of Frame N-1
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
 
             VectorArray(SourceBlkNoV,SourceBlkNoH,:)= [I1(I2)-Nshift-1,I2-Nshift-1]; %without matching error elimination
        
            
        end
        
        iBlkNo = ((SourceBlkNoV-1)*iHorizontalrun  + SourceBlkNoH) ;
        waitbar(iBlkNo/iTotalBlk,hWB,sprintf('Matching... Shift:%d, Frame:%d/%d, Blk:%d/%d ',Nshift,m,StopFrame,iBlkNo,iTotalBlk)); 
       end
   %%
   
   for SourceBlkNoV = 1:1:iVerticalrun
        for SourceBlkNoH  = 1:1:iHorizontalrun
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
   
   FrameDiffPreviousi16 = int16(OrigFrameCurrentui8)-int16(OrigFramePreviousui8); % get array of difference
   FrameDiffPreviousd = double(FrameDiffPreviousi16);
   FDPreviousMSE = sum(sum(FrameDiffPreviousd.^2),2) /NPixelsFrame;
   MSEprofile(m-1,1) =  FDPreviousMSE;
   
   FrameDFDPreviousi16 = int16(OrigFrameCurrentui8)-int16( FrameDisplacedFramePreviousui8); % get array of difference
   FrameDFDPreviousd= double(FrameDFDPreviousi16 );
   DFDMSEPrevious = sum(sum(FrameDFDPreviousd.^2),2) /NPixelsFrame;
   MSEprofile(m-1,2) =  DFDMSEPrevious;
   
   
    AvgMSEprofile = mean(MSEprofile(StartFrame:m-1,:),1)
    figure(1);
    plot(StartFrame+1:m,MSEprofile(StartFrame:m-1,1),'-o',...
        StartFrame+1:m,MSEprofile(StartFrame:m-1,2),'-x');
    legend('FD','DFD','Location','NorthWest');
    title('MSE Graph');
    ylabel('MSE')
    xlabel('Frame Index');
    hold on 
    plot(StartFrame+1:m,ones(1,m - (StartFrame+1)+1)*AvgMSEprofile(1),'--',...
        StartFrame+1:m,ones(1,m - (StartFrame+1)+1)*AvgMSEprofile(2),'--');
    hold off
    FrameCaptureStruct = getframe;  
    CapturedFrame = imresize(FrameCaptureStruct.cdata,[iFrameHeight,iFrameWidth]);
   
   
   %GridOut = [OrigFrameCurrentui8,FrameDisplacedFramePreviousui8,uint8(FrameDiffPreviousi16*DFDFactor+127),uint8(FrameDFDPreviousi16*DFDFactor+127)];
   %strText= strcat(,', );
   
   strText= strcat('Original');
   RGBOrig = insertText( OrigFrameCurrentui8, [10,10],strText, 'FontSize', 18, 'FontSize', 18);
    strText= strcat('Displaced Frame');
   RGBDF = insertText( FrameDisplacedFramePreviousui8, [10,10],strText, 'FontSize', 18, 'FontSize', 18);
    strText= strcat('FD: ',num2str(FDPreviousMSE,'%3.2f'));
   RGBFD = insertText( uint8(FrameDiffPreviousi16*DFDFactor+127), [10,10],strText, 'FontSize', 18, 'FontSize', 18);
    strText= strcat('DFD: ', num2str(DFDMSEPrevious,'%3.2f'));
   RGBDFD = insertText( uint8(FrameDFDPreviousi16*DFDFactor+127), [10,10],strText, 'FontSize', 18, 'FontSize', 18);
   
   
   
   RGB = cat(2,RGBOrig,RGBDF,RGBFD,RGBDFD,CapturedFrame);
   
   
   
   figure(2);
   imshow(RGB);
   
   imwrite(RGB,strcat('MC_',RNAME,num2str(iBlkSize),'_',num2str(Nshift),'_',num2str(m),'.bmp'),'bmp');
%    imwrite(OrigFrameCurrentui8,strcat('140123 -','N','.bmp'),'bmp');
%    imwrite(OrigFramePreviousui8,strcat('140123 -','N-1','.bmp'),'bmp');
%    imwrite(FrameDisplacedFramePreviousui8,strcat('140123 -','DF','.bmp'),'bmp');
%    imwrite(uint8(FrameDiffPreviousi16*DFDFactor+127),strcat('140123 -','FD','.bmp'),'bmp');
%    imwrite(uint8(FrameDFDPreviousi16*DFDFactor+127),strcat('140123 -','DFD','.bmp'),'bmp');
   
   [img,idx] = rgb2ind(RGB,255);
   if m == StartFrame+1;
		imwrite(img,idx,strcat('MC_',RNAME,num2str(iBlkSize),'_',num2str(Nshift),'.gif'),'gif','LoopCount',Inf,'DelayTime',0.13);
   else
		imwrite(img,idx,strcat('MC_',RNAME,num2str(iBlkSize),'_',num2str(Nshift),'.gif'),'gif','WriteMode','append','DelayTime',0.13);
   end
   
   
    
    
   
  
    
    
    
    
    
   
end
close(hWB);
