clc;
clear all;
close all;
%%
%=========User-Defined Varibles============
%C:\\Users\\Chayut\\Dropbox\\Academic File 3.2\\Intern (img)

%FNAME = 'C:\\Users\\Chayut\\Dropbox\\Academic File 3.2\\Intern (img)\\Resource\\Foreman Sequences\\CIF BMP\\foreman_cif_%d.bmp';
FNAME = '..\\..\\..\\Resource\\Foreman Sequences\\CIF BMP\\foreman_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Waterfall Sequence\\waterfall_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Tempete Sequence\\tempete_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Stefan Sequence\\stefan_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Foreman Sequences\\CIF BMP\\foreman_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Akiyo Sequence\\akiyo_qcif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Coastguard Sequence\\cif\\coastguard_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Flower Sequence\\flower_cif_%d.bmp';
RNAME= 'Stefan_';

StartFrame = 148; %Frame No. of the sequence to start
StopFrame = 161;
iBlkSize = 8;
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


NRows =2*Nshift+1;
iPixPerBlk = iBlkSize^2;
iVerticalrun =iFrameHeight/iBlkSize;
iHorizontalrun =iFrameWidth/iBlkSize;
iTotalBlk =iHorizontalrun*iVerticalrun;
NPixelsFrame = iFrameHeight*iFrameWidth;


%%


%======================================================

%%
hWB = waitbar(0,'Loading Source...');
%Load Source
for i=StartFrame:StopFrame
    OrigLumui8{i} = rgb2gray(imread(sprintf(FNAME,i))); % read file with conversion to RGB
end





%============================================

for  m=StartFrame+1:StopFrame-1 % loop for all frame-pair
    
    
        
        
    %%
    %============================================
    VectorArrayPrevious= zeros(iVerticalrun,iHorizontalrun,2);
    VectorArrayNext= zeros(iVerticalrun,iHorizontalrun,2);
    OrigFrameCurrentui8 = OrigLumui8{m};
    OrigFramePrevisouui8 = OrigLumui8{m-1};
    OrigFrameNextui8 = OrigLumui8{m+1};

    FramePreviousPadui8 = padarray(OrigFramePrevisouui8,[Nshift Nshift],'replicate','both'); % Padding of Frame N-1
    FrameNextPadui8 = padarray(OrigFrameNextui8,[Nshift Nshift],'replicate','both'); % Padding of Frame N+1
   
   %--------------------------------------------
   waitbar(0,hWB,'Matching...');   
   for SourceBlkNoV = 1:1:iVerticalrun
        for SourceBlkNoH  = 1:1:iHorizontalrun
            Vsourceblkoffset = (SourceBlkNoV-1) * iBlkSize;
            Hsourceblkoffset = (SourceBlkNoH-1) * iBlkSize;

            BlkOrig = OrigFrameCurrentui8(...
                Vsourceblkoffset+1:Vsourceblkoffset+iBlkSize,...
                 Hsourceblkoffset+1:Hsourceblkoffset+iBlkSize);

             
             
             %===Matching Against Previous Frame===
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

             VectorArrayPrevious(SourceBlkNoV,SourceBlkNoH,:)= [I1(I2)-Nshift-1,I2-Nshift-1]; %without matching error elimination

             %===Matching Against Next Frame
             aMADBlkOrig = zeros(NRows,NRows);
             for i =-Nshift:1:Nshift 
                for j =-Nshift:1:Nshift

                    BlkOrigSearch =  FrameNextPadui8(...
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

             VectorArrayNext(SourceBlkNoV,SourceBlkNoH,:)= [I1(I2)-Nshift-1,I2-Nshift-1]; %without matching error elimination
             
             %=============
             aQuiverX (SourceBlkNoV,SourceBlkNoH) = Vsourceblkoffset +1;%for quiver plotting
             aQuiverY (SourceBlkNoV,SourceBlkNoH) = Hsourceblkoffset +1;%for quiver plotting
        end 
        iBlkNo = ((SourceBlkNoV-1)*iHorizontalrun  + SourceBlkNoH) ;
        waitbar(iBlkNo/iTotalBlk,hWB,sprintf('Matching... Shift:%d, Frame:%d/%d, Blk:%d/%d ',Nshift,m,StopFrame,iBlkNo,iTotalBlk)); 
       end
   %%

   for SourceBlkNoV = 1:1:iVerticalrun
        for SourceBlkNoH  = 1:1:iHorizontalrun

            Vsourceblkoffset = (SourceBlkNoV-1) * iBlkSize;
            Hsourceblkoffset = (SourceBlkNoH-1) * iBlkSize;
            
            %Retrieve Vector map of Previous Frame
            iVerticalVectorPrevious = VectorArrayPrevious(SourceBlkNoV,SourceBlkNoH,1);
            iHorizontalVectorPrevous = VectorArrayPrevious(SourceBlkNoV,SourceBlkNoH,2);
            
            %Retrieve Vector map of Next Frame
            iVerticalVectorNext = VectorArrayNext(SourceBlkNoV,SourceBlkNoH,1);
            iHorizontalVectorNext = VectorArrayNext(SourceBlkNoV,SourceBlkNoH,2);

            %===Obtaining Blocks===
            %obtain Previous blocks from vector map
            BlkFramePrevious =  FramePreviousPadui8(...
                         Nshift+iVerticalVectorPrevious +Vsourceblkoffset +1:...
                         Nshift+iVerticalVectorPrevious +Vsourceblkoffset +iBlkSize,...
                         Nshift+iHorizontalVectorPrevous+Hsourceblkoffset +1:...
                         Nshift+iHorizontalVectorPrevous+Hsourceblkoffset +iBlkSize);
            %obtain Next blocks from vector map
            BlkFrameNext =  FrameNextPadui8(...
                         Nshift+iVerticalVectorNext +Vsourceblkoffset +1:...
                         Nshift+iVerticalVectorNext +Vsourceblkoffset +iBlkSize,...
                         Nshift+iHorizontalVectorNext+Hsourceblkoffset +1:...
                         Nshift+iHorizontalVectorNext+Hsourceblkoffset +iBlkSize);      
                     
           BlkAvg = uint8((uint16(BlkFrameNext)+uint16(BlkFramePrevious))/2);       
           
           
           %===Adaptive Weighting===
           
           BlkOrig = OrigFrameCurrentui8(...
                Vsourceblkoffset+1:Vsourceblkoffset+iBlkSize,...
                 Hsourceblkoffset+1:Hsourceblkoffset+iBlkSize);
           
           BlkStack = cat(3,BlkFramePrevious,BlkFrameNext,BlkAvg);
           BlkOrigStack = repmat(double(BlkOrig),[1 1 3]);
           
           aMADResult = sum(sum((double(BlkStack)-BlkOrigStack).^2,1),2);
           [R1,I1] = min(aMADResult);
           DecisionArray(SourceBlkNoV,SourceBlkNoH) = I1;
           BlkAdt= BlkStack(:,:,I1); %choose Block according to the result
           
           
           %===Quiver Plotting===
           locMVPrevious = DecisionArray==1;
           locMVNext = DecisionArray ==2;
           locMVAVG = DecisionArray ==3;
           
         
            %===Placing Block===
            %-place DF N-1
            FrameDisplacedFramePreviousui8(...
                         Vsourceblkoffset  +1:...
                         Vsourceblkoffset  +iBlkSize,...
                         Hsourceblkoffset  +1:...
                         Hsourceblkoffset  +iBlkSize)...
                =  BlkFramePrevious ;
            FrameDisplacedAvgui8(...
                         Vsourceblkoffset  +1:...
                         Vsourceblkoffset  +iBlkSize,...
                         Hsourceblkoffset  +1:...
                         Hsourceblkoffset  +iBlkSize)...
                =  BlkAvg ;
            FrameDisplacedAdtui8(...
                         Vsourceblkoffset  +1:...
                         Vsourceblkoffset  +iBlkSize,...
                         Hsourceblkoffset  +1:...
                         Hsourceblkoffset  +iBlkSize)...
                =  BlkAdt ;
            
            

        end
        iBlkNo = ((SourceBlkNoV-1)*iHorizontalrun  + SourceBlkNoH) ;
        waitbar(iBlkNo/iTotalBlk,hWB,sprintf('Reconstructing... Shift:%d, Frame:%d/%d, Blk:%d/%d ',Nshift,m,StopFrame,iBlkNo,iTotalBlk)); 
       
   end

   waitbar(1,hWB,'Plotting...');    
   FrameDiffPreviousi16 = int16(OrigFrameCurrentui8)-int16(OrigFramePrevisouui8); % get array of difference
   FrameDiffPreviousd = double(FrameDiffPreviousi16);
   FDPreviousMSE = sum(sum(FrameDiffPreviousd.^2),2) /NPixelsFrame;
   MSEprofile(m-1,1) =  FDPreviousMSE;

   FrameDFDPreviousi16 = int16(OrigFrameCurrentui8)-int16( FrameDisplacedFramePreviousui8); % get array of difference
   FrameDFDPreviousd= double(FrameDFDPreviousi16 );
   DFDMSEPrevious = sum(sum(FrameDFDPreviousd.^2),2) /NPixelsFrame;
   MSEprofile(m-1,2) =  DFDMSEPrevious;

   FrameDFDAvgi16 = int16(OrigFrameCurrentui8)-int16( FrameDisplacedAvgui8); % get array of difference
   FrameDFDAvgd= double(FrameDFDAvgi16 );
   DFDMSEAvg= sum(sum(FrameDFDAvgd.^2),2) /NPixelsFrame;
   MSEprofile(m-1,3) =  DFDMSEAvg;
   
   FrameDFDAdti16 = int16(OrigFrameCurrentui8)-int16( FrameDisplacedAdtui8); % get array of difference
   FrameDFDAdtd= double(FrameDFDAdti16 );
   DFDMSEAdt= sum(sum(FrameDFDAdtd.^2),2) /NPixelsFrame;
   MSEprofile(m-1,4) =  DFDMSEAdt;
   
    
   
    AvgMSEprofile = mean(MSEprofile(StartFrame:m-1,:),1);
    figure(1);
    plot(StartFrame+1:m,MSEprofile(StartFrame:m-1,1),'-o',...
        StartFrame+1:m,MSEprofile(StartFrame:m-1,2),'-x',...
        StartFrame+1:m,MSEprofile(StartFrame:m-1,3),'-*',...
        StartFrame+1:m,MSEprofile(StartFrame:m-1,4),'-+');
    legend('FD','DFD(Previous)','DFD(Avg)','DFD(Adt)','Location','NorthWest');
    title('MSE Graph');
    ylabel('MSE')
    xlabel('Frame Index');
    hold on 
    plot(StartFrame+1:m,ones(1,m - (StartFrame+1)+1)*AvgMSEprofile(1),'--',...
        StartFrame+1:m,ones(1,m - (StartFrame+1)+1)*AvgMSEprofile(2),'--',...
        StartFrame+1:m,ones(1,m - (StartFrame+1)+1)*AvgMSEprofile(3),'--',...
        StartFrame+1:m,ones(1,m - (StartFrame+1)+1)*AvgMSEprofile(4),'--');
    hold off
    
    figure(2);
    plot(StartFrame+1:m,MSEprofile(StartFrame:m-1,2),'-x',...
        StartFrame+1:m,MSEprofile(StartFrame:m-1,3),'-*',...
        StartFrame+1:m,MSEprofile(StartFrame:m-1,4),'-+');
    legend('DFD(Previous)','DFD(Avg)','DFD(Adt)','Location','NorthWest');
    title('MSE Graph');
    ylabel('MSE')
    xlabel('Frame Index');
    hold on 
    plot(StartFrame+1:m,ones(1,m - (StartFrame+1)+1)*AvgMSEprofile(2),'--',...
        StartFrame+1:m,ones(1,m - (StartFrame+1)+1)*AvgMSEprofile(3),'--',...
        StartFrame+1:m,ones(1,m - (StartFrame+1)+1)*AvgMSEprofile(4),'--');
    hold off
    
    figure(3);
    plot(StartFrame+1:m,MSEprofile(StartFrame:m-1,3),'-*',...
        StartFrame+1:m,MSEprofile(StartFrame:m-1,4),'-+');
    legend('DFD(Avg)','DFD(Adt)','Location','NorthWest');
    title('MSE Graph');
    ylabel('MSE')
    xlabel('Frame Index');
    hold on 
    plot(StartFrame+1:m,ones(1,m - (StartFrame+1)+1)*AvgMSEprofile(3),'--',...
        StartFrame+1:m,ones(1,m - (StartFrame+1)+1)*AvgMSEprofile(4),'--');
    hold off
   
   GridOut1 = [FrameDisplacedFramePreviousui8,FrameDisplacedAvgui8,FrameDisplacedAdtui8];
   GridOut2 = [uint8(FrameDFDPreviousi16+127),uint8(FrameDFDAvgi16+127),uint8(FrameDFDAdti16+127)];
   GridOut = cat(1,GridOut1,GridOut2);
   strText = 'Previous, Average, Adaptive';
   RGB = insertText( GridOut, [10,10],strText, 'FontSize', 18, 'FontSize', 18);
   figure(4);
   imshow(RGB);
   
   %===Average Quiver plot overlay===
   FigHandle = figure(5);
   set(FigHandle, 'Position', [100, 100, 1049, 895]);
   imshow(OrigFrameCurrentui8);
   hold on
   quiver(aQuiverY,aQuiverX,-VectorArrayPrevious(:,:,2),-VectorArrayPrevious(:,:,1),'color',[1 0 0]);
   hold on
   quiver(aQuiverY,aQuiverX,-VectorArrayNext(:,:,2),-VectorArrayNext(:,:,1),'color',[0 0 1]);
   
   hold off
   
    FrameCaptureStruct = getframe;  
    CapturedFrameAvg = imresize(FrameCaptureStruct.cdata,[iFrameHeight*2,iFrameWidth*2]);
    
   %===Adaptive Quiver Plot Overlay====
   FigHandle = figure(6);
   set(FigHandle, 'Position', [100, 100, 1049, 895]);
   imshow(OrigFrameCurrentui8);
   hold on
   quiver(aQuiverY.*locMVPrevious,aQuiverX.*locMVPrevious,-VectorArrayPrevious(:,:,2).*locMVPrevious,-VectorArrayPrevious(:,:,1).*locMVPrevious,'color',[1 0 0]);
   hold on
   quiver(aQuiverY.*locMVNext,aQuiverX.*locMVNext,-VectorArrayNext(:,:,2).*locMVNext,-VectorArrayNext(:,:,1).*locMVNext,'color',[0 0 1]);
   hold on
   quiver(aQuiverY.*locMVAVG,aQuiverX.*locMVAVG,-VectorArrayPrevious(:,:,2).*locMVAVG,-VectorArrayPrevious(:,:,1).*locMVAVG,'color',[1 1 0]);
   hold on
   quiver(aQuiverY.*locMVAVG,aQuiverX.*locMVAVG,-VectorArrayNext(:,:,2).*locMVAVG,-VectorArrayNext(:,:,1).*locMVAVG,'color',[0 1 1]);
   hold off

    FrameCaptureStruct = getframe;  
    CapturedFrameAdt = imresize(FrameCaptureStruct.cdata,[iFrameHeight*2,iFrameWidth*2]);
    
   %%===imwrite===
   imwrite(CapturedFrameAdt,strcat('Adt_',RNAME,num2str(Nshift),'_',num2str(m),'.bmp'),'bmp');
  
   [img,idx] = rgb2ind(CapturedFrameAvg,255);
   if m == StartFrame+1;
		imwrite(img,idx,strcat('Adv_',RNAME,'.gif'),'gif','LoopCount',Inf,'DelayTime',0.13);
   else
		imwrite(img,idx,strcat('Adv_',RNAME,'.gif'),'gif','WriteMode','append','DelayTime',0.13);
   end
   
   imwrite(CapturedFrameAdt,strcat('Adv_',RNAME,num2str(iBlkSize),'_',num2str(m),'.bmp'),'bmp');
  
   [img,idx] = rgb2ind(CapturedFrameAdt,255);
   if m == StartFrame+1;
		imwrite(img,idx,strcat('Adt_',RNAME,'.gif'),'gif','LoopCount',Inf,'DelayTime',0.25);
   else
		imwrite(img,idx,strcat('Adt_',RNAME,'.gif'),'gif','WriteMode','append','DelayTime',0.25);
   end
   
 end

close(hWB);
