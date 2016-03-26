% =====================================================
% Author:Chayut Orapinpatipat 2013 
% Created: 10:25am 24.10.13
% IBPBI Sequence II
% File: IBPBI2_131024_1.m
% generate output
% B frame costruct take from both original frame
%%
%======================================================
close all;
clear all;
clc;
tic;
%============User-defined Varible=====================
%FNAME = '..\\..\\..\\Resource\\Coastguard Sequence\\coastguard_qcif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Foreman Sequences\\QCIF BMP\\foreman_qcif_%d.bmp';
FNAME = '..\\..\\..\\Resource\\Foreman Sequences\\CIF BMP\\foreman_cif_%d.bmp';
%FNAME = '..\\..\\..\\Resource\\Akiyo Sequence\\akiyo_qcif_%d.bmp';

StartFrame = 1; %Frame No. of the sequence to start
StopFrame = 40;
iFrameHeight =288; 
iFrameLength =  352;
Nshift = 11;
iBlkSize =8;


%=====================================================
NFrames =StopFrame - StartFrame; %No of Frames to be loaded
iVertivalrun =iFrameHeight/iBlkSize;
iHorizontalrun =iFrameLength/iBlkSize;
NRows =2*Nshift+1;
MADResult = zeros(NRows,NRows);%results in each run
ResultVectorArrayFrameP = zeros(iVertivalrun,iHorizontalrun,2);%array of vector for Frame P
ResultVectorArrayFrameB1Previous = zeros(iVertivalrun,iHorizontalrun,2);%array of vectors for Frame B2
ResultVectorArrayFrameB1Next = zeros(iVertivalrun,iHorizontalrun,2);%array of vectors for Frame B1
ResultVectorArrayFrameB2Previous = zeros(iVertivalrun,iHorizontalrun,2);%array of vectors for Frame B2
ResultVectorArrayFrameB2Next = zeros(iVertivalrun,iHorizontalrun,2);%array of vectors for Frame B1

NPixelsFrame = iFrameHeight * iFrameLength ;
OrigLumui8 = cell(1,NFrames);
CompareLumui8 = cell(1,NFrames);

MSEprofile = zeros(NFrames-1,4);
%======================================================





%======================================================

%%
%Load Source
for i=StartFrame:StopFrame
    OrigLumui8{i} = rgb2gray(imread(sprintf(FNAME,i))); % read file with conversion to RGB
end


%%    

    for m=StartFrame:StopFrame-4 % loop for all frame-pair

        %----for plotting quiver (with clearing result from previous
        %pass
        QuiverArrayFrameP = zeros(iVertivalrun,iHorizontalrun,4);
        QuiverArrayFrameB1Previous = zeros(iVertivalrun,iHorizontalrun,4);
        QuiverArrayFrameB1Next = zeros(iVertivalrun,iHorizontalrun,4);
        QuiverArrayFrameB1BlendPrevious = zeros(iVertivalrun,iHorizontalrun,4);
        QuiverArrayFrameB1BlendNext = zeros(iVertivalrun,iHorizontalrun,4);
        QuiverArrayFrameB2Previous = zeros(iVertivalrun,iHorizontalrun,4);
        QuiverArrayFrameB2Next = zeros(iVertivalrun,iHorizontalrun,4);
        QuiverArrayFrameB2BlendPrevious = zeros(iVertivalrun,iHorizontalrun,4);
        QuiverArrayFrameB2BlendNext = zeros(iVertivalrun,iHorizontalrun,4);
        
        %-------------------------------------------------------------
        FrameSourceI1ui8 = OrigLumui8{m};
        FrameSourceB1ui8 = OrigLumui8{m+1};
        FrameSourcePui8 = OrigLumui8{m+2};
        FrameSourceB2ui8 = OrigLumui8{m+3};
        FrameSourceI2ui8 = OrigLumui8{m+4};

        
        FrameSourceI1Padui8 = padarray(FrameSourceI1ui8,[Nshift Nshift],'replicate','both'); % Padding of Frame I1
        FrameSourceI2Padui8 = padarray(FrameSourceI2ui8,[Nshift Nshift],'replicate','both'); % Padding of Frame I2
            
        
        %%    
        %==================Vector-map FrameP ====================
        for SourceBlkNoV = 1:1:iVertivalrun
            for SourceBlkNoH  = 1:1:iHorizontalrun
                Vsourceblkoffset = (SourceBlkNoV-1) * iBlkSize;
                Hsourceblkoffset = (SourceBlkNoH-1) * iBlkSize;
                %Obtain FrameN
                BlkSourceP=   FrameSourcePui8(...
                    Vsourceblkoffset+1:...
                    Vsourceblkoffset+iBlkSize,...
                    Hsourceblkoffset+1:...
                    Hsourceblkoffset+iBlkSize);


                %Obtain I1 and P (with shift offset)
                for i =-Nshift:1:Nshift 
                    for j =-Nshift:1:Nshift
 
                        BlkSourceI1 = FrameSourceI1Padui8(...
                             Nshift+i +Vsourceblkoffset  +1:...
                             Nshift+i +Vsourceblkoffset  +iBlkSize,...
                             Nshift+j +Hsourceblkoffset  +1:...
                             Nshift+j +Hsourceblkoffset + +iBlkSize); 

                        %Pixel-wise matching criterion operation 

                         iMADi16 =  abs(int16(BlkSourceP)-int16(BlkSourceI1));

                        MADResult(i+Nshift+1,j+Nshift+1) = sum(sum(iMADi16),2);%sumation of all pixel
                    end
                end

                 [r,c]=find(MADResult(:,:) == min(min( MADResult(:,:) ))); %find location of lowest value of resule array

                 r1=r(1,1);
                 c1=c(1,1);

                ResultVectorArrayFrameP(SourceBlkNoV,SourceBlkNoH,:) = [r1-Nshift-1,c1-Nshift-1]; %store vector on array
                
                %update quiverarray for P Frame         
                QuiverArrayFrameP(SourceBlkNoV,SourceBlkNoH,1) = SourceBlkNoH*iBlkSize;
                QuiverArrayFrameP(SourceBlkNoV,SourceBlkNoH,2) = SourceBlkNoV*iBlkSize;
                QuiverArrayFrameP(SourceBlkNoV,SourceBlkNoH,3) = ResultVectorArrayFrameP(SourceBlkNoV,SourceBlkNoH,1);
                QuiverArrayFrameP(SourceBlkNoV,SourceBlkNoH,4) = ResultVectorArrayFrameP(SourceBlkNoV,SourceBlkNoH,2);
                   
                
             
                %-place DF Frame P
                FrameDisplacedFramePui8(...
                             Vsourceblkoffset  +1:...
                             Vsourceblkoffset  +iBlkSize,...
                             Hsourceblkoffset  +1:...
                             Hsourceblkoffset  +iBlkSize)...
                    =  FrameSourceI1Padui8(...
                             Nshift+ResultVectorArrayFrameP(SourceBlkNoV,SourceBlkNoH,1) +Vsourceblkoffset +1:...
                             Nshift+ResultVectorArrayFrameP(SourceBlkNoV,SourceBlkNoH,1)+Vsourceblkoffset +iBlkSize,...
                             Nshift+ResultVectorArrayFrameP(SourceBlkNoV,SourceBlkNoH,2)+Hsourceblkoffset +1:...
                             Nshift+ResultVectorArrayFrameP(SourceBlkNoV,SourceBlkNoH,2)+Hsourceblkoffset +iBlkSize);
       
             
            end
        end
        %%
        %=================Frame B1========================
        SourceFramePPadui8 = padarray(FrameSourcePui8,[Nshift Nshift],'replicate','both'); % Padding of Frame I2
        
        
        for SourceBlkNoV = 1:1:iVertivalrun
            for SourceBlkNoH  = 1:1:iHorizontalrun
                Vsourceblkoffset = (SourceBlkNoV-1) * iBlkSize;
                Hsourceblkoffset = (SourceBlkNoH-1) * iBlkSize;
                %Obtain FrameN
                BlkB1=   FrameSourceB1ui8(...
                    Vsourceblkoffset+1:...
                    Vsourceblkoffset+iBlkSize,...
                    Hsourceblkoffset+1:...
                    Hsourceblkoffset+iBlkSize);


                %Obtain FrameI1 and B1 (with shift offset)
                for i =-Nshift:1:Nshift 
                    for j =-Nshift:1:Nshift

                        BlkSourceI1 =   FrameSourceI1Padui8(...
                             Nshift+i +Vsourceblkoffset  +1:...
                             Nshift+i +Vsourceblkoffset  +iBlkSize,...
                             Nshift+j +Hsourceblkoffset  +1:...
                             Nshift+j +Hsourceblkoffset + +iBlkSize); 

                        %Pixel-wise matching criterion operation 
                         iMADi16 =  abs(int16(BlkB1)-int16(BlkSourceI1));

                        MADResult(i+Nshift+1,j+Nshift+1) = sum(sum(iMADi16),2);%sumation of all pixel
                    end
                end

                 [r,c]=find(MADResult(:,:) == min(min( MADResult(:,:) ))); %find location of lowest value of resule array

                 r1=r(1,1);
                 c1=c(1,1);

                 ResultVectorArrayFrameB1Previous(SourceBlkNoV,SourceBlkNoH,:) = [r1-Nshift-1,c1-Nshift-1];





                %====Obtain FrameP and B1(with shift offset)

                for i =-Nshift:1:Nshift 
                    for j =-Nshift:1:Nshift

                        BlkSourceFrameP = SourceFramePPadui8(...
                             Nshift+i +Vsourceblkoffset  +1:...
                             Nshift+i +Vsourceblkoffset  +iBlkSize,...
                             Nshift+j +Hsourceblkoffset  +1:...
                             Nshift+j +Hsourceblkoffset + +iBlkSize); 

                        %Pixel-wise matching criterion operation 
                          iMADi16 =  abs(int16(BlkB1)-int16(BlkSourceFrameP));

                        MADResult(i+Nshift+1,j+Nshift+1) = sum(sum(iMADi16),2);%sumation of all pixel
                    end
                end

                 [r,c]=find(MADResult(:,:) == min(min( MADResult(:,:) ))); %find location of lowest value of resule array

                 r1=r(1,1);
                 c1=c(1,1);

                 ResultVectorArrayFrameB1Next(SourceBlkNoV,SourceBlkNoH,:) = [r1-Nshift-1,c1-Nshift-1];
                 
                 %------interpolate from vector map
                    iVerticalVectorPrevious = ResultVectorArrayFrameB1Previous(SourceBlkNoV,SourceBlkNoH,1);
                    iHorizontalVectorPrevous = ResultVectorArrayFrameB1Previous(SourceBlkNoV,SourceBlkNoH,2);
                    iVerticalVectorNext = ResultVectorArrayFrameB1Next(SourceBlkNoV,SourceBlkNoH,1);
                    iHorizontalVectorNext = ResultVectorArrayFrameB1Next(SourceBlkNoV,SourceBlkNoH,2);
                    
                    %obtain blocks from vector map
                      BlkFramePrevious =   FrameSourceI1Padui8(...
                                 Nshift+iVerticalVectorPrevious +Vsourceblkoffset +1:...
                                 Nshift+iVerticalVectorPrevious +Vsourceblkoffset +iBlkSize,...
                                 Nshift+iHorizontalVectorPrevous+Hsourceblkoffset +1:...
                                 Nshift+iHorizontalVectorPrevous+Hsourceblkoffset +iBlkSize);
                      BlkFrameNext =  SourceFramePPadui8(...
                                 Nshift+iVerticalVectorNext+Vsourceblkoffset +1:...
                                 Nshift+iVerticalVectorNext+Vsourceblkoffset +iBlkSize,...
                                 Nshift+iHorizontalVectorNext+Hsourceblkoffset +1:...
                                 Nshift+iHorizontalVectorNext+Hsourceblkoffset +iBlkSize);
                             
                      BlkFrameBlend = (BlkFrameNext+ BlkFramePrevious)/2;
                      BlkFrameOri =    FrameSourceB1ui8(...
                        Vsourceblkoffset+1:...
                        Vsourceblkoffset+iBlkSize,...
                        Hsourceblkoffset+1:...
                        Hsourceblkoffset+iBlkSize);

                    iMADBlkPreviousi16 =  sum(sum(abs(int16(BlkFramePrevious)-int16(BlkFrameOri))));
                    iMADBlkNexti16  =  sum(sum(abs(int16(BlkFrameNext)-int16(BlkFrameOri))));      
                    iMADBlkBlendi16 =  sum(sum(abs(int16(BlkFrameBlend)-int16(BlkFrameOri))));
                    
                    
                    if((iMADBlkPreviousi16<=iMADBlkBlendi16) & (iMADBlkPreviousi16<=iMADBlkNexti16)) %if N-1 is lowest
                        BlkInter = BlkFramePrevious;
                        QuiverArrayFrameB1Previous(SourceBlkNoV,SourceBlkNoH,1) = SourceBlkNoH*iBlkSize;
                        QuiverArrayFrameB1Previous(SourceBlkNoV,SourceBlkNoH,2) = SourceBlkNoV*iBlkSize;
                        QuiverArrayFrameB1Previous(SourceBlkNoV,SourceBlkNoH,3) = ResultVectorArrayFrameB1Previous(SourceBlkNoV,SourceBlkNoH,1);
                        QuiverArrayFrameB1Previous(SourceBlkNoV,SourceBlkNoH,4) = ResultVectorArrayFrameB1Previous(SourceBlkNoV,SourceBlkNoH,2);
                    elseif (( iMADBlkNexti16  <= iMADBlkBlendi16) & (iMADBlkNexti16  <= iMADBlkPreviousi16)) %if N+1 is lowest
                        BlkInter = BlkFrameNext;
                        QuiverArrayFrameB1Next(SourceBlkNoV,SourceBlkNoH,1) = SourceBlkNoH*iBlkSize;
                        QuiverArrayFrameB1Next(SourceBlkNoV,SourceBlkNoH,2) = SourceBlkNoV*iBlkSize;
                        QuiverArrayFrameB1Next(SourceBlkNoV,SourceBlkNoH,3) = ResultVectorArrayFrameB1Next(SourceBlkNoV,SourceBlkNoH,1);
                        QuiverArrayFrameB1Next(SourceBlkNoV,SourceBlkNoH,4) = ResultVectorArrayFrameB1Next(SourceBlkNoV,SourceBlkNoH,2);            
                    else
                        BlkInter = BlkFrameBlend;
                        QuiverArrayFrameB1BlendPrevious(SourceBlkNoV,SourceBlkNoH,1) = SourceBlkNoH*iBlkSize;
                        QuiverArrayFrameB1BlendPrevious(SourceBlkNoV,SourceBlkNoH,2) = SourceBlkNoV*iBlkSize;
                        QuiverArrayFrameB1BlendPrevious(SourceBlkNoV,SourceBlkNoH,3) = ResultVectorArrayFrameB1Previous(SourceBlkNoV,SourceBlkNoH,1);
                        QuiverArrayFrameB1BlendPrevious(SourceBlkNoV,SourceBlkNoH,4) = ResultVectorArrayFrameB1Previous(SourceBlkNoV,SourceBlkNoH,2);
                        QuiverArrayFrameB1BlendNext(SourceBlkNoV,SourceBlkNoH,1) = SourceBlkNoH*iBlkSize;
                        QuiverArrayFrameB1BlendNext(SourceBlkNoV,SourceBlkNoH,2) = SourceBlkNoV*iBlkSize;
                        QuiverArrayFrameB1BlendNext(SourceBlkNoV,SourceBlkNoH,3) = ResultVectorArrayFrameB1Next(SourceBlkNoV,SourceBlkNoH,1);
                        QuiverArrayFrameB1BlendNext(SourceBlkNoV,SourceBlkNoH,4) =ResultVectorArrayFrameB1Next(SourceBlkNoV,SourceBlkNoH,2);
                        
                    end

                    
                     FrameDisplacedFrameB1ui8(...
                                 Vsourceblkoffset  +1:...
                                 Vsourceblkoffset  +iBlkSize,...
                                 Hsourceblkoffset  +1:...
                                 Hsourceblkoffset  +iBlkSize)...  
                                 =BlkInter;   
                    
            end
        end 
        %%
        %=======Frame B2======================
        for SourceBlkNoV = 1:1:iVertivalrun
            for SourceBlkNoH  = 1:1:iHorizontalrun
                Vsourceblkoffset = (SourceBlkNoV-1) * iBlkSize;
                Hsourceblkoffset = (SourceBlkNoH-1) * iBlkSize;
                %Obtain FrameN
                BlkB2=   FrameSourceB2ui8(...
                    Vsourceblkoffset+1:...
                    Vsourceblkoffset+iBlkSize,...
                    Hsourceblkoffset+1:...
                    Hsourceblkoffset+iBlkSize);


                %Obtain FrameP and B2 (with shift offset)
                 for i =-Nshift:1:Nshift 
                    for j =-Nshift:1:Nshift

                        BlkSourceFrameP = SourceFramePPadui8(...
                             Nshift+i +Vsourceblkoffset  +1:...
                             Nshift+i +Vsourceblkoffset  +iBlkSize,...
                             Nshift+j +Hsourceblkoffset  +1:...
                             Nshift+j +Hsourceblkoffset + +iBlkSize); 

                        %Pixel-wise matching criterion operation 
                          iMADi16 =  abs(int16(BlkB2)-int16(BlkSourceFrameP));

                        MADResult(i+Nshift+1,j+Nshift+1) = sum(sum(iMADi16),2);%sumation of all pixel
                    end
                end
                 [r,c]=find(MADResult(:,:) == min(min( MADResult(:,:) ))); %find location of lowest value of resule array

                 r1=r(1,1);
                 c1=c(1,1);

                 ResultVectorArrayFrameB2Previous(SourceBlkNoV,SourceBlkNoH,:) = [r1-Nshift-1,c1-Nshift-1];





                %====Obtain FrameI2 and B1(with shift offset)
                for i =-Nshift:1:Nshift 
                    for j =-Nshift:1:Nshift

                        BlkSourceI2 =   FrameSourceI2Padui8(...
                             Nshift+i +Vsourceblkoffset  +1:...
                             Nshift+i +Vsourceblkoffset  +iBlkSize,...
                             Nshift+j +Hsourceblkoffset  +1:...
                             Nshift+j +Hsourceblkoffset + +iBlkSize); 

                        %Pixel-wise matching criterion operation 
                         iMADi16 =  abs(int16(BlkB2)-int16(BlkSourceI2));

                        MADResult(i+Nshift+1,j+Nshift+1) = sum(sum(iMADi16),2);%sumation of all pixel
                    end
                end
               

                 [r,c]=find(MADResult(:,:) == min(min( MADResult(:,:) ))); %find location of lowest value of resule array

                 r1=r(1,1);
                 c1=c(1,1);

                 ResultVectorArrayFrameB2Next(SourceBlkNoV,SourceBlkNoH,:) = [r1-Nshift-1,c1-Nshift-1];
                 
                 %------interpolate from vector map
                    iVerticalVectorPrevious = ResultVectorArrayFrameB2Previous(SourceBlkNoV,SourceBlkNoH,1);
                    iHorizontalVectorPrevous = ResultVectorArrayFrameB2Previous(SourceBlkNoV,SourceBlkNoH,2);
                    iVerticalVectorNext = ResultVectorArrayFrameB2Next(SourceBlkNoV,SourceBlkNoH,1);
                    iHorizontalVectorNext = ResultVectorArrayFrameB2Next(SourceBlkNoV,SourceBlkNoH,2);
                    
                    %obtain blocks from vector map
                      BlkFramePrevious =  SourceFramePPadui8(...
                                 Nshift+iVerticalVectorPrevious +Vsourceblkoffset +1:...
                                 Nshift+iVerticalVectorPrevious +Vsourceblkoffset +iBlkSize,...
                                 Nshift+iHorizontalVectorPrevous+Hsourceblkoffset +1:...
                                 Nshift+iHorizontalVectorPrevous+Hsourceblkoffset +iBlkSize);
                      BlkFrameNext =  FrameSourceI2Padui8(...
                                 Nshift+iVerticalVectorNext+Vsourceblkoffset +1:...
                                 Nshift+iVerticalVectorNext+Vsourceblkoffset +iBlkSize,...
                                 Nshift+iHorizontalVectorNext+Hsourceblkoffset +1:...
                                 Nshift+iHorizontalVectorNext+Hsourceblkoffset +iBlkSize);
                             
                      BlkFrameBlend = (BlkFrameNext+ BlkFramePrevious)/2;
                      BlkFrameOri =    FrameSourceB2ui8(...
                        Vsourceblkoffset+1:...
                        Vsourceblkoffset+iBlkSize,...
                        Hsourceblkoffset+1:...
                        Hsourceblkoffset+iBlkSize);

                    iMADBlkPreviousi16 =  sum(sum(abs(int16(BlkFramePrevious)-int16(BlkFrameOri))));
                    iMADBlkNexti16  =  sum(sum(abs(int16(BlkFrameNext)-int16(BlkFrameOri))));      
                    iMADBlkBlendi16 =  sum(sum(abs(int16(BlkFrameBlend)-int16(BlkFrameOri))));
                    
                    
                    if((iMADBlkPreviousi16<=iMADBlkBlendi16) & (iMADBlkPreviousi16<=iMADBlkNexti16)) %if N-1 is lowest
                        BlkInter = BlkFramePrevious;
                        QuiverArrayFrameB2Previous(SourceBlkNoV,SourceBlkNoH,1) = SourceBlkNoH*iBlkSize;
                        QuiverArrayFrameB2Previous(SourceBlkNoV,SourceBlkNoH,2) = SourceBlkNoV*iBlkSize;
                        QuiverArrayFrameB2Previous(SourceBlkNoV,SourceBlkNoH,3) = ResultVectorArrayFrameB2Previous(SourceBlkNoV,SourceBlkNoH,1);
                        QuiverArrayFrameB2Previous(SourceBlkNoV,SourceBlkNoH,4) = ResultVectorArrayFrameB2Previous(SourceBlkNoV,SourceBlkNoH,2);
                    elseif (( iMADBlkNexti16  <= iMADBlkBlendi16) & (iMADBlkNexti16  <= iMADBlkPreviousi16)) %if N+1 is lowest
                        BlkInter = BlkFrameNext;
                        QuiverArrayFrameB2Next(SourceBlkNoV,SourceBlkNoH,1) = SourceBlkNoH*iBlkSize;
                        QuiverArrayFrameB2Next(SourceBlkNoV,SourceBlkNoH,2) = SourceBlkNoV*iBlkSize;
                        QuiverArrayFrameB2Next(SourceBlkNoV,SourceBlkNoH,3) = ResultVectorArrayFrameB2Next(SourceBlkNoV,SourceBlkNoH,1);
                        QuiverArrayFrameB2Next(SourceBlkNoV,SourceBlkNoH,4) = ResultVectorArrayFrameB2Next(SourceBlkNoV,SourceBlkNoH,2);            
                    else
                        BlkInter = BlkFrameBlend;
                        QuiverArrayFrameB2BlendPrevious(SourceBlkNoV,SourceBlkNoH,1) = SourceBlkNoH*iBlkSize;
                        QuiverArrayFrameB2BlendPrevious(SourceBlkNoV,SourceBlkNoH,2) = SourceBlkNoV*iBlkSize;
                        QuiverArrayFrameB2BlendPrevious(SourceBlkNoV,SourceBlkNoH,3) = ResultVectorArrayFrameB2Previous(SourceBlkNoV,SourceBlkNoH,1);
                        QuiverArrayFrameB2BlendPrevious(SourceBlkNoV,SourceBlkNoH,4) = ResultVectorArrayFrameB2Previous(SourceBlkNoV,SourceBlkNoH,2);
                        QuiverArrayFrameB2BlendNext(SourceBlkNoV,SourceBlkNoH,1) = SourceBlkNoH*iBlkSize;
                        QuiverArrayFrameB2BlendNext(SourceBlkNoV,SourceBlkNoH,2) = SourceBlkNoV*iBlkSize;
                        QuiverArrayFrameB2BlendNext(SourceBlkNoV,SourceBlkNoH,3) = ResultVectorArrayFrameB2Next(SourceBlkNoV,SourceBlkNoH,1);
                        QuiverArrayFrameB2BlendNext(SourceBlkNoV,SourceBlkNoH,4) =ResultVectorArrayFrameB2Next(SourceBlkNoV,SourceBlkNoH,2);
                        
                    end
                    
                     FrameDisplacedFrameB2ui8(...
                                 Vsourceblkoffset  +1:...
                                 Vsourceblkoffset  +iBlkSize,...
                                 Hsourceblkoffset  +1:...
                                 Hsourceblkoffset  +iBlkSize)...  
                                 =BlkInter;   
                    
            end
        end 
       
        %%
        %==============calculate Frame Difference DFD=================
        FrameDiffB1i16  = int16(FrameSourceB1ui8)-int16(FrameDisplacedFrameB1ui8); % get array of difference
        FrameDiffPi16  = int16(FrameSourcePui8)-int16(FrameDisplacedFramePui8); % get array of difference
        FrameDiffB2i16  = int16(FrameSourceB2ui8)-int16(FrameDisplacedFrameB2ui8); % get array of difference
        
        FrameDiffB1d = sum(sum(FrameDiffB1i16.^2),2) /NPixelsFrame;
        FrameDiffPd = sum(sum(FrameDiffPi16.^2),2) /NPixelsFrame;
        FrameDiffB2d = sum(sum(FrameDiffB2i16.^2),2) /NPixelsFrame;
        
        MSEprofile(m,1) =  FrameDiffB1d;
        MSEprofile(m,2) =  FrameDiffPd;
        MSEprofile(m,3) =  FrameDiffB2d;
        
        figure(1);
        plot(MSEprofile)
        %%
        %=====================Grid Display===========================
        
        Grid = cat(2, FrameSourceI1ui8,...
            uint8(FrameDiffB1i16+127),...
            uint8(FrameDiffPi16+127),...
            uint8(FrameDiffB2i16+127),...
            FrameSourceI2ui8);
        
        
        %%
        %---write output file

        prefix_image='GridOut_';
        fileformat='.bmp';
        image_index=m;

        X= strcat(prefix_image,num2str(image_index),fileformat);
        imwrite(Grid,X,'bmp');
      save(strcat('MSE_',num2str(StartFrame),'_',num2str(StopFrame)),'MSEprofile'); 
  
    end
toc



