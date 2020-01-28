%   Clustering über die beiden Dimensionen der Daten unter verwendung der
%   R-Funktion hcluster aus dem Paket amap
% 
% ClusterLarge(E,file,normieren,meth,chipAnnotation,doplot,linkageArrays, linkGenes, clusterdims)
% 
%       file         Beginn der Dateinamen 
%       normieren   0: normale Daten (Default) 
%                   1: Jedes Gen erhält Mittelwert 0 und Standardabweichung 1 
%                   2: Jedes Gen erhält Mittelwert 0 
%                   3: Skaliere global alles zwischen -1 und 1 via
%                   range(y(:))
%                   4: ranks from -1 to 1 
%       meth        the distance measure to be used. This must be one of
%           'euclidean', 'maximum', 'manhattan', 'canberra'
%           'binary' 'pearson', 'correlation' or  'spearman'. Any
%           unambiguous substring can be given.
% 
%       chipAnnotation  Annotation fuer die Arrays, die in der Heatmap mit
%                       geplottet wird
%                       Die Farbe wird äquidistant gewaehlt, die Reihenfolge
%                       nach der Sortierung gewaehlt.
%       doplot      1: Default
%                   0: keine Abbildungen
%                   2: Paper, keine Strings am Dendrogramm, mit colorbar
% 
%       linkageArrays "ward", "single", "complete", "average" (DEFAULT), "mcquitty",
%                       "median" or "centroid","centroid2"
% 
%       clusterdims     Default [1,2]  = both (genes & samples)
%                   1 = genes will be clustered
%                   2 = samples will be clustered
%         
%   global onlyHeatmap  Default: 0
% 
function erg = ClusterLarge(E,file,normieren, meth,chipAnnotation,doplot, linkageArrays, linkGenes, clusterdims)
global onlyHeatmap
global colortype
if(isempty(onlyHeatmap))
    onlyHeatmap = 0;
end
if(~exist('doplot','var') | isempty(doplot))
    doplot = 1;
end
if(~exist('file','var') | isempty(file))
    file = 'Cluster';
end
if(~exist('normieren','var') | isempty(normieren))
    normieren = 0;
end
if(~exist('linkageArrays','var') | isempty(linkageArrays))
    linkageArrays = 'average';
end
if(~exist('linkGenes','var') | isempty(linkGenes))
    linkGenes = 'average';
end
if(~exist('meth','var') | isempty(meth))
    meth = 'correlation';
else
    if(strcmp(meth,'cor')==1)
        meth = 'correlation';
    end
end
if(~exist('clusterdims','var') || isempty(clusterdims))
    clusterdims = 1:2;
end

if contains(file,' ')
    warning('ClusterLarge.m: Blanks are not allowed in the filename. They are now removed.');
    file = deblank(file);
end

erg = [];

if(~exist('chipAnnotation','var') | isempty(chipAnnotation))
    chipAnnotation = [];
    fn_anno = cell(0);
elseif(isstruct(chipAnnotation))
    tmp = chipAnnotation;
    fn_anno = fieldnames(tmp);
    
    chipAnnotation = NaN*ones(length(fn_anno),get(E,'na'));
    for i=1:length(fn_anno)
        if(~iscell(tmp.(fn_anno{i})) & length(tmp.(fn_anno{i})) == get(E,'na'))
            chipAnnotation(i,:) = tmp.(fn_anno{i});
        elseif(iscell(tmp.(fn_anno{i})))
            chipAnnotation(i,:) = cell2levels(tmp.(fn_anno{i}));
%             for j=1:length(tmp.(fn_anno{i}))
%                 chipAnnotation(i,tmp.(fn_anno{i}){j}) = j-1;            
%             end
        elseif(isnumeric(tmp.(fn_anno{i})));
            chipAnnotation(i,:) = tmp.(fn_anno{i});
        end
        chipAnnotation(i,:) = cell2levels(num2strArray(chipAnnotation(i,:)));
%         chipAnnotation(i,:) = 5*(chipAnnotation(i,:)/max(chipAnnotation(i,:))-.5);
    end
elseif(isnumeric(chipAnnotation))
    if(length(chipAnnotation)~=get(E,'na'))
        save error
        error('length of chipAnnotation ~= na')
    end
    fn_anno = {'anno'};
else
    error('chipAnnotation has to be a struct')
end
chipAnnoIn = chipAnnotation;
    
pw = strrep([pwd,filesep],filesep,'/');

sd = get(E,'sd');
if(sum(sd==0)>0)
    warning(['ClusterLarge.m: ',num2str(sum(sd==0)),' genes with SD==0 in data.'])
end
dat = get(E,'data');
if(normieren == 1)
    dat = Normiere(dat);
elseif(normieren == 2)
    dat = Normiere(dat,1);
elseif normieren==3 % Skaliere global alles zwischen -1 und 1
    dat = dat-nanmin(dat(:));
    dat = dat./range(dat(:)) *2 -1;
elseif normieren==4
    for i=1:size(dat,1)
        notnan = find(~isnan(dat(i,:)));
        r = Rankasgn(dat(i,:));
        r = (r-1)/(length(notnan)-1)*2-1; % scaling to -1,1 for ~isnan
        dat(i,notnan) = r(notnan);
    end
end

clusterGenesOK=0;
if(get(E,'ngene')>2)
	openR;
	
	putRdata('SavePath',pw);
	putRdata('DataPath',pw);
	
	putRdata('HybrName',get(E,'hybrname'));
    gn = get(E,'genenames');
    if(isempty(gn))
        gn = get(E,'IDs');        
    end
    putRdata('IDs',gn);
	putRdata('ChipData',dat);
		
	evalR('setwd(SavePath)');
    evalR('require(amap)')

% 	evalR('datumLang <- unlist(strsplit(date()," "))');
% 	evalR('datum     <- paste(datumLang[2],datumLang[3],"-",datumLang[5],sep="",collapse="")');
    evalR('datum     <- format(Sys.time(), "%b%d-%Y")');
	% evalR('memory.limit(size=10000)')
	try
        evalR('IDs <- t(IDs)')
%     evalR('save.image("test.RData")')

%         evalR('ChipDataAll <- ChipData')
        evalR('	anzbreaks <- 100	')
        
        evalR('	potenz <- 1/3	')
    	evalR('minChipData <- min(ChipData,na.rm=TRUE)')
    	evalR('maxChipData <- max(ChipData,na.rm=TRUE)')
    	evalR('plotGrenze <- ceiling(max(  c(abs(minChipData),abs(maxChipData))  )) ')
		evalR('breaks <- seq(-plotGrenze,plotGrenze,len=2*anzbreaks+2) ')
% 		evalR('reds <- rgb(((1:anzbreaks)/anzbreaks)^potenz, g=0,b=0, names=paste("red",1:anzbreaks,sep=".")) ')
% 		evalR('blacks <- rgb(r=0, g=0,b=0, names="black.0") ')
% 		evalR('greens <- rgb(r=0,g=((anzbreaks:1)/anzbreaks)^potenz,b=0, names=paste("red",1:anzbreaks,sep=".")) ')
%         evalR('cols <- c(greens,blacks,reds)	 ')
        R_colors(colortype)
        
    	evalR('ps.options(append=T)')
		evalR('xlim <- c(-0.12*length(ChipData[,1]),length(ChipData[,1])+1) ')
		evalR('ylim <- c(-0.12*length(ChipData[1,]),length(ChipData[1,])+1) ')

       	if(get(E,'narrays')>30)
		    evalR('Hoehe <- 10') %Hoehe der graph. Ausgabe in inches
    	else
		    evalR('Hoehe <- 5');
        end

%         evalR('save.image("test.RData")')
%         evalR(['h <- hcluster(t(ChipData), method = "euclid", link = "average")']);
        if(sum(clusterdims==2)>0)
            evalR(['h <- hcluster(t(ChipData), method = "',meth,'", link = "',linkageArrays,'")']);
            clusterSamplesOK = 1;
        else
            evalR('h<-list(order=1:dim(ChipData)[2])')
            clusterSamplesOK = 0;
        end
        evalR('orderChips <- h$order');
%         evalR('save.image("test.RData")')
        evalR(' memory.limit(size=4000)')
        
        try            
            if(sum(clusterdims==1)>0)
                evalR(['hgene <- hcluster(ChipData, method = "',meth,'", link = "',linkGenes,'")']);
                clusterGenesOK = 1;
            else
                evalR('hgene<-list(order=1:dim(ChipData)[1])')
            end
        catch
            if(get(E,'ngene')>1000)
                warning('Clustering genes failes')
                clusterGenesOK = 0;
                evalR('hgene<-list(order=1:dim(ChipData)[1])')
            else
                error(lasterr)
            end
        end
        evalR('orderGene <- hgene$order')

        if(doplot~=0)
            if(onlyHeatmap~=1)
                if(clusterSamplesOK==1)
                    evalR(['postscript("',file,'_Dendrogram_Arrays.ps", horizontal=TRUE, onefile=TRUE,height=Hoehe, width=10,pointsize=5)'])
                    if(doplot==2)
                        evalR('plot(h,labels=F,main="", xlab="", ylab="", sub="")')
                    else
                        evalR('plot(h,label=HybrName,main="", xlab="", ylab="", sub="")')
                    end
                    evalR('dev.off()')
                end
                
                if(clusterGenesOK==1)
                    evalR(['postscript("',file,'_Dendrogram_Genes.ps", horizontal=TRUE, onefile=TRUE,height=Hoehe, width=10,pointsize=5)'])
                    if(doplot==2 & clusterGenesOK==1)
                        evalR(' plot(hgene,labels=F,main="", xlab="", ylab="", sub="") ')
                    elseif(clusterGenesOK==1)
                        evalR(' plot(hgene,label=IDs,main="", xlab="", ylab="", sub="") ')
                    end
                    evalR('dev.off() ')
                end

                if(doplot==2) % Farbskala
                    evalR('anzSteps <- 100')
                    evalR(['postscript("',file,'_Farbskala.ps", horizontal=TRUE, onefile=TRUE,height=10, width=2)'])
                    evalR('par(lab=c(1,5,1))')
                    if(normieren==1)
                        ylab = '[ 1/SD ]';                
                    elseif normieren == 3
                        ylab = 'a.u. after scaling to [-1,1]';
                    elseif normieren == 4
                        ylab = 'Ranks scaled to [-1,1]';
                    else
                        ylab = 'expression [log2]';
                    end
                    evalR(['image(1,seq(-3,3,len=anzSteps),t(matrix(seq(-3,3,len=anzSteps))),breaks=breaks,col=cols,xlab="",ylab="',ylab,'",add=FALSE,xaxt="n")'])
                    %                 evalR('image(1,c(1:anzSteps)/10-5,t(farben),col=cols,xlab="",ylab="",add=FALSE)')
                    evalR('dev.off()')

                    evalR(['postscript("',file,'_Farbskala2.ps", horizontal=TRUE, onefile=TRUE,height=10, width=2)'])
                    evalR('par(lab=c(1,5,1))')
                    if(normieren==1)
                        ylab = '[ 1/SD ]';
                    else
                        ylab = 'expression [log2]';
                    end
                    evalR(['image(1,seq(-plotGrenze,plotGrenze,len=anzSteps),t(matrix(seq(-plotGrenze,plotGrenze,len=anzSteps))),breaks=breaks,col=cols,xlab="",ylab="',ylab,'",add=FALSE,xaxt="n")'])
                    %                 evalR('image(1,c(1:anzSteps)/10-5,t(farben),col=cols,xlab="",ylab="",add=FALSE)')
                    evalR('dev.off()')
                 setwdR
                saveRimage
               end
            end

            evalR(['postscript("',file,'_Heatmap.ps", horizontal=TRUE, onefile=TRUE,height=Hoehe, width=10,pointsize=5)'])
            if(isempty(chipAnnotation))
                evalR('image(1:length(ChipData[,1]),1:length(ChipData[1,]),ChipData[hgene$order,h$order],col=cols,ylab="",xlab="",breaks=breaks,xlim=xlim,ylim=ylim) ')
%                 saveRimage
                if(get(E,'narrays')<51)
                    evalR('for(i in 1:length(HybrName)){text(-0.1*length(ChipData[,1]),i,  HybrName[h$order[i]],adj=0)}')
                elseif(get(E,'narrays')<101)
                    evalR('for(i in 1:length(HybrName)){text(-0.1*length(ChipData[,1]),i,  HybrName[h$order[i]],adj=0,cex=.5)}')
                elseif(get(E,'narrays')<151)
                    evalR('for(i in 1:length(HybrName)){text(-0.04*length(ChipData[,1]),i,  HybrName[h$order[i]],adj=0,cex=.35)}')
                elseif(get(E,'narrays')<201)
                    evalR('for(i in 1:length(HybrName)){text(-0.04*length(ChipData[,1]),i,  HybrName[h$order[i]],adj=0,cex=.25)}')
                elseif(get(E,'narrays')<401)
                    evalR('for(i in 1:length(HybrName)){text(-0.01*length(ChipData[,1]),i,  HybrName[h$order[i]],adj=0,cex=.1)}')
                end
                if(get(E,'ngene')<101)
                    evalR('for(i in 1:length(IDs)) {text(i,-0.1*length(ChipData[1,]),  IDs[hgene$order[i]], srt=90,adj=0)}')
                elseif(get(E,'ngene')<151)
                    evalR('for(i in 1:length(IDs)) {text(i,-0.1*length(ChipData[1,]),  IDs[hgene$order[i]], srt=90,adj=0,cex=0.75)}')
                elseif(get(E,'ngene')<201)
                    evalR('for(i in 1:length(IDs)) {text(i,-0.1*length(ChipData[1,]),  IDs[hgene$order[i]], srt=90,adj=0,cex=0.55)}')
                end
            else % mit Annotation
                anzlev = NaN(1,size(chipAnnotation,1));
                for i=1:size(chipAnnotation,1)
                    levs = levels(chipAnnotation(i,:));
                    levs = levs(~isnan(levs));
                    anzlev(i) = length(levs);
                end
                maxlev = max(anzlev);
                if(maxlev<5)
                    dc = [0,0,1
                        1,0,1
                        0    1.0000    0.5000
                        1.0000    0.7812    0.4975];
                    dc = dc(1:maxlev,:);
                else
                    dc = defaultcolors(maxlev+4);
                    dc = dc(5:end,:);
                end
                for j=1:maxlev
                    putRdata('momcol',dc(j,:));
                    evalR('cols <- c(cols, rgb(momcol[1],momcol[2],momcol[3]))')
                    evalR('breaks <- c(breaks,breaks[length(breaks)]+1)')
                end

                putRdata('anno',chipAnnotation)
                putRdata('stepSizeX' ,size(dat,1)/500*length(fn_anno));
                evalR('xvals <- 1:length(ChipData[,1])')
            
                putRdata('fn_anno',fn_anno)
                evalR('newdat     <- NULL')
                evalR('anno_xtick <- NULL')
                evalR('anno_xticklabel<- NULL')
                evalR('xvals  <- c(xvals,max(xvals)+1)')
                evalR('newdat <- rbind(newdat,anno[1,]*NaN+plotGrenze-0.5)')
                for i=1:length(fn_anno)
                    lastbreak = rterm('breaks[length(breaks)]');
                    chipAnnotation(i,:) = chipAnnotation(i,:) + 100*eps + lastbreak;

%                     evalR('orderGene <- c(orderGene, length(orderGene)+1)')
                    evalR('xvals <- c(xvals,max(xvals)+stepSizeX)')
                    evalR('xvals <- c(xvals,max(xvals)+stepSizeX)')
                    evalR('anno_xtick <- c(anno_xtick,xvals[length(xvals)])');
                    evalR(['anno_xticklabel <- c(anno_xticklabel,"',fn_anno{i},'")']);
                    evalR(['newdat <- rbind(newdat,anno[',num2str(i),',]*NaN+plotGrenze-0.5)'])
%                     evalR(['newdat <- rbind(newdat,anno[',num2str(i),',]+plotGrenze-0.5)'])
%                     evalR(['tmp <- 4*(rank(anno[',num2str(i),',])/length(anno[',num2str(i),',])-.5)'])
                    evalR(['tmp <- 0.5*plotGrenze*(anno[',num2str(i),',]/nlevels(as.factor(anno[',num2str(i),',]))-.5)'])
                    evalR(['newdat <- rbind(newdat,sign(tmp)*tmp^2)'])
                end
                putRdata('anno',chipAnnotation)
                if(~isempty(chipAnnotation))
%                     evalR('newdat <- rbind(matrix(data=NaN,nrow=ceiling(stepSizeX/2),ncol=length(ChipData[1,])),t(matrix(rep(t(anno),ceiling(stepSizeX/2)),nrow=length(ChipData[1,]))))')
                    evalR('ChipDataAnno <- rbind(ChipData,newdat)')
                    evalR('neworder <- dim(ChipData)[1]+(1:dim(newdat)[1])')
                else
                    evalR('ChipDataAnno <- rbind(ChipData')
                    evalR('neworder <- NULL')
                end

                evalR('xlimAnno <- c(-0.12*length(ChipDataAnno[,1]),xvals[length(xvals)]+stepSizeX) ')
%                 setwdR
%                 saveRimage
                evalR(' image(xvals,1:length(ChipDataAnno[1,]),ChipDataAnno[c(orderGene,neworder),h$order],col=cols, ylab="",xlab="",breaks=breaks ,xlim=xlimAnno ,ylim=ylim) ')
                evalR('IDs_anno <- c(IDs,"",fn_anno)')
                if(get(E,'ngene')<201)
                    evalR('for(i in 1:length(xvals))     {text(xvals[i],-0.1*length(ChipDataAnno[1,]), IDs_anno[orderGene[i]], srt=90,adj=0)}')
                else
                    evalR('for(i in length(ChipData[,1])+2:length(xvals))     {text(xvals[i],-0.1*length(ChipDataAnno[1,]), IDs_anno[orderGene[i]], srt=90,adj=0)}')
                end
                evalR('for(i in 1:length(anno_xtick))  {text(anno_xtick[i],-0.1*length(ChipDataAnno[1,]), anno_xticklabel[i], srt=90,adj=0)}')

            
                if(get(E,'narrays')<51)
                    evalR('for(i in 1:length(HybrName)){text(-0.1*length(ChipData[,1]),i,  HybrName[h$order[i]],adj=0)}')
                elseif(get(E,'narrays')<76)
                    evalR('for(i in 1:length(HybrName)){text(-0.09*length(ChipData[,1]),i,  HybrName[h$order[i]],adj=0,cex=.8)}')
                elseif(get(E,'narrays')<101)
                    evalR('for(i in 1:length(HybrName)){text(-0.08*length(ChipData[,1]),i,  HybrName[h$order[i]],adj=0,cex=.5)}')
                elseif(get(E,'narrays')<201)
                    evalR('for(i in 1:length(HybrName)){text(-0.04*length(ChipData[,1]),i,  HybrName[h$order[i]],adj=0,cex=.25)}')
                elseif(get(E,'narrays')<401)
                    evalR('for(i in 1:length(HybrName)){text(-0.01*length(ChipData[,1]),i,  HybrName[h$order[i]],adj=0,cex=.1)}')
                end
                
            end
            evalR(' dev.off() ')
%         evalR('graphics.off()')
            evalR('ps.options(append=F)')
        end
        
        erg.order = double(getRdata('orderGene'));
        erg.orderChips  = double(getRdata('orderChips'));
        erg.IDs         = get(E,'IDs');
        erg.gene        = get(E,'gene');
        erg.hybrname    = get(E,'hybrname');
        erg.data        = get(E,'data');
        
        erg.ChipAnnotation = chipAnnoIn;
        erg.ChipAnnotation_Name = fn_anno;

        if(clusterGenesOK==1)
            if(onlyHeatmap~=1)
                try
                    erg.anz = cell(1,50);
                    for nk = [2:2:20,25:5:50]
                        evalR(['cut',num2str(nk),' <- cutree(hgene,k=',num2str(nk),')'])
                    end
                    for nk = [2:2:20,25:5:50]
                        eval(['erg.anz{',num2str(nk),'} = double(getRdata(''cut',num2str(nk),'''));'])
                    end
                    erg.anz10 = erg.anz{10};
                    erg.anz20 = erg.anz{20};
                    erg.anz30 = erg.anz{30};
                    erg.anz50 = erg.anz{50};
                catch
                    disp(lasterr)
                end

                try
                    na = get(E,'na');
                    erg.cutArrays = cell(1,na);
                    
%                     for nk = 2:na
%                         evalR(['cutArrays',num2str(nk),' <- cutree(h,k=',num2str(nk),')'])
%                     end
%                     for nk = 2:na
%                         eval(['erg.cutArrays{',num2str(nk),'} = double(getRdata(''cutArrays',num2str(nk),'''));'])
%                     end
                catch
                end
            end
        end
        
        evalR('orderGene <- hgene$order')
        evalR('orderChips <- h$order')

    catch
%         saveRimage('ClusterLarge_catch.rdata')
        disp(lasterr)
%         erg = [];
	end
	closeR;

    if doplot>0
        try
            if ~system(['ps2pdf ',file,'_Heatmap.ps']);
                system(['rm ',file,'_Heatmap.ps']);
            end
            if onlyHeatmap~=1
                if ~system(['ps2pdf ',file,'_Dendrogram_Arrays.ps']);
                    system(['rm ',file,'_Dendrogram_Arrays.ps']);
                end
            end
        end
        if doplot==2
            try
                if ~system(['ps2pdf ',file,'_Farbskala.ps']);
                    system(['rm ',file,'_Farbskala.ps']);
                end
                if ~system(['ps2pdf ',file,'_Farbskala2.ps']);
                    system(['rm ',file,'_Farbskala2.ps']);
                end
            catch
                file
            end
        end
    end
    
    
    if(doplot==1 && clusterGenesOK==1 && onlyHeatmap~=1)
        try
%             WriteResultClustering(erg,[file,'.xls']);
            WriteData(FilterArrays(FilterGene(E,erg.order),erg.orderChips),[file,'.xls'],',');
        catch
            disp(lasterr)
        end
    end
    
else  %% genuegend Gene
    erg = [];
    warning('ClusterLarge.m: Weniger als 3 Gene.')
end


