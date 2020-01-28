% Erzeugt im R-workspace die Variable cols
% Die Variablen anzbreaks und potenz müssen vorhanden sein
% z.B.  anzbreaks <- 100
%       potenz <- 1/3
% 
% Usage Example:
%   global colortype
%   colortype = 'mm2';

function R_colors(colortype)
if(~exist('colortype','var') | isempty(colortype))
    colortype = '';
end

try
switch lower(colortype)
    case ''
        evalR('reds <- rgb(((1:anzbreaks)/anzbreaks)^potenz, g=0,b=0, names=paste("red",1:anzbreaks,sep=".")) ')
        evalR('blacks <- rgb(r=0, g=0,b=0, names="black.0") ')
        evalR('greens <- rgb(r=0,g=((anzbreaks:1)/anzbreaks)^potenz,b=0, names=paste("green",1:anzbreaks,sep=".")) ')
        evalR('cols <- c(greens,blacks,reds)	 ')
        
    case {'rb','redblue','jet'}
        %     evalR('reds <- rgb(((1:anzbreaks)/anzbreaks)^potenz, g=0,b=0, names=paste("red",1:anzbreaks,sep=".")) ')
        %     evalR('blacks <- rgb(r=0, g=0,b=0, names="black.0") ')
        %     evalR('blues <- rgb(r=0,g=0,b=((anzbreaks:1)/anzbreaks)^potenz, names=paste("green",1:anzbreaks,sep=".")) ')
        %     evalR('cols <- c(blues,blacks,reds)	 ')
        
        evalR('require(GA)')
        evalR('maincols <- col2rgb(jet.colors(11))/2^8')
        %     evalR('maincols <- col2rgb(heat.colors(n=11))/2^8');
        %     evalR('maincols <- rbind(seq(maincols[1,1],maincols[1,2],len=11),seq(maincols[2,1],maincols[2,2],len=11),seq(maincols[3,1],maincols[3,2],len=11))')
        evalR('colx <- seq(0,1,len=6)^(1/potenz)');
        evalR('colx <- c(-colx[rev(2:6)],colx)');
        evalR('coly <- seq(-1,1,len=2*anzbreaks+1)');
        evalR('cols <- rgb(r=approx(colx,maincols[1,],coly)$y,g=approx(colx,maincols[2,],coly)$y,b=approx(colx,maincols[3,],coly)$y)');
        
        
    case {'mm'} % colorbrewer & sigmoidal
        evalR('require(RColorBrewer)');
        evalR('maincols <- cbind(col2rgb(brewer.pal(9,"Reds"))[,9:1],rbind(255,255,255),col2rgb(brewer.pal(9,"Blues")))/255');
%         evalR('maincols <- rbind(approx(1:18,maincols18[1,],seq(1,18,len=91))$y,approx(1:18,maincols18[2,],seq(1,18,len=91))$y,approx(1:18,maincols18[3,],seq(1,18,len=91))$y)');
%         evalR('colx <- seq(0,1,len=46)');
%         evalR('colx <- (colx^3) / ((0.01+colx)^3) * ((0.01+1)^3)') % MM
%         evalR(['maincols <- col2rgb(brewer.pal(11,"RdBu"))/2^8']);
%         evalR('addcol1 <- rbind(approx(c(0,1),c(maincols[1,5],maincols[1,6]),0.25)$y,approx(c(0,1),c(maincols[2,5],maincols[2,6]),0.25)$y,approx(c(0,1),c(maincols[3,5],maincols[3,6]),0.25)$y)');
%         evalR('addcol2 <- rbind(approx(c(0,1),c(maincols[1,6],maincols[1,7]),0.25)$y,approx(c(0,1),c(maincols[2,6],maincols[2,7]),0.25)$y,approx(c(0,1),c(maincols[3,6],maincols[3,7]),0.25)$y)');
%         evalR('maincols <- cbind(maincols[,1:5],addcol1,maincols[,6],addcol2,maincols[,7:11])')


        evalR('colx <- c(0,10^seq(-2,0,len=9))#c(0,0.02,0.05,0.1,0.2,.5,1)');
        evalR('colx <- c(-colx[rev(2:length(colx))],colx)');
        evalR('coly <- rev(seq(-1,1,len=2*anzbreaks+1))');        
        evalR('cols <- rgb(r=approx(colx,maincols[1,],coly)$y,g=approx(colx,maincols[2,],coly)$y,b=approx(colx,maincols[3,],coly)$y)');
        
    case {'mm2'} % colorbrewer & sigmoidal
        evalR('require(RColorBrewer)');
        evalR('maincols <- cbind(col2rgb(brewer.pal(9,"Reds"))[,9:1],rbind(255,255,255),col2rgb(brewer.pal(9,"Blues")))/255');

        evalR('colx <-  c(0,10^seq(-2,0,len=9))^.75');
        evalR('colx <- c(-colx[rev(2:length(colx))],colx)');
        evalR('coly <- rev(seq(-1,1,len=2*anzbreaks+1))');        
        evalR('cols <- rgb(r=approx(colx,maincols[1,],coly)$y,g=approx(colx,maincols[2,],coly)$y,b=approx(colx,maincols[3,],coly)$y)');
        
    case {'mm3'} % colorbrewer & sigmoidal
        evalR('require(RColorBrewer)');
        evalR('maincols <- cbind(col2rgb(brewer.pal(9,"Reds"))[,9:1],rbind(255,255,255),col2rgb(brewer.pal(9,"Blues")))/255');

        global mm3_string 
        if(~isempty(mm3_string))
            evalR(['colx <- ',mm3_string]);
        else
            evalR('colx <-  seq(0,1,len=10)^2#c(0,0.02,0.05,0.1,0.2,.5,1)');
        end
        evalR('colx <- c(-colx[rev(2:length(colx))],colx)');
        evalR('coly <- rev(seq(-1,1,len=2*anzbreaks+1))');        
        evalR('cols <- rgb(r=approx(colx,maincols[1,],coly)$y, g=approx(colx,maincols[2,],coly)$y, b=approx(colx,maincols[3,],coly)$y)');
        
    otherwise % colorbrewer
        evalR('require(RColorBrewer)');
        evalR(['maincols <- col2rgb(brewer.pal(11,"',colortype,'"))/2^8']);
        evalR('colx <- seq(0,1,len=6)^(1/(0.5*potenz))');
        evalR('colx <- c(-colx[rev(2:6)],colx)');
        evalR('coly <- rev(seq(-1,1,len=2*anzbreaks+1))');
        evalR('cols <- rgb(r=approx(colx,maincols[1,],coly)$y,g=approx(colx,maincols[2,],coly)$y,b=approx(colx,maincols[3,],coly)$y)');
end

catch
    evalR('load("cols.Rdata"')
end
