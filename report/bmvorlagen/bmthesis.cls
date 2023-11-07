%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kenndaten
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{bmvorlagen/bmthesis}[2011/11/11 BM Thesis Class]
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisierung
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage{ifthen} \newboolean{bmcolorlinks} %
\setboolean{bmcolorlinks}{false} %
\newboolean{bmshowl}\setboolean{bmshowl}{false} %
\newboolean{bmgerman}\setboolean{bmgerman}{false} %
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deklarieren von Optionen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\DeclareOption{bmcolorlinks}{\setboolean{bmcolorlinks}{true}}
\DeclareOption{bmshowlabels}{\setboolean{bmshowl}{true}}
\DeclareOption{german}{\setboolean{bmgerman}{true}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{scrbook}}  
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ausf�hren von Optionen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\ProcessOptions\relax   
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Laden von Paketen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\LoadClass{scrbook}%[2001/07/30]   % scrbook-Version mind. 2.8j von 2001/07/30
\usepackage{scrhack}

\usepackage[
	format=hang,
	font={footnotesize, sf},
	labelfont={bf},
	margin=0.5cm,
	aboveskip=5pt,
	belowskip=5pt
	]{caption,subcaption}      										% mehrzeilige Captions ausrichten; subcaption: Untergrafiken 

\usepackage[centertags]{amsmath} 								% AMS-Mathematik, centertags zentriert Nummer bei split
\usepackage{amssymb}             								% zus�tzliche Symbole
\usepackage{lmodern}
\usepackage{icomma}
\usepackage{bm}
\usepackage{graphicx}            								% zum Einbinden von Grafiken
\usepackage[svgnames,table,hyperref]{xcolor}
\usepackage{float}               								% u.a. genaue Plazierung von Gleitobjekten mit H
\usepackage[pdftex,pstarrows]{pict2e}
\usepackage{array}
\usepackage{listings}														% Code-Einbindungsumgebung
%\usepackage{luximono}
\usepackage{courier}       											%verwende Courier statt cmtt als monospace-schrift

\usepackage{setspace}            								% Zeilenabstand einstellbar
\usepackage{rotating}
\onehalfspacing                  								% eineinhalbzeilig einstellen
\usepackage{longtable}													% Erm�glicht Tabellen die �ber den Seitenumbruch gehen (s. Symbolverzeichnis)
\usepackage{scrlayer-scrpage}            				% Kopf und Fusszeilen-Layout passt besser zur Dokumentklasse KOMA-Skript (scrbook) als das Pake fancyhdr, sonst ziemlich gleichwertig
%%% Literatur und sonstige Referenzen 
% \usepackage{natbib}															%braucht man f�r natdin.bst (Literaturzitate mit Namen und Jahr im Text)
% \setcitestyle{square,aysep={}}             			%Eckige Klammern um Zitat

\typearea[current]{current}        							% Neuberechnung des Satzspiegels mit alten Werten nach �nderung von Zeilenabstand,etc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aufrufen anderer bm.sty Dateien
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\usepackage{bmvorlagen/bm_core}

%-------------------------------------------------------------------------
% Aufrufen selber definierter Kommandos
%-------------------------------------------------------------------------
\input{bmvorlagen/bm_commands}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hauptteil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% showlabels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\ifthenelse{\boolean{bmshowl}}{%
\RequirePackage{showkeys}          % gibt die labels mit aus.
\renewcommand{\showkeyslabelformat}[1]{\begin{rotate}{70}\fbox{\normalfont\tiny\ttfamily#1}\end{rotate}}}{}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New- and Renew-Commands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand{\headfont}{\normalfont\sffamily}    % Kolumnentitel serifenlos
\renewcommand{\pnumfont}{\normalfont\ttfamily}    % Seitennummern typewriter und fett 
\pagestyle{scrheadings}
\ihead[]{\headmark}              % Kolumnentitel immer oben innen
\ohead[\pagemark]{\pagemark}     % Seitennummern immer oben aussen
\lefoot[]{} 
\rofoot[]{}                      % Seitennummern in der Fusszeile loeschen
%\pdfminorversion=6


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% German - English
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\ifthenelse{\boolean{bmgerman}}{%

\renewcommand{\bibname}{Literatur} % Literaturverzeichnis wird zu Literatur
\renewcommand{\figurename}{Bild}   % Abbildung wird zu Bild
\renewcommand{\listfigurename}{Bildverzeichnis}%
\usepackage{ngerman}             % neue Rechtschreibung
% \bibliographystyle{bmvorlagen/dinat}
\bibliographystyle{unsrt}
}%
% {\bibliographystyle{bmvorlagen/dinat_e}}
{\bibliographystyle{unsrt}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pdf/hyperref/colorlinks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\pdfoutput=1                  % PDF-Ausgabe

\AtBeginDocument{
  \hypersetup{
    pdftitle = {\bmtitle},%
    pdfauthor = {\bmauthor},%
    pdfkeywords={\bmkeywords},%
		pdfstartview={FitH},%
    pdfstartpage={\bmstartpage}%
  }
}

\ifthenelse{\boolean{bmcolorlinks}}{
\usepackage[pdftex,  % 
     bookmarks,%
     raiselinks,%
     pageanchor,%
     hyperindex,%
     colorlinks,%
     citecolor=green!60!black,%
     linkcolor=red!70!black,%
     urlcolor=magenta!70!black,%
     filecolor=magenta!70!black,%
     menucolor=orange!70!black,%
    ]{hyperref} % 
}{
\usepackage[pdftex,  % 
     bookmarks,%
     raiselinks,%
     pageanchor,%
     hyperindex,%
     hidelinks%
    ]{hyperref} % 
		}



