%---------- Kenndaten
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{bmvorlagen/bmbericht}[2011/11/11 BM Bericht Klasse TUM]
%
%---------- Initialisierung
\RequirePackage{ifthen} \newboolean{bmcolorlinks} %
\setboolean{bmcolorlinks}{false} %
\newboolean{bmshowl}\setboolean{bmshowl}{false} %
\newboolean{bmgerman}\setboolean{bmgerman}{false} %
%
%---------- Deklarieren von Optionen
\DeclareOption{bmcolorlinks}{\setboolean{bmcolorlinks}{true}}
\DeclareOption{bmshowlabels}{\setboolean{bmshowl}{true}}
\DeclareOption{german}{\setboolean{bmgerman}{true}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{scrreprt}}  
%
%---------- Ausführen von Optionen
\ProcessOptions\relax   
%
%---------- Laden von Paketen
\LoadClass{scrreprt}[2001/07/30]   % scrbook-Version mind. 2.8j von 2001/07/30

\usepackage[
	format=hang,
	font={footnotesize, sf},
	labelfont={bf},
	margin=0.5cm,
	aboveskip=5pt,
	belowskip=5pt
	]{caption,subfig}      % mehrzeilige Captions ausrichten; subfig: untergrafiken 

\usepackage[centertags]{amsmath} % AMS-Mathematik, centertags zentriert Nummer bei split
\usepackage{amssymb}             % zus�tzliche Symbole
\usepackage{icomma}
\usepackage{bm}
\usepackage{graphicx}            % zum Einbinden von Grafiken
\usepackage[svgnames,table,hyperref]{xcolor}
\usepackage{float}               % u.a. genaue Plazierung von Gleitobjekten mit H
\usepackage[pdftex,pstarrows]{pict2e}
\usepackage{array}
\usepackage{listings}						% Code-Einbindungsumgebung
%\usepackage{luximono}
\usepackage{courier}       %verwende Courier statt cmtt als monospace-schrift

\usepackage{setspace}            % Zeilenabstand einstellbar
\usepackage{rotating}
\onehalfspacing                  % eineinhalbzeilig einstellen
\usepackage{longtable}					% Erm�glicht Tabellen die �ber den Seitenumbruch gehen (s. Symbolverzeichnis)
\usepackage{scrpage2}            % Kopf und Fusszeilen-Layout passt besser zur Dokumentklasse KOMA-Skript (scrbook) als das Pake fancyhdr, sonst ziemlich gleichwertig

%%% Literatur und sonstige Referenzen 
\usepackage{natbib}								%braucht man f�r natdin.bst (Literaturzitate mit Namen und Jahr im Text)
%\setcitestyle{square,aysep={}}             %Eckige Klammern um Zitat
%
%\typearea[current]{current}        % Neuberechnung des Satzspiegels mit alten Werten nach �nderung von Zeilenabstand,etc


%% LADEN VON BMCORE
\usepackage{bmvorlagen/bm_core}

%-------------------------------------------------------------------------
% Aufrufen selber definierter Kommandos
%-------------------------------------------------------------------------
\include{bmvorlagen/bm_commands}


\setcitestyle{square,aysep={}}             %Eckige Klammern um Zitat

\typearea[current]{current}        % Neuberechnung des Satzspiegels mit alten Werten nach �nderung von Zeilenabstand,etc
\headheight = 1.3cm
\setlength{\mathindent}{0.5cm}
%\addtolength {\textheight}{-1.5cm}
%\addtolength {\footskip}{1.5cm}
%
%---------- Hauptteil

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% showlabels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\ifthenelse{\boolean{bmshowl}}{%
\RequirePackage{showkeys}          % gibt die labels mit aus.
\renewcommand{\showkeyslabelformat}[1]{\begin{rotate}{70}\fbox{\normalfont\tiny\ttfamily#1}\end{rotate}}}{}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New- and Renew-Commands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand{\headfont}{\scriptsize\normalfont\sffamily}    % Kolumnentitel serifenlos
\renewcommand{\pnumfont}{\footnotesize\normalfont\ttfamily\bfseries}    % Seitennummern typewriter und fett 
\pagestyle{scrheadings}
\lehead[{\includegraphics[width=0.85cm, keepaspectratio=true]{bvlogo_sw} \hspace{3mm} \includegraphics[width=1.05cm, keepaspectratio=true]{bmlogo_sw} }]{{\includegraphics[width=0.85cm, keepaspectratio=true]{bvlogo_sw} \hspace{3mm} \includegraphics[width=1.05cm, keepaspectratio=true]{bmlogo_sw} }}
\lohead[{\includegraphics[width=0.85cm, keepaspectratio=true]{bvlogo_sw} \hspace{3mm} \includegraphics[width=1.05cm, keepaspectratio=true]{bmlogo_sw} }]{{\includegraphics[width=0.85cm, keepaspectratio=true]{bvlogo_sw} \hspace{3mm} \includegraphics[width=1.05cm, keepaspectratio=true]{bmlogo_sw}} }
%
\rehead[{\includegraphics[height=0.75cm, keepaspectratio=true]{tumlogo_out_sw} }]{{\includegraphics[height=0.75cm, keepaspectratio=true]{tumlogo_out_sw} } }
\rohead[{\includegraphics[height=0.75cm, keepaspectratio=true]{tumlogo_out_sw} }]{{\includegraphics[height=0.75cm, keepaspectratio=true]{tumlogo_out_sw} }}

\chead{\begin{spacing}{1.0}TUM Lehrstuhl f�r Baumechanik \\ \bmtitle \end{spacing}}

%\pdfminorversion=6


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% German - English
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\ifthenelse{\boolean{bmgerman}}{%

\renewcommand{\bibname}{Literatur} % Literaturverzeichnis wird zu Literatur
\renewcommand{\figurename}{Bild}   % Abbildung wird zu Bild
\renewcommand{\listfigurename}{Bildverzeichnis}%
\usepackage{ngerman}             % neue Rechtschreibung
\bibliographystyle{bmvorlagen/dinat}
}%
{\bibliographystyle{bmvorlagen/dinat_e}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pdf/hyperref/colorlinks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\pdfoutput=1                  % PDF-Ausgabe
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

\AtBeginDocument{
  \hypersetup{
    pdftitle = {\bmtitle},%
    pdfauthor = {\bmauthor},%
    pdfkeywords={\bmkeywords},%
    pdfstartview={FitH},%
    pdfstartpage={\bmstartpage},%
  }
}

