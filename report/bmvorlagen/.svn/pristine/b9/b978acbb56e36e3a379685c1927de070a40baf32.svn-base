%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kenndaten
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\NeedsTeXFormat{LaTeX2e}

\ProvidesClass{bmvorlagen/bmskript}[2012/04/27 BM Scriptum Class]
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisierung
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage{ifthen}
\newboolean{bmcolorlinks} %
\setboolean{bmcolorlinks}{false} %
\newboolean{bmshowl}\setboolean{bmshowl}{false} %
\newboolean{bmgerman}\setboolean{bmgerman}{false} %
\newboolean{bmkaros}\setboolean{bmkaros}{true}
\newcommand{\smodus}{student}
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deklarieren von Optionen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\DeclareOption{bmcolorlinks}{\setboolean{bmcolorlinks}{true}}
\DeclareOption{bmnokaros}{\setboolean{bmkaros}{false}}
\DeclareOption{bmsmodus=dozent}{\renewcommand{\smodus}{dozent}\setboolean{bmshowl}{false}}
\DeclareOption{bmsmodus=dozentnolabel}{\renewcommand{\smodus}{dozent}\setboolean{bmshowl}{false}}
\DeclareOption{bmsmodus=buch}{\renewcommand{\smodus}{buch}\setboolean{bmshowl}{false}}
\DeclareOption{bmsmodus=student}{\renewcommand{\smodus}{student}\setboolean{bmshowl}{false}}
\DeclareOption{bmshowlabels}{\setboolean{bmshowl}{true}}
\DeclareOption{german}{\setboolean{bmgerman}{true}}
%\DeclareOption{fleqn}{\relax}
%\PassOptionsToClass{fleqn}{scrbook}
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
\LoadClass{scrbook}
\usepackage{scrhack}

%-------------------------------------------------------------------------
% Packages
%-------------------------------------------------------------------------
%

\usepackage[
	format=hang,
	font={footnotesize, sf},
	labelfont={bf},
	margin=1cm,
	aboveskip=5pt,
	belowskip=5pt
	]{caption,subfig}      													% mehrzeilige Captions ausrichten; subfig: untergrafiken 


\usepackage[centertags]{amsmath} 									% AMS-Mathematik, centertags zentriert Nummer bei split
\usepackage{amssymb}             								  % zusaetzliche Symbole
\usepackage{amsfonts} 														% zus�tzlich noch weitere fonts 
\usepackage{icomma}
\usepackage{geometry}
\usepackage{bm}
\usepackage{graphicx}            									% zum Einbinden von Grafiken
\usepackage[svgnames,table,hyperref]{xcolor}
\usepackage{float}              			 						% u.a. genaue Plazierung von Gleitobjekten mit H
\renewcommand{\topfraction}{0.75}
\renewcommand{\floatpagefraction}{0.75}
\usepackage[pdftex,pstarrows]{pict2e}
\usepackage{array}
\usepackage{listings}															% Code-Einbindungsumgebung
\usepackage{courier}       												% verwende Courier statt cmtt als monospace-schrift
\usepackage{setspace}            									% Zeilenabstand einstellbar
\usepackage{longtable}														% Ermoeglicht Tabellen die �ber den Seitenumbruch gehen (s.Symbolverzeichnis) 
\usepackage{xr}																		% zum Verkn�pfung zweier Dokumente, wie z.B. in der SoilDynamics �bung auf das SoilDynamics Skript verwiesen wird.
\usepackage{lmodern}

%%%																							
% Packages aus TM2 Skript
%%%
\usepackage{paralist}
\usepackage{mdwlist}
%\usepackage{needspace}														% verschoben nach bm_core um auch in den Aufgabenbl�ttern etc. zur Verf�gung zu stehen
\usepackage{empheq}
\usepackage{textcomp}
\usepackage{tabto}
\usepackage{listliketab}
\usepackage{tabularx}
\usepackage{mathtools}
\usepackage{calc}
\usepackage{multirow}
\usepackage{booktabs}
\usepackage{extarrows}
\usepackage{xr}
\usepackage{scrlayer-scrpage}          									% Kopf und Fusszeilen-Layout passt besser zur Dokumentklasse 		
																									% KOMA-Skript (scrbook) als das Paket fancyhdr, sonst fast gleichwertig
%\usepackage{wasysym}															% Spezielle Symbole er�brigt sich bei Verwendung von lmodern
\usepackage{siunitx}                              % Einheiten
\sisetup{locale = DE}
\sisetup{quotient-mode=fraction}
\sisetup{per-mode=fraction}
                                                  %Er macht dann aus dem \per einen Bruch und auch wenn du 3/4 eingibst einen Bruch draus.

                                                  %Eigenkreisfrequenz: $\omega_E = \sqrt{\frac{k}{m}} = \SI{5}{\radian\per\second}$\\[1em]
                                                  %Eigenfrequenz: $f_E = \frac{\omega_E}{2\pi} = \SI{5/2\pi}{\hertz} = \SI{0,796}{\hertz}$\\[1em]
                                                  %Periodendauer: $T = \frac{1}{f_E} = \SI{1,26}{\second}$

\usepackage{mathtools}                            % Fallunterscheidung mit geschweiften Klammern links und rechts

%%%																							
% Packages aus SD Skript
%%%
\usepackage{physics}
\usepackage{nicefrac}
\usepackage{xspace}

%%%
% Literatur und sonstige Referenzen
%%%

\usepackage{natbib}																% braucht man f�r natdin.bst (Literaturzitate mit Namen und Jahr im Text)
\setcitestyle{square,aysep={}}             				% Eckige Klammern um Zitat
\typearea[current]{current}        								% Neuberechnung des Satzspiegels mit alten Werten nach �nderung von Zeilenabstand

%-------------------------------------------------------------------------
% Aufrufen anderer bm.sty Dateien
%-------------------------------------------------------------------------
\usepackage{bmvorlagen/bm_core}
\usepackage{bmvorlagen/bm_karos}
\setlength{\bmkarowidth}{17.0cm}
\setlength{\bmkaroheight}{23.0cm}
\usepackage{bmvorlagen/bm_skriptfeatures}

%-------------------------------------------------------------------------
% Aufrufen selber definierter Kommandos
%-------------------------------------------------------------------------
\include{bmvorlagen/bm_commands}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hauptteil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% showlabels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\ifthenelse{\boolean{bmshowl}}{%
\RequirePackage[notcite,notref]{showkeys}          % gibt die labels mit aus.
\renewcommand{\showkeyslabelformat}[1]{\begin{rotate}{70}\fbox{\normalfont\tiny\ttfamily#1}\end{rotate}}}{}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New- and Renew-Commands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seitenlayout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\typearea[current]{current}
\headheight = 1.3cm
\mathindent = 0.5cm
\renewcommand{\headfont}{\scriptsize\normalfont\sffamily}    %  						% Kolumnentitel serifenlos
\renewcommand{\pnumfont}{\footnotesize\normalfont\ttfamily}  %ttfamily\bfseries}    % Seitennummern typewriter und fett 
\pagestyle{scrheadings}

% BM Logos
%\lehead[{\includegraphics[height=0.85cm, keepaspectratio=true]{bmlogo_sr_sw} }]{{\includegraphics[height=0.85cm, keepaspectratio=true]{bmlogo_sr_sw} }}
%\lohead[{\includegraphics[height=0.85cm, keepaspectratio=true]{bmlogo_sr_sw} }]{{\includegraphics[height=0.85cm, keepaspectratio=true]{bmlogo_sr_sw} }}
\lehead[{}]{{}}
\lohead[{}]{{}}
%

% TUM Logos
\rehead[{\includegraphics[height=0.75cm, keepaspectratio=true]{TUM_Logo_schwarz_rgb_2016} }]{{\includegraphics[height=0.75cm, keepaspectratio=true]{TUM_Logo_schwarz_rgb_2016} }}
\rohead[{\includegraphics[height=0.75cm, keepaspectratio=true]{TUM_Logo_schwarz_rgb_2016} }]{{\includegraphics[height=0.75cm, keepaspectratio=true]{TUM_Logo_schwarz_rgb_2016} }}
%
\lofoot[\begin{spacing}{1.0}~ \\ \veranstaltung \\ \semester \end{spacing}]{\begin{spacing}{1.0}~ \\ \veranstaltung \\ \semester \end{spacing}}
\lefoot[\begin{spacing}{1.0}~ \\ \veranstaltung \\ \semester \end{spacing}]{\begin{spacing}{1.0}~ \\ \veranstaltung \\ \semester \end{spacing}}
%Einschub f�r pers�nliche Kopie unten mittig
%\cfoot[\begin{spacing}{1.0}~ \\ ~\\ ~\\ ~\\ ~\\ ~\\ ~\\ \centering{pers�nliche Kopie f�r }\end{spacing}]{\begin{spacing}{1.0}~ \\ ~\\ ~\\~\\ ~\\ ~\\ ~\\ \centering{pers�nliche Kopie f�r }\end{spacing}}
%
\rofoot[~ \\ \pagemark]{\begin{spacing}{1.0} ~ \\  \leftmark \\ \pagemark\end{spacing}}
\refoot[~ \\ \pagemark]{\begin{spacing}{1.0} ~ \\  \leftmark \\ \pagemark\end{spacing}}
%\pdfminorversion=6


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% German - English
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\ifthenelse{\boolean{bmgerman}}{%
\renewcommand{\bibname}{Literatur} 									% Literaturverzeichnis wird zu Literatur
\renewcommand{\figurename}{Bild}   									% Abbildung wird zu Bild
\renewcommand{\listfigurename}{Bildverzeichnis}%
\usepackage{ngerman}             										% neue Rechtschreibung
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
    pdftitle = {\veranstaltung \semester},%
    pdfauthor = {\lehrer},%
%    pdfkeywords={},%
%    pdfstartview={FitH},%
%    pdfstartpage={1},%
  }
}

