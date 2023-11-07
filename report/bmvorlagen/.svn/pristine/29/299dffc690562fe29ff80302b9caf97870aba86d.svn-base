%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kenndaten
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\NeedsTeXFormat{LaTeX2e}

\ProvidesClass{bmvorlagen/bmblatt}[2012/04/27 BM Blatt Class]
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
\DeclareOption{bmsmodus=dozent}{\renewcommand{\smodus}{dozent}\setboolean{bmshowl}{true}}
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

\usepackage[
	format=hang,
	font={footnotesize, sf},
	labelfont={bf},
	margin=1cm,
	aboveskip=5pt,
	belowskip=5pt
	]{caption,subfig}      												% mehrzeilige Captions ausrichten; subfig: untergrafiken 

\usepackage[centertags]{amsmath} 								% AMS-Mathematik, centertags zentriert Nummer bei split
\usepackage{amssymb}             							  % zusaetzliche Symbole
\usepackage{amsfonts}														% ersetzt bzw inkludiert amssymb
\usepackage{icomma}
\usepackage{bm}
\usepackage{graphicx}            								% zum Einbinden von Grafiken
\usepackage[svgnames,table,hyperref]{xcolor}
\usepackage{float}               								% u.a. genaue Plazierung von Gleitobjekten mit H
\renewcommand{\topfraction}{0.75}
\renewcommand{\floatpagefraction}{0.75}
\usepackage[pdftex,pstarrows]{pict2e}
\usepackage{array}
\usepackage{listings}														% Code-Einbindungsumgebung
\usepackage{courier}       											% verwende Courier statt cmtt als monospace-schrift
\usepackage[rightbars]{changebar}								% Grauer Balken am Seitenrand
\setcounter{changebargrey}{45}
\usepackage{setspace}            								% Zeilenabstand einstellbar                							
\usepackage{longtable}													% Ermoeglicht Tabellen die �ber den Seitenumbruch gehen (s. Symbolverzeichnis)
\usepackage{fancyhdr}
\usepackage{totpages}
\usepackage[rightbars]{changebar}
\setcounter{changebargrey}{45}
%\usepackage{scrpage2}            							% Kopf und Fusszeilen-Layout passt besser zur Dokumentklasse KOMA-Skript (scrbook) als das Pake fancyhdr, sonst ziemlich gleichwertig
%\onehalfspacing 																% eineinhalbzeilig einstellen

%%%																							
% Packages aus SD Skript
%%%
\usepackage{physics}
\usepackage{nicefrac}
\usepackage{xspace}

% ---Literatur und sonstige Referenzen ---
\usepackage{natbib}															% f�r natdin.bst (Literaturzitate mit Namen und Jahr im Text)
\setcitestyle{square,aysep={}}             			% Eckige Klammern um Zitat

\typearea[current]{current}       	 						% Neuberechnung des Satzspiegels mit alten Werten nach �nderung von Zeilenabstand,etc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aufrufen anderer bm.sty Dateien
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\usepackage{bmvorlagen/bm_core}

\usepackage{bmvorlagen/bm_karos}
\setlength{\bmkarowidth}{16.0cm}
\setlength{\bmkaroheight}{25.0cm}

%\usepackage{bmvorlagen/bm_skriptfeatures}

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

%Aufz�hlungsstile
\makeatletter
\ifdefined\theblattnummer{\renewcommand{\theenumi}{\theblattnummer.\arabic{enumi}}}\\\fi
\renewcommand{\labelenumi}{\textbf{\theenumi)}}
\renewcommand{\p@enumi}{}
\renewcommand{\theenumii}{\alph{enumii}}
\renewcommand{\labelenumii}{\textbf{\theenumii)}}
\renewcommand{\p@enumi}{}
\makeatother

% Sprachunterscheidungen

\newcommand{\vorname}{\ifthenelse{\boolean{bmgerman}}{Vorname}{First Name}}
\newcommand{\nachname}{\ifthenelse{\boolean{bmgerman}}{Nachname}{Last Name}}
\newcommand{\matrikelnummer}{\ifthenelse{\boolean{bmgerman}}{Matrikelnummer}{Matr. Number}}
\newcommand{\studiengang}{\ifthenelse{\boolean{bmgerman}}{Studiengang}{fields of study}}
\newcommand{\sitzplatz}{\ifthenelse{\boolean{bmgerman}}{Sitzplatz}{seat}}
\newcommand{\aufgabewort}{\ifthenelse{\boolean{bmgerman}}{Aufgabe}{Problem}}
\newcommand{\hinweisarbeitszeit}{\ifthenelse{\boolean{bmgerman}}{Bearbeitungszeit \ref{bearbeitungszeit}~Min.}{Examination time \ref{bearbeitungszeit}~min.}}
\newcommand{\arbeitszeit}[1]{\ifthenelse{\boolean{bmgerman}}{(ca. #1~Min.)}{(approx #1~min.)}}
\newcommand{\seitevon}[2]{\ifthenelse{\boolean{bmgerman}}{Seite #1 von #2}{page #1 of #2}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seitenlayout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Seitenstile
\pagestyle{fancy}
\fancyhf{}
\fancyfoot[C]{\sffamily \footnotesize \seitevon{\thepage}{\ref{TotPages}}}
\renewcommand{\headrulewidth}{0pt}
\parindent0pt
\setlength{\mathindent}{0.5cm}
\raggedright

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% German - English
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\ifthenelse{\boolean{bmgerman}}{%
\renewcommand{\bibname}{Literatur} 								% Literaturverzeichnis wird zu Literatur
\renewcommand{\figurename}{Bild}   								% Abbildung wird zu Bild
\renewcommand{\listfigurename}{Bildverzeichnis}%
\usepackage{ngerman}             									% neue Rechtschreibung
\bibliographystyle{bmvorlagen/dinat}
}%

%{\bibliographystyle{bmvorlagen/dinat_e}}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pdf/hyperref/colorlinks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\pdfoutput=1                  										% PDF-Ausgabe
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
























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALT ALT ALT ALT ALT ALT ALT ALT ALT ALT ALT ALT ALT ALT ALT ALT ALT ALT ALT ALTALT ALT ALT ALT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seitenlayout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%\typearea[current]{current}
%\headheight = 1.3cm
%\mathindent = 0.5cm
%\renewcommand{\headfont}{\scriptsize\normalfont\sffamily}    % Kolumnentitel serifenlos
%\renewcommand{\pnumfont}{\footnotesize\normalfont\ttfamily\bfseries}    % Seitennummern typewriter und fett 
%\pagestyle{scrheadings}
%\lehead[{\includegraphics[height=0.85cm, keepaspectratio=true]{bmlogo_sr_sw} }]{{\includegraphics[height=0.85cm, keepaspectratio=true]{bmlogo_sr_sw} }}
%\lohead[{\includegraphics[height=0.85cm, keepaspectratio=true]{bmlogo_sr_sw} }]{{\includegraphics[height=0.85cm, keepaspectratio=true]{bmlogo_sr_sw} }}
%%
%\rehead[{\includegraphics[height=0.75cm, keepaspectratio=true]{tumlogo_out_sw} }]{{\includegraphics[height=0.75cm, keepaspectratio=true]{tumlogo_out_sw} } }
%\rohead[{\includegraphics[height=0.75cm, keepaspectratio=true]{tumlogo_out_sw} }]{{\includegraphics[height=0.75cm, keepaspectratio=true]{tumlogo_out_sw} }}
%%
%\lofoot[\begin{spacing}{1.0}~ \\ \veranstaltung \\ \semester \end{spacing}]{\begin{spacing}{1.0}~ \\ \veranstaltung \\ \semester \end{spacing}}
%\lefoot[\begin{spacing}{1.0}~ \\ \veranstaltung \\ \semester \end{spacing}]{\begin{spacing}{1.0}~ \\ \veranstaltung \\ \semester \end{spacing}}
%%
%\rofoot[~ \\ \pagemark]{\begin{spacing}{1.0} ~ \\  \leftmark \\ \pagemark\end{spacing}}
%\refoot[~ \\ \pagemark]{\begin{spacing}{1.0} ~ \\  \leftmark \\ \pagemark\end{spacing}}
%\pdfminorversion=6



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pdf/hyperref/colorlinks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%\AtBeginDocument{
  %\hypersetup{
    %pdftitle = {\veranstaltung \semester},%
    %pdfauthor = {\lehrer},%
%    pdfkeywords={},%
%    pdfstartview={FitH},%
%    pdfstartpage={1},%
 % }
%}