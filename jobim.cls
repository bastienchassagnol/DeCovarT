%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                           %%
%%   Copyright © 2009-2010 -- LIRMM/CNRS                                     %%
%%                           (Laboratoire d'Informatique, de Robotique et de %%
%%                            Microélectronique de Montpellier /             %%
%%                            Centre National de la Recherche Scientifique). %%
%%                                                                           %%
%% Auteur/Author: DoccY <alban.mancheron@lirmm.fr>                           %%
%%                                                                           %%
%% ------------------------------------------------------------------------- %%
%%                                                                           %%
%% Ce fichier fait partie des sources des actes de Jobim 2010.               %%
%%                                                                           %%
%% Jobim est la conférence annuelle francophone de bioinformatique.          %%
%%                                                                           %%
%% Ce logiciel est régi  par la licence CeCILL  soumise au droit français et %%
%% respectant les principes  de diffusion des logiciels libres.  Vous pouvez %%
%% utiliser, modifier et/ou redistribuer ce programme sous les conditions de %%
%% la licence CeCILL  telle que diffusée par le CEA,  le CNRS et l'INRIA sur %%
%% le site "http://www.cecill.info".                                         %%
%%                                                                           %%
%% En contrepartie de l'accessibilité au code source et des droits de copie, %%
%% de modification et de redistribution accordés par cette licence, il n'est %%
%% offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons, %%
%% seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le %%
%% titulaire des droits patrimoniaux et les concédants successifs.           %%
%%                                                                           %%
%% À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques %%
%% associés  au chargement,  à  l'utilisation,  à  la modification  et/ou au %%
%% développement  et à la reproduction du  logiciel par  l'utilisateur étant %%
%% donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à %%
%% manipuler et qui le réserve donc à des développeurs et des professionnels %%
%% avertis  possédant  des  connaissances  informatiques  approfondies.  Les %%
%% utilisateurs  sont donc  invités  à  charger  et  tester  l'adéquation du %%
%% logiciel  à leurs besoins  dans des conditions  permettant  d'assurer  la %%
%% sécurité de leurs systêmes et ou de leurs données et,  plus généralement, %%
%% à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.        %%
%%                                                                           %%
%% Le fait  que vous puissiez accéder  à cet en-tête signifie  que vous avez %%
%% pris connaissance  de la licence CeCILL,  et que vous en avez accepté les %%
%% termes.                                                                   %%
%%                                                                           %%
%% ------------------------------------------------------------------------- %%
%%                                                                           %%
%% This file is part of the proceedings source of Jobim 2010.                %%
%%                                                                           %%
%% Jobim is the yearly french-speaking conference on bioinformatics.         %%
%%                                                                           %%
%% This software is governed by the CeCILL license under French law and      %%
%% abiding by the rules of distribution of free software. You can use,       %%
%% modify and/ or redistribute the software under the terms of the CeCILL    %%
%% license as circulated by CEA, CNRS and INRIA at the following URL         %%
%% "http://www.cecill.info".                                                 %%
%%                                                                           %%
%% As a counterpart to the access to the source code and rights to copy,     %%
%% modify and redistribute granted by the license, users are provided only   %%
%% with a limited warranty and the software's author, the holder of the      %%
%% economic rights, and the successive licensors have only limited           %%
%% liability.                                                                %%
%%                                                                           %%
%% In this respect, the user's attention is drawn to the risks associated    %%
%% with loading, using, modifying and/or developing or reproducing the       %%
%% software by the user in light of its specific status of free software,    %%
%% that may mean that it is complicated to manipulate, and that also         %%
%% therefore means that it is reserved for developers and experienced        %%
%% professionals having in-depth computer knowledge. Users are therefore     %%
%% encouraged to load and test the software's suitability as regards their   %%
%% requirements in conditions enabling the security of their systems and/or  %%
%% data to be ensured and, more generally, to use and operate it in the same %%
%% conditions as regards security.                                           %%
%%                                                                           %%
%% The fact that you are presently reading this means that you have had      %%
%% knowledge of the CeCILL license and that you accept its terms.            %%
%%                                                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% $Id: jobim2010.cls,v 0.8 2010-06-28 13:33:10 doccy Exp $
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% $Log: jobim2010.cls,v $
%% Revision 0.8  2010-06-28 13:33:10  doccy
%% Update of class version.
%%
%% Revision 0.7  2010-06-28 13:29:40  doccy
%% Bug Fix.
%%
%% Revision 0.6  2010-03-19 20:14:24  doccy
%% Finally, only abstract is optional for short papers.
%% Keywords are mandatory for both long and short contributions.
%%
%% Revision 0.5  2010-03-19 08:55:36  doccy
%% Abstract and keywords are now mandatory only for long papers.
%% Reporting this modification on the templates.
%% More details about:
%% * LaTeX Class:
%%   - add long/short class option and some funny macros.
%%
%% * OpenOffice/RTF Templates:
%%   - cleaning styles.
%%
%% Revision 0.4  2010-03-15 10:34:09  doccy
%% Removing the shadow for french title/subtitle/abstract and keywords (since it seems to be to difficult to be done efficiently under MS Word, and can't be done with RTF).
%%
%% Revision 0.3  2010-03-05 19:11:05  doccy
%% Corrigendum about file encoding.
%%
%% Revision 0.2  2010-03-05 19:05:19  doccy
%% Update LaTeX & RTF templates according to the
%% recommandations of O. GASCUEL.
%%
%% Revision 0.1  2010-03-03 11:44:21  doccy
%% Adding LaTeX Class for authors.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K-B\pounds\M\N-A\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h-C\mathonesuperior(<*(B\k-B\maththreesuperior\m\n-A\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) (:#(B
%%   Dollar        -=\"a     Percent       (<&(B     Ampersand     \^O-A
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon   <  \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    (:_(B
%%   Grave accent  \`     Left brace    (<<(B     Vertical bar  (<r(B
%%   Right brace   (<>(B     Tilde         \~}
%%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{jobim}[2010/06/28 v0.8 JOBIM 2010 Proceedings class]

\newif\if@jobim@longpaper\@jobim@longpapertrue
\newif\if@jobim@keywords\@jobim@keywordstrue
\newif\if@jobim@referencepaper\@jobim@referencepaperfalse
\def\jobim@longpaper@page@limit{8}
\def\jobim@posterpaper@page@limit{1}
\def\jobim@highlight@page@limit{2}
\DeclareOption{long}{%
  \@jobim@longpapertrue%
  \@jobim@referencepaperfalse%
  \let\jobim@page@limit\jobim@longpaper@page@limit%
}
\DeclareOption{poster}{%
  \@jobim@longpaperfalse%
  \@jobim@referencepaperfalse%
  \@jobim@keywordsfalse%
  \let\jobim@page@limit\jobim@posterpaper@page@limit%
}
\DeclareOption{highlight}{%
  \@jobim@longpaperfalse%
  \@jobim@referencepapertrue%
  \let\jobim@page@limit\jobim@highlight@page@limit%
}
\newif\if@showframe\@showframefalse
\DeclareOption{showframe}{\@showframetrue}
\DeclareOption{draft}{%
  \PassOptionsToClass{\CurrentOption}{article}%
}
\DeclareOption{final}{%
  \PassOptionsToClass{\CurrentOption}{article}%
}
\ProcessOptions\relax

\LoadClass[twoside, onecolumn, a4paper, 11pt]{article}
\RequirePackage[%
  pdfstartview=Fit, pdfpagelayout=SinglePage, plainpages=false,%
  pdfpagelabels, naturalnames, frenchlinks, pdfpagemode=UseOutlines,%
  bookmarks, bookmarksnumbered, bookmarksopenlevel=1,%
  colorlinks, urlcolor=blue, linkcolor=black, hyperindex
]{hyperref}
%\RequirePackage{times}
\RequirePackage{color}
\RequirePackage{xspace}
\RequirePackage{graphicx}
\RequirePackage{fancyhdr}
\RequirePackage{fancybox}
\RequirePackage{indentfirst}
\RequirePackage[%
  margin=20mm,%
  bindingoffset=0mm,%
  a4paper,
  footskip=10mm,%
  headheight=5mm,%
  headsep=5mm,%
  marginpar=0mm%
  \if@showframe, showframe\fi%
]{geometry}
\def\geometry#1{
  \ClassWarning{jobim}{%
    You are not allowed to use \string\geometry\space with this class%
  }%
}
\RequirePackage[english]{babel}

% languages
\let\switcht@@therlang\relax
\def\ds@english{\def\switcht@@therlang{\switcht@english}}
\def\ds@francais{\def\switcht@@therlang{\switcht@francais}}
\AtBeginDocument{%
  \@ifundefined{url}{\def\url#1{#1}}{}%
  \@ifpackageloaded{babel}{%
    \@ifundefined{extrasenglish}{}{\addto\extrasenglish{\switcht@english}}%
    \@ifundefined{extrasfrench}{}{\addto\extrasfrench{\switcht@francais}}%
    \@ifundefined{extrasfrenchb}{}{\addto\extrasfrenchb{\switcht@francais}}%
    \@ifundefined{extrasfrancais}{}{\addto\extrasfrancais{\switcht@francais}}%
  }{%
    \switcht@@therlang%
  }%
}

\AtEndDocument{%
  \advance\c@page by -1
  \ifnum\the\c@page>\jobim@page@limit%
    \ClassWarningNoLine{jobim}{%
      Your paper has too many pages. You must not exceed \jobim@page@limit\space
      pages%
    }%
  \fi%
  \advance\c@page by 1
}

\def\JobimPaperMaxPages{\jobim@page@limit\xspace}
\def\JobimLongPaperMaxPages{\jobim@longpaper@page@limit\xspace}
\def\JobimPosterPaperMaxPages{\jobim@posterpaper@page@limit\xspace}
\def\JobimHighlightMaxPages{\jobim@highlight@page@limit\xspace}

%% Fonts
\def\jobim@title@font{\fontsize{14pt}{14pt}\rm\bfseries\bf\selectfont}
\def\jobim@subtitle@font{\fontsize{14pt}{14pt}\rm\selectfont}
\def\jobim@corrauthor@font{\tt\fontsize{9pt}{9pt}\selectfont}
\def\jobim@refpaper@font{\rm\fontsize{10pt}{10pt}\selectfont}
\def\jobim@institute@font{\tt\fontsize{9pt}{9pt}\selectfont}
\def\jobim@acknowledgements@font{\fontsize{10pt}{10pt}\selectfont}
\def\jobim@references@font{\fontsize{10pt}{12pt} \selectfont}

%%

\fancypagestyle{jobim@title}{
  \fancyhf{}
  \fancyhead[LO,RE]{\print@session}%
  \fancyfoot[RO,LE]{\thepage}%
  \renewcommand{\headrulewidth}{0pt}
  \renewcommand{\footrulewidth}{0pt}
}

\fancypagestyle{jobim}{
  \fancyhf{}
  \fancyhead[RO]{\footnotesize\foreignlanguage{english}{\@title}}%
  \fancyhead[LE]{\footnotesize\list@authors*}%
  \fancyhead[LO,RE]{\print@session}%
  \fancyfoot[RO,LE]{\thepage}%
  \renewcommand{\headrulewidth}{0pt}
  \renewcommand{\footrulewidth}{0pt}
}

\pagestyle{jobim}

%\frenchspacing
\widowpenalty=10000
\clubpenalty=10000

\renewcommand\@pnumwidth{2em}
\renewcommand\@tocrmarg{3.5em}

\def\switcht@english{%
  % \typeout{-> english}
  \gdef\abstractname{Abstract}
  \gdef\keywordsname{Keywords}
  \gdef\ackname{Acknowledgements}
  \gdef\andname{and}
  \gdef\lastandname{and}
  \gdef\appendixname{Appendix}
  \gdef\chaptername{Chapter}
  \gdef\claimname{Claim}
  \gdef\conjecturename{Conjecture}
  \gdef\contentsname{Table of Contents}
  \gdef\corollaryname{Corollary}
  \gdef\definitionname{Definition}
  \gdef\examplename{Example}
  \gdef\exercisename{Exercise}
  \gdef\paperreferencename{Reference paper}
  \gdef\corrauthorname{Corresponding author}
  \gdef\figurename{Fig.}
  \gdef\sectionname{Sec.}
  \gdef\indexname{Index}
  \gdef\equationname{Eq.}
  \gdef\lemmaname{Lemma}
  \gdef\listfigurename{List of Figures}
  \gdef\listtablename{List of Tables}
  \gdef\mailname{e-mail}
  \gdef\noteaddname{Note added in proof}
  \gdef\notename{Note}
  \gdef\partname{Part}
  \gdef\problemname{Problem}
  \gdef\proofname{Proof}
  \gdef\propertyname{Property}
  \gdef\propositionname{Proposition}
  \gdef\questionname{Question}
  \gdef\remarkname{Remark}
  \gdef\seename{see}
  \gdef\solutionname{Solution}
  \gdef\tablename{Tab.}
  \gdef\theoremname{Theorem}
}

% French section
\def\switcht@francais{%
  % \typeout{-> francais}
  \gdef\abstractname{R\'esum\'e}
  \gdef\ackname{Remerciements}
  \gdef\keywordsname{Mots-cl\'es}
  \gdef\andname{et}%
  \gdef\lastandname{et}%
  \gdef\appendixname{Appendice}
  \gdef\chaptername{Chapitre}%
  \gdef\claimname{Pr\'etention}%
  \gdef\conjecturename{Hypoth\`ese}%
  \gdef\contentsname{Table des mati\`eres}%
  \gdef\corollaryname{Corollaire}%
  \gdef\definitionname{D\'efinition}%
  \gdef\examplename{Exemple}%
  \gdef\exercisename{Exercice}%
  \gdef\paperreferencename{Article pr\'esent\'e}
  \gdef\corrauthorname{Auteur r\'ef\'erent}
  \gdef\figurename{Fig.}%
  \gdef\indexname{Index}
  \gdef\lemmaname{Lemme}%
  \gdef\listfigurename{Liste des figures}%
  \gdef\listtablename{Liste des tables}%
  \gdef\mailname{e-mail}
  \gdef\noteaddname{Note ajout\'ee \`a l'\'epreuve}%
  \gdef\notename{Remarque}%
  \gdef\partname{Partie}%
  \gdef\problemname{Probl\`eme}%
  \gdef\proofname{Preuve}%
  \gdef\propertyname{Propri\'et\'e}%
  \gdef\propositionname{Proposition}%
  \gdef\questionname{Question}%
  \gdef\remarkname{Remarque}%
  \gdef\seename{voir}
  \gdef\solutionname{Solution}%
  \gdef\tablename{Tab.}%
  \gdef\theoremname{Th\'eor\`eme}%
}
%

% Ragged bottom for the actual page
\def\thisbottomragged{%
  \def\@textbottom{%
    \vskip\z@ plus.0001fil%
    \global\let\@textbottom\relax%
  }%
}


%% Little trick to disable hyperref warning
\let\jobim@addcontentsline=\addcontentsline
\def\addcontentsline{\phantomsection\jobim@addcontentsline}

\renewcommand\section{\@startsection{section}{1}{\z@}%
                     {0.5ex}%
                     {0.1ex}%
                     {\normalfont\selectfont\large\bfseries\boldmath
                      \rightskip=\z@ \@plus 8em\pretolerance=10000 }}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                     {0.3ex}%
                     {0.1ex}%
                     {\normalfont\selectfont\normalsize\bfseries\boldmath
                      \rightskip=\z@ \@plus 8em\pretolerance=10000 }}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                     {-12\p@ \@plus -4\p@ \@minus -4\p@}%
                     {-0.5em \@plus -0.22em \@minus -0.1em}%
                     {\normalfont\selectfont\normalsize\bfseries\boldmath}}
\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                     {-12\p@ \@plus -4\p@ \@minus -4\p@}%
                     {-0.5em \@plus -0.22em \@minus -0.1em}%
                     {\normalfont\normalsize\itshape}}
\renewcommand\subparagraph[1]{%
  \ClassWarning{jobim}{%
  You should not use
  \string\subparagraph\space with this class%
}%
\vskip0.5cm
You should not use \verb|\subparagraph| with this class.\vskip0.5cm%
}

\setlength\leftmargini  {17\p@}
\setlength\leftmargin    {\leftmargini}
\setlength\leftmarginii  {\leftmargini}
\setlength\leftmarginiii {\leftmargini}
\setlength\leftmarginiv  {\leftmargini}
\setlength  \labelsep  {.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}

\def\@listI{\leftmargin\leftmargini
            \parsep 0\p@ \@plus1\p@ \@minus\p@
            \topsep 8\p@ \@plus2\p@ \@minus4\p@
            \itemsep0\p@}
\let\@listi\@listI
\@listi
\def\@listii {\leftmargin\leftmarginii
              \labelwidth\leftmarginii
              \advance\labelwidth-\labelsep
              \topsep    0\p@ \@plus2\p@ \@minus\p@}
\def\@listiii{\leftmargin\leftmarginiii
              \labelwidth\leftmarginiii
              \advance\labelwidth-\labelsep
              \topsep    0\p@ \@plus\p@\@minus\p@
              \parsep    \z@
              \partopsep \p@ \@plus\z@ \@minus\p@}

\renewcommand\labelitemi{\normalfont\bfseries --}
\renewcommand\labelitemii{$\m@th\bullet$}

\setlength\arraycolsep{1.4\p@}
\setlength\tabcolsep{1.4\p@}


\def\@dotsep{2}


\renewcommand\refname{R\'ef\'erences}
%\ifx\oribibl\undefined
%\ifx\citeauthoryear\undefined
\renewenvironment{thebibliography}[1]
     {\section*{\refname}
      \def\@biblabel##1{[##1]}
      \jobim@references@font
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
        \renewcommand\newblock{\hskip .11em \@plus.33em \@minus.07em}%
      \sloppy\clubpenalty4000\widowpenalty4000%
      \sfcode`\.=\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\def\@lbibitem[#1]#2{\item[{[#1]}\hfill]\if@filesw
     {\let\protect\noexpand\immediate
     \write\@auxout{\string\bibcite{#2}{#1}}}\fi\ignorespaces}
\newcount\@tempcntc
\def\@citex[#1]#2{\if@filesw\immediate\write\@auxout{\string\citation{#2}}\fi
  \@tempcnta\z@\@tempcntb\m@ne\def\@citea{}\@cite{\@for\@citeb:=#2\do
    {\@ifundefined
       {b@\@citeb}{\@citeo\@tempcntb\m@ne\@citea\def\@citea{,}{\bfseries
        ?}\@warning
       {Citation `\@citeb' on page \thepage \space undefined}}%
    {\setbox\z@\hbox{\global\@tempcntc0\csname b@\@citeb\endcsname\relax}%
     \ifnum\@tempcntc=\z@ \@citeo\@tempcntb\m@ne
       \@citea\def\@citea{,}\hbox{\csname b@\@citeb\endcsname}%
     \else
      \advance\@tempcntb\@ne
      \ifnum\@tempcntb=\@tempcntc
      \else\advance\@tempcntb\m@ne\@citeo
      \@tempcnta\@tempcntc\@tempcntb\@tempcntc\fi\fi}}\@citeo}{#1}}
\def\@citeo{\ifnum\@tempcnta>\@tempcntb\else
               \@citea\def\@citea{,\,\hskip\z@skip}%
               \ifnum\@tempcnta=\@tempcntb\the\@tempcnta\else
               {\advance\@tempcnta\@ne\ifnum\@tempcnta=\@tempcntb \else
                \def\@citea{--}\fi
      \advance\@tempcnta\m@ne\the\@tempcnta\@citea\the\@tempcntb}\fi\fi}

\def\idxquad{\hskip 10\p@}% space that divides entry from number

\def\@idxitem{\par\hangindent 10\p@}

\def\subitem{\par\setbox0=\hbox{--\enspace}% second order
                \noindent\hangindent\wd0\box0}% index entry

\def\subsubitem{\par\setbox0=\hbox{--\,--\enspace}% third
                \noindent\hangindent\wd0\box0}% order index entry

\def\indexspace{\par \vskip 10\p@ plus5\p@ minus3\p@\relax}


\renewcommand\footnoterule{%
  \kern-3\p@
  \hrule\@width 2truecm
  \kern2.6\p@}
  \newdimen\fnindent
  \fnindent1em
\long\def\@makefntext#1{%
    \parindent \fnindent%
    \leftskip \fnindent%
    \noindent
    \llap{\hb@xt@1em{\hss\@makefnmark\ }}\ignorespaces#1}

\long\def\jobim@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{{\bfseries #1.} #2}%
  \ifdim \wd\@tempboxa >\hsize
    {\bfseries #1.} #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}

\def \@floatboxreset {%
        \reset@font
        \fontsize{10pt}{12pt} \selectfont
        %\small
        \@setnobreak
        \@setminipage
}

\def\fps@figure{htbp}
\def\fnum@figure{{\figurename\thinspace\thefigure}}
\renewenvironment{figure}
               {\setlength\abovecaptionskip{0\p@}%
                \setlength\belowcaptionskip{0\p@}%
                \@float{figure}}
               {\end@float}
\renewenvironment{figure*}
               {\setlength\abovecaptionskip{0\p@}%
                \setlength\belowcaptionskip{0\p@}%
                \@dblfloat{figure}}
               {\end@dblfloat}

\def\fps@table{htbp}
\def\fnum@table{{\tablename\thinspace\thetable}}
\renewenvironment{table}
               {\setlength\abovecaptionskip{0\p@}%
                \setlength\belowcaptionskip{0\p@}%
                \@float{table}}
               {\end@float}
\renewenvironment{table*}
               {\setlength\abovecaptionskip{0\p@}%
                \setlength\belowcaptionskip{0\p@}%
                \@dblfloat{table}}
               {\end@dblfloat}

\long\def\@caption#1[#2]#3{%
  \par%
  \addcontentsline{\csname ext@#1\endcsname}{#1}%
     {\protect\numberline{\csname the#1\endcsname}{\ignorespaces #2}}%
  \begingroup%
  \@parboxrestore%
  \@makecaption{\csname fnum@#1\endcsname}{\ignorespaces #3}\par%
  \endgroup%
}

\newenvironment{acknowledgements}{\section*{\ackname}\jobim@acknowledgements@font}{}

% LaTeX does not provide a command to enter the authors institute
% addresses. The \institute command is defined here.

\newcounter{@inst}
\newcounter{@auth}
\newdimen\instindent
\newbox\authrun
\newtoks\authorrunning
\newtoks\tocauthor
\newbox\titrun
\newtoks\titlerunning
\newtoks\toctitle

\def\clearheadinfo{%
  \author{%
    No Author Given (please use \texttt{\string\author} in the preamble of your
    document)%
  }%
  \title{%
    No Title Given (please use \texttt{\string\title} in the preamble of your
    document)%
  }%
  \subtitle{}%
  \session{}%
  \if@jobim@longpaper%
    \abstract{%
      No Abstract Given (please use \texttt{\string\abstract} in the preamble of your
      document)%
    }%
  \else%
    \abstract{}%
  \fi%
  \keywords{%
    No Keywords Given (please use \texttt{\string\keywords} in the preamble of your
    document)%
  }%
  \motscles{%
    Pas de mots-cl\'es fourni (veuillez utiliser \texttt{\string\motscles} dans le pr\'eambule
    de votre document)%
  }%
  \resume{%
    Pas de r\'esum\'e fourni (veuillez utiliser \texttt{\string\resume} dans le pr\'eambule
    de votre document)%
  }%
  \institute{%
    No Institute Given (please use \texttt{\string\institute} in the preamble of your
    document)%
  }%
  \gdef\@thanks{}%
  \global\titlerunning={}%
  \global\authorrunning={}%
  \global\toctitle={}%
  \global\tocauthor={}%
}

\def\institute#1{\gdef\@institute{#1}}

\def\institutename{\par
 \begingroup
 \parskip=\z@
 \parindent=\z@
 \setcounter{@inst}{1}%
 \def\and{\par\stepcounter{@inst}%
 \noindent$^{\the@inst}$\enspace\ignorespaces}%
 \setbox0=\vbox{\def\thanks##1{}\@institute}%
 \ifnum\c@@inst=1\relax
   \gdef\fnnstart{0}%
 \else
   \xdef\fnnstart{\c@@inst}%
   \setcounter{@inst}{1}%
   \noindent$^{\the@inst}$\enspace
 \fi
 \ignorespaces
 \@institute\par
 \endgroup}

\def\@fnsymbol#1{%
  \ensuremath{%
    \ifcase#1%
    \or\star%
    \or{\star\star}%
    \or{\star\star\star}%
    \or\dagger%
    \or\ddagger%
    \or\mathchar "278%
    \or\mathchar "27B%
    \or\|%
    \or **%
    \or\dagger\dagger%
    \or\ddagger\ddagger%
    \else\@ctrerr%
    \fi%
  }%
}

\def\papername{\gdef\@papername}
\def\corresponding{\gdef\@corresponding}

\def\inst#1{\unskip$^{#1}$}
\def\fnmsep{\unskip$^,$}
\def\email#1{{\tt#1}}

\def\subtitle{\gdef\@subtitle}

\def\list@authors{%
  \@ifstar{%
    \begingroup%
    \let\inst=\@gobble%
    \list@@authors%
    \endgroup
  }{%
    \list@@authors%
  }%
}
\def\list@@authors{%
  \begingroup%
  \def\lastand{%
    \ifnum\value{@inst}=2\relax%
      \andname%
    \else%
      \lastandname%
    \fi%
  }%
  \def\and{%
    \stepcounter{@auth}\relax%
    \ifnum\value{@auth}=\value{@inst}%
      \lastand%
    \else%
      \unskip,%
    \fi%
    \xspace%
  }%
  \setcounter{@auth}{1}%
  \@author%
  \endgroup%
}

\def\maketitle{%
  \newpage
  \refstepcounter{section}%
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  \setcounter{figure}{0}
  \setcounter{table}{0}
  \setcounter{equation}{0}
  \setcounter{footnote}{0}%
  \begingroup
  \parindent=\z@
  \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
  \if@twocolumn
    \ClassError{jobim}{%
      You must typeset your paper with two columns%
    }%
  \else
    \@twocolumnfalse
    \newpage
    \global\@topnum\z@   % Prevents figures from going at top of page.
    \@maketitle
  \fi
  % \this{empty}
  \@thanks
  \endgroup
  \setcounter{footnote}{\fnnstart}%
}

\def\jobim@@maketitle{%
  \newpage%
  \markboth{}{}%
  \begin{center}%
    \let\newline\\%
    \begingroup
    \jobim@title@font%
    \pretolerance=10000%
    \@title\par%
    \endgroup
    \if!\@subtitle!%
    \else
      \begingroup
      \jobim@substitle@font
      \pretolerance=10000%
      \@subtitle\par%
      \endgroup
    \fi%
    \vskip .2cm%
    \setbox0=\vbox{%
      \setcounter{@auth}{1}%
      \def\and{\stepcounter{@auth}}%
      \def\thanks##1{}%
      \@author
      \par%
    }%
    \global\value{@inst}=\value{@auth}%
    \begingroup
    \lineskip .5em%
    \noindent\ignorespaces%
    \list@authors%\vskip.2cm%
    \endgroup
    \begingroup
    \jobim@institute@font
    \institutename%
    \global\value{@inst}=\value{@auth}%
    \endgroup
  \end{center}%

\jobim@corrauthor@font
\corrauthorname: \@corresponding
\selectfont\\

\if@jobim@referencepaper
  \jobim@refpaper@font
  \textbf{\emph{\paperreferencename:} \@papername}
\fi
\normalfont

  \begingroup%
  \let\inst\@gobble%
  \def\and{| }%
  \hypersetup{%
    pdftitle={\@title\ifx\@subtitle\@empty\else\space--\space\@subtitle\fi},%
    pdfauthor={\@author},%
    pdfsubject={\ifx\@session\@empty\else\@session\space--\space\jobim@session\fi},%
    pdfcreator={JOBIM class for \LaTeX},%
    pdfkeywords={\@keywords}%
  }%
  \endgroup%
  \newif\if@abstract@provided
  \if@jobim@longpaper
    \@abstract@providedtrue
  \else%
    \ifx\@abstract\@empty%
      \@abstract@providedfalse
    \else%
      \@abstract@providedtrue
    \fi%
  \fi%
  \if@abstract@provided
    \begin{list}{}{%
        \advance\topsep by 2.5mm\relax%
        \normalsize%
        \leftmargin=1cm%
        \labelwidth=\z@%
        \listparindent=\z@%
        \itemindent\listparindent%
        \rightmargin\leftmargin%
      }%
        \item[\hskip\labelsep\bfseries\small\abstractname]%
        \begingroup%
        \itshape%
        \@abstract%
        \endgroup%
        \vskip 2.5mm%
      \if@jobim@keywords
        \item[\hskip\labelsep\bfseries\small\keywordsname]%
          \@keywords%
      \fi
    \end{list}%
  \fi%
}

\AtBeginDocument{
  \let\@makecaption=\jobim@makecaption%
  \let\@maketitle=\jobim@@maketitle%
}

\def\abstract{\gdef\@abstract}
\def\keywords{\gdef\@keywords}
\def\resume{\gdef\@resume}
\def\motscles{\gdef\@motscles}
\def\session{\gdef\@session}

\def\print@session{%
  \if\@session\@empty%
  \else%
    [Session \jobim@session]%
  \fi%
}

\def\jobim@session{%
  \ifnum\@session>\c@@session%
    ???%
  \else%
    \ifnum\@session=0\c@@session%
      ???%
    \else%
      \expandafter\csname jobim@session\@session name\endcsname%
    \fi%
  \fi%
}

\newcounter{@session}
\def\jobim@AddSession#1#2{%
  \stepcounter{@session}%
  \expandafter\gdef\csname jobim@session\the\c@@session name\endcsname{%
    \iflanguage{french}{#1}{#2}%
  }%
}

\jobim@AddSession{%
  Algorithmique et statistiques pour l'analyse de s\'equences%
}{%
  Algorithmics and Statistics for Sequence Analysis%
}
\jobim@AddSession{%
  \'Evolution, g\'enomique comparative, phylog\'enie%
}{%
  Evolution, Comparative Genomics, Phylogeny%
}
\jobim@AddSession{%
  Structures macromol\'eculaires et leurs interactions%
}{%
  Macromolecular Structures and their Interactions%
}
\jobim@AddSession{%
  Transcriptomique, prot\'eomique,  biologie syst\'emique et r\'eseaux biologique%
}{%
  Transcriptomics, Proteomics, Systemic Biology and Biological Networks%
}
\jobim@AddSession{%
  Extraction et gestion de connaissances biologiques%
}{%
  Biological Knowledge Discovery and Management%
}

\newdimen\headlineindent             % dimension for space between
\headlineindent=1.166cm              % number and text of headings.

\setlength\arraycolsep{1.4\p@}
\setlength\tabcolsep{1.4\p@}

\def\jobim{\textsc{JOBIM}\xspace}
\def\JOBIM#1#2#3#4{\jobim\oldstylenums{#1#2#3#4}}

\bibliographystyle{unsrt}

%\newcommand{\email}[1]{{\tt #1}}

\def\@begintheorem#1#2{
   \par\noindent\bgroup{\sc #1\ #2. }\it\ignorespaces}
\def\@opargbegintheorem#1#2#3{
   \par\bgroup{\sc #1\ #2\ (#3). }\it\ignorespaces}
\def\@endtheorem{\egroup}
\def\proof{\par{\it \proofname}. \ignorespaces}
\def\endproof{{\ \vbox{\hrule\hbox{%
   \vrule height1.3ex\hskip0.8ex\vrule}\hrule
  }}\par}
\newtheorem{theorem}{\theoremname}[section]
\newtheorem{lemma}[theorem]{\lemmaname}
\newtheorem{definition}[theorem]{\definitionname}
\newtheorem{conjecture}[theorem]{\conjecturename}
\newtheorem{claim}[theorem]{\claimname}
\newtheorem{example}[theorem]{\examplename}
\newtheorem{exercise}[theorem]{\exercisename}
\newtheorem{proposition}[theorem]{\propositionname}
\newtheorem{property}[theorem]{\propertyname}
\newtheorem{question}[theorem]{\questionname}
\newtheorem{problem}[theorem]{\problemname}
\newtheorem{remark}[theorem]{\remarkname}
\newtheorem{solution}[theorem]{\solutionname}
\clearheadinfo
\endinput
