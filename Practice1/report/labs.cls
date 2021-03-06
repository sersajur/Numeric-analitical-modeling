% ***************** ���� ��� ��������� ���� ******************
% ************************ ��������� �.�. **********************
% **************** �� ����� ����� ������� �.�. ****************
% ************************* 25.05.14 ***************************
% ************************* 27.05.14 ***************************
% ************************* 29.05.14 ***************************
% *********** ����������� ������� �.�. �� 02.06.14 ************
% *********** ����������� ������� �.�. �� 02.06.14 ************
% ************************* 05.06.14 ***************************
% ************************* 08.06.14 ***************************
% ************************* 09.06.14 ***************************
% *********** ����������� ������� �.�. �� 10.06.14 ************
% *********** ����������� ������� �.�. �� 11.06.14 ************
% *********** ����������� �������� �.�. �� 12.06.14 ***********
% ************************* 17.02.15 ***************************

% �������� �������:
% \maketitlepage[�����]{�������}{����}{�������}{ϲ�}{�����}{������'������ ����}
% \annotation{���������}
% \total{��������}, �� �������� -- �� ���� � ��������� �������:
%	page -- ������� �������
%	pages -- ������� ������� ��� �������
%	equations -- ������� ������
%	figures -- ������� �������
%	tables -- ������� �������
%	bibitemcount -- ������� �������� ������ ��������
%	appendnum --  ������� �������
% \tableofcontents
% \shortings
% \intro
% \conclusion
% \begin{thebibliography} .. \end{thebibliography}
% \append{�����}

% 1. ������������ ����
\ProvidesClass{labs}[2014/06/11 v2.2 ���� ��� ��������� ���� �� ���� 3008-95]
\NeedsTeXFormat{LaTeX2e} % �������� ������

% 2. ����������� "����������� ����"
\LoadClass[a4paper,14pt,ukrainian,oneside]{extreport}


% 3. ����������� ������
% ����� ������
\RequirePackage{cmap}% ��� ��������� ������ � pdf
\RequirePackage[T2A]{fontenc} % ����� ��� �������� ��������� ������
\RequirePackage[cp1251]{inputenc} % ����� ��� �������� Windows'�������� ��������� cp1251
\RequirePackage[russian,ukrainian]{babel} % ����� ��� �������� ��������� �� �������� ���
%\RequirePackage{pscyr} % ����� �� ������� ���������� ��������
\RequirePackage{mathtext} % ������� ������� � ����� ������������ ������
\RequirePackage{amssymb,amsfonts,amsmath,amsthm} % ������� ������� � ����� ������������ ������
\RequirePackage{textcomp} % ������� ������� � ����� ������������ ������
\RequirePackage{upgreek} % ����� ��� ������ �������� ����

% �������� ������
\RequirePackage{graphicx} % �������� ����� ��� ������� �������
\RequirePackage{subfig} % ����� ��� ������� ���������
\RequirePackage{flafter} % ����� �������, �� ������� �� ������� ������ �� �'�������� ������ ��������� �� ���� ����
\RequirePackage{caption} % ����� ��� ������������� ��������� ������� �� �������
\RequirePackage{xtab} % ����� ��� ������ �������
\RequirePackage{longtable} % ����� ��� ������ �������
\RequirePackage{multirow} % ����� ��� ��'������� ����� �������

% ������ ��� ������������ ������� �� ����
\RequirePackage[includeheadfoot, headsep=0pt, top=2cm, bottom=2cm, left=3cm, right=2cm]{geometry} % ������ ���� �������
\RequirePackage{indentfirst} % ������ ��� ������� ����� ������
\RequirePackage{hanging} % ������ ��� ��� ����� ����� ������
\RequirePackage{enumitem}	% ����� ��� ������������ ������
\RequirePackage{setspace}	% ����� ��� ����, ��� ����� ���� ������ ���������� ��������
\RequirePackage{fancyhdr} % ����� �������� ��� ����������� ������� ������ �������
\RequirePackage{etoolbox} % ����� ��� ���������� �������� ������� � ����������

% ���� ��������� �� �����
\RequirePackage[rm,small,raggedright]{titlesec} % ����� ��� ������� ��������� ������ ���������
\RequirePackage[titles,subfigure]{tocloft} % ����� ��� ������� ��������� ������ �����

% ����� ��� ���������� ������������� ��������-�������� � "��������"
\RequirePackage{url}

% ����� ��� ���������� ������ ���� ��� � ���������� ������
\RequirePackage{icomma}

% ����� ��� ����������� ������������ �������� ����
\RequirePackage{listings}

% ����� ��� ������������� ������
\RequirePackage{xkeyval}

% ����� ��� ���������
\RequirePackage{totcount}

% ����� ��� ��������� ����� ����� ��� ����������� ������� � ��������
\RequirePackage{array}

% 4. ������� �������
\linespread{1.5} % �������� "1,5"
\renewcommand{\rmdefault}{ftm} % ������ ����� �� ������������ --- "Times New Roman" �� ������ pscyr

\newcommand{\setfontsize}[1]{\fontsize{#1pt}{#1pt}\selectfont} % ������� ��� ��������� ������ ������
\newcommand{\invcommas}[1]{\guillemotleft #1\guillemotright} % ������� ��� ����������� ������� �����

% ����� ������������ ���� ���� � ������, �� � �������� �����
\urlstyle{same}

% ����� ������� ��� ������� ������ �������
\fancypagestyle{plain}{
    \fancyhf{}
    \rhead{\thepage}
		\renewcommand{\headrulewidth}{0pt}
		\renewcommand{\footrulewidth}{0pt}
}

% ������� ��� �������� ������� ������ � ���������� ����
\binoppenalty=10000
\relpenalty=10000

% ������������� ���������
\@ifpackageloaded{totcount}{
	%\regtotcounter{part}
	%\regtotcounter{chapter}
	%\regtotcounter{section}
	%\regtotcounter{subsection}
	%\regtotcounter{subsubsection}
	%\regtotcounter{paragraph}
	%\regtotcounter{subparagraph}
	%\regtotcounter{footnote}
	%\regtotcounter{mpfootnote}
	\regtotcounter{page}
	\newtotcounter{equations}\AtBeginEnvironment{equation}{\stepcounter{equations}}
	\newtotcounter{figures}\AtBeginEnvironment{figure}{\stepcounter{figures}}
	\newtotcounter{tables}\AtBeginEnvironment{table}{\stepcounter{tables}}
	\newtotcounter{bibitemcount}\pretocmd\bibitem{\stepcounter{bibitemcount}}{}{}
}{}

% ������� ��� ��������� �������������� ������, ���� ��� �� ����� ���������� � ����
\newcommand{\uchapter}[1]{\chapter*{#1}\addcontentsline{toc}{chapter}{#1}} 

% ��� ��������� ������ ����������� ������� \uchapter
\newcommand{\intro}{\uchapter{�����}}

% ��� ��������� �������� ����������� ������� \uchapter
\newcommand{\conclusion}{\uchapter{��������}}

% ��� ��������� ������ ������� ��������� ����������� ������� \shortings
\newcommand{\shortings}{\uchapter{������ ������� ���������, ��������� � �������}}

% ����� ������� ��� ��������� �������
\newcounter{pages}\@ifpackageloaded{totcount}{\regtotcounter{pages}}{}
\newcounter{appendnum}\@ifpackageloaded{totcount}{\regtotcounter{appendnum}}{}
\newcommand{\append}[1]{
\titleformat{\chapter}{\centering\normalsize}{\hspace{1cm}\thechapter}{0.5em}{}
\titlespacing{\chapter}{0pt}{0pt}{0pt}
\ifnum\value{appendnum}=1 \setcounter{pages}{\value{page}}\fi
\chapter*{������� \Asbuk{appendnum}}
\begin{center}
#1
\end{center}
\titlespacing{\chapter}{0pt}{0pt}{49pt}
\titleformat{\chapter}{\centering\normalsize\MakeUppercase}{\hspace{1cm}\thechapter}{0.5em}{}
\addcontentsline{toc}{chapter}{������� \Asbuk{appendnum}\ #1}
\setcounter{figure}{0}
\setcounter{table}{0}
\renewcommand{\thefigure}{\Asbuk{appendnum}.\arabic{figure}}
\renewcommand{\thetable}{\Asbuk{appendnum}.\arabic{table}}
}
% ��������� ���������� ����� � ���������� �������
\pretocmd\append{
	\stepcounter{appendnum}
	\ifnum\Asbuk{appendnum}=� \stepcounter{appendnum}\fi
	\ifnum\Asbuk{appendnum}=� \stepcounter{appendnum}\fi
	\ifnum\Asbuk{appendnum}=� \stepcounter{appendnum}\fi
	\ifnum\Asbuk{appendnum}=� \stepcounter{appendnum}\fi
	\ifnum\Asbuk{appendnum}=� \stepcounter{appendnum}\fi
	\ifnum\Asbuk{appendnum}=� \stepcounter{appendnum}\fi
	\ifnum\Asbuk{appendnum}=� \stepcounter{appendnum}\fi
	\ifnum\Asbuk{appendnum}=� \stepcounter{appendnum}\fi
	\ifnum\Asbuk{appendnum}=� \stepcounter{appendnum}\fi
}{}{}

% ������ �������� ���� ���������
% �������� --- ��� ������
\titleformat{\chapter} % �������, �� ���������� ���� �����
			{\centering\normalsize\MakeUppercase} % ������� ������ ����� (������, "�� �����")
			{\hspace{1cm}\thechapter} % ������� ������ ������ ������: �� ���� ������ �����, ��� ������
			{0.5em} % ������� �� ������� � ������
			{} % �����, �� ������ ����
	
% ����� --- ��� ��������
\titleformat{\section}[block] % �������, �� ���������� ���� �������
			{\normalsize} % ������� ������ ����� (������)
			{\hspace{1.25cm}\thesection} % ������� ������ ������ ������: �� ���� ������ �����, ��� ������
			{0.5em} % ������� �� ������� � ������
			{} % �����, �� ������ ����

% ����� --- ��� �������� �� ������� ����
\titleformat{\subsection}[block] % �������, �� ���������� ���� ������� �������� ����
			{\normalsize} % ������� ������ ����� (������)
			{\hspace{1.25cm}\thesubsection} % ������� ������ ������ ������: �� ���� ������ �����, ��� ������
			{0.5em} % ������� �� ������� � ������
			{} % �����, �� ������ ����
			
% ����� ������ ������� ��� ��������� ������
\titlespacing{\chapter} % �������, �� ���������� ���� �����
			{0pt} % ������ ����
			{0pt} % ������ ������ (����)
			{42pt} % ������ ����� (2 �����)

% ����� ������ ������� ��� ��������� ��������			
\titlespacing{\section} % �������, �� ���������� ���� �������
			{0pt} % ������ ����
			{42pt} % ������ ������ (2 �����)
			{42pt} % ������ ����� (2 �����)
			
% ����� ������ ������� ��� ��������� �������� �� ������� ����
\titlespacing{\subsection} % �������, �� ���������� ���� ������� �������� ����
			{0pt} % ������ ����
			{42pt} % ������ ������ (2 �����)
			{42pt} % ������ ����� (2 �����)

% ����� ������� ���, ��� �� ����� ����������� ���� �������� �������
\preto\section{%
	\ifnum\value{section}=0\addvspace{-42pt}\fi%
}

\preto\subsection{%
	\ifnum\value{subsection}=0\addvspace{-42pt}\fi%
}
			
% ����������� ����
\makeatletter
\renewcommand{\tableofcontents}{
	\clearpage
	\thispagestyle{empty}
	\linespread{1.48}

	\begin{center}
	  ����
	\end{center}
	
	\@starttoc{toc}

	\linespread{1.5}
}
\makeatother

% ��������� ������ ����� � �������
\renewcommand{\arraystretch}{1.5}

% ������� ������ ��������� ������
\renewcommand{\theequation}{\arabic{chapter}.\arabic{equation}}

% ������ �������� ������, � ����������� ���������
% ����������� ������� ����
\AddEnumerateCounter{\asbuk}{\@asbuk}{\cyrm}
\setlist[itemize,1]{label=-, partopsep=0pt, parsep=0pt, topsep=0pt, itemsep=0pt, itemindent=1.7cm, leftmargin=0pt}
\setlist[enumerate,1]{label=\asbuk*), partopsep=0pt, parsep=0pt, topsep=0pt, itemsep=0pt, itemindent=1.9cm, leftmargin=0pt}
\setlist[itemize,2]{label=-, partopsep=0pt, parsep=0pt, topsep=0pt, itemsep=0pt, itemindent=2.6cm, leftmargin=0pt}
\setlist[enumerate,2]{label=\arabic*), partopsep=0pt, parsep=0pt, topsep=0pt, itemsep=0pt, itemindent=2.6cm, leftmargin=0pt}

% ����� ��������� ������
% ������ ������� �� �������
\renewcommand{\cftchapdotsep}{1}
\renewcommand{\cftsecdotsep}{1}
\renewcommand{\cftsubsecdotsep}{1}

% ������ ����� ����, �� ���������� �� ����� �������
\cftsetpnumwidth{1em}

% ������ ������� �� ������� � ����
\setlength\cftbeforechapskip{0pt}

% ��������� ��������
\renewcommand{\cftchapaftersnum}{} % �� ���� �������� ���� ������ ������ (������)
\renewcommand{\cftchapindent}{0cm} % ������ ������ ������ �� ����� ����
\renewcommand{\cftchapnumwidth}{1em} % ����� ����, �� ���������� �� ����� ������
\renewcommand{\cftchappagefont}{\mdseries} % ����� ������ ������� (���������, 14pt)
\renewcommand{\cftchapfont}{\mdseries} % ����� ��������� (���������)
\renewcommand{\cftchapleader}{\cftdotfill{\cftchapdotsep}} % ��� ����������� ������� �� ��������� �� ������ ������� (�������)

% ����� ��������� ����������
\renewcommand{\cftsecaftersnum}{} % �� ���� �������� ���� ������ �������� (������)
\renewcommand{\cftsecindent}{1em} % ������ ������ �������� �� ����� ����
\renewcommand{\cftsecnumwidth}{2em} % ����� ����, �� ���������� �� ����� ��������
\renewcommand{\cftsecpagefont}{\mdseries} % ����� ������ ������� (���������, 14pt)
\renewcommand{\cftsecfont}{\mdseries} % ����� ��������� (���������)
\renewcommand{\cftsecleader}{\cftdotfill{\cftsecdotsep}} % ��� ����������� ������� �� ��������� �� ������ ������� (�������)

% ����� ��������� ���������� �� ������� ����
\renewcommand{\cftsubsecaftersnum}{} % �� ���� �������� ���� ������ �������� �������� ���� (������)
\renewcommand{\cftsubsecindent}{3em} % ������ ������ �������� �������� ���� �� ����� ����
\renewcommand{\cftsubsecnumwidth}{2.5em} % ����� ����, �� ���������� �� ����� �������� �������� ����
\renewcommand{\cftsubsecpagefont}{\mdseries} % ����� ������ ������� (���������, 14pt)
\renewcommand{\cftsubsecfont}{\mdseries} % ����� ��������� (���������)
\renewcommand{\cftsubsecleader}{\cftdotfill{\cftsubsecdotsep}} % ��� ����������� ������� �� ��������� �� ������ ������� (�������)

% ��������� ������� � �������
\DeclareCaptionLabelFormat{tablecap}{������� #2} % ������ �� ������������ "����." �� ��������
\DeclareCaptionLabelFormat{figurecap}{������� #2} % ������ �� ������������ "���." �� ��������
\DeclareCaptionFormat{listing}{˳����� ����� \lstname\ --- #3} % ������ �� ������������ "Listing" �� ��������

\captionsetup[table]{margin=0pt,singlelinecheck=false, justification=justified, labelsep=endash, labelformat=tablecap, font = {stretch = 1.5}} % ������ ��� �������
\captionsetup[figure]{justification=centering, labelsep=endash, labelformat=figurecap, font = {stretch = 1.5}} % ������ ��� �������
\captionsetup[lstlisting]{margin=0pt, format=listing, singlelinecheck=false, font = {stretch = 1.5}} % ������ ��� �������

% �������� ���� ��� ����������� ������� � ��������
\makeatletter
\@ifpackageloaded{array}{
	\newcolumntype{L}[1]{>{\raggedright\let\newline\\\arraybackslash}p{#1}}
	\newcolumntype{C}[1]{>{\centering\let\newline\\\arraybackslash}p{#1}}
	\newcolumntype{R}[1]{>{\raggedleft\let\newline\\\arraybackslash}p{#1}}
}{}
\makeatother

% �������, �� ����� �������� ����� �� �������� � ����������
\renewcommand{\thesubfigure}{\asbuk{subfigure}}

% ���� �������� ��� ������ �������
\renewenvironment{table}[4]{
	\begin{center}
	\begin{longtable}{#1}
	\kill
	\caption{#2}\label{#3}\\
	#4
	\endfirsthead
	\kill
	\caption*{����������� ������� \thetable}\\
	#4
	\endhead
	\hline
	\endfoot
}
	{%
	\end{longtable}
	\end{center}	
}
\setlength{\LTpre}{-14pt}
\setlength{\LTpost}{-14pt}
\setlength{\LTcapwidth}{0pt}

% �������� ������� �� ��������� � ����������
\patchcmd{\thebibliography}
  {\advance\leftmargin\labelsep}
  {\leftmargin=0pt\itemindent=1.9cm\advance\itemindent\labelsep}
  {}{}

% ������� ������ ������ ��� ������ ������ (����� + ������)
\makeatletter
\renewcommand{\@biblabel}[1]{#1.}
\makeatother

% ������� ���������� ��������, ������� ������� ������ ������ ������ � ����
\let\tempbib\thebibliography
\let\tempbibend\endthebibliography
\renewenvironment{thebibliography}{
	\renewcommand{\bibname}{����˲� ��������}% ������� ���������
	\tempbib{99.}% ����������� �������� ����� ������� --- 99
	% ��� �������� ������� �� ���������
	\setlength\partopsep{0pt}
	\setlength\parskip{0pt}
	\setlength\parsep{0pt}
	\setlength\topsep{0pt}
	\setlength\itemsep{0pt}
	\addcontentsline{toc}{chapter}{������ ��������}% ������ � ����
	}
	{
	\tempbibend
}

% ���� ������� ��� ���������� ������
% �������, ��������� ������������
\newcommand{\maketitlepage}[6][]{
	  \thispagestyle{empty}
	  \setlength{\parindent}{0cm}
		\begin{center} % �����������
			\MakeUppercase{̳���������� ����� � ����� ������}

			\MakeUppercase{������������ ��������� ����������� ������}
			
			\MakeUppercase{\invcommas{�������� ������������ ��������}}

			\addvspace{6pt} % ��������� ������
			
			������� ��������� ����������

			\vfill 
			
			\setfontsize{14pt}

			\MakeUppercase{��������� ������  \ifstrempty{#1}{}{\No~{#1}}}

			�� ���������� \invcommas{#2}

			\MakeUppercase{#3}
			
		\end{center}

		\vfill % ��������� ������
		
		{\renewcommand{\arraystretch}{1.0}
		\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lr}
			& ��������:\\
			& �������� ����� #6\\
			#5
			&\\
			& �������:\\
			& ������ ������� ��������� ����������,\\
			& �������� ��������� ����\\
			& {#4}\\
			&\\
 			& ������ �� ������ \_\_\_\_\_\_\_\_\_\_\_\_\_\\
		\end{tabular*}
		}
		
		\vfill % ����������� ���������� �������

		\begin{center}
			��� --- \the\year\
		\end{center}
		
		\clearpage % �������� � ���� �������

		\setcounter{page}{2} % ����� ������� ���� 2
		\pagestyle{plain}
		
		\setlength{\parindent}{1.25cm} % ������ ������ ������� ����� ������ (5 �����, ��� 1,25 ��)
		
		% ������� ����� �� �� ���� �������
		\setlength{\belowdisplayskip}{14pt}
		\setlength{\abovedisplayskip}{14pt}
		\setlength{\belowdisplayshortskip}{14pt}
		\setlength{\abovedisplayshortskip}{14pt}
}

% ���� ������� ��� ��������
\newcommand{\annotation}[1]{\chapter*{#1}\thispagestyle{empty}\renewcommand{\baselinestretch}{1.5}\setfontsize{14pt}}

%����������� ������� � ����������� � ��������
\makeatletter
\def\latex@@rquote{^\bgroup\prim@s} 
\let\hang@@h@ngrquote\h@ngrquote
\renewcommand{\h@ngrquote}{
	\relax\ifmmode
		\expandafter\latex@@rquote
	\else
		\expandafter\hang@@h@ngrquote
	\fi}
\makeatother