lcmd=latex -interaction=batchmode -halt-on-error
bcmd=bibtex
pspdf=ps2pdf -dMaxSubsetPct=100 -dCompatibilityLevel=1.4 -dSubsetFonts=true -dEmbedAllFonts=true -dAutoFilterColorImages=false -dAutoFilterGrayImages=false -dColorImageFilter=/FlateEncode -dGrayImageFilter=/FlateEncode -dMonoImageFilter=/FlateEncode 
#pspdf=ps2pdf

all: main

main:
	@echo
	@echo
	@echo "[1;32m################################### main #######################################[0m"
	$(lcmd)  main
	bibtex   main
	$(lcmd)  main
	dvips    main.dvi
	$(pspdf) main.ps

cl:
	rm *.aux *.log *.bbl *.blg *.out *.toc *.tit *.spl *.ps *.dvi
