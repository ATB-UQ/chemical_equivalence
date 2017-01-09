make_pdfs:
	for file in $$(ls *.svg); do \
	echo $$file; \
	case $$file in \
	*_corrected.svg) \
	inkscape --file=$$file --export-pdf=$${file/.svg/.pdf} --export-latex ; \
	;; \
	*) \
	inkscape --file=$$file --export-pdf=$${file/.svg/.pdf}; \
	;;
	esac
.PHONY: make_pdfs

clean:
	-@rm *_corrected.svg
.PHONY: clean
