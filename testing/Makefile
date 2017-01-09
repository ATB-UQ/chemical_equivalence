make_pdfs:
	for file in $$(ls *.svg); do \
	echo $$file; \
	inkscape --file=$$file --export-pdf=$${file/.svg/.pdf} --export-latex ; \
	done
.PHONY: make_pdfs

clean:
	-@rm *_corrected.svg
.PHONY: clean
