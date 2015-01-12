DREADNAUT_SCRIPT = dreadnaut_script.txt

test:
	cd testing && \
	cat $(DREADNAUT_SCRIPT) | dreadnaut | tail -n 1
