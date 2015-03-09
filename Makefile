DREADNAUT_SCRIPT = dreadnaut_script.txt

python:
	python calcChemEquivalency.py 

test:
	cd testing && \
	cat $(DREADNAUT_SCRIPT) | dreadnaut
