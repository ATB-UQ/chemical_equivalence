DREADNAUT_SCRIPT = dreadnaut_script.txt

python:
	python NautyInterface.py 

test:
	cd testing && \
	cat $(DREADNAUT_SCRIPT) | dreadnaut
