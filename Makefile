PYTHON_BIN_DIR = /usr/local/python35/bin

PYTHON_EXEC = PYTHONPATH=$(PYTHONPATH) $(PYTHON_BIN_DIR)/python3

vimdiff: refactor.log
	vimdiff reference.log $<

refactor.log:
	make test > $@
.PHONY: refactor.log

test:
	$(PYTHON_EXEC) test.py
.PHONY : test

errors:
	PYTHONPATH=$(PYTHONPATH) $(PYTHON_BIN_DIR)/pylint -E $$(find . -name '*.py')
.PHONY: errors

mypy: $(PYTHON_BIN_DIR)/mypy
	MYPYPATH=$(PYTHONPATH) $(PYTHON_BIN_DIR)/mypy calcChemEquivalency.py
.PHONY: mypy

timing_table_rows:
	make test | egrep '&' 2>&1
.PHONY: timing_table_rows
