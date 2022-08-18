## Makefile
#
# Development scripts for cell2sentence
#
# @author Rahul Dhodapkar <rahul.dhodapkar@yale.edu>
#

.PHONY: install test

install:
	echo "Installing locally using setup.py"
	python -m pip install -e .

test:
	pytest src/cell2sentence/tests
