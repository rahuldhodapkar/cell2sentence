## Makefile
#
# Development scripts for cell2sentence
#
# @author Rahul Dhodapkar <rahul.dhodapkar@yale.edu>
#

.PHONY: install test lint auto-fix-lint

install:
	echo "Installing locally using setup.py"
	python -m pip install -e .

test:
	pytest src/cell2sentence/tests

lint:
	pylint cell2sentence

auto-fix-lint:
	autopep8 --recursive --in-place --aggressive --aggressive src