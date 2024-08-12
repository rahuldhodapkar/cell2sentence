## Makefile
#
# Development scripts for cell2sentence
#
# @author Rahul Dhodapkar <rahul.dhodapkar@yale.edu>
#

.PHONY: install test lint autopep build distribute-test distribute

install:
	echo "Installing locally using setup.py"
	python -m pip install -e .

test:
	pytest src/cell2sentence/tests

lint:
	pylint cell2sentence

autopep:
	autopep8 --recursive --in-place --aggressive --aggressive src

build:
	python -m build

distribute-test:
	python -m twine upload --repository testpypi dist/*

distribute:
	python -m twine upload dist/*
