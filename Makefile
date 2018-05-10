SHELL=bash

all: MIDAS

MIDAS: 
	-rm -rf build/*
	python setup_no_fort.py install --home=$(INSTALL_PATH)


