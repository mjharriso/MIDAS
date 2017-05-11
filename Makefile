SHELL=bash

all: MIDAS

MIDAS: 
	-rm -rf build/*
	python setup.py install --home=$(INSTALL_PATH)


