SHELL=bash
INSTALL_PATH = "/home/$(USER)/local"

all: MIDAS

MIDAS: 
	-rm -rf build/*
	python setup.py install --home=$(INSTALL_PATH)


