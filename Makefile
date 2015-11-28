all: 	py2

py2:
	rm -rf build/ *.so
	python2 setup.py build_ext --inplace

py3:
	rm -rf build/ *.so
	python3 setup.py build_ext --inplace

