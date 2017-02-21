# Build wrapper for testing.
# Python module dependencies for testing: numpy, scipy, matplotlib, unittest2

# distutils is br0ken; contents of builddir will ALL be removed by our clean target 

builddir=./tmp
libdir=./install

#

PYTHON=python
SETUPOPTS=build --build-base=${builddir} install_lib --install-dir=${libdir}

NOSETESTS=nosetests

.PHONY: all steps clean test test-unit test-validation pythonpath

all: steps

steps: 
	mkdir -p '${libdir}'
	env PYTHONPATH="${abspath ${libdir}}:$$PYTHONPATH" ${PYTHON} setup.py ${SETUPOPTS}

pythonpath:
	@echo ${abspath ${libdir}}:$$PYTHONPATH

test: test-unit test-validation
	$(MAKE) test-reports/test-results.xml

test-unit: steps
	env PYTHONPATH="${abspath ${libdir}}:$$PYTHONPATH" ${NOSETESTS} -v --with-xunit --xunit-file=${abspath test-reports}/nose-unit-tests.xml --all-modules -w test/unit_test model_test

test-validation: steps
	env PYTHONPATH="${abspath ${libdir}}:$$PYTHONPATH" ${NOSETESTS} -v --with-xunitmp --xunitmp-file=${abspath test-reports}/nose-validation-tests.xml --all-modules --processes=-1 --process-timeout=600 -w test validation

# Allow merging of test results independently of running tests

test-reports/test-results.xml: test/test-results-in.xml \
                               test/collate-junit.xml \
                               test-reports/nose-unit-tests.xml test-reports/nose-validation-tests.xml
	xsltproc --path test-reports test/collate-junit.xml $< > $@


clean:
	rm -rf '${builddir}'

testclean:
	rm -f test-reports/*.xml

realclean: clean testclean
	rm -rf '${libdir}/steps'
