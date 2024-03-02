DOCKER_IMAGE := ydmt/dreem
VERSION := $(shell git describe --always --dirty --long)
PYPI_PASSWORD := $(shell cat pypi_pass.txt)

default: 
	python setup.py install


push_to_pypi:
	rm -fr dist
	python setup.py sdist
	twine upload -r pypi dist/* --user yvesmartindestaillades --password $(PYPI_PASSWORD)
	rm -fr dist

test:
	pytest dms_ci/test.py -v