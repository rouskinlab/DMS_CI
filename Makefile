DOCKER_IMAGE := ydmt/dreem
VERSION := $(shell git describe --always --dirty --long)
PYPI_TOKEN := $(shell cat ~/.pypi_token.txt)

default: 
	python setup.py install

push_to_pypi:
	rm -fr dist
	python setup.py sdist
	twine upload -r pypi dist/* --user __token__ --password $(PYPI_TOKEN)
	rm -fr dist

test:
	pytest dms_ci/test.py -v