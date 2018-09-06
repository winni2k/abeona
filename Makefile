ENV = pipenv run
PYTHON = $(ENV) python
PYTEST = $(ENV) pytest
TEST_ARGS?=''

test:
	tox $(TEST_ARGS)

setup:
	scripts/setup

benchmark:
	$(PYTEST) --benchmark-autosave benchmark


check_git_dirty:
	git status --porcelain
	test -z "$$(git status --porcelain)"

clean:
	rm -rf dist

dist:
	$(PYTHON) setup.py sdist
	$(PYTHON) setup.py bdist_wheel

build: clean
	$(MAKE) dist

upload:
	$(ENV) twine upload dist/*

deploy: check_git_dirty
	$(MAKE) test TEST_ARGS='-r'
	$(MAKE) build
	$(MAKE) upload

.PHONY: test dist benchmark
