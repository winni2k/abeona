ENV = pipenv run
PYTHON = $(ENV) python
PYTEST = $(ENV) pytest

test:
	tox

setup:
	scripts/setup

benchmark:
	$(PYTEST) --benchmark-autosave benchmark

dist:
	$(PYTHON) setup.py sdist
	$(PYTHON) setup.py bdist_wheel

.PHONY: test dist benchmark
