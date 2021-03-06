.PHONY: docs

pypi_package:
	git clean -xdf
	python3 setup.py sdist bdist_wheel

pypi_upload_test: pypi_package
	python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*

pypi_upload: pypi_package
	python3 -m twine upload dist/*

docs:
	rm -rf docs
	sphinx-build -c sphinx -b html . docs
	git add docs

lint:
	python3 -m pylint mpf
