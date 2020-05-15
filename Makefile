profile:
	python3 -m cProfile -s 'time' clipkit-runner.py  tests/integration/samples/EOG092C0CZK_aa_aln.fasta -o output/test > output/large_profile.txt

run:
	python3 -m clipkit-runner tests/integration/samples/12_YIL115C_Anc_2.253_aa_aln.fasta -o output/test --log

run.simple:
	python3 -m clipkit-runner tests/integration/samples/simple.fa -o output/simple --log	

install:
	# install so clipkit command is available in terminal
	python setup.py install

develop:
	# https://setuptools.readthedocs.io/en/latest/setuptools.html#development-mode
	python setup.py develop

test: test.unit test.integration

test.unit:
	python -m pytest -m "not integration"

test.integration:
	rm -rf output/
	mkdir output/
	python -m pytest --basetemp=output -m "integration"

test.fast:
	python -m pytest -m "not (integration or slow)"
	rm -rf output/
	mkdir output/
	python -m pytest --basetemp=output -m "integration and not slow"

# used by GitHub actions during CI workflow
test.coverage: coverage.unit coverage.integration

coverage.unit:
	python -m pytest --cov=./ -m "not integration" --cov-report=xml:unit.coverage.xml

coverage.integration:
	rm -rf output/
	mkdir output/
	python -m pytest --basetemp=output --cov=./ -m "integration" --cov-report=xml:integration.coverage.xml
